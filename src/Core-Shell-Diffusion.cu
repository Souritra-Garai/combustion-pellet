#include "species/Aluminium.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.cuh"
#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "pde-problems/Core-Shell-Diffusion.cuh"

#include "utilities/Keyboard_Interrupt.hpp"
#include "utilities/File_Generator.hpp"

#include <iostream>

#define MAX_ITER 1E6

size_t particle_num_grid_points = 1001;
double time_step = 1E-6;	// s

double temperature = 1500;	// K
double core_radius = 32.5E-6;	// m
double overall_radius = 39.5E-6;	// m

__device__ CoreShellDIffusion::Diffusion *diffusion_problem;
__device__ ArrheniusDiffusivityModel *Alawieh_diffusivity;
__device__ bool combustion_complete;
__device__ double mass_A;
__device__ double mass_B;
__device__ double Y_A;
__device__ double Y_B;
__device__ double Y_AB;

__global__ void allocateMemory(
	double time_step,
	double num_grid_points,
	double core_radius, 
	double overall_radius,
	double *concentration_array_A,
	double *concentration_array_B
);

__global__ void deallocateMemory();
__global__ void initializeCoreShellParticle();
__global__ void initIteration(double temperature);
__global__ void setUpEquations();
__global__ void solve();
__global__ void updateMassFractions();

__host__ void printMass(std::ostream &output_stream);
__host__ void printMassFractions(std::ostream &output_stream);
__host__ void printConcentrationArray(std::ostream &output_stream, double concentration_array[], size_t n);

int main(int argc, char const *argv[])
{
	double step = 0.005 / time_step;

	double *concentration_array_A, *concentration_array_B;
	cudaMalloc(&concentration_array_A, particle_num_grid_points * sizeof(double));
	cudaMalloc(&concentration_array_B, particle_num_grid_points * sizeof(double));

	double concentration_array_A_host[particle_num_grid_points], concentration_array_B_host[particle_num_grid_points];
	bool combustion_complete = false;

	allocateMemory<<<1,1>>>(
		time_step,
		particle_num_grid_points,
		core_radius,
		overall_radius,
		concentration_array_A,
		concentration_array_B
	);
	cudaDeviceSynchronize();

	FileGenerator folder;
	
	std::ofstream concentration_array_A_file = folder.getCSVFile("concentration_A");
	std::ofstream concentration_array_B_file = folder.getCSVFile("concentration_B");

	std::ofstream mass_file = folder.getCSVFile("mass");
	std::ofstream mass_fractions_file = folder.getCSVFile("mass_fractions");

	initializeCoreShellParticle<<<1, particle_num_grid_points>>>();
	cudaDeviceSynchronize();

	initIteration<<<1,1>>>(temperature);
	cudaDeviceSynchronize();

	std::cout << "Initialized Diffusion Problems\nStarting Iterations\n(Press Ctrl+C to break)\n";

	setUpKeyboardInterrupt();

	try
	{
		for (size_t i = 0; i < MAX_ITER; i++)
		{
			size_t i_top = i + step;
			for (; i < i_top; i++)
			{
				setUpEquations<<<1, particle_num_grid_points>>>();
				solve<<<1,2>>>();
				updateMassFractions<<<1,1>>>();

				cudaMemcpyFromSymbol(&combustion_complete, ::combustion_complete, sizeof(bool));
				if (combustion_complete) break;
			}

			cudaDeviceSynchronize();

			std::cout << "Iteration Completed : " << i << "\n";

			cudaMemcpy(concentration_array_A_host, concentration_array_A, particle_num_grid_points * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(concentration_array_B_host, concentration_array_B, particle_num_grid_points * sizeof(double), cudaMemcpyDeviceToHost);

			printConcentrationArray(concentration_array_A_file, concentration_array_A_host, particle_num_grid_points);
			printConcentrationArray(concentration_array_B_file, concentration_array_B_host, particle_num_grid_points);

			printMass(mass_file);
			printMassFractions(mass_fractions_file);

			if (combustion_complete) break;
		}
	}

	catch (InterruptException& e)
    {
        std::cout << "\nCaught signal " << e.S << std::endl;

		cudaDeviceSynchronize();

		cudaMemcpy(concentration_array_A_host, concentration_array_A, particle_num_grid_points * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(concentration_array_B_host, concentration_array_B, particle_num_grid_points * sizeof(double), cudaMemcpyDeviceToHost);

		printConcentrationArray(concentration_array_A_file, concentration_array_A_host, particle_num_grid_points);
		printConcentrationArray(concentration_array_B_file, concentration_array_B_host, particle_num_grid_points);

		printMass(mass_file);
		printMassFractions(mass_fractions_file);
    }

	concentration_array_A_file.close();
	concentration_array_B_file.close();
	mass_fractions_file.close();
	mass_file.close();

	cudaFree(concentration_array_A);
	cudaFree(concentration_array_B);

	deallocateMemory<<<1,1>>>();
	cudaDeviceSynchronize();

    return 0;
}

__global__ void allocateMemory(
	double time_step,
	double num_grid_points,
	double core_radius,
	double overall_radius,
	double concentration_array_A[],
	double concentration_array_B[]
)
{
	loadAluminium();
	loadNickel();
	loadNickelAlumnide();

	CoreShellParticle::initialize(
		aluminium,
		nickel,
		nickel_aluminide,
		core_radius,
		overall_radius
	);

	CoreShellDIffusion::setGridSize(num_grid_points);
	CoreShellDIffusion::setTimeStep(time_step);

	diffusion_problem = new CoreShellDIffusion::Diffusion();
	Alawieh_diffusivity = new ArrheniusDiffusivityModel(2.56E-6, 102.191E3);

	diffusion_problem->setArrayAddresses(concentration_array_A, concentration_array_B);
}

__global__ void deallocateMemory()
{
	delete Alawieh_diffusivity;
	delete diffusion_problem;

	CoreShellDIffusion::deallocate();

	unloadAluminium();
	unloadNickel();
	unloadNickelAluminide();
}

__global__ void initializeCoreShellParticle()
{
	size_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < CoreShellDIffusion::n) 
	{
		diffusion_problem->setInitialState(i);
	}
}

__global__ void initIteration(double temperature)
{
	diffusion_problem->setDiffusivity(Alawieh_diffusivity->getDiffusivity(temperature));
	combustion_complete = false;
}

__global__ void setUpEquations()
{
	size_t i = blockDim.x * blockIdx.x + threadIdx.x;

	diffusion_problem->setUpEquations(i);
}

__global__ void solve()
{
	size_t i = blockDim.x * blockIdx.x + threadIdx.x;

	diffusion_problem->solveEquations(i);
}

__global__ void updateMassFractions()
{
	diffusion_problem->updateMassFractions();

	mass_A = diffusion_problem->getAtmMassA();
	mass_B = diffusion_problem->getAtmMassB();

	Y_A  = diffusion_problem->getMassFractionsCoreMaterial();
	Y_B  = diffusion_problem->getMassFractionsShellMaterial();
	Y_AB = diffusion_problem->getMassFractionsProductMaterial();

	combustion_complete = diffusion_problem->isReactionComplete(0.01);
}

__host__ void printMass(std::ostream &output_stream)
{
	double mass_A, mass_B;

	cudaMemcpyFromSymbol(&mass_A, ::mass_A, sizeof(double));
	cudaMemcpyFromSymbol(&mass_B, ::mass_B, sizeof(double));
	
	output_stream << mass_A << ',' << mass_B << ',' << mass_A + mass_B << '\n';
}

__host__ void printMassFractions(std::ostream &output_stream)
{
	double Y_A, Y_B, Y_AB;

	cudaMemcpyFromSymbol(&Y_A,  ::Y_A,  sizeof(double));
	cudaMemcpyFromSymbol(&Y_B,  ::Y_B,  sizeof(double));
	cudaMemcpyFromSymbol(&Y_AB, ::Y_AB, sizeof(double));

	output_stream << Y_A << ',' << Y_B << ',' << Y_AB << ',' << Y_A + Y_B + Y_AB << '\n';
}

__host__ void printConcentrationArray(std::ostream &output_stream, double conc_array[], size_t n)
{
	for (size_t i = 0; i < n-1; i++) output_stream << conc_array[i] << ',';
	output_stream << conc_array[n-1] << '\n';
}