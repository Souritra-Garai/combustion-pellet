#include "species/Aluminium.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.cuh"
#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "pde-problems/Core-Shell-Diffusion.cuh"

#include "utilities/Keyboard_Interrupt.hpp"
#include "utilities/File_Generator.hpp"

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
__global__ void initIteration();
__global__ void setUpEquations();
__global__ void solve();
__global__ void printMass();

__host__ void printMass(std::ostream &output_stream);
__host__ void printMassFractions(std::ostream &output_stream);
__host__ void printConcentrationArrays(std::ostream &output_stream);

int main(int argc, char const *argv[])
{
	double step = 0.005 / time_step;
	double *concentration_array_A, *concentration_array_B;
	cudaMalloc(&concentration_array_A, particle_num_grid_points * sizeof(double));
	cudaMalloc(&concentration_array_B, particle_num_grid_points * sizeof(double));

	allocateMemory<<<1,1>>>(
		time_step,
		particle_num_grid_points,
		core_radius,
		overall_radius,
		concentration_array_A,
		concentration_array_B
	);

	initializeCoreShellParticle<<<1, particle_num_grid_points>>>();
	cudaDeviceSynchronize();

	initIteration<<<1,1>>>();
	cudaDeviceSynchronize();

	for (size_t i = 0; i < MAX_ITER; i++)
	{
		size_t i_top = i + 
		for ()
		setUpEquations<<<1, particle_num_grid_points>>>();
		solve<<<1,2>>>();
		updateMassFractions<<<1,1>>>();
	}

	printMass<<<1,1>>>();
	cudaDeviceSynchronize();

	deallocateMemory<<<1,1>>>();
	cudaDeviceSynchronize();

    return 0;
}

__global__ void allocateMemory(double time_step, double num_grid_points, double core_radius, double overall_radius)
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

__global__ void initIteration(double temeperature)
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
}

__global__ void printMass()
{
	mass_A = diffusion_problem->getAtmMassA();
	mass_B = diffusion_problem->getAtmMassB();

	Y_A  = diffusion_problem->getMassFractionsCoreMaterial();
	Y_B  = diffusion_problem->getMassFractionsShellMaterial();
	Y_AB = diffusion_problem->getMassFractionsProductMaterial();
}