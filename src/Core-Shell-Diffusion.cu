#include "species/Aluminium.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.cuh"
#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "pde-problems/Core-Shell-Diffusion.cuh"

#include "utilities/Keyboard_Interrupt.hpp"
#include "utilities/File_Generator.hpp"

#include <iostream>
#include <boost/program_options.hpp>

#define MAX_ITER 1E6

size_t particle_num_grid_points = 1001;
double time_step = 1E-6;	// s

double temperature = 1500;	// K
double core_radius = 32.5E-6;	// m
double overall_radius = 39.5E-6;	// m
double activation_energy = 102.191E3;	// J / mol.
double pre_exponential_factor = 2.56E-6;	// m2 / s

__device__ CoreShellDIffusion::Diffusion *diffusion_problem;
__device__ ArrheniusDiffusivityModel *diffusivity_model;
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
	double activation_energy,
	double pre_exponential_factor,
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

__host__ int setConfiguration(int argc, char const *argv[]);

int main(int argc, char const *argv[])
{
	int flag = setConfiguration(argc, argv);

	if (flag != 2) return flag;

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
		activation_energy,
		pre_exponential_factor,
		concentration_array_A,
		concentration_array_B
	);
	cudaDeviceSynchronize();

	FileGenerator folder;
	
	std::ofstream concentration_array_A_file = folder.getCSVFile("concentration_A");
	std::ofstream concentration_array_B_file = folder.getCSVFile("concentration_B");

	std::ofstream mass_file = folder.getCSVFile("mass");
	std::ofstream mass_fractions_file = folder.getCSVFile("mass_fractions");
	
	std::ofstream config_file = folder.getTXTFile("configuration");
	
	config_file << "Temperature :\t" << temperature << "\tK\n\n\n";
	config_file << "Arrhenius Diffusivity Model Parameters\n\n";
	config_file << "Pre-exponential Factor :\t" << pre_exponential_factor << "\tm2 / s\n";
	config_file << "Activation Energy :\t" << activation_energy << "\tJ / mol. - K\n\n\n";

	CoreShellDIffusion::printConfiguration(config_file);
	config_file.close();

	CoreShellDIffusion::printGridPoints(concentration_array_A_file, ',');
	CoreShellDIffusion::printGridPoints(concentration_array_B_file, ',');

	initializeCoreShellParticle<<<1, particle_num_grid_points>>>();
	cudaDeviceSynchronize();

	cudaMemcpy(concentration_array_A_host, concentration_array_A, particle_num_grid_points * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(concentration_array_B_host, concentration_array_B, particle_num_grid_points * sizeof(double), cudaMemcpyDeviceToHost);

	CoreShellDIffusion::printConcentrationArray(concentration_array_A_file, concentration_array_A_host, 0.0, ',');
	CoreShellDIffusion::printConcentrationArray(concentration_array_B_file, concentration_array_B_host, 0.0, ',');

	initIteration<<<1,1>>>(temperature);
	cudaDeviceSynchronize();

	std::cout << "Initialized Diffusion Problems\nStarting Iterations\n(Press Ctrl+C to break)\n";

	setUpKeyboardInterrupt();

	try
	{
		double time = 0.0;

		for (size_t i = 0; i < MAX_ITER;)
		{
			size_t i_top = i + step;
			for (; i < i_top; i++)
			{
				setUpEquations<<<1, particle_num_grid_points>>>();
				
				solve<<<1,2>>>();
				time += time_step;

				updateMassFractions<<<1,1>>>();

				cudaMemcpyFromSymbol(&combustion_complete, ::combustion_complete, sizeof(bool));
				if (combustion_complete) break;
			}

			cudaDeviceSynchronize();

			std::cout << "Iteration Completed : " << i << "\n";

			cudaMemcpy(concentration_array_A_host, concentration_array_A, particle_num_grid_points * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(concentration_array_B_host, concentration_array_B, particle_num_grid_points * sizeof(double), cudaMemcpyDeviceToHost);

			CoreShellDIffusion::printConcentrationArray(concentration_array_A_file, concentration_array_A_host, time, ',');
			CoreShellDIffusion::printConcentrationArray(concentration_array_B_file, concentration_array_B_host, time, ',');

			printMass(mass_file);
			printMassFractions(mass_fractions_file);

			if (combustion_complete) break;
		}
	}

	catch (InterruptException& e)
    {
		std::cout << "\nQuitting...\n";

		cudaDeviceSynchronize();

		cudaMemcpy(concentration_array_A_host, concentration_array_A, particle_num_grid_points * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(concentration_array_B_host, concentration_array_B, particle_num_grid_points * sizeof(double), cudaMemcpyDeviceToHost);

		CoreShellDIffusion::printConcentrationArray(concentration_array_A_file, concentration_array_A_host, particle_num_grid_points, ',');
		CoreShellDIffusion::printConcentrationArray(concentration_array_B_file, concentration_array_B_host, particle_num_grid_points, ',');

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
	double activation_energy,
	double pre_exponential_factor,
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
	diffusivity_model = new ArrheniusDiffusivityModel(pre_exponential_factor, activation_energy);

	diffusion_problem->setArrayAddresses(concentration_array_A, concentration_array_B);
}

__global__ void deallocateMemory()
{
	delete diffusivity_model;
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
	diffusion_problem->setDiffusivity(diffusivity_model->getDiffusivity(temperature));
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

int setConfiguration(int argc, const char * argv[])
{
	try
	{
		boost::program_options::options_description options_description_obj("Usage: Core-Shell-Diffusion [options]\nOptions");

		options_description_obj.add_options()
			("help",	"Produce help message")
			("Du",		"Set diffusivity parameters to Du et al's model. Default is Alawieh et al's model.")
			("num",		boost::program_options::value<size_t>(),	"Set number of grid points for particle diffusion pde-solver to arg. Default is 1001.")
			("Dt",		boost::program_options::value<double>(),	"Set time step in seconds to arg. Default is 1E-6 seconds.")
			("temp",	boost::program_options::value<double>(),	"Set particle temperature in kelvin to arg. Default is 1500 K.")
			("radius",	boost::program_options::value<double>(),	"Set particle radius in metre to arg. Default is 39.5E-6 m.")
			("core",	boost::program_options::value<double>(),	"Set particle core radius in metre to arg. Default is 32.5E-6 m");

		boost::program_options::variables_map variables_map_obj;
		boost::program_options::store(
			boost::program_options::parse_command_line(argc, argv, options_description_obj),
			variables_map_obj
		);
		boost::program_options::notify(variables_map_obj);

		if (variables_map_obj.count("help"))
		{
			std::cout << options_description_obj << std::endl;
			return 0;
		}

		if (variables_map_obj.count("Du"))
		{
			activation_energy = 26E3;
			pre_exponential_factor = 9.54E-8;
			std::cout << "Using Du et al's Diffusivity parameters." << std::endl;
		}

		if (variables_map_obj.count("num"))
		{
			particle_num_grid_points = variables_map_obj["num"].as<size_t>();
		}

		if (variables_map_obj.count("Dt"))
		{
			time_step = variables_map_obj["Dt"].as<double>();
			
			if (time_step <= 0.0)
			{
				std::cerr << "Dt cannot be negative. Given value : " << time_step << std::endl;
				return 1;
			}
		}

		if (variables_map_obj.count("temp"))
		{
			temperature = variables_map_obj["temp"].as<double>();

			if (temperature <= 0.0)
			{
				std::cerr << "Temperature cannot be negative. Given value : " << temperature << std::endl;
				return 1;
			}
		}

		if (variables_map_obj.count("radius"))
		{
			overall_radius = variables_map_obj["radius"].as<double>();

			if (overall_radius <= 0.0)
			{
				std::cerr << "Overall radius cannot be negative. Given value : " << overall_radius << std::endl;
				return 1;
			}
		}

		if (variables_map_obj.count("core"))
		{
			core_radius = variables_map_obj["core"].as<double>();

			if (core_radius <= 0.0)
			{
				std::cerr << "Core radius cannot be negative. Given value : " << core_radius << std::endl;
				return 1;
			}

			else if (core_radius >= overall_radius)
			{
				std::cerr << "Core radius cannot be greater than overall radius. Given value : " << core_radius << std::endl;
				std::cerr << "Overall radius : " << overall_radius << std::endl;
				return 1;
			}
		}
	}

	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	return 2;
}