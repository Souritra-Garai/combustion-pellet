#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>

#include "species/Aluminium.cuh"
#include "species/Argon.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "pde-problems/Pellet-Flame-Propagation.cuh"

#include "utilities/File_Generator.hpp"
#include "utilities/Keyboard_Interrupt.hpp"

#define MAX_ITER 1E8

double delta_t = 1E-6;
size_t num_grid_points_particle = 101;
size_t num_grid_points_pellet = 101;

__device__ PelletFlamePropagation::FlamePropagation *flame_propagation_problem;
__device__ ArrheniusDiffusivityModel *diffusivity_model;
__device__ bool combustion_complete;

__device__ double core_radius = 32.5E-6;
__device__ double overall_radius = 39.5E-6;

__device__ double length = 6.35E-3;
__device__ double diameter = 6.35E-3;

__device__ double delta_T = 0.001;

__device__ double implicitness_diffusion_term = 0.5;
__device__ double implicitness_source_term = 0.5;

__device__ double initial_ignition_fraction = 0.1;
__device__ double initial_ignition_temperature = 1500;

__device__ double particle_volume_fractions = 0.7;

__device__ double pre_exponential_factor = 2.56E-6;
__device__ double activation_energy = 102.191E3;

__global__ void allocateMemory(double delta_t, size_t num_grid_points_pellet, size_t num_grid_points_particle);
__global__ void allocateParticles();
__global__ void deallocateParticles();
__global__ void deallocateMemory();

__global__ void setArrayAddresses(double *temperature, double *concentration_A, double *concentration_B);
__global__ void setInitialConditions();

__host__ __forceinline__ void evolveMainParticles();
__host__ __forceinline__ void evolveAuxiliaryParticles();
__host__ __forceinline__ void iterate();

int main(int argc, char const *argv[])
{
	cudaDeviceSetLimit(cudaLimitMallocHeapSize, 1024 * 1024 * 500);

	allocateMemory<<<1,1>>>(delta_t, num_grid_points_pellet, num_grid_points_particle);
	allocateParticles<<<num_grid_points_pellet, 1>>>();

	double *temperature_array_device;
	double *concentration_array_A_device;
	double *concentration_array_B_device;

	cudaMalloc(&temperature_array_device, num_grid_points_pellet * sizeof(double));
	cudaMalloc(&concentration_array_A_device, 3 * num_grid_points_particle * num_grid_points_pellet * sizeof(double));
	cudaMalloc(&concentration_array_B_device, 3 * num_grid_points_particle * num_grid_points_pellet * sizeof(double));

	setArrayAddresses<<<1,1>>>(temperature_array_device, concentration_array_A_device, concentration_array_B_device);
	setInitialConditions<<<num_grid_points_pellet, num_grid_points_particle>>>();

	cudaDeviceSynchronize();

	std::cout << "Initialization complete.\n";

	FileGenerator folder;
	std::ofstream config_file = folder.getTXTFile("configurations");
	std::ofstream temperature_file = folder.getCSVFile("temperature");

	PelletFlamePropagation::printConfiguration(config_file);
	config_file.close();

	double time = 0.0;
	double temperature_array_host[num_grid_points_pellet];
	cudaMemcpy(temperature_array_host, temperature_array_device, num_grid_points_pellet * sizeof(double), cudaMemcpyDeviceToHost);

	PelletFlamePropagation::printGridPoints(temperature_file, ',');
	PelletFlamePropagation::printTemperatureArray(temperature_file, temperature_array_host, time, ',');

	std::cout << "Starting Iterations\n(Press Ctrl+C to break)\n" ;

	setUpKeyboardInterrupt();

	try
	{
		size_t step = 0.001 / delta_t;

		bool combustion_complete = false;

		evolveMainParticles();

		for (size_t i = 0; i < MAX_ITER && !combustion_complete;)
		{
			size_t i_steps = i + step;

			for (; i < i_steps && !combustion_complete; i++)
			{
				iterate();
				time += delta_t;

				cudaMemcpyFromSymbol(&combustion_complete, ::combustion_complete, sizeof(bool));
			}

			std::cout << "Iterations completed : " << i << "\n";

			cudaDeviceSynchronize();

			cudaMemcpy(temperature_array_host, temperature_array_device, num_grid_points_pellet * sizeof(double), cudaMemcpyDeviceToHost);
			PelletFlamePropagation::printTemperatureArray(temperature_file, temperature_array_host, time, ',');
		}
	}

	catch (InterruptException& e)
    {
		std::cout << "\nQuitting...\n";

		cudaDeviceSynchronize();

		cudaMemcpy(temperature_array_host, temperature_array_device, num_grid_points_pellet * sizeof(double), cudaMemcpyDeviceToHost);
		PelletFlamePropagation::printTemperatureArray(temperature_file, temperature_array_host, time, ',');
    }

	temperature_file.close();

	deallocateParticles<<<num_grid_points_pellet, 1>>>();
	deallocateMemory<<<1,1>>>();

	cudaFree(temperature_array_device);
	cudaFree(concentration_array_A_device);
	cudaFree(concentration_array_B_device);

	return 0;
}

__global__ void allocateMemory(
	double delta_t,
	size_t num_grid_points_pellet,
	size_t num_grid_points_particle
) {
	loadAluminium(1);
	loadArgon();
	loadNickel(1);
	loadNickelAlumnide(1);

	CoreShellParticle::initialize(
		aluminium,
		nickel,
		nickel_aluminide,
		core_radius,
		overall_radius
	);

	CoreShellDIffusion::setGridSize(num_grid_points_particle);
	CoreShellDIffusion::setTimeStep(delta_t);

	PackedPellet::setPelletDimensions(length, diameter);
	PackedPellet::setDegassingFluid(argon);
	PackedPellet::setHeatLossParameters(0, 0, 0);
	PackedPellet::setTemperatureParameters(298.15, 900.0);

	PelletFlamePropagation::setNumGridPoints(num_grid_points_pellet);
	PelletFlamePropagation::setImplicitness(implicitness_source_term, implicitness_diffusion_term);
	PelletFlamePropagation::setInitialConditions(initial_ignition_temperature, initial_ignition_fraction * length);
	PelletFlamePropagation::setInfinitesimalChangeInTemperature(delta_T);

	flame_propagation_problem = new PelletFlamePropagation::FlamePropagation(particle_volume_fractions);
	diffusivity_model = new ArrheniusDiffusivityModel(pre_exponential_factor, activation_energy);

	flame_propagation_problem->setDiffusivityModel(diffusivity_model);
}

__global__ void allocateParticles()
{
	flame_propagation_problem->allocateParticleMemory(blockIdx.x);
}

__global__ void deallocateParticles()
{
	flame_propagation_problem->deallocateParticleMemory(blockIdx.x);
}

__global__ void deallocateMemory()
{
	delete flame_propagation_problem;
	delete diffusivity_model;

	unloadAluminium();
	unloadArgon();
	unloadNickel();
	unloadNickelAluminide();
}

__global__ void setArrayAddresses(double *temperature, double *concentration_A, double *concentration_B)
{
	flame_propagation_problem->setArrayAddresses(temperature, concentration_A, concentration_B);
}

__global__ void setInitialConditions()
{
	flame_propagation_problem->initializePellet(blockIdx.x, threadIdx.x);
}

__global__ void setMainParticleDiffusivity()
{
	flame_propagation_problem->setParticleDiffusivity(0, blockIdx.x);
}

__global__ void setMainParticleEquations()
{
	flame_propagation_problem->setParticleEquation(0, blockIdx.x, threadIdx.x);
}

__global__ void solveMainParticleEquations()
{
	flame_propagation_problem->solveParticleEquation(0, blockIdx.x, threadIdx.x);
}

__global__ void updateMainParticles()
{
	flame_propagation_problem->collectEnthalpy(blockIdx.x);
	flame_propagation_problem->updateParticleMassFraction(0, blockIdx.x);
	flame_propagation_problem->updateThermalConductivity(blockIdx.x);
}

__host__ __forceinline__ void evolveMainParticles()
{
	setMainParticleDiffusivity<<<num_grid_points_pellet, 1>>>();
	setMainParticleEquations<<<num_grid_points_pellet, num_grid_points_particle>>>();
	solveMainParticleEquations<<<num_grid_points_pellet, 2>>>();
	updateMainParticles<<<num_grid_points_pellet, 1>>>();
}

__global__ void setAuxiliaryParticleDiffusivity()
{
	flame_propagation_problem->setParticleDiffusivity(1 + threadIdx.x, blockIdx.x);
}

__global__ void setAuxiliaryParticleEquations()
{
	flame_propagation_problem->setParticleEquation(1, blockIdx.x, threadIdx.x);
	flame_propagation_problem->setParticleEquation(2, blockIdx.x, threadIdx.x);
}

__global__ void solveAuxiliaryParticleEquations()
{
	flame_propagation_problem->solveParticleEquation(
		1 + (threadIdx.x > 1) * 1,
		blockIdx.x,
		threadIdx.x - (threadIdx.x > 1) * 2);
}

__global__ void updateAuxiliaryParticles()
{
	flame_propagation_problem->updateParticleMassFraction(1 + threadIdx.x, blockIdx.x);
}

__host__ __forceinline__ void evolveAuxiliaryParticles()
{
	setAuxiliaryParticleDiffusivity<<<num_grid_points_pellet, 2>>>();
	setAuxiliaryParticleEquations<<<num_grid_points_pellet, num_grid_points_particle>>>();
	solveAuxiliaryParticleEquations<<<num_grid_points_pellet, 4>>>();
	updateAuxiliaryParticles<<<num_grid_points_pellet, 2>>>();
}

__global__ void setTemperatureEquations()
{
	flame_propagation_problem->setUpEquation(blockIdx.x);
}

__global__ void solveTemperatureEquations()
{
	flame_propagation_problem->solveEquation();

	combustion_complete = flame_propagation_problem->isCombustionComplete();
}

__host__ __forceinline__ void iterate()
{
	evolveAuxiliaryParticles();

	setTemperatureEquations<<<num_grid_points_pellet, 1>>>();
	solveTemperatureEquations<<<1, 1>>>();

	evolveMainParticles();
}

int setConfiguration(int argc, const char * argv[])
{
	try
	{
		boost::program_options::options_description options_description_obj("Usage: Core-Shell-Diffusion [options]\nOptions");

		options_description_obj.add_options()
			("help",		"Produce help message")
			("Du",			"Set diffusivity parameters to Du et al's model. Default is Alawieh et al's model.")
			("N",			boost::program_options::value<size_t>(),	"Set number of grid points for particle diffusion pde-solver to arg. Default is 101.")
			("M",			boost::program_options::value<size_t>(),	"Set number of grid points for pellet enthalpy equation pde-solver to arg. Default is 101.")
			("Dt",			boost::program_options::value<double>(),	"Set time step in seconds to arg. Default is 1E-6 seconds.")
			("DT",			boost::program_options::value<double>(),	"Set infinitesimal change in temperature (in kelvin) for temperature derivative. Default is 0.001 seconds.")
			("r_P",			boost::program_options::value<double>(),	"Set particle radius in metre to arg. Default is 39.5E-6 m.")
			("r_C",			boost::program_options::value<double>(),	"Set particle core radius in metre to arg. Default is 32.5E-6 m")
			("phi",			boost::program_options::value<double>(),	"Set particle volume fractions in the pellet to arg. Default is 0.7.")
			("gamma",		boost::program_options::value<double>(),	"Set implicitness of source terms to arg. Default is 0.5.")
			("kappa",		boost::program_options::value<double>(),	"Set implicitness of diffusion terms to arg. Default is 0.5")
			("length",		boost::program_options::value<double>(),	"Set length of pellet in metre to arg. Default is 6.35E-3 m.")
			("diameter",	boost::program_options::value<double>(),	"Set diameter of pellet in metre to arg. Default is 6.35E-3 m.");

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
			double activation_energy = 26E3;
			double pre_exponential_factor = 9.54E-8;

			cudaMemcpyToSymbol(::activation_energy, &activation_energy, sizeof(double));
			cudaMemcpyToSymbol(::pre_exponential_factor, &pre_exponential_factor, sizeof(double));

			std::cout << "Using Du et al's Diffusivity parameters." << std::endl;
		}

		if (variables_map_obj.count("N"))
		{
			num_grid_points_particle = variables_map_obj["N"].as<size_t>();
		}

		if (variables_map_obj.count("M"))
		{
			num_grid_points_pellet = variables_map_obj["M"].as<size_t>();
		}

		if (variables_map_obj.count("Dt"))
		{
			delta_t = variables_map_obj["Dt"].as<double>();
			
			if (delta_t <= 0.0)
			{
				std::cerr << "Dt cannot be negative. Given value : " << delta_t << std::endl;
				return 1;
			}
		}

		if (variables_map_obj.count("DT"))
		{
			double delta_T = variables_map_obj["DT"].as<double>();

			cudaMemcpyToSymbol(::delta_T, &delta_T, sizeof(double));
			
			if (delta_T <= 0.0)
			{
				std::cerr << "DT cannot be negative. Given value : " << delta_T << std::endl;
				return 1;
			}
		}

		if (variables_map_obj.count("r_P"))
		{
			double overall_radius = variables_map_obj["r_P"].as<double>();

			cudaMemcpyToSymbol(::overall_radius, &overall_radius, sizeof(double));

			if (overall_radius <= 0.0)
			{
				std::cerr << "Overall radius cannot be negative. Given value : " << overall_radius << std::endl;
				return 1;
			}
		}

		if (variables_map_obj.count("r_C"))
		{
			double core_radius = variables_map_obj["r_C"].as<double>();

			cudaMemcpyToSymbol(::core_radius, &core_radius, sizeof(double));

			if (core_radius <= 0.0)
			{
				std::cerr << "Core radius cannot be negative. Given value : " << core_radius << std::endl;
				return 1;
			}

			double overall_radius;
			cudaMemcpyFromSymbol(&overall_radius, ::overall_radius, sizeof(double));

			if (core_radius >= overall_radius)
			{
				std::cerr << "Core radius cannot be greater than overall radius. Given value : " << core_radius << std::endl;
				std::cerr << "Overall radius : " << overall_radius << std::endl;
				return 1;
			}
		}

		if (variables_map_obj.count("phi"))
		{
			double phi = variables_map_obj["phi"].as<double>();

			cudaMemcpyToSymbol(::particle_volume_fractions, &phi, sizeof(double));
			
			if (phi <= 0.0 || phi >= 1.0)
			{
				std::cerr << "phi must belong to (0, 1). Given value : " << phi << std::endl;
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