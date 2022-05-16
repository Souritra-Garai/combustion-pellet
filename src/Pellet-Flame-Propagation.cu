#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>

#include "species/Aluminium.cuh"
#include "species/Argon.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "pde-problems/Pellet-Flame-Propagation.cuh"

#include "utilities/File_Generator.hpp"
#include "utilities/Program_Options.cuh"
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

__device__ double sharpness_coefficient = 0.1;

__device__ double2 diffusivity_parameters = {2.56E-6, 102.191E3};
__device__ double2 diffusivity_parameters_low = {2.08E-7, 92.586E3};

__global__ void allocateMemory(double delta_t, size_t num_grid_points_pellet, size_t num_grid_points_particle);
__global__ void allocateParticles();
__global__ void deallocateParticles();
__global__ void deallocateMemory();

__global__ void setArrayAddresses(double *temperature, double *concentration_A, double *concentration_B);
__global__ void setInitialConditions();

__host__ __forceinline__ void evolveMainParticles();
__host__ __forceinline__ void evolveAuxiliaryParticles();
__host__ __forceinline__ void iterate();

__host__ void parseProgramOptions(int argc, char const *argv[]);

int main(int argc, char const *argv[])
{
	parseProgramOptions(argc, argv);

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

	double phi;
	cudaMemcpyFromSymbol(&phi, ::particle_volume_fractions, sizeof(double));
	double2 diffusivity_parameters;
	cudaMemcpyFromSymbol(&diffusivity_parameters, ::diffusivity_parameters, sizeof(double2));

	PelletFlamePropagation::printConfiguration(config_file, phi, diffusivity_parameters);
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

	std::cout << "Combustion completed.\n";

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
	loadAluminium(sharpness_coefficient);
	loadArgon();
	loadNickel(sharpness_coefficient);
	loadNickelAlumnide(sharpness_coefficient);

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
	PackedPellet::setHeatLossParameters(15.681, 15.715, 0.25);
	PackedPellet::setTemperatureParameters(298.15, 900.0);

	PelletFlamePropagation::setNumGridPoints(num_grid_points_pellet);
	PelletFlamePropagation::setImplicitness(implicitness_source_term, implicitness_diffusion_term);
	PelletFlamePropagation::setInitialConditions(initial_ignition_temperature, initial_ignition_fraction * length);
	PelletFlamePropagation::setInfinitesimalChangeInTemperature(delta_T);

	flame_propagation_problem = new PelletFlamePropagation::FlamePropagation(particle_volume_fractions);
	diffusivity_model = new ArrheniusDiffusivityModel(diffusivity_parameters.x, diffusivity_parameters.y);
	diffusivity_model->setParametersLow(diffusivity_parameters_low.x, diffusivity_parameters_low.y);

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

__host__ void parseProgramOptions(int argc, char const *argv[])
{
	ez::ezOptionParser opt;

	opt.overview	= "Flame Propagation in a pellet packed with energetic intermetallic core-shell particles.";
	opt.syntax		= "Pellet-Flame-Propagation [OPTIONS]";
	opt.example		= "Pellet-Flame-Propagation --phi 0.68\n\n";
	opt.footer		= "Developed by Souritra Garai, 2021-22.\n";

	setHelpOption(opt);
	setSharpnessCoefficientOption(opt);
	setParticleGridSizeOption(opt);
	setPelletGridSizeOption(opt);
	setTimeStepOption(opt);
	setTemperatureStepOption(opt);
	setDiffusivityModelOption(opt);
	setCoreRadiusOption(opt);
	setOverallRadiusOption(opt);
	setLengthOption(opt);
	setDiameterOption(opt);
	setGammaOption(opt);
	setKappaOption(opt);
	setPhiOption(opt);

	opt.parse(argc, argv);

	displayHelpOption(opt);

	double sharpness_coefficient = 0.1;
	sharpness_coefficient = getSharpnessCoefficientOption(opt, sharpness_coefficient);
	cudaMemcpyToSymbol(::sharpness_coefficient, &sharpness_coefficient, sizeof(double));
	
	double delta_T = 0.001;
	delta_T = getTemperatureStepOption(opt, delta_T);
	cudaMemcpyToSymbol(::delta_T, &delta_T, sizeof(delta_T));

	num_grid_points_particle = getParticleGridSizeOption(opt, num_grid_points_particle);
	num_grid_points_pellet   = getPelletGridSizeOption(opt, num_grid_points_pellet);

	delta_t = getTimeStepOption(opt, delta_t);

	double2 parameters = {2.56E-6, 102.191E3};
	double2 parameters_low = {2.08E-7, 92.586E3};
	parameters = getDiffusivityModelOption(opt, parameters, parameters_low);
	cudaMemcpyToSymbol(::diffusivity_parameters, &parameters, sizeof(double2));
	cudaMemcpyToSymbol(::diffusivity_parameters_low, &parameters_low, sizeof(double2));

	double core_radius = 32.5E-6;
	core_radius = getCoreRadiusOption(opt, core_radius);
	cudaMemcpyToSymbol(::core_radius, &core_radius, sizeof(double));
	double overall_radius = 39.5E-6;
	overall_radius = getOverallRadiusOption(opt, overall_radius);
	cudaMemcpyToSymbol(::overall_radius, &overall_radius, sizeof(double));

	double pellet_length = 6.35E-3;
	pellet_length = getLengthOption(opt, pellet_length);
	cudaMemcpyToSymbol(::length, &pellet_length, sizeof(double));
	double pellet_diameter = 6.35E-3;
	pellet_diameter = getDiameterOption(opt, pellet_diameter);
	cudaMemcpyToSymbol(::diameter, &pellet_diameter, sizeof(double));

	double diffusion_term_implicitness = 0.5;
	diffusion_term_implicitness = getKappaOption(opt, diffusion_term_implicitness);
	cudaMemcpyToSymbol(::implicitness_diffusion_term, &diffusion_term_implicitness, sizeof(double));
	double source_term_implicitness = 0.5;
	source_term_implicitness = getGammaOption(opt, source_term_implicitness);
	cudaMemcpyToSymbol(::implicitness_source_term, &source_term_implicitness, sizeof(double));

	double phi = 0.68;
	phi = getPhiOption(opt, phi);
	cudaMemcpyToSymbol(::particle_volume_fractions, &phi, sizeof(double));
}