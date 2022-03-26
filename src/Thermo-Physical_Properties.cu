#include "species/Aluminium.cuh"
#include "species/Argon.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "utilities/File_Generator.hpp"

#include <boost/program_options.hpp>
#include <iostream>

#define TEMPERATURE_UPPER_BOUND 2500.0	// K
#define TEMPERATURE_LOWER_BOUND 273.15	// K

#define GET_TEMPERATURE(i, n) (TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i) / ((double) n - 1.0))

const unsigned int n = 10000;
double sharpness_coefficient = 100.0;

__host__ void initializeSpecies();
__host__ void deinitializeSpecies();
__host__ void saveSpecies(std::ofstream &output_file);

Species *aluminium_dev_ptr, *nickel_dev_ptr, *nickel_aluminide_dev_ptr;
IdealGas *argon_dev_ptr;

// Arrays for storing temperature and corresponding thermophysical properties on host memory
double heat_capacities[n], enthalpies[n], thermal_conductivities[n];

// Pointer to arrays for storing temperature and corresponding thermophysical properties
// on device memory
double *heat_capacities_dev_ptr, *enthalpies_dev_ptr, *thermal_conductivities_dev_ptr;

__global__ void getThermalConductivities(unsigned int n, double *thermal_conductivities, Species *species_ptr);
__global__ void getHeatCapacities(unsigned int n, double *heat_capacities, Species *species_ptr);
__global__ void getEnthalpies(unsigned int n, double *enthalpies, Species *species_ptr);

__global__ void getThermalConductivities(unsigned int n, double *thermal_conductivities, IdealGas *gas);
__global__ void getHeatCapacities(unsigned int n, double *heat_capacities, IdealGas *gas);
__global__ void getEnthalpies(unsigned int n, double *enthalpies, IdealGas *gas);

int main(int argc, char const *argv[])
{
	try
	{
		boost::program_options::options_description options_description_obj("Usage: Thermo-Physical_Properties [options]\nOptions");

		options_description_obj.add_options()
			("help",	"produce help message")
			("alpha",	boost::program_options::value<double>(),		"set sharpness coefficient to arg, default is 100.0");

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

		if (variables_map_obj.count("alpha"))
		{
			sharpness_coefficient = variables_map_obj["alpha"].as<double>();
		}
	}

	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return(1);
	}
	
	cudaMalloc(&thermal_conductivities_dev_ptr, n * sizeof(double));
	cudaMalloc(&heat_capacities_dev_ptr, n * sizeof(double));
	cudaMalloc(&enthalpies_dev_ptr, n * sizeof(double));

	initializeSpecies();

	FileGenerator solution_folder;

	// Aluminium
	getThermalConductivities<<<(n + 255) / 256, 256>>>(n, thermal_conductivities_dev_ptr, aluminium_dev_ptr);
	getHeatCapacities<<<(n + 255) / 256, 256>>>(n, heat_capacities_dev_ptr, aluminium_dev_ptr);
	getEnthalpies<<<(n + 255) / 256, 256>>>(n, enthalpies_dev_ptr, aluminium_dev_ptr);

	cudaDeviceSynchronize();

	cudaMemcpy(thermal_conductivities, thermal_conductivities_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(heat_capacities, heat_capacities_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(enthalpies, enthalpies_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	
	std::ofstream aluminium_file = solution_folder.getCSVFile("aluminium");
	saveSpecies(aluminium_file);
	aluminium_file.close();

	// Argon
	getThermalConductivities<<<(n + 255) / 256, 256>>>(n, thermal_conductivities_dev_ptr, argon_dev_ptr);
	getHeatCapacities<<<(n + 255) / 256, 256>>>(n, heat_capacities_dev_ptr, argon_dev_ptr);
	getEnthalpies<<<(n + 255) / 256, 256>>>(n, enthalpies_dev_ptr, argon_dev_ptr);

	cudaDeviceSynchronize();

	cudaMemcpy(thermal_conductivities, thermal_conductivities_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(heat_capacities, heat_capacities_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(enthalpies, enthalpies_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	
	std::ofstream argon_file = solution_folder.getCSVFile("argon");
	saveSpecies(argon_file);
	argon_file.close();

	// Nickel
	getThermalConductivities<<<(n + 255) / 256, 256>>>(n, thermal_conductivities_dev_ptr, nickel_dev_ptr);
	getHeatCapacities<<<(n + 255) / 256, 256>>>(n, heat_capacities_dev_ptr, nickel_dev_ptr);
	getEnthalpies<<<(n + 255) / 256, 256>>>(n, enthalpies_dev_ptr, nickel_dev_ptr);

	cudaDeviceSynchronize();

	cudaMemcpy(thermal_conductivities, thermal_conductivities_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(heat_capacities, heat_capacities_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(enthalpies, enthalpies_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	
	std::ofstream nickel_file = solution_folder.getCSVFile("nickel");
	saveSpecies(nickel_file);
	nickel_file.close();

	// Nickel Aluminide
	getThermalConductivities<<<(n + 255) / 256, 256>>>(n, thermal_conductivities_dev_ptr, nickel_aluminide_dev_ptr);
	getHeatCapacities<<<(n + 255) / 256, 256>>>(n, heat_capacities_dev_ptr, nickel_aluminide_dev_ptr);
	getEnthalpies<<<(n + 255) / 256, 256>>>(n, enthalpies_dev_ptr, nickel_aluminide_dev_ptr);

	cudaDeviceSynchronize();

	cudaMemcpy(thermal_conductivities, thermal_conductivities_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(heat_capacities, heat_capacities_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(enthalpies, enthalpies_dev_ptr, n * sizeof(double), cudaMemcpyDeviceToHost);
	
	std::ofstream nickel_aluminide_file = solution_folder.getCSVFile("nickel_aluminide");
	saveSpecies(nickel_aluminide_file);
	nickel_aluminide_file.close();

	deinitializeSpecies();

	cudaFree(thermal_conductivities_dev_ptr);
	cudaFree(heat_capacities_dev_ptr);
	cudaFree(enthalpies_dev_ptr);

	return 0;
}

__host__ void initializeSpecies()
{
	loadAluminium(sharpness_coefficient);
	loadArgon();
	loadNickel(sharpness_coefficient);
	loadNickelAlumnide(sharpness_coefficient);

	cudaMalloc(&aluminium_dev_ptr, sizeof(Species));
	cudaMemcpy(aluminium_dev_ptr, &aluminium, sizeof(Species), cudaMemcpyHostToDevice);

	cudaMalloc(&nickel_dev_ptr, sizeof(Species));
	cudaMemcpy(nickel_dev_ptr, &nickel, sizeof(Species), cudaMemcpyHostToDevice);

	cudaMalloc(&nickel_aluminide_dev_ptr, sizeof(Species));
	cudaMemcpy(nickel_aluminide_dev_ptr, &nickel_aluminide, sizeof(Species), cudaMemcpyHostToDevice);

	cudaMalloc(&argon_dev_ptr, sizeof(IdealGas));
	cudaMemcpy(argon_dev_ptr, &argon, sizeof(IdealGas), cudaMemcpyHostToDevice);
}

__host__ void deinitializeSpecies()
{
	aluminium.deallocateMemory();
	nickel.deallocateMemory();
	nickel_aluminide.deallocateMemory();

	cudaFree(aluminium_dev_ptr);
	cudaFree(argon_dev_ptr);
	cudaFree(nickel_dev_ptr);
	cudaFree(nickel_aluminide_dev_ptr);
}

__host__ void saveSpecies(std::ofstream &file)
{
	file << "Temperature (K)," << "Heat Capacity (J / mol. - K)," << "Enthalpy (J / mol.)," << "Thermal Conductivity (W / m - K)\n";

	for (unsigned int i = 0; i < n; i++)
	{
		file << GET_TEMPERATURE(i, n) << ',' << heat_capacities[i] << ',' << enthalpies[i] << ',' <<
		thermal_conductivities[i] << '\n';
	}
}

__global__ void getThermalConductivities(unsigned int n, double *thermal_conductivities, Species *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n) thermal_conductivities[i] = species_ptr->getThermalConductivity(GET_TEMPERATURE(i, n));
}

__global__ void getHeatCapacities(unsigned int n, double *heat_capacities, Species *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n) heat_capacities[i] = species_ptr->getHeatCapacity(GET_TEMPERATURE(i, n));
}

__global__ void getEnthalpies(unsigned int n, double *enthalpies, Species *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n) enthalpies[i] = species_ptr->getEnthalpy(GET_TEMPERATURE(i, n));
}

__global__ void getThermalConductivities(unsigned int n, double *thermal_conductivities, IdealGas *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n) thermal_conductivities[i] = species_ptr->getThermalConductivity(GET_TEMPERATURE(i, n));
}

__global__ void getHeatCapacities(unsigned int n, double *heat_capacities, IdealGas *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n) heat_capacities[i] = species_ptr->getHeatCapacity(GET_TEMPERATURE(i, n));
}

__global__ void getEnthalpies(unsigned int n, double *enthalpies, IdealGas *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n) enthalpies[i] = species_ptr->getEnthalpy(GET_TEMPERATURE(i, n));
}