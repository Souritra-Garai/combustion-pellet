#include "species/Aluminium.cuh"
#include "species/Argon.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "utilities/File_Generator.hpp"

#include <boost/program_options.hpp>
#include <iostream>

#define N 10000
#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 2000.0	// K
#define GET_TEMPERATURE(i) (TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i) / ((double) N - 1.0))

double sharpness_coefficient = 0.1;

__global__ void allocateMemory(double sharpness_coefficient)
{
	loadAluminium(sharpness_coefficient);
	loadArgon();
	loadNickel(sharpness_coefficient);
	loadNickelAlumnide(sharpness_coefficient);
}

__global__ void deallocateMemory()
{
	unloadAluminium();
	unloadArgon();
	unloadNickel();
	unloadNickelAluminide();
}

__global__ void getThermalConductivities(double *thermal_conductivities, Species *species_ptr);
__global__ void getHeatCapacities(double *heat_capacities, Species *species_ptr);
__global__ void getEnthalpies(double *enthalpies, Species *species_ptr);

__global__ void getThermalConductivities(double *thermal_conductivities, IdealGas *gas);
__global__ void getHeatCapacities(double *heat_capacities, IdealGas *gas);
__global__ void getEnthalpies(double *enthalpies, IdealGas *gas);

__host__ void saveSpecies(std::ofstream &file);

double enthalpy_array_h[N], heat_capacity_array_h[N], thermal_conductivity_array_h[N];

double *enthalpy_array_d, *heat_capacity_array_d, *thermal_conductivity_array_d;

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

	cudaMalloc(&enthalpy_array_d, N*sizeof(double));
	cudaMalloc(&heat_capacity_array_d, N*sizeof(double));
	cudaMalloc(&thermal_conductivity_array_d, N*sizeof(double));

	allocateMemory<<<1,1>>>(sharpness_coefficient);

	Species *aluminium, *nickel, *nickel_aluminide;
	IdealGas *argon;

	cudaMemcpyFromSymbol(&aluminium, ::aluminium, sizeof(Species*));
	cudaMemcpyFromSymbol(&argon, ::argon, sizeof(IdealGas*));
	cudaMemcpyFromSymbol(&nickel, ::nickel, sizeof(Species*));
	cudaMemcpyFromSymbol(&nickel_aluminide, ::nickel_aluminide, sizeof(Species*));
	cudaDeviceSynchronize();

	FileGenerator solution_folder;

	// Aluminium
	getThermalConductivities<<<(N + 255) / 256, 256>>>(thermal_conductivity_array_d, aluminium);
	getHeatCapacities<<<(N + 255) / 256, 256>>>(heat_capacity_array_d, aluminium);
	getEnthalpies<<<(N + 255) / 256, 256>>>(enthalpy_array_d, aluminium);

	cudaDeviceSynchronize();

	cudaMemcpy(thermal_conductivity_array_h, thermal_conductivity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(heat_capacity_array_h, heat_capacity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(enthalpy_array_h, enthalpy_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	
	std::ofstream aluminium_file = solution_folder.getCSVFile("aluminium");
	saveSpecies(aluminium_file);
	aluminium_file.close();

	// Argon
	getThermalConductivities<<<(N + 255) / 256, 256>>>(thermal_conductivity_array_d, argon);
	getHeatCapacities<<<(N + 255) / 256, 256>>>(heat_capacity_array_d, argon);
	getEnthalpies<<<(N + 255) / 256, 256>>>(enthalpy_array_d, argon);

	cudaDeviceSynchronize();

	cudaMemcpy(thermal_conductivity_array_h, thermal_conductivity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(heat_capacity_array_h, heat_capacity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(enthalpy_array_h, enthalpy_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	
	std::ofstream argon_file = solution_folder.getCSVFile("argon");
	saveSpecies(argon_file);
	argon_file.close();

	// Nickel
	getThermalConductivities<<<(N + 255) / 256, 256>>>(thermal_conductivity_array_d, nickel);
	getHeatCapacities<<<(N + 255) / 256, 256>>>(heat_capacity_array_d, nickel);
	getEnthalpies<<<(N + 255) / 256, 256>>>(enthalpy_array_d, nickel);

	cudaDeviceSynchronize();

	cudaMemcpy(thermal_conductivity_array_h, thermal_conductivity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(heat_capacity_array_h, heat_capacity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(enthalpy_array_h, enthalpy_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	
	std::ofstream nickel_file = solution_folder.getCSVFile("nickel");
	saveSpecies(nickel_file);
	aluminium_file.close();

	// Nickel Aluminide
	getThermalConductivities<<<(N + 255) / 256, 256>>>(thermal_conductivity_array_d, nickel_aluminide);
	getHeatCapacities<<<(N + 255) / 256, 256>>>(heat_capacity_array_d, nickel_aluminide);
	getEnthalpies<<<(N + 255) / 256, 256>>>(enthalpy_array_d, nickel_aluminide);

	cudaDeviceSynchronize();

	cudaMemcpy(thermal_conductivity_array_h, thermal_conductivity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(heat_capacity_array_h, heat_capacity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(enthalpy_array_h, enthalpy_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	
	std::ofstream nickel_aluminide_file = solution_folder.getCSVFile("nickel_aluminide");
	saveSpecies(nickel_aluminide_file);
	nickel_aluminide_file.close();

	deallocateMemory<<<1,1>>>();

	cudaFree(thermal_conductivity_array_d);
	cudaFree(heat_capacity_array_d);
	cudaFree(enthalpy_array_d);

	return 0;
}

__host__ void saveSpecies(std::ofstream &file)
{
	file << "Temperature (K)," << "Heat Capacity (J / mol. - K)," << "Enthalpy (J / mol.)," << "Thermal Conductivity (W / m - K)\n";

	for (unsigned int i = 0; i < N; i++)
	{
		file << GET_TEMPERATURE(i) << ',' << heat_capacity_array_h[i] << ',' << enthalpy_array_h[i] << ',' <<
		thermal_conductivity_array_h[i] << '\n';
	}
}

__global__ void getThermalConductivities(double *thermal_conductivities, Species *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < N) thermal_conductivities[i] = species_ptr->getThermalConductivity(GET_TEMPERATURE(i));
}

__global__ void getHeatCapacities(double *heat_capacities, Species *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < N) heat_capacities[i] = species_ptr->getHeatCapacity(GET_TEMPERATURE(i)) * species_ptr->getMolarMass();
}

__global__ void getEnthalpies(double *enthalpies, Species *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < N) enthalpies[i] = species_ptr->getEnthalpy(GET_TEMPERATURE(i)) * species_ptr->getMolarMass();
}

__global__ void getThermalConductivities(double *thermal_conductivities, IdealGas *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < N) thermal_conductivities[i] = species_ptr->getThermalConductivity(GET_TEMPERATURE(i));
}

__global__ void getHeatCapacities(double *heat_capacities, IdealGas *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < N) heat_capacities[i] = species_ptr->getHeatCapacity(GET_TEMPERATURE(i)) * species_ptr->getMolarMass();
}

__global__ void getEnthalpies(double *enthalpies, IdealGas *species_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < N) enthalpies[i] = species_ptr->getEnthalpy(GET_TEMPERATURE(i)) * species_ptr->getMolarMass();
}