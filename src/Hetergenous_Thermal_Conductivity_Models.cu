#include "species/Aluminium.cuh"
#include "species/Argon.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "thermo-physical-properties/Thermal_Conductivity_Pellet.cuh"

#include "utilities/File_Generator.hpp"

#include <iostream>
#include <boost/program_options.hpp>

#define N 1000
#define GET_PARTICLE_VOLUME_FRACTIONS(i) (((double) i) / (double) N)

double temperature = 298.15;	// K
double core_radius = 32.5E-6;	// m
double overall_radius = 39.5E-6;	// m

__device__ double particle_thermal_conductivity;
__device__ double gas_thermal_conductivity;

__device__ double *thermal_conductivity_ME1;
__device__ double *thermal_conductivity_EMT;
__device__ double *thermal_conductivity_MEB;
__device__ double *thermal_conductivity_CC;
__device__ double *thermal_conductivity_ME2;

__global__ void allocateMemory(
	double core_radius, 
	double overall_radius,
	double temperature,
	double *thermal_conductivity_ME1,
	double *thermal_conductivity_EMT,
	double *thermal_conductivity_MEB,
	double *thermal_conductivity_CC,
	double *thermal_conductivity_ME2
);

__global__ void deallocateMemory();
__global__ void calculateThermalConductivities();

__host__ int setConfiguration(int argc, char const *argv[]);

int main(int argc, char const *argv[])
{
	int flag = setConfiguration(argc, argv);

	if (flag != 2) return flag;
	
	double *thermal_conductivity_ME1;
	double *thermal_conductivity_EMT;
	double *thermal_conductivity_MEB;
	double *thermal_conductivity_CC;
	double *thermal_conductivity_ME2;
	
	cudaMalloc(&thermal_conductivity_ME1, N * sizeof(double));
	cudaMalloc(&thermal_conductivity_EMT, N * sizeof(double));
	cudaMalloc(&thermal_conductivity_MEB, N * sizeof(double));
	cudaMalloc(&thermal_conductivity_CC,  N * sizeof(double));
	cudaMalloc(&thermal_conductivity_ME2, N * sizeof(double));

	double thermal_conductivity_ME1_host[N];
	double thermal_conductivity_EMT_host[N];
	double thermal_conductivity_MEB_host[N];
	double thermal_conductivity_CC_host[N];
	double thermal_conductivity_ME2_host[N];
	
	allocateMemory<<<1,1>>>(
		core_radius,
		overall_radius,
		temperature,
		thermal_conductivity_ME1,
		thermal_conductivity_EMT,
		thermal_conductivity_MEB,
		thermal_conductivity_CC,
		thermal_conductivity_ME2
	);

	cudaDeviceSynchronize();

	FileGenerator folder;
	
	std::ofstream file = folder.getCSVFile("thermal_conductivity_pellet");
	std::ofstream config_file = folder.getTXTFile("configuration");

	config_file << "Temperature :\t" << temperature << "\tK\n\n";

	double particle_thermal_conductivity, gas_thermal_conductivity;
	cudaMemcpyFromSymbol(&particle_thermal_conductivity, ::particle_thermal_conductivity, sizeof(double));
	cudaMemcpyFromSymbol(&gas_thermal_conductivity, ::gas_thermal_conductivity, sizeof(double));

	config_file << "Particle Thermal Conductivity :\t" << particle_thermal_conductivity << "\tW / m - K\n";
	config_file << "Gas Thermal Conductivity :\t" << gas_thermal_conductivity << "\tW / m - K\n\n\n";

	CoreShellParticle::printConfiguration(config_file);
	config_file.close();

	file << "phi,ME1,EMT,MEB,CC,ME2\n";

	calculateThermalConductivities<<<5, N>>>();
	cudaDeviceSynchronize();

	cudaMemcpy(thermal_conductivity_ME1_host, thermal_conductivity_ME1, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(thermal_conductivity_EMT_host, thermal_conductivity_EMT, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(thermal_conductivity_MEB_host, thermal_conductivity_MEB, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(thermal_conductivity_CC_host,  thermal_conductivity_CC,  N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(thermal_conductivity_ME2_host, thermal_conductivity_ME2, N * sizeof(double), cudaMemcpyDeviceToHost);

	for (size_t i = 0; i < N; i++)

		file << GET_PARTICLE_VOLUME_FRACTIONS(i+1) << ',' << 
		thermal_conductivity_ME1_host[i] << ',' <<
		thermal_conductivity_EMT_host[i] << ',' <<
		thermal_conductivity_MEB_host[i] << ',' <<
		thermal_conductivity_CC_host[i]  << ',' <<
		thermal_conductivity_ME2_host[i] << '\n';

	file.close();

	cudaFree(thermal_conductivity_ME1);
	cudaFree(thermal_conductivity_EMT);
	cudaFree(thermal_conductivity_MEB);
	cudaFree(thermal_conductivity_CC);
	cudaFree(thermal_conductivity_ME2);

	deallocateMemory<<<1,1>>>();
	cudaDeviceSynchronize();

    return 0;
}

__global__ void allocateMemory(
	double core_radius,
	double overall_radius,
	double temperature,
	double *thermal_conductivity_ME1,
	double *thermal_conductivity_EMT,
	double *thermal_conductivity_MEB,
	double *thermal_conductivity_CC,
	double *thermal_conductivity_ME2
)
{
	loadAluminium();
	loadArgon();
	loadNickel();
	loadNickelAlumnide();

	CoreShellParticle::initialize(
		aluminium,
		nickel,
		nickel_aluminide,
		core_radius,
		overall_radius
	);
	
	particle_thermal_conductivity = CoreShellParticle::Particle().getThermalConductivity(temperature);
	gas_thermal_conductivity = argon->getThermalConductivity(temperature);

	::thermal_conductivity_ME1 = thermal_conductivity_ME1;
	::thermal_conductivity_EMT = thermal_conductivity_EMT;
	::thermal_conductivity_MEB = thermal_conductivity_MEB;
	::thermal_conductivity_CC  = thermal_conductivity_CC;
	::thermal_conductivity_ME2 = thermal_conductivity_ME2;
}

__global__ void deallocateMemory()
{
	unloadAluminium();
	unloadNickel();
	unloadArgon();
	unloadNickelAluminide();
}

__global__ void calculateThermalConductivities()
{
	size_t i = threadIdx.x;

	if (i < N)
	{
		switch (blockIdx.x)
		{
			case 0:
				
				thermal_conductivity_ME1[i] = PelletThermalConductivity::getThermalConductivityME1(
					1.0 - GET_PARTICLE_VOLUME_FRACTIONS(i+1),
					gas_thermal_conductivity,
					particle_thermal_conductivity
				);
				
				break;

			case 1:
				
				thermal_conductivity_EMT[i] = PelletThermalConductivity::getThermalConductivityEMT(
					1.0 - GET_PARTICLE_VOLUME_FRACTIONS(i+1),
					gas_thermal_conductivity,
					particle_thermal_conductivity
				);
				
				break;

			case 2:
			
				thermal_conductivity_MEB[i] = PelletThermalConductivity::getThermalConductivityMEB(
					1.0 - GET_PARTICLE_VOLUME_FRACTIONS(i+1),
					gas_thermal_conductivity,
					particle_thermal_conductivity
				);
				
				break;

			case 3:
			
				thermal_conductivity_CC[i] = PelletThermalConductivity::getThermalConductivityCC(
					1.0 - GET_PARTICLE_VOLUME_FRACTIONS(i+1),
					gas_thermal_conductivity,
					particle_thermal_conductivity
				);
				
				break;

			case 4:
			
				thermal_conductivity_ME2[i] = PelletThermalConductivity::getThermalConductivityME2(
					1.0 - GET_PARTICLE_VOLUME_FRACTIONS(i+1),
					gas_thermal_conductivity,
					particle_thermal_conductivity
				);
				
				break;
			
			default:

				break;
		}
	}
}

int setConfiguration(int argc, const char * argv[])
{
	try
	{
		boost::program_options::options_description options_description_obj("Usage: Core-Shell-Diffusion [options]\nOptions");

		options_description_obj.add_options()
			("help",	"Produce help message")
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