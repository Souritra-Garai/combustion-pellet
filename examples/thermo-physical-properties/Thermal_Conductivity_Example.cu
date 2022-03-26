#include "thermo-physical-properties/Thermal_Conductivity.cuh"

#include <iostream>

#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 1000.0	// K

__global__ void calcThermalConductivities(unsigned int n, double *temperatures, double *thermal_conductivities, ThermalConductivity *material_conductivity)
{
	unsigned int i = threadIdx.x;

	if (i < n) 
	{
		thermal_conductivities[i] = material_conductivity->getThermalConductivity(temperatures[i]);
	}
}

int main(int argc, char const *argv[])
{
	unsigned int n = 1000;
	
	double h_temperatures[n], h_thermal_conductivities[n];
	double *d_temperatures, *d_thermal_conductivities;

	cudaMalloc(&d_temperatures, n * sizeof(double));
	cudaMalloc(&d_thermal_conductivities, n * sizeof(double));

	for (unsigned int i = 0; i < n; i++)
	
		h_temperatures[i] = TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i / ((double) n - 1));

	cudaMemcpy(d_temperatures, h_temperatures, n * sizeof(double), cudaMemcpyHostToDevice);

	ThermalConductivity *d_material_conductivity;
	
	{
		ThermalConductivity material_conductivity;
		material_conductivity.assignCoefficients(248.0, -0.067, 0.0);

		cudaMalloc(&d_material_conductivity, sizeof(material_conductivity));
		cudaMemcpy(d_material_conductivity, &material_conductivity, sizeof(material_conductivity), cudaMemcpyHostToDevice);
	}

	calcThermalConductivities<<<1,n>>>(n, d_temperatures, d_thermal_conductivities, d_material_conductivity);	

	cudaDeviceSynchronize();

	cudaMemcpy(h_thermal_conductivities, d_thermal_conductivities, n * sizeof(double), cudaMemcpyDeviceToHost);

	for (int i = 0; i < n; i++)

		std::cout << "Temperature : " << h_temperatures[i] << " K\t" << "Thermal Conductivity : " << h_thermal_conductivities[i] << " W / m - K\n";

	cudaFree(d_temperatures);
	cudaFree(d_thermal_conductivities);

	cudaFree(d_material_conductivity);

	return 0;
}
