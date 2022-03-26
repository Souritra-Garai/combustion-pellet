#include "thermo-physical-properties/Enthalpy.cuh"

#include <iostream>

#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 1000.0	// K

__global__ void calcEnthalpies(unsigned int n, double *temperatures, double *enthalpies, double *heat_capacities, Enthalpy *material_enthalpy)
{
	unsigned int i = threadIdx.x;

	if (i < n) 
	{
		enthalpies[i] = material_enthalpy->getEnthalpy(temperatures[i]);
		heat_capacities[i] = material_enthalpy->getHeatCapacity(temperatures[i]);
	}
}

int main(int argc, char const *argv[])
{
	unsigned int n = 1001;
	
	double h_temperatures[n], h_enthalpies[n], h_heat_capacities[n];
	double *d_temperatures, *d_enthalpies, *d_heat_capacities;

	cudaMalloc(&d_temperatures, n * sizeof(double));
	cudaMalloc(&d_enthalpies, n * sizeof(double));
	cudaMalloc(&d_heat_capacities, n * sizeof(double));

	for (unsigned int i = 0; i < n; i++)
	
		h_temperatures[i] = TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i / ((double) n - 1));

	cudaMemcpy(d_temperatures, h_temperatures, n * sizeof(double), cudaMemcpyHostToDevice);

	Enthalpy *d_material_enthalpy;

	{
		Enthalpy material_enthalpy;
		material_enthalpy.assignCoefficients(
			28.08920,
			-5.414849,
			8.560423,
			3.427370,
			-0.277375,
			-9.147187
		);

		cudaMalloc(&d_material_enthalpy, sizeof(material_enthalpy));
		cudaMemcpy(d_material_enthalpy, &material_enthalpy, sizeof(material_enthalpy), cudaMemcpyHostToDevice);
	}

	calcEnthalpies<<<1,n>>>(n, d_temperatures, d_enthalpies, d_heat_capacities, d_material_enthalpy);	

	cudaDeviceSynchronize();

	cudaMemcpy(h_enthalpies, d_enthalpies, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_heat_capacities, d_heat_capacities, n * sizeof(double), cudaMemcpyDeviceToHost);

	for (int i = 0; i < n; i++)

		std::cout << "Temperature : " << h_temperatures[i] << " K\t" << "Heat Capacity : " << h_heat_capacities[i] << " J / mol. - K\t" <<
		"Enthalpy : " << h_enthalpies[i] << " J / mol.\n";

	cudaFree(d_temperatures);
	cudaFree(d_heat_capacities);
	cudaFree(d_enthalpies);

	cudaFree(d_material_enthalpy);

	return 0;
}
