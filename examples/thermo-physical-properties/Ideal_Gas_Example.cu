#include <iostream>

#include "thermo-physical-properties/Ideal_Gas.cuh"

#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 1000.0	// K

__global__ void calcThermoPhysicalProperties(
	unsigned int n,
	double *temperatures,
	double *densities,
	double *enthalpies,
	double *heat_capacities,
	double *thermal_conductivities,
	IdealGas *gas
) {
	unsigned int i = threadIdx.x;

	if (i < n)
	{
		densities[i] = gas->getDensity(temperatures[i]);
		enthalpies[i] = gas->getEnthalpy(temperatures[i]);
		heat_capacities[i] = gas->getHeatCapacity(temperatures[i]);
		thermal_conductivities[i] = gas->getThermalConductivity(temperatures[i]);
	}
}

int main(int argc, char const *argv[])
{
	unsigned int n = 1001;
	
	double h_temperatures[n], h_densities[n], h_enthalpies[n], h_heat_capacities[n], h_thermal_conductivities[n];
	double *d_temperatures, *d_densities, *d_enthalpies, *d_heat_capacities, *d_thermal_conductivities;

	cudaMalloc(&d_temperatures, n * sizeof(double));
	cudaMalloc(&d_densities, n * sizeof(double));
	cudaMalloc(&d_enthalpies, n * sizeof(double));
	cudaMalloc(&d_heat_capacities, n * sizeof(double));
	cudaMalloc(&d_thermal_conductivities, n * sizeof(double));

	for (unsigned int i = 0; i < n; i++)
	
		h_temperatures[i] = TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i / ((double) n - 1.0));

	cudaMemcpy(d_temperatures, h_temperatures, n * sizeof(double), cudaMemcpyHostToDevice);

	IdealGas *d_argon;

	{
		IdealGas argon;

		Enthalpy enthalpy_Ar;
		enthalpy_Ar.assignCoefficients(
			20.78600,
			2.825911E-7,
			-1.464191E-7,
			1.092131E-8,
			-3.661371E-8,
			-6.197350
		);

		ThermalConductivity conductivity_Ar;
		conductivity_Ar.assignCoefficients(1.49E-3, 5.98E-5, -1.92E-8);
		
		argon.initialize(
			39.948E-3,
			enthalpy_Ar,
			conductivity_Ar
		);
		
		cudaMalloc(&d_argon, sizeof(argon));
		cudaMemcpy(d_argon, &argon, sizeof(argon), cudaMemcpyHostToDevice);
	}

	calcThermoPhysicalProperties<<<1,n>>>(n, d_temperatures, d_densities, d_enthalpies, d_heat_capacities, d_thermal_conductivities, d_argon);	

	cudaDeviceSynchronize();

	cudaMemcpy(h_densities, d_densities, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_enthalpies, d_enthalpies, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_heat_capacities, d_heat_capacities, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_thermal_conductivities, d_thermal_conductivities, n * sizeof(double), cudaMemcpyDeviceToHost);

	for (int i = 0; i < n; i++)

		std::cout << "Temperature : " << h_temperatures[i] << " K\t" << "Density : " << h_densities[i] << " kg / m3\t" <<
		"Heat Capacity : " << h_heat_capacities[i] << " J / mol. - K\t" << "Enthalpy : " << h_enthalpies[i] << " J / mol.\n" <<
		"Thermal Conductivity : " << h_thermal_conductivities[i] << " W / m - K\n";

	cudaFree(d_temperatures);
	cudaFree(d_densities);
	cudaFree(d_heat_capacities);
	cudaFree(d_enthalpies);
	cudaFree(d_thermal_conductivities);
	cudaFree(d_argon);

	return 0;
}
