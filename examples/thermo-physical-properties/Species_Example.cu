#include <iostream>

#include "thermo-physical-properties/Species.cuh"

#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 1000.0	// K

__global__ void calcThermoPhysicalProperties(
	int n,
	double *temperatures,
	double *enthalpies,
	double *heat_capacities,
	Species *species
) {
	unsigned int i = threadIdx.x;

	if (i < n)
	{
		// densities[i] = species->getDensity(temperatures[i]);
		enthalpies[i] = species->getEnthalpy(temperatures[i]);
		heat_capacities[i] = species->getHeatCapacity(temperatures[i]);
		// thermal_conductivities[i] = species->getThermalConductivity(temperatures[i]);
	}
}

int main(int argc, char const *argv[])
{
	unsigned int n = 512;
	
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

	Species aluminium, *aluminium_ptr;
	{
		Enthalpy enthalpy_solid_Al, enthalpy_liquid_Al;

		enthalpy_solid_Al.assignCoefficients(
			28.08920,
			-5.414849,
			8.560423,
			3.427370,
			-0.277375,
			-9.147187
		);

		enthalpy_liquid_Al.assignCoefficients(
			31.75104,
			3.935826E-8,
			-1.786515E-8,
			2.694171E-9,
			5.480037E-9,
			-0.945684
		);

		ThermalConductivity conductivity_solid_Al, conductivity_liquid_Al;

		conductivity_solid_Al.assignCoefficients(248.0, -0.067, 0.0);
		conductivity_liquid_Al.assignCoefficients(33.9, 7.892E-2, -2.099E-5);
		
		Phase phases[2];
		
		phases[0].initialize(
			2700,
			enthalpy_solid_Al,
			conductivity_solid_Al,
			273,
			933,
			10
		);

		phases[1].initialize(
			2375.0,
			enthalpy_liquid_Al,
			conductivity_liquid_Al,
			933.47,
			INFINITY,
			10
		);
		
		aluminium.initialize(2, phases, 26.9815386E-3);

		cudaMalloc(&aluminium_ptr, sizeof(Species));
		cudaMemcpy(aluminium_ptr, &aluminium, sizeof(Species), cudaMemcpyHostToDevice);
	}

	// std::cout << "Size of species pointer " << sizeof(aluminium_ptr) << "\n";

	calcThermoPhysicalProperties<<<1,n>>>(n, d_temperatures, d_enthalpies, d_heat_capacities, aluminium_ptr);	

	cudaDeviceSynchronize();

	cudaMemcpy(h_densities, d_densities, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_enthalpies, d_enthalpies, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_heat_capacities, d_heat_capacities, n * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_thermal_conductivities, d_thermal_conductivities, n * sizeof(double), cudaMemcpyDeviceToHost);

	// for (int i = 0; i < n; i++)

		// std::cout << "Temperature : " << h_temperatures[i] << " K\t" << "Density : " << h_densities[i] << " kg / m3\t" <<
		// "Heat Capacity : " << h_heat_capacities[i] << " J / mol. - K\t" << "Enthalpy : " << h_enthalpies[i] << " J / mol.\n" <<
		// "Thermal Conductivity : " << h_thermal_conductivities[i] << " W / m - K\n";

	aluminium.deallocateMemory();

	cudaFree(d_temperatures);
	cudaFree(d_densities);
	cudaFree(d_heat_capacities);
	cudaFree(d_enthalpies);
	cudaFree(d_thermal_conductivities);
	cudaFree(aluminium_ptr);

	return 0;
}