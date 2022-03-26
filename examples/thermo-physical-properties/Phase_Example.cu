#include <iostream>

#include <thermo-physical-properties/Phase.cuh>

#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 1000.0	// K

__global__ void calcThermoPhysicalProperties(
	unsigned int n,
	double *temperatures,
	double *densities,
	double *enthalpies,
	double *heat_capacities,
	double *thermal_conductivities,
	Phase *species
) {
	unsigned int i = threadIdx.x;

	if (i < n)
	{
		densities[i] = species->getDensity(temperatures[i]);
		enthalpies[i] = species->getEnthalpy(temperatures[i]);
		heat_capacities[i] = species->getHeatCapacity(temperatures[i]);
		thermal_conductivities[i] = species->getThermalConductivity(temperatures[i]);
	}
}

int main(int argc, char const *argv[])
{
	unsigned int n = 500;
	
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

	Phase *d_solid_Al;

	{
		Enthalpy enthalpy_solid_Al;
		enthalpy_solid_Al.assignCoefficients(
			28.08920,
			-5.414849,
			8.560423,
			3.427370,
			-0.277375,
			-9.147187
		);

		ThermalConductivity conductivity_solid_Al;
		conductivity_solid_Al.assignCoefficients(248.0, -0.067, 0.0);
		
		Phase solid_Al;
		solid_Al.initialize(
			2700,
			enthalpy_solid_Al,
			conductivity_solid_Al,
			273,
			933,
			10
		);

		cudaMalloc(&d_solid_Al, sizeof(solid_Al));
		cudaMemcpy(d_solid_Al, &solid_Al, sizeof(solid_Al), cudaMemcpyHostToDevice);
	}

	calcThermoPhysicalProperties<<<1,n>>>(n, d_temperatures, d_densities, d_enthalpies, d_heat_capacities, d_thermal_conductivities, d_solid_Al);	

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
	cudaFree(d_solid_Al);

	return 0;
}