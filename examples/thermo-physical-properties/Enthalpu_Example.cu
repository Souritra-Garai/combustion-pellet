#include "thermo-physical-properties/Enthalpy.cuh"

#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 1000.0	// K

__global__ void calcEnthalpies(unsigned int n, double *temperatures, double *enthalpies, double *heat_capacities, double coefficients[6])
{
	unsigned int i = threadIdx.x;

	enthalpies[i] = getStandardEnthalpy(temperatures[i], coefficients);
	heat_capacities[i] = getHeatCapacity(temperatures[i], coefficients);
}

int main(int argc, char const *argv[])
{
	unsigned int n = 1000;
	
	double h_temperatures[n], h_enthalpies[n], h_heat_capacities[n];
	double d_temperatures[n], d_enthalpies[n], d_heat_capacities[n];

	for (unsigned int i = 0; i < n; i++)
	
		h_temperatures[i] = TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i / ((double) n));

	

	return 0;
}
