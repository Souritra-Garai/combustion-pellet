#include "thermo-physical-properties/Species.cuh"

#include <iostream>

#define N 1000
#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 1000.0	// K
#define GET_TEMPERATURE(i) (TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i) / ((double) N - 1.0))

__device__ Species *species;
__device__ Phase *species_phases;

__global__ void allocateMemory()
{
	Enthalpy species_enthalpy_1(
		28.08920,
		-5.414849,
		8.560423,
		3.427370,
		-0.277375,
		-9.147187
	);

	Enthalpy species_enthalpy_2(
		31.75104,
		3.935826E-8,
		-1.786515E-8,
		2.694171E-9,
		5.480037E-9,
		-0.945684
	);

	ThermalConductivity species_thermal_conductivity_1(248.0, -0.067, 0.0);
	ThermalConductivity species_thermal_conductivity_2(33.9, 7.89E-2, -2.099E-5);

	Phase phase_1(
		2700,
		species_enthalpy_1,
		species_thermal_conductivity_1,
		273.15,
		933.47,
		10
	);

	Phase phase_2(
		2700,
		species_enthalpy_2,
		species_thermal_conductivity_2,
		933.47,
		INFINITY,
		10
	);

	species_phases = new Phase[2]{phase_1, phase_2};

	species = new Species(
		2,
		species_phases,
		26.9815386E-3
	);
}

__global__ void deallocateMemory()
{
	delete species;
	delete [] species_phases;
}

__global__ void calcEnthalpies(double *enthalpy_array, double *heat_capacity_array, double *thermal_conductivity_array)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < N) 
	{
		double T = GET_TEMPERATURE(i);
		enthalpy_array[i] = species->getEnthalpy(T);
		heat_capacity_array[i] = species->getHeatCapacity(T);
		thermal_conductivity_array[i] = species->getThermalConductivity(T);
	}
}

int main(int argc, char const *argv[])
{
	double enthalpy_array_h[N], heat_capacity_array_h[N], thermal_conductivity_array_h[N];

	double *enthalpy_array_d, *heat_capacity_array_d, *thermal_conductivity_array_d;
	cudaMalloc(&enthalpy_array_d, N*sizeof(double));
	cudaMalloc(&heat_capacity_array_d, N*sizeof(double));
	cudaMalloc(&thermal_conductivity_array_d, N*sizeof(double));

	allocateMemory<<<1,1>>>();

	cudaDeviceSynchronize();

	calcEnthalpies<<<(N+255)/256, 256>>>(enthalpy_array_d, heat_capacity_array_d, thermal_conductivity_array_d);

	cudaDeviceSynchronize();
	
	cudaMemcpy(enthalpy_array_h, enthalpy_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(heat_capacity_array_h, heat_capacity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(thermal_conductivity_array_h, thermal_conductivity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);

	for (int i = 0; i < N; i++)

		std::cout << "Temperature : " << GET_TEMPERATURE(i) << " K\t" << "Heat Capacity : " << heat_capacity_array_h[i] << " J / kg - K\t" <<
		"Enthalpy : " << enthalpy_array_h[i] << " J / kg\t" << "Thermal Conductivity : " << thermal_conductivity_array_h[i] << " W / m - K\n";

	deallocateMemory<<<1,1>>>();

	cudaFree(thermal_conductivity_array_d);
	cudaFree(heat_capacity_array_d);
	cudaFree(enthalpy_array_d);

	return 0;
}
