#include "thermo-physical-properties/Thermal_Conductivity.cuh"

#include <iostream>

#define N 1000
#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 1000.0	// K
#define GET_TEMPERATURE(i) (TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i) / ((double) N - 1.0))

__device__ ThermalConductivity *species_thermal_conductivity;

__global__ void allocateMemory()
{
	species_thermal_conductivity = new ThermalConductivity(248.0, -0.067, 0.0);
}

__global__ void deallocateMemory()
{
	delete species_thermal_conductivity;
}

__global__ void calcEnthalpies(double *thermal_conductivity_array)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < N) 
	{
		double T = GET_TEMPERATURE(i);
		thermal_conductivity_array[i] = species_thermal_conductivity->getThermalConductivity(T);
	}
}

int main(int argc, char const *argv[])
{
	double thermal_conductivity_array_h[N];

	double *thermal_conductivity_array_d;
	cudaMalloc(&thermal_conductivity_array_d, N*sizeof(double));

	allocateMemory<<<1,1>>>();

	cudaDeviceSynchronize();

	calcEnthalpies<<<(N+255)/256, 256>>>(thermal_conductivity_array_d);

	cudaDeviceSynchronize();
	
	cudaMemcpy(thermal_conductivity_array_h, thermal_conductivity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);

	for (int i = 0; i < N; i++)

		std::cout << "Temperature : " << GET_TEMPERATURE(i) << " K\t" << "Thermal Conductivity : " << thermal_conductivity_array_h[i] << " W / m - K\n";

	deallocateMemory<<<1,1>>>();

	cudaFree(thermal_conductivity_array_d);

	return 0;
}
