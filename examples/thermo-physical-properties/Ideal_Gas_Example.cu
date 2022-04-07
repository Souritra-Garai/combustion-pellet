#include "thermo-physical-properties/Ideal_Gas.cuh"

#include <iostream>

#define N 1000
#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 1000.0	// K
#define GET_TEMPERATURE(i) (TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i) / ((double) N - 1.0))

__device__ IdealGas *gas;

__global__ void allocateMemory()
{
	Enthalpy gas_enthalpy(
		20.78600,
		2.825911E-7,
		-1.464191E-7,
		1.092131E-8,
		-3.661371E-8,
		-6.197350
	);

	ThermalConductivity gas_thermal_conductivity(1.49E-3, 5.98E-5, -1.92E-8);

	gas = new IdealGas(
		39.948E-3,
		gas_enthalpy,
		gas_thermal_conductivity
	);
}

__global__ void deallocateMemory()
{
	delete gas;
}

__global__ void calcEnthalpies(double *enthalpy_array, double *heat_capacity_array, double *thermal_conductivity_array)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < N) 
	{
		double T = GET_TEMPERATURE(i);
		enthalpy_array[i] = gas->getEnthalpy(T);
		heat_capacity_array[i] = gas->getHeatCapacity(T);
		thermal_conductivity_array[i] = gas->getThermalConductivity(T);
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
