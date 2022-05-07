#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "thermo-physical-properties/Packed-Pellet.cuh"

#include "species/Aluminium.cuh"
#include "species/Argon.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include <iostream>

#define N 1000
#define TEMPERATURE_LOWER_BOUND	250.0	// K
#define TEMPERATURE_UPPER_BOUND 1000.0	// K
#define GET_TEMPERATURE(i) (TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i) / ((double) N - 1.0))

__device__ PackedPellet::Pellet *pellet;
__device__ CoreShellParticle::Particle *Ni_clad_Al_particle;

__global__ void allocateMemory()
{
	loadAluminium();
	loadArgon();
	loadNickel();
	loadNickelAlumnide();

	CoreShellParticle::initialize(
		aluminium,
		nickel,
		nickel_aluminide,
		32.5E-6,
		39.5E-6
	);

	PackedPellet::setPelletDimensions(6.35E-3, 6.35E-3);
	PackedPellet::setDegassingFluid(argon);
	PackedPellet::setHeatLossParameters(0, 0, 0);
	PackedPellet::setTemperatureParameters(298, 933);

	Ni_clad_Al_particle = new CoreShellParticle::Particle();
	Ni_clad_Al_particle->setMassFractions(0.3, 0.3, 0.4);
	pellet = new PackedPellet::Pellet(0.7);
}

__global__ void deallocateMemory()
{
	delete pellet;
	delete Ni_clad_Al_particle;

	unloadAluminium();
	unloadArgon();
	unloadNickel();
	unloadNickelAluminide();
}

__global__ void calcThermalConductivity(double *thermal_conductivity_array)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < N) 
	{
		double T = GET_TEMPERATURE(i);
		thermal_conductivity_array[i] = pellet->getThermalConductivity(Ni_clad_Al_particle, T);
	}
}

int main(int argc, char const *argv[])
{
	double thermal_conductivity_array_h[N];

	double *thermal_conductivity_array_d;
	cudaMalloc(&thermal_conductivity_array_d, N*sizeof(double));

	allocateMemory<<<1,1>>>();

	cudaDeviceSynchronize();

	calcThermalConductivity<<<(N+255)/256, 256>>>(thermal_conductivity_array_d);

	cudaDeviceSynchronize();
	
	cudaMemcpy(thermal_conductivity_array_h, thermal_conductivity_array_d, N * sizeof(double), cudaMemcpyDeviceToHost);

	for (int i = 0; i < N; i++)

		std::cout << "Temperature : " << GET_TEMPERATURE(i) << " K\t" << "Thermal Conductivity : " << thermal_conductivity_array_h[i] << " W / m - K\n";

	deallocateMemory<<<1,1>>>();

	cudaFree(thermal_conductivity_array_d);
	
	return 0;
}
