#include <iostream>

#include "species/Aluminium.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "thermo-physical-properties/Core-Shell-Particle.cuh"

#define TEMPERATURE_UPPER_BOUND 2500.0	// K
#define TEMPERATURE_LOWER_BOUND 273.15	// K

#define GET_TEMPERATURE(i, n) (TEMPERATURE_LOWER_BOUND + (TEMPERATURE_UPPER_BOUND - TEMPERATURE_LOWER_BOUND) * ((double) i) / ((double) n - 1.0))

double core_radius = 32.5E-6;
double overall_radius = 39.5E-6;

__global__ void initialize(CoreShellParticle::Particle *particle_ptr)
{
	particle_ptr->initialize();
}

__global__ void getEnthalpies(unsigned int n, double *enthalpies, CoreShellParticle::Particle *particle_ptr)
{
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n) enthalpies[i] = particle_ptr->getEnthalpy(GET_TEMPERATURE(i, n));
}

int main(int argc, char const *argv[])
{
	loadAluminium();
	loadNickel();
	loadNickelAlumnide();

	CoreShellParticle::initialize(
		&aluminium, &nickel, &nickel_aluminide,
		core_radius, overall_radius
	);

    CoreShellParticle::Particle *Ni_clad_Al_particle;
	cudaMalloc(&Ni_clad_Al_particle, sizeof(CoreShellParticle::Particle));

	initialize<<<1,1>>>(Ni_clad_Al_particle);

	unsigned int n = 10;

    double *enthalpies;
	cudaMalloc(&enthalpies, n * sizeof(double));

	getEnthalpies<<<(n+255)/256, 256>>>(n, enthalpies, Ni_clad_Al_particle);

	cudaDeviceSynchronize();

	double enthalpies_host[n];

	cudaMemcpy(enthalpies_host, enthalpies, n * sizeof(double), cudaMemcpyDeviceToHost);

	for (int i = 0; i < n; i++)

		std::cout << "Temperature : " << GET_TEMPERATURE(i, n) << " K\t" << "Enthalpy : " << enthalpies_host[i] << " J / mol.\n";

	cudaFree(enthalpies);

	cudaFree(Ni_clad_Al_particle);

	aluminium.deallocateMemory();
	nickel.deallocateMemory();
	nickel_aluminide.deallocateMemory();

	CoreShellParticle::deallocate();

    return 0;
}