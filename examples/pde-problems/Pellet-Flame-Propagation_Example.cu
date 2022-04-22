#include "species/Aluminium.cuh"
#include "species/Argon.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.cuh"
#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "thermo-physical-properties/Packed-Pellet.cuh"
#include "pde-problems/Pellet-Flame-Propagation.cuh"
#include "pde-problems/Core-Shell-Diffusion.cuh"

#include <iostream>

#define MAX_ITER 1E3
#define Dt 0.000001

#define N 1001
#define M 1001

__device__ PelletFlamePropagation::FlamePropagation *flame_propagation;
__device__ ArrheniusDiffusivityModel *diffusivity_model;

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

	CoreShellDIffusion::setGridSize(N);
	CoreShellDIffusion::setTimeStep(Dt);

	PackedPellet::setPelletDimensions(6.35E-3, 6.35E-3);
	PackedPellet::setDegassingFluid(argon);
	PackedPellet::setHeatLossParameters(0, 0, 0);
	PackedPellet::setTemperatureParameters(298, 933);

	PelletFlamePropagation::setNumGridPoints(M);
	PelletFlamePropagation::setImplicitness(1, 1);
	PelletFlamePropagation::setInitialConditions(1500, 0.1 * 6.35E-3);
	PelletFlamePropagation::setInfinitesimalChangeInTemperature(0.001);
	
	flame_propagation = new PelletFlamePropagation::FlamePropagation(0.7);
	diffusivity_model = new ArrheniusDiffusivityModel(2.56E-6, 102.191E3);

	flame_propagation->setDiffusivityModel(diffusivity_model);
}

__global__ void allocateParticleMemory()
{
	flame_propagation->allocateParticleMemory(blockIdx.x);
}

__global__ void deallocateParticleMemory()
{
	flame_propagation->deallocateParticleMemory(blockIdx.x);
}

__global__ void deallocateMemory()
{
	delete diffusivity_model;
	delete flame_propagation;

	unloadAluminium();
	unloadArgon();
	unloadNickel();
	unloadNickelAluminide();
}

__global__ void setArrayAddresses(double *temperature_array, double *concentration_array_A, double *concentration_array_B)
{
	flame_propagation->setArrayAddresses(temperature_array, concentration_array_A, concentration_array_B);
}

__global__ void setInitialConditions()
{
	flame_propagation->initializePellet(blockIdx.x, threadIdx.x);
}

__global__ void setMainParticleDiffusivity()
{
	flame_propagation->setParticleDiffusivity(0, blockIdx.x);
}

__global__ void setMainParticleEquations()
{
	flame_propagation->setParticleEquation(0, blockIdx.x, threadIdx.x);
}

__global__ void solveMainParticleEquations()
{
	flame_propagation->solveParticleEquation(0, blockIdx.x, threadIdx.x);
}

__global__ void updateMainParticles()
{
	flame_propagation->collectEnthalpy(blockIdx.x);
	flame_propagation->updateParticleMassFraction(0, blockIdx.x);
	flame_propagation->updateThermalConductivity(blockIdx.x);
}

__host__ __forceinline__ void evolveMainParticles()
{
	setMainParticleDiffusivity<<<M, 1>>>();
	setMainParticleEquations<<<M, N>>>();
	solveMainParticleEquations<<<M, 2>>>();
	updateMainParticles<<<M, 1>>>();
}

__global__ void setAuxiliaryParticleDiffusivity()
{
	flame_propagation->setParticleDiffusivity(1 + threadIdx.x, blockIdx.x);
}

__global__ void setAuxiliaryParticleEquations()
{
	flame_propagation->setParticleEquation(1, blockIdx.x, threadIdx.x);
	flame_propagation->setParticleEquation(2, blockIdx.x, threadIdx.x);
}

__global__ void solveAuxiliaryParticleEquations()
{
	flame_propagation->solveParticleEquation(
		1 + (threadIdx.x > 1) * 1,
		blockIdx.x,
		threadIdx.x - (threadIdx.x > 1) * 2);
}

__global__ void updateAuxiliaryParticles()
{
	flame_propagation->updateParticleMassFraction(1 + threadIdx.x, blockIdx.x);
}

__host__ __forceinline__ void evolveAuxiliaryParticles()
{
	setAuxiliaryParticleDiffusivity<<<M, 2>>>();
	setAuxiliaryParticleEquations<<<M, N>>>();
	solveAuxiliaryParticleEquations<<<M, 4>>>();
	updateAuxiliaryParticles<<<M, 2>>>();
}

__global__ void setTemperatureEquations()
{
	flame_propagation->setUpEquation(blockIdx.x);
}

__global__ void solveTemperatureEquations()
{
	flame_propagation->solveEquation();
}

__host__ __forceinline__ void iterate()
{
	evolveAuxiliaryParticles();

	setTemperatureEquations<<<M, 1>>>();
	solveTemperatureEquations<<<1, 1>>>();

	evolveMainParticles();
}

int main(int argc, char const *argv[])
{
	allocateMemory<<<1,1>>>();

	allocateParticleMemory<<<M, 1>>>();

	double *temperature_array;
	double *concentration_array_A;
	double *concentration_array_B;

	cudaMalloc(&temperature_array, M * sizeof(double));
	cudaMalloc(&concentration_array_A, 3 * M * N * sizeof(double));
	cudaMalloc(&concentration_array_B, 3 * M * N * sizeof(double));

	setArrayAddresses<<<1,1>>>(temperature_array, concentration_array_A, concentration_array_B);

	cudaDeviceSynchronize();

	setInitialConditions<<<M,N>>>();

	evolveMainParticles();

	cudaDeviceSynchronize();

	std::cout << "Initialized\n";

	for (size_t i = 0; i < MAX_ITER; i++) 
	{
		iterate();
		if (i % 1000 == 0) 
		{
			cudaDeviceSynchronize();
			std::cout << "Iterations completed : " << i << "\n";
		}
	}

	cudaDeviceSynchronize();

	cudaFree(temperature_array);
	cudaFree(concentration_array_A);
	cudaFree(concentration_array_B);

	deallocateParticleMemory<<<M, 1>>>();
	deallocateMemory<<<1,1>>>();

	return 0;
}
