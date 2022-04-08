#include <iostream>

#include "species/Aluminium.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "thermo-physical-properties/Core-Shell-Particle.cuh"
#include "pde-problems/Core-Shell-Diffusion.cuh"

// #include "utilities/Keyboard_Interrupt.hpp"
// #include "utilities/File_Generator.hpp"

#define MAX_ITER 1E6
#define Dt 0.000001

#define N 1001

// ArrheniusDiffusivityModel<long double> Alawieh_diffusivity(2.56E-6, 102.191E3);
// ArrheniusDiffusivityModel<long double> Du_diffusivity(9.54E-8, 26E3);

// void printState(size_t iteration_number, CoreShellDiffusion<long double> &particle);

__device__ CoreShellDIffusion::Diffusion *diffusion_problem;
// __device__ double diffusivity = 1E-10;

__global__ void allocateMemory()
{
	loadAluminium();
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

	diffusion_problem = new CoreShellDIffusion::Diffusion();
}

__global__ void deallocateMemory()
{
	delete diffusion_problem;

	CoreShellDIffusion::deallocate();

	unloadAluminium();
	unloadNickel();
	unloadNickelAluminide();
}

__global__ void initializeCoreShellParticle()
{
	size_t i = blockDim.x * blockIdx.x + threadIdx.x;

	if (i < CoreShellDIffusion::n) 
	{
		diffusion_problem->setInitialState(i);
	}
}

__global__ void initIteration()
{
	diffusion_problem->setCoefficient_1(1E-9);
}

__global__ void iterate()
{
	size_t i = blockDim.x * blockIdx.x + threadIdx.x;

	diffusion_problem->setUpEquations(i);
}

__global__ void solve()
{
	size_t i = blockDim.x * blockIdx.x + threadIdx.x;

	diffusion_problem->solveEquations(i);
}

__global__ void printMass()
{
	double mass_A = diffusion_problem->getAtmMassA();
	double mass_B = diffusion_problem->getAtmMassB();

	printf("Mass A : %e kg\n", mass_A);
	printf("Mass B : %e kg\n", mass_B);
	printf("Mass : %e kg\n\n", mass_A + mass_B);
}

int main(int argc, char const *argv[])
{
	allocateMemory<<<1,1>>>();

	initializeCoreShellParticle<<<1, N>>>();

	double mass;
	cudaMemcpyFromSymbol(&mass, CoreShellParticle::core_mass, sizeof(double));
	printf("Core mass : %e kg\n", mass);
	cudaMemcpyFromSymbol(&mass, CoreShellParticle::shell_mass, sizeof(double));
	printf("Shell mass : %e kg\n", mass);
	cudaMemcpyFromSymbol(&mass, CoreShellParticle::mass, sizeof(double));
	printf("Overall mass : %e kg\n\n", mass);

	printMass<<<1,1>>>();

	cudaDeviceSynchronize();

	initIteration<<<1,1>>>();

	cudaDeviceSynchronize();

	for (size_t i = 0; i < 1000; i++)
	{
		iterate<<<1, N>>>();

		cudaDeviceSynchronize();

		solve<<<1,2>>>();

		cudaDeviceSynchronize();
	}

	printMass<<<1,1>>>();
	cudaDeviceSynchronize();

	deallocateMemory<<<1,1>>>();
	cudaDeviceSynchronize();

    return 0;
}