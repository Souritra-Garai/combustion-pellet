#include "species/Aluminium.cuh"
#include "species/Nickel.cuh"
#include "species/NickelAluminide.cuh"

#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.cuh"
#include "pde-problems/Core-Shell-Diffusion.cuh"

#define MAX_ITER 1E6
#define Dt 0.000001

#define N 1001

__device__ CoreShellDIffusion::Diffusion *diffusion_problem;
__device__ ArrheniusDiffusivityModel *Alawieh_diffusivity;

__global__ void allocateMemory(double *concentration_array_A, double *concentration_array_B)
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
	Alawieh_diffusivity = new ArrheniusDiffusivityModel(2.56E-6, 102.191E3);

	diffusion_problem->setArrayAddresses(concentration_array_A, concentration_array_B);
}

__global__ void deallocateMemory()
{
	delete Alawieh_diffusivity;
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
	diffusion_problem->setDiffusivity(Alawieh_diffusivity->getDiffusivity(1500));
}

__global__ void setUpEquations()
{
	size_t i = blockDim.x * blockIdx.x + threadIdx.x;

	diffusion_problem->setUpEquations(i);
}

__global__ void solve()
{
	size_t i = blockDim.x * blockIdx.x + threadIdx.x;

	diffusion_problem->solveEquations(i);
}

__global__ void updateMassFractions()
{
	diffusion_problem->updateMassFractions();
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
	double *concentration_array_A, *concentration_array_B;
	cudaMalloc(&concentration_array_A, N * sizeof(double));
	cudaMalloc(&concentration_array_B, N * sizeof(double));

	allocateMemory<<<1,1>>>(concentration_array_A, concentration_array_B);

	initializeCoreShellParticle<<<1, N>>>();

	cudaDeviceSynchronize();

	double mass;
	cudaMemcpyFromSymbol(&mass, CoreShellParticle::core_mass, sizeof(double));
	printf("Core mass : %e kg\n", mass);
	cudaMemcpyFromSymbol(&mass, CoreShellParticle::shell_mass, sizeof(double));
	printf("Shell mass : %e kg\n", mass);
	cudaMemcpyFromSymbol(&mass, CoreShellParticle::mass, sizeof(double));
	printf("Overall mass : %e kg\n\n", mass);

	printMass<<<1,1>>>();

	double conc_array[N];

	initIteration<<<1,1>>>();

	for (size_t i = 0; i < 10000; i++)
	{
		setUpEquations<<<1, N>>>();
		solve<<<1,2>>>();
		updateMassFractions<<<1,1>>>();
		// cudaMemcpy(conc_array, concentration_array_A, N * sizeof(double), cudaMemcpyDeviceToHost);
	}

	// for(size_t i = 0; i < N; i++) printf("%e\n", conc_array[i]);

	printMass<<<1,1>>>();
	cudaDeviceSynchronize();

	deallocateMemory<<<1,1>>>();
	cudaDeviceSynchronize();

	cudaFree(concentration_array_A);
	cudaFree(concentration_array_B);

    return 0;
}