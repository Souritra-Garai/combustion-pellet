#include "lusolver/LU_Solver.cuh"

#include <iostream>
#include <time.h>

#include <curand_kernel.h>

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* a, double b) { return b; }
#endif

__global__ void initializeCurandState(unsigned int n, double seed, curandState *curand_state_ptr)
{
	unsigned int i = threadIdx.x % n;

	curand_init(seed, i, 0, curand_state_ptr + i);
}

__global__ void initializeMatrix(TridiagonalMatrix::Matrix A, double *vector_ptr, curandState *curand_state_ptr)
{
	unsigned int i = threadIdx.x;

	if (i == 0)
	{
		A.setElement(0, 0, curand(curand_state_ptr));
		A.setElement(0, 1, curand(curand_state_ptr));

		vector_ptr[0] = curand(curand_state_ptr);
	}

	else if (i < A.getDim() - 1)
	{
		A.setElement(i, i-1, curand(curand_state_ptr + i));
		A.setElement(i, i,   curand(curand_state_ptr + i));
		A.setElement(i, i+1, curand(curand_state_ptr + i));

		vector_ptr[i] = curand(curand_state_ptr + i);
	}
	
	else if (i == A.getDim() - 1)
	{
		A.setElement(i, i-1, curand(curand_state_ptr + i));
		A.setElement(i, i,   curand(curand_state_ptr + i));

		vector_ptr[i] = curand(curand_state_ptr + i);
	}
}

__global__ void initializeSolver(TridiagonalMatrix::Matrix A, double *b, LUSolver solver)
{
	unsigned int i = threadIdx.x;

	if (i==0)

		solver.setEquationFirstRow(A.getElement(0, 0), A.getElement(0, 1), b[0]);

	else if (i < solver.getDim() - 1)

		solver.setEquation(i, A.getElement(i, i-1), A.getElement(i, i), A.getElement(i, i+1), b[i]);

	else if (i == solver.getDim() - 1)

		solver.setEquationLastRow(A.getElement(i, i-1), A.getElement(i, i), b[i]);
}

__global__ void solve(LUSolver solver, double *x)
{
	solver.getSolution(x);
}

__global__ void calcError(unsigned int n, double *x, double *y, double *sum)
{
	unsigned int i = threadIdx.x;
	
	double error = (x[i] - y[i]) * (x[i] - y[i]);

	// printf("%f\n", error);

	atomicAdd(sum, error);
}

int main(int argc, char const *argv[])
{
	const unsigned int n = 1024;

	double *x, *b;

	curandState *curand_states;

	cudaMalloc(&x, n * sizeof(double));
	cudaMalloc(&b, n * sizeof(double));

	cudaMalloc(&curand_states, n * sizeof(curandState));

	TridiagonalMatrix::Matrix matrix;
	matrix.allocateMemoryFromHost(n);

	initializeCurandState<<<1,n>>>(n, time(0), curand_states);

	cudaDeviceSynchronize();

	initializeMatrix<<<1,n>>>(matrix, x, curand_states);

	cudaDeviceSynchronize();

	TridiagonalMatrix::multiply<<<1,n>>>(matrix, x, b);

	cudaDeviceSynchronize();

    LUSolver my_solver;
	my_solver.allocateMemoryFromHost(n);

    initializeSolver<<<1,n>>>(matrix, b, my_solver);

	cudaDeviceSynchronize();

	// my_solver.printMatrixEquation(std::cout);

    double *x_soln;

	cudaMalloc(&x_soln, n * sizeof(double));

    solve<<<1,1>>>(my_solver, x_soln);

	cudaDeviceSynchronize();

    double *MSE, *MSE_h;

	cudaMalloc(&MSE, sizeof(double));

	MSE_h = new double;
	*MSE_h = 0;

	cudaMemcpy(MSE, MSE_h, sizeof(double), cudaMemcpyHostToDevice);

	calcError<<<1,n>>>(n, x, x_soln, MSE);

	cudaDeviceSynchronize();

	cudaMemcpy(MSE_h, MSE, sizeof(double), cudaMemcpyDeviceToHost);

    std::cout << "\nMSE : " << *MSE_h << std::endl;

	cudaFree(x_soln);
	cudaFree(MSE);
	delete MSE_h;

	cudaFree(x);
	cudaFree(b);
	cudaFree(curand_states);

	matrix.deallocateMemoryFromHost();
	my_solver.deallocateMemoryFromHost();

    return 0;
}