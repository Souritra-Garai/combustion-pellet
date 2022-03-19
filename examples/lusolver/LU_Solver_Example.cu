#include "lusolver/LU_Solver.hpp"

#include <iostream>
#include <time.h>

#include <curand_kernel.h>

__global__ void initializeCurandState(unsigned int n, double seed, curandState *curand_state_ptr)
{
	unsigned int i = threadIdx.x % n;

	curand_init(seed, i, 0, curand_state_ptr + i);
}

__global__ void initializeMatrix(unsigned int n, double *matrix_ptr, double *vector_ptr, curandState *curand_state_ptr)
{
	unsigned int i = threadIdx.x;

	if (i == 0)
	{
		matrix_ptr[getTridiagonalMatrixIndex(0, 0)] = curand(curand_state_ptr);
		matrix_ptr[getTridiagonalMatrixIndex(0, 1)] = curand(curand_state_ptr);

		vector_ptr[0] = curand(curand_state_ptr);
	}

	else if (i < n-1)
	{
		matrix_ptr[getTridiagonalMatrixIndex(i, i-1)] = curand(curand_state_ptr + i);
		matrix_ptr[getTridiagonalMatrixIndex(i, i)]   = curand(curand_state_ptr + i);
		matrix_ptr[getTridiagonalMatrixIndex(i, i+1)] = curand(curand_state_ptr + i);

		vector_ptr[i] = curand(curand_state_ptr + i);
	}
	
	else if (i == n-1)
	{
		matrix_ptr[getTridiagonalMatrixIndex(n-1, n-2)] = curand(curand_state_ptr + n-1);
		matrix_ptr[getTridiagonalMatrixIndex(n-1, n-1)] = curand(curand_state_ptr + n-1);

		vector_ptr[n-1] = curand(curand_state_ptr + n-1);
	}
}

__global__ void initializeSolver(unsigned int n, double *A, double *b, double *solver_A, double *solver_b)
{
	unsigned int i = threadIdx.x;

	if (i==0)

		setEquationFirstRow(solver_A, solver_b, A[getTridiagonalMatrixIndex(0,0)], A[getTridiagonalMatrixIndex(0,1)], b[0]);

	else if (i < n-1)

		setEquation(i, solver_A, solver_b, A[getTridiagonalMatrixIndex(i, i-1)], A[getTridiagonalMatrixIndex(i, i)], A[getTridiagonalMatrixIndex(i, i+1)], b[i]);

	else if (i == n-1)

		setEquationLastRow(n, solver_A, solver_b, A[getTridiagonalMatrixIndex(n-1, n-2)], A[getTridiagonalMatrixIndex(n-1, n-1)], b[n-1]);
}

__global__ void calcError(unsigned int n, double *x, double *y, float *sum)
{
	unsigned int i = threadIdx.x;
	
	float error = (float) (x[i] - y[i]) * (x[i] - y[i]);

	atomicAdd(sum, error);
}

int main(int argc, char const *argv[])
{
    const unsigned int n = 10000;

	double *tridiagonal_matrix, *x, *b;

	curandState *curand_states;

	cudaMalloc(&x, n * sizeof(double));
	cudaMalloc(&b, n * sizeof(double));
	cudaMalloc(&tridiagonal_matrix, getTridiagonalMatrixSize(n) * sizeof(double));

	cudaMalloc(&curand_states, n * sizeof(curandState));

	initializeCurandState<<<1,n>>>(n, time(0), curand_states);

	cudaDeviceSynchronize();

	initializeMatrix<<<1,n>>>(n, tridiagonal_matrix, x, curand_states);

	cudaDeviceSynchronize();

	multiplyTridiagonalMatrix<<<1,n>>>(n, tridiagonal_matrix, x, b);

    cudaDeviceSynchronize();

    LUSolver my_solver(n);

    initializeSolver<<<1,n>>>(n, tridiagonal_matrix, b, my_solver.getMatrixPtr(), my_solver.getVectorPtr());

	cudaDeviceSynchronize();

	// my_solver.printMatrixEquation(std::cout);

    double *x_soln;

	cudaMalloc(&x_soln, n * sizeof(double));

    my_solver.getSolution(x_soln);

    float *MSE, *MSE_h;

	cudaMalloc(&MSE, sizeof(float));

	MSE_h = new float;
	*MSE_h = 0;

	cudaMemcpy(MSE, MSE_h, sizeof(float), cudaMemcpyHostToDevice);

	calcError<<<1,n>>>(n, x, x_soln, MSE);

	cudaDeviceSynchronize();
    // for (int i = 0; i < N; i++) MSE += (x[i] - x_soln[i]) * (x[i] - x_soln[i]);

	cudaMemcpy(MSE_h, MSE, sizeof(float), cudaMemcpyDeviceToHost);

    std::cout << "\nMSE : " << *MSE_h << std::endl;

	cudaFree(MSE);
	delete MSE_h;

    return 0;
}