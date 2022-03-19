#include "lusolver/Tridiagonal_Matrix.hpp"

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

int main(int argc, char const *argv[])
{
	const unsigned int n = 5;

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

	printTridiagonalMatrix(std::cout, n, tridiagonal_matrix);

	cudaFree(x);
	cudaFree(b);
	cudaFree(tridiagonal_matrix);
	cudaFree(curand_states);

	return 0;
}
