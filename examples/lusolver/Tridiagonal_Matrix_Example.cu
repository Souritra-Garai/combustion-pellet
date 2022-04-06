#include "lusolver/Tridiagonal_Matrix.cuh"

#include <iostream>
#include <time.h>

#include <curand_kernel.h>

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

int main(int argc, char const *argv[])
{
	const unsigned int n = 5;

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

	matrix.print(std::cout);

	cudaFree(x);
	cudaFree(b);
	cudaFree(curand_states);

	matrix.deallocateMemoryFromHost();

	return 0;
}
