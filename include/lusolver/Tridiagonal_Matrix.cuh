#ifndef __TRIDIAGONAL_MATRIX__
#define __TRIDIAGONAL_MATRIX__

#include <ostream>

class TridiagonalMatrix
{
	public:

		__host__ void allocateMemory(unsigned int n)
		{
			_n = n;

			cudaMalloc(&_array, getSize(n) * sizeof(double));
		}

		__host__ void deallocateMemory()
		{
			cudaFree(_array);

			_n = 0;
		}

		__device__ __forceinline__ void setElement(
			unsigned int row_index,
			unsigned int column_index,
			double value
		) {
			_array[getIndex(row_index, column_index)] = value;
		}

		__device__ __forceinline__ double getElement(
			unsigned int row_index,
			unsigned int column_index
		) {
			return _array[getIndex(row_index, column_index)];
		}

		__host__ __device__ __forceinline__ unsigned int getDim()
		{
			return _n;
		}

		__host__ void printMatrix(std::ostream &output_stream);

	private:

		unsigned int _n;

		double *_array;

		__host__ __device__ __forceinline__ unsigned int getIndex(unsigned int i, unsigned int j)
		{
			return 2U * i + j;
		}

		__host__ __device__ __forceinline__ unsigned int getSize(unsigned int n)
		{
			return 3U * n - 2U;
		}
};

__global__ void multiplyTridiagonalMatrix(
	TridiagonalMatrix A,
	double *x,
	double *b
) {
	unsigned int i = threadIdx.x;

	if (i==0)

		b[0] = A.getElement(0, 0) * x[0] + A.getElement(0, 1) * x[1];

	else if (i < A.getDim() - 1)

		b[i] = A.getElement(i, i-1) * x[i-1] + A.getElement(i, i) * x[i] + A.getElement(i, i+1) * x[i+1];

	else if (i == A.getDim() - 1)

		b[i] = A.getElement(i, i-1) * x[i-1] + A.getElement(i, i) * x[i];
}

__host__ void TridiagonalMatrix::printMatrix(std::ostream &output_stream) 
{
	unsigned int i, j;

	double matrix[getSize(_n)];

	cudaMemcpy(matrix, _array, getSize(_n) * sizeof(double), cudaMemcpyDeviceToHost);

	output_stream << matrix[getIndex(0,0)] << '\t' << matrix[getIndex(0,1)];

	for (i = 2; i < _n; i++) output_stream << "\t0.0";

	output_stream << '\n';

	
	for (i = 1; i < _n-1; i++)
	{
		for (j = 0; j < i-1; j++) output_stream << "0.0\t";

		output_stream << matrix[getIndex(i, i-1)] << '\t' << matrix[getIndex(i, i)] << '\t' << matrix[getIndex(i, i+1)];

		for (j = i+2; j < _n; j++) output_stream << "\t0.0";

		output_stream << '\n';
	}


	for (i = 0; i < _n-2; i++) output_stream << "0.0\t";

	output_stream << matrix[getIndex(_n-1, _n-2)] << '\t' << matrix[getIndex(_n-1, _n-1)];

	output_stream << '\n';
}

#endif