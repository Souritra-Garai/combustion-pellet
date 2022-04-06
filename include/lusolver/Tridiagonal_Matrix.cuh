#ifndef __TRIDIAGONAL_MATRIX__
#define __TRIDIAGONAL_MATRIX__

#include <ostream>

namespace TridiagonalMatrix
{
	__host__ __device__ __forceinline__ size_t getIndex(size_t i, size_t j)
	{
		return 2U * i + j;
	}

	__host__ __device__ __forceinline__ size_t getSize(size_t n)
	{
		return 3U * n - 2U;
	}

	class Matrix
	{
		public:

			__host__ void allocateMemoryFromHost(size_t n)
			{
				_n = n;
				cudaMalloc(&_array, getSize(n) * sizeof(double));
			}

			__device__ void allocateMemoryFromDevice(size_t n)
			{
				_n = n;
				_array = new double[getSize(n)];
			}

			__host__ void deallocateMemoryFromHost()
			{
				cudaFree(_array);
				_n = 0;
			}

			__device__ void deallocateMemoryFromDevice()
			{
				delete [] _array;
				_n = 0;
			}

			__device__ __forceinline__ void setElement(
				size_t row_index,
				size_t column_index,
				double value
			) {
				_array[getIndex(row_index, column_index)] = value;
			}

			__device__ __forceinline__ double getElement(
				size_t row_index,
				size_t column_index
			) {
				return _array[getIndex(row_index, column_index)];
			}

			__host__ __device__ __forceinline__ size_t getDim()
			{
				return _n;
			}

			__host__ void print(std::ostream &output_stream);

		private:

			size_t _n;

			double *_array;
	};

	__global__ void multiply(
		Matrix A,
		double *x,
		double *b
	) {
		size_t i = threadIdx.x;

		if (i==0)

			b[0] = A.getElement(0, 0) * x[0] + A.getElement(0, 1) * x[1];

		else if (i < A.getDim() - 1)

			b[i] = A.getElement(i, i-1) * x[i-1] + A.getElement(i, i) * x[i] + A.getElement(i, i+1) * x[i+1];

		else if (i == A.getDim() - 1)

			b[i] = A.getElement(i, i-1) * x[i-1] + A.getElement(i, i) * x[i];
	}

	__host__ void Matrix::print(std::ostream &output_stream)
	{
		size_t i, j;

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

} // namespace TridiagonalMatrix

#endif