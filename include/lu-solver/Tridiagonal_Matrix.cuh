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
		private:

			size_t _n;

			double *_array;

		public:

			__host__ __device__ Matrix(size_t n)
			{
				_n = n;
				_array = new double[getSize(n)];
			}

			__host__ __device__ ~Matrix()
			{
				delete [] _array;
			}

			__host__ __device__ __forceinline__ void setElement(
				size_t row_index,
				size_t column_index,
				double value
			) {
				_array[getIndex(row_index, column_index)] = value;
			}

			__host__ __device__ __forceinline__ double getElement(
				size_t row_index,
				size_t column_index
			) {
				return _array[getIndex(row_index, column_index)];
			}

			__host__ __device__ __forceinline__ size_t getDim()
			{
				return _n;
			}

			__host__ __device__ void print();
	};

	__global__ void multiply(
		Matrix *A,
		double *x,
		double *b
	) {
		size_t i = threadIdx.x;

		if (i==0)

			b[0] = A->getElement(0, 0) * x[0] + A->getElement(0, 1) * x[1];

		else if (i < A->getDim() - 1)

			b[i] = A->getElement(i, i-1) * x[i-1] + A->getElement(i, i) * x[i] + A->getElement(i, i+1) * x[i+1];

		else if (i == A->getDim() - 1)

			b[i] = A->getElement(i, i-1) * x[i-1] + A->getElement(i, i) * x[i];
	}

	__host__ __device__ void Matrix::print()
	{
		size_t i, j;

		printf("%f\t%f", getElement(0,0), getElement(0,1));

		for (i = 2; i < getDim(); i++) printf("\t0.0");

		printf("\n");

		
		for (i = 1; i < getDim()-1; i++)
		{
			for (j = 0; j < i-1; j++) printf("0.0\t");

			printf("%f\t%f\t%f", getElement(i, i-1), getElement(i, i), getElement(i, i+1));

			for (j = i+2; j < getDim(); j++) printf("\t0.0");

			printf("\n");
		}


		for (i = 0; i < getDim()-2; i++) printf("0.0\t");

		printf("%f\t%f\n", getElement(getDim()-1, getDim()-2), getElement(getDim()-1, getDim()-1));
	}
	
} // namespace TridiagonalMatrix

#endif