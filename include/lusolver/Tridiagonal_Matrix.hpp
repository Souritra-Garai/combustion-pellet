#ifndef __TRIDIAGONAL_MATRIX__
#define __TRIDIAGONAL_MATRIX__

#include <ostream>

__host__ unsigned int getTridiagonalMatrixSize(
	unsigned int square_matrix_dimension
) {
	return 3U * square_matrix_dimension - 2U;
}

__host__ __device__ inline unsigned int getTridiagonalMatrixIndex(
	unsigned int row_index,
	unsigned int column_index
) {
	return 2U * row_index + column_index;
}

__global__ void multiplyTridiagonalMatrix(
	unsigned int n,
	double *A,
	double *x,
	double *b
) {
	unsigned int i = threadIdx.x;

	if (i==0)

		b[0] = A[getTridiagonalMatrixIndex(0, 0)] * x[0] + A[getTridiagonalMatrixIndex(0, 1)] * x[1];

	else if (i < n-1)

		b[i] = A[getTridiagonalMatrixIndex(i, i-1)] * x[i-1] + A[getTridiagonalMatrixIndex(i, i)] * x[i] + A[getTridiagonalMatrixIndex(i, i+1)] * x[i+1];

	else if (i == n-1)

		b[n-1] = A[getTridiagonalMatrixIndex(n-1, n-2)] * x[n-2] + A[getTridiagonalMatrixIndex(n-1, n-1)] * x[n-1];
}

__host__ void printTridiagonalMatrix(
	std::ostream &output_stream,
	unsigned int n,
	double *matrix_ptr
) {
	unsigned int i, j;

	double matrix[getTridiagonalMatrixSize(n)];

	cudaMemcpy(matrix, matrix_ptr, getTridiagonalMatrixSize(n) * sizeof(double), cudaMemcpyDeviceToHost);

	output_stream << matrix[getTridiagonalMatrixIndex(0,0)] << '\t' << matrix[getTridiagonalMatrixIndex(0,1)];

	for (i = 2; i < n; i++) output_stream << "\t0.0";

	output_stream << '\n';

	
	for (i = 1; i < n-1; i++)
	{
		for (j = 0; j < i-1; j++) output_stream << "0.0\t";

		output_stream << matrix[getTridiagonalMatrixIndex(i, i-1)] << '\t' << matrix[getTridiagonalMatrixIndex(i, i)] << '\t' << matrix[getTridiagonalMatrixIndex(i, i+1)];

		for (j = i+2; j < n; j++) output_stream << "\t0.0";

		output_stream << '\n';
	}


	for (i = 0; i < n-2; i++) output_stream << "0.0\t";

	output_stream << matrix[getTridiagonalMatrixIndex(n-1, n-2)] << '\t' << matrix[getTridiagonalMatrixIndex(n-1, n-1)];

	output_stream << '\n';
}

#endif