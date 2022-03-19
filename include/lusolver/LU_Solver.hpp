#ifndef __LU_SOLVER__
#define __LU_SOLVER__

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <ostream>
#include "lusolver/Tridiagonal_Matrix.hpp"

__global__ void LU_DecompositionAndForwardSubstitution(unsigned int n, double *A, double *b)
{
	double Lower_Matrix_Diagonal_less_1;
	
	for (unsigned int i=1; i<n; i++)
	{
		// Set l_{i,i-1} = a_{i,i-1} / u_{i-1,i-1}
		Lower_Matrix_Diagonal_less_1 = A[getTridiagonalMatrixIndex(i, i-1)] / A[getTridiagonalMatrixIndex(i-1, i-1)];

		// Set u_{i,i} = a_{i,i} - a_{i-1,i} * l_(i,i-1)
		A[getTridiagonalMatrixIndex(i, i)] = A[getTridiagonalMatrixIndex(i, i)] - A[getTridiagonalMatrixIndex(i-1, i)] * Lower_Matrix_Diagonal_less_1;
		
		// Set d_{i} = b_{i} - d_{i-1} * l_{i,i-1}
		b[i] -= b[i-1] * Lower_Matrix_Diagonal_less_1;
	}
}

__global__ void backwardSubstitution(unsigned int n, double *A, double *b, double *x)
{
	// Backward substitution
	// Last element is simply x_{n-1} = d_{n-1} / U_{n-1, n-1}
	x[n-1] = b[n-1] / A[getTridiagonalMatrixIndex(n-1, n-1)];

	for (unsigned int i=n-2; i>0; i--)

		// Set x_{i} = ( d_{i} - U_{i,i+1} * x_{i+1} ) / U_{i,i}
		x[i] = ( b[i] - A[getTridiagonalMatrixIndex(i, i+1)] * x[i+1] ) / A[getTridiagonalMatrixIndex(i, i)];
		
	x[0] = ( b[0] - A[getTridiagonalMatrixIndex(0, 1)] * x[1] ) / A[getTridiagonalMatrixIndex(0, 0)];
}

class LUSolver
{
	private :

		const unsigned int _n;

		double *_tridiagonal_matrix;

		double *_b_vector;

	public :

		LUSolver(unsigned int n) : _n(n)
		{
			cudaMalloc(&_tridiagonal_matrix, getTridiagonalMatrixSize(n) * sizeof(double));
			cudaMalloc(&_b_vector, n * sizeof(double));
		}

		~LUSolver()
		{
			cudaFree(_tridiagonal_matrix);
			cudaFree(_b_vector);
		}

		__host__ inline double* getMatrixPtr() { return _tridiagonal_matrix; }
		__host__ inline double* getVectorPtr() { return _b_vector; }

        __host__ void printMatrixEquation(std::ostream &output_stream)
		{
			output_stream << "Matrix Eqn A.x = b\n\n";
    
			output_stream << "A\t" << _n << " x " << _n << '\n';

			printTridiagonalMatrix(output_stream, _n, _tridiagonal_matrix);

			output_stream << '\n';

			output_stream << "b\t" << _n << " x 1\n";

			double b[_n];

			cudaMemcpy(b, _b_vector, _n * sizeof(double), cudaMemcpyDeviceToHost);

			for (unsigned int i = 0; i < _n; i++) output_stream << b[i] << '\n';
		}

        __host__ void getSolution(double *x)
		{
			LU_DecompositionAndForwardSubstitution<<<1,1>>>(_n, _tridiagonal_matrix, _b_vector);

			cudaDeviceSynchronize();

			backwardSubstitution<<<1,1>>>(_n, _tridiagonal_matrix, _b_vector, x);

			cudaDeviceSynchronize();
		}
};

__device__ inline void setEquation(
	unsigned int i,
	double *tridiagonal_matrix,
	double *b_vector,
	double e,
	double f,
	double g,
	double b
) {
	// \f$ e x_{i-1} + f x_i + g x_{i+1} = b \f$

	// Set \f$ a_{i,i-1} \f$ to e_val
	tridiagonal_matrix[getTridiagonalMatrixIndex(i, i-1)] = e;
	// Set \f$ a_{i,i} \f$ to f_val
	tridiagonal_matrix[getTridiagonalMatrixIndex(i, i)]   = f;
	// Set \f$ a_{i,i+1} \f$ to g_val
	tridiagonal_matrix[getTridiagonalMatrixIndex(i, i+1)] = g;

	// Set the constant
	b_vector[i] = b;
}

__device__ void setEquationFirstRow(
	double *tridiagonal_matrix,
	double *b_vector,
	double f,
	double g,
	double b
) {
	// \f$ f x_i + g x_{i+1} = b \f$

	// Set \f$ a_{i,i} \f$ to f_val
	tridiagonal_matrix[getTridiagonalMatrixIndex(0, 0)] = f;
	// Set \f$ a_{i,i+1} \f$ to g_val
	tridiagonal_matrix[getTridiagonalMatrixIndex(0, 1)] = g;

	// Set the constant
	b_vector[0] = b;
}

__device__ void setEquationLastRow(
	unsigned int n,
	double *tridiagonal_matrix,
	double *b_vector,
	double e,
	double f,
	double b
) {
	// \f$ e x_{i-1} f x_i = b \f$

	// Set \f$ a_{i,i-1} \f$ to e_val
	tridiagonal_matrix[getTridiagonalMatrixIndex(n-1, n-2)] = e;
	// Set \f$ a_{i,i} \f$ to f_val
	tridiagonal_matrix[getTridiagonalMatrixIndex(n-1, n-1)] = f;

	// Set the constant
	b_vector[n-1] = b;
}

#endif