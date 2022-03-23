#ifndef __LU_SOLVER__
#define __LU_SOLVER__

#include <ostream>
#include "lusolver/Tridiagonal_Matrix.cuh"

class LUSolver
{
	public :

		__host__ void allocateMemory(unsigned int n)
		{
			_A_matrix.allocateMemory(n);

			cudaMalloc(&_b_vector, n * sizeof(double));
		}

		__host__ void deallocateMemory()
		{
			_A_matrix.deallocateMemory();

			cudaFree(_b_vector);
		}

		__device__ __forceinline__ void setEquation(
			unsigned int i,
			double e,
			double f,
			double g,
			double b
		) {
			// e x_{i-1} + f x_i + g x_{i+1} = b

			// Set a_{i,i-1) to e
			_A_matrix.setElement(i, i-1, e);
			// Set a_{i,i) to f
			_A_matrix.setElement(i, i, f);
			// Set a_{i,i+1} to g
			_A_matrix.setElement(i, i+1, g);

			// Set the constant
			_b_vector[i] = b;
		}

		__device__ __forceinline__ void setEquationFirstRow(
			double f,
			double g,
			double b
		) {
			// f x_0 + g x_1 = b

			// Set a_{0,0} to f
			_A_matrix.setElement(0, 0, f); 
			// Set a_{0,1} to g
			_A_matrix.setElement(0, 1, g);

			// Set the constant
			_b_vector[0] = b;
		}
		
		__device__ __forceinline__ void setEquationLastRow(
			double e,
			double f,
			double b
		) {
			// e x_{n-1} + f x_n = b

			// Set a_{n,n-1} to e
			_A_matrix.setElement(getDim() - 1, getDim() - 2, e);
			// Set a_{n,n} to f
			_A_matrix.setElement(getDim() - 1, getDim() - 1, f);

			// Set the constant
			_b_vector[getDim() - 1] = b;
		}

        __device__ void getSolution(double *x)
		{
			LU_DecompositionAndForwardSubstitution();

			backwardSubstitution(x);
		}

		__host__ __device__ __forceinline__ unsigned int getDim()
		{
			return _A_matrix.getDim();
		}

        __host__ void printMatrixEquation(std::ostream &output_stream);

	private :

		unsigned int _n;

		TridiagonalMatrix _A_matrix;

		double *_b_vector;

		__device__ void LU_DecompositionAndForwardSubstitution()
		{
			double lower_matrix_diagonal_less_1;

			for (unsigned int i = 1U; i < _A_matrix.getDim(); i++)
			{
				// Set l_{i,i-1} = a_{i,i-1} / u_{i-1,i-1}
				lower_matrix_diagonal_less_1 = _A_matrix.getElement(i, i-1) / _A_matrix.getElement(i-1, i-1);

				// Set u_{i,i} = a_{i,i} - a_{i-1,i} * l_(i,i-1)
				_A_matrix.setElement(i, i,
					_A_matrix.getElement(i, i) - _A_matrix.getElement(i-1, i) * lower_matrix_diagonal_less_1
				);
				
				// Set d_{i} = b_{i} - d_{i-1} * l_{i,i-1}
				_b_vector[i] -= _b_vector[i-1] * lower_matrix_diagonal_less_1;
			}
		}

		__device__ void backwardSubstitution(double *x_vector)
		{
			// Last element is simply x_{n-1} = d_{n-1} / U_{n-1, n-1}
			x_vector[getDim() - 1] = _b_vector[getDim() - 1] / _A_matrix.getElement(getDim() - 1, getDim() - 1);

			for (unsigned int i = getDim() - 2; i > 0; i--)

				// Set x_{i} = ( d_{i} - U_{i,i+1} * x_{i+1} ) / U_{i,i}
				x_vector[i] = ( _b_vector[i] - _A_matrix.getElement(i, i+1) * x_vector[i+1] ) / _A_matrix.getElement(i, i);
				
			x_vector[0] = ( _b_vector[0] - _A_matrix.getElement(0, 1) * x_vector[1] ) / _A_matrix.getElement(0, 0);
		}
};

__host__ void LUSolver::printMatrixEquation(std::ostream &output_stream)
{
	output_stream << "Matrix Eqn A.x = b\n\n";
    
	output_stream << "A\t" << _n << " x " << _n << '\n';

	_A_matrix.printMatrix(output_stream);

	output_stream << '\n';

	output_stream << "b\t" << _n << " x 1\n";

	double b[_n];

	cudaMemcpy(b, _b_vector, _n * sizeof(double), cudaMemcpyDeviceToHost);

	for (unsigned int i = 0; i < _n; i++) output_stream << b[i] << '\n';
}

#endif