#ifndef __LU_SOLVER__
#define __LU_SOLVER__

#include "math/Data-Type.hpp"
#include "lusolver/Tridiagonal-Matrix.hpp"

class LUSolver
{
	private:

		// size of x vector for the linear algebra problem A.x = b
		const unsigned int _n;

		TridiagonalMatrix _A;

		real_t *_b;

		// Required during LU decomposition
		real_t _lower_matrix_off_diagonal_element;

		// Aliases
		real_t * &_d = _b;
		
		// Upper Diagonal matrix
		// U[i][j], j > i-1
		inline real_t& _U(unsigned int i, unsigned int j)
		{
			// if (j < i) throw std::out_of_range("Indices cannot be j < i for upper diagonal matrices.");
			return _A(i, j);
		}

		// Function to decompose the matrix A into Upper and Lower matrices
		// and simultaneously perform forward substitution
		inline void LU_DecompositionAndForwardSubstitution()
		{
			for (int i=1; i<_n; i++)
			{
				// L[i][i-1] = A[i][i-1] / U[i-1][i-1]
				_lower_matrix_off_diagonal_element = _A(i, i-1) / _U(i-1, i-1);

				// U[i][i] = A[i][i] - A[i-1][i] * L[i][i-1]
				// _U(i, i) = _A(i, i) - _A(i-1, i) * _lower_matrix_off_diagonal_element;
				_U(i, i) -= _A(i-1, i) * _lower_matrix_off_diagonal_element;
				
				// Set d[i] = b[i] - d[i-1] * L[i][i-1]
				// _d[i] = _b[i] - _d[i-1] * _lower_matrix_off_diagonal_element;
				_d[i] -= _d[i-1] * _lower_matrix_off_diagonal_element;
			}
		}

	public:

		LUSolver(unsigned int N);

		~LUSolver();

		// Set up equation represented by ith row of the matrix equation, i.e.,
		// e * x[i-1] + f * x[i] + g * x[i+1] = b
		inline void setEquation(
			unsigned int i,
			real_t e,
			real_t f,
			real_t g,
			real_t b
		) {
			_A(i, i-1) = e;
			_A(i, i)   = f;
			_A(i, i+1) = g;
			
			_b[i] = b;
		}

		// Set up equation represented by ith row of the matrix equation, i.e.,
		// e * x[i-1] + f * x[i] + g * x[i+1] = b
		// Assuming rows of the matrix are populated chronologically by row index i
		// Decomposition and forward substitution are done simultaneously
		inline void setEquationSerially(
			unsigned int i,
			real_t e,
			real_t f,
			real_t g,
			real_t b
		) {
			_A(i, i+1) = g;

			// L[i][i-1] = A[i][i-1] / U[i-1][i-1]
			_lower_matrix_off_diagonal_element = e / _U(i-1, i-1);

			// U[i][i] = A[i][i] - A[i-1][i] * L[i][i-1]
			_U(i, i) = f - _A(i-1, i) * _lower_matrix_off_diagonal_element;

			// Set d[i] = b[i] - d[i-1] * L[i][i-1]
			_d[i] = b - _d[i-1] * _lower_matrix_off_diagonal_element;
		}

		// Set up equation represented by the first row of the matrix equation
		// f * x[i] + g * x[i+1] = b
		inline void setEquationFirstRow(
			real_t f,
			real_t g,
			real_t b
		) {
			_A(0, 0) = f;
			_A(0, 1) = g;
			
			_b[0] = b;
		}

		// Set up equation represented by the last row of the matrix equation
		// e * x[n-2] + f * x[n-1] = b
		inline void setEquationLastRow(
			real_t e,
			real_t f,
			real_t b
		) {
			_A(_n-1, _n-2) = e;
			_A(_n-1, _n-1) = f;
			
			_b[_n-1] = b;
		}

		// Set up equation represented by the last row of the matrix equation
		// e * x[n-2] + f * x[n-1] = b
		// Assuming rows of the matrix are populated chronologically by row index i
		// Decomposition and forward substitution are done simultaneously
		inline void setEquationLastRowSerially(
			real_t e,
			real_t f,
			real_t b
		) {
			// L[i][i-1] = A[i][i-1] / U[i-1][i-1]
			_lower_matrix_off_diagonal_element = e / _U(_n-2, _n-2);

			// U[i][i] = A[i][i] - A[i-1][i] * L[i][i-1]
			_U(_n-1, _n-1) = f - _A(_n-2, _n-1) * _lower_matrix_off_diagonal_element;

			// Set d[i] = b[i] - d[i-1] * L[i][i-1]
			_d[_n-1] = b - _d[_n-2] * _lower_matrix_off_diagonal_element;
		}

		// Prints the matrix equation A.x = b
		void printMatrixEquation();

		// Solve the matrix equation and store it to array x
		inline void getSolution(real_t *x)
		{
			LU_DecompositionAndForwardSubstitution();

			// Backward substitution
			x[_n-1] = _d[_n-1] / _U(_n-1, _n-1);

			for (int i=_n-2; i>=0; i--)

				x[i] = ( _d[i] - _U(i, i+1) * x[i+1] ) / _U(i, i);
		}

		// Solve the matrix equation and store it to array x
		// Assuming rows of the matrix were populated chronologically by row index i
		// Decomposition and forward substitution were done simultaneously
		inline void getSolutionSerially(real_t *x)
		{
			// Backward substitution
			x[_n-1] = _d[_n-1] / _U(_n-1, _n-1);

			for (int i=_n-2; i>=0; i--)

				x[i] = ( _d[i] - _U(i, i+1) * x[i+1] ) / _U(i, i);
		}
};

#endif