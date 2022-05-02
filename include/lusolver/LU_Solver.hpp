#ifndef __LU_SOLVER__
#define __LU_SOLVER__

#include "lusolver/Tridiagonal_Matrix.hpp"

template<typename real_t>
class LUSolver
{
	private:

		// size of x vector for the linear algebra problem A.x = b
		const unsigned int _n;

		TridiagonalMatrix<real_t> _A;

		// pointer to array of real type numbers representing
		// constant vector b in the matrix equation
		real_t *b;

		// helper function to decompose the A matrix in Upper and Lower matrices
		inline void LU_DecompositionAndForwardSubstitution()
		{
			real_t Lower_Matrix_Diagonal_less_1;
			
			for (int i=1; i<_n; i++)
			{
				// Set l_{i,i-1} = a_{i,i-1} / u_{i-1,i-1}
				Lower_Matrix_Diagonal_less_1 = _A.getElement(i, i-1) / _A.getElement(i-1, i-1);

				// Set u_{i,i} = a_{i,i} - a_{i-1,i} * l_(i,i-1)
				_A.setElement(i, i,
					_A.getElement(i, i) - _A.getElement(i-1, i) * Lower_Matrix_Diagonal_less_1
				);
				
				// Set d_{i} = b_{i} - d_{i-1} * l_{i,i-1}
				b[i] -= b[i-1] * Lower_Matrix_Diagonal_less_1;
			}
		}

	public:

		LUSolver(unsigned int N);

		~LUSolver();

        /**
         * @brief Set up equation represented by ith
         * row of the matrix equation
         * \f$ e x_{i-1} + f x_i + g x_{i+1} = b \f$
         * 
         * @param index i
         * @param e Coefficient to \f$ x_{i-1} \f$
         * @param f Coefficient to \f$ x_{i} \f$
         * @param g Coefficient to \f$ x_{i+1} \f$
         * @param b Constant
         */
        inline void setEquation(
            unsigned int i,
            real_t e_val,
            real_t f_val,
            real_t g_val,
            real_t b_val
        ) {
			// \f$ e x_{i-1} + f x_i + g x_{i+1} = b \f$

			// Set \f$ a_{i,i-1} \f$ to e_val
			_A.setElement(i, i-1, e_val);
			// Set \f$ a_{i,i} \f$ to f_val
			_A.setElement(i, i,   f_val);
			// Set \f$ a_{i,i+1} \f$ to g_val
			_A.setElement(i, i+1, g_val);

			// Set the constant
			b[i] = b_val;
		}

        /**
         * @brief Set up equation represented by the first row
         * of the matrix equation
         * \f$ f x_i + g x_{i+1} = b \f$
         * 
         * @param f Coefficient to \f$ x_{i} \f$
         * @param g Coefficient to \f$ x_{i+1} \f$
         * @param b Constant
         */
        inline void setEquationFirstRow(
            real_t f_val,
            real_t g_val,
            real_t b_val
        ) {
			// \f$ f x_i + g x_{i+1} = b \f$

			// Set \f$ a_{i,i} \f$ to f_val
			_A.setElement(0, 0, f_val);
			// Set \f$ a_{i,i+1} \f$ to g_val
			_A.setElement(0, 1, g_val);

			// Set the constant
			b[0] = b_val;
		}

        /**
         * @brief up equation represented by the last row
         * of the matrix equation
         * \f$ e x_{i-1} f x_i = b \f$
         * 
         * @param e Coefficient to \f$ x_{i-1} \f$
         * @param f Coefficient to \f$ x_{i} \f$
         * @param b Constant
         */
        inline void setEquationLastRow(
            real_t e_val,
            real_t f_val,
            real_t b_val
        ) {
			// \f$ e x_{i-1} f x_i = b \f$

			// Set \f$ a_{i,i-1} \f$ to e_val
			_A.setElement(_n-1, _n-2, e_val);
			// Set \f$ a_{i,i} \f$ to f_val
			_A.setElement(_n-1, _n-1, f_val);

			// Set the constant
			b[_n-1] = b_val;
		}

        /**
         * @brief Prints the matrix \f$ A \f$ and vector \f$ b \f$
         */
        void printMatrixEquation();

        /**
         * @brief Finds the solution to matrix equation and saves
         * it to array x
         * 
         * @param x Array to store the solution of the matrix equation
         */
        inline void getSolution(real_t *x)
		{
			LU_DecompositionAndForwardSubstitution();

			// Backward substitution
			// Last element is simply x_{n-1} = d_{n-1} / U_{n-1, n-1}
			x[_n-1] = b[_n-1] / _A.getElement(_n-1, _n-1);

			for (int i=_n-2; i>=0; i--)

				// Set x_{i} = ( d_{i} - U_{i,i+1} * x_{i+1} ) / U_{i,i}
				x[i] = ( b[i] - _A.getElement(i, i+1) * x[i+1] ) / _A.getElement(i, i);
		}
};

#endif