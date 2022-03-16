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
		void LU_DecompositionAndForwardSubstitution();

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
        void setEquation(
            unsigned int index,
            real_t e,
            real_t f,
            real_t g,
            real_t b
        );

        /**
         * @brief Set up equation represented by the first row
         * of the matrix equation
         * \f$ f x_i + g x_{i+1} = b \f$
         * 
         * @param f Coefficient to \f$ x_{i} \f$
         * @param g Coefficient to \f$ x_{i+1} \f$
         * @param b Constant
         */
        void setEquationFirstRow(
            real_t f,
            real_t g,
            real_t b
        );

        /**
         * @brief up equation represented by the last row
         * of the matrix equation
         * \f$ e x_{i-1} f x_i = b \f$
         * 
         * @param e Coefficient to \f$ x_{i-1} \f$
         * @param f Coefficient to \f$ x_{i} \f$
         * @param b Constant
         */
        void setEquationLastRow(
            real_t e,
            real_t f,
            real_t b
        );

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
        void getSolution(real_t *x);
};

#endif