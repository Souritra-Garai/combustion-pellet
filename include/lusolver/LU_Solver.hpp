#ifndef __LU_SOLVER__
#define __LU_SOLVER__

template<typename real_t>
class LUSolver
{
	private:

		// size of x vector for the linear algebra problem A.x = b
		const unsigned int N;

		real_t *A_matrix_diagonal;
		real_t *A_matrix_diagonal_plus_1;
		real_t *A_matrix_diagonal_less_1;

		// vectors / arrays carrying values of diagonals 
		// of decomposed upper and lower matrices
		
		// // Upper matrix
		// // Main Diagonal elements - total size n
		// real_t *Upper_Matrix_Diagonal;
		// // Diagonal immediately to right of main diagonal
		// // total size n; Upper_Matrix_Diagonal_plus_1[n-1] is ignored
		// real_t *Upper_Matrix_Diagonal_plus_1;

		real_t *b;
		// real_t *d;

		// helper function to decompose the A matrix in Upper and Lower matrices
		void LU_Decomposition();

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