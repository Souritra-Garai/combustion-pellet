#include "lusolver/LU_Solver.hpp"

// Required for the printf function
#include <stdio.h>

template<typename real_t>
LUSolver<real_t>::LUSolver(unsigned int n) : N(n)
{
	A_matrix_diagonal = new real_t[N];
	A_matrix_diagonal_plus_1 = new real_t[N];
	A_matrix_diagonal_less_1 = new real_t[N];

	// Upper_Matrix_Diagonal = new real_t[N];
	// Upper_Matrix_Diagonal_plus_1 = new real_t[N];

	b = new real_t[N];
	// d = new real_t[N];
}

template<typename real_t>
LUSolver<real_t>::~LUSolver()
{
	delete [] A_matrix_diagonal;
	delete [] A_matrix_diagonal_less_1;
	delete [] A_matrix_diagonal_plus_1;

	// delete [] Upper_Matrix_Diagonal;
	// delete [] Upper_Matrix_Diagonal_plus_1;

	delete [] b;
	// delete [] d;
}

template<typename real_t>
void LUSolver<real_t>::setEquation(
	unsigned int i,
	real_t e_val,
    real_t f_val,
    real_t g_val,
    real_t b_val
) {
	// \f$ e x_{i-1} + f x_i + g x_{i+1} = b \f$

	// Set \f$ a_{i,i-1} \f$ to e_val
    A_matrix_diagonal_less_1[i] = e_val;
    // Set \f$ a_{i,i} \f$ to f_val
    A_matrix_diagonal[i] = f_val;
    // Set \f$ a_{i,i+1} \f$ to g_val
    A_matrix_diagonal_plus_1[i] = g_val;

    // Set the constant
    b[i] = b_val;
}

template<typename real_t>
void LUSolver<real_t>::setEquationFirstRow(
    real_t f_val,
    real_t g_val,
    real_t b_val
) {
    // \f$ f x_i + g x_{i+1} = b \f$

    // Set \f$ a_{i,i} \f$ to f_val
    A_matrix_diagonal[0] = f_val;
    // Set \f$ a_{i,i+1} \f$ to g_val
    A_matrix_diagonal_plus_1[0] = g_val;

    // Set the constant
    b[0] = b_val;
}

template<typename real_t>
void LUSolver<real_t>::setEquationLastRow(
    real_t e_val,
    real_t f_val,
    real_t b_val
) {
    // \f$ e x_{i-1} f x_i = b \f$

    // Set \f$ a_{i,i-1} \f$ to e_val
    A_matrix_diagonal_less_1[N-1] = e_val;
    // Set \f$ a_{i,i} \f$ to f_val
    A_matrix_diagonal[N-1] = f_val;

    // Set the constant
    b[N-1] = b_val;
}

template<typename real_t>
void LUSolver<real_t>::printMatrixEquation()
{
    printf("Matrix Eqn A.x = b\n");
    printf("\nA_{%d \\times %d}\n", N, N);
    // Print \f$ A \f$ matrix
    
	printf("%LE\t%LE\t",
		(long double) A_matrix_diagonal[0],
		(long double) A_matrix_diagonal_plus_1[0]);
	for (int i = 2; i < N; i++) printf("0\t");
	printf("\n");
	
	for (int i = 1; i < N-1; i++)
	{
		for (int j = 0; j < i-1; j++) printf("0\t");

		printf("%LE\t%LE\t%LE\t",
			(long double) A_matrix_diagonal_less_1[i],
			(long double) A_matrix_diagonal[i],
			(long double) A_matrix_diagonal_plus_1[i]
		);

		for (int j = i+2; j < N; j++) printf("0\t");

		printf("\n");
	}

	for (int i = 0; i < N-2; i++) printf("0\t");
	printf("%LE\t%LE\t",
		(long double) A_matrix_diagonal_less_1[N-1],
		(long double) A_matrix_diagonal[N-1]);
	printf("\n");

    printf("\nb_{%d \\times 1}\n", N);
    // Print \f$ b \f$ vector
    for (int k = 0; k < N; k++) printf("%LE\n", (long double) b[k]);
}

template<typename real_t>
void LUSolver<real_t>::LU_Decomposition()
{
	// Upper_Matrix_Diagonal[0] = A_matrix_diagonal[0];
    // Upper_Matrix_Diagonal_plus_1[0] = A_matrix_diagonal_plus_1[0];

	// d[0] = b[0];

	real_t Lower_Matrix_Diagonal_less_1;
    
    for (int i=1; i<N; i++)
    {
        Lower_Matrix_Diagonal_less_1 = A_matrix_diagonal_less_1[i] / A_matrix_diagonal[i-1];

        A_matrix_diagonal[i] -= A_matrix_diagonal_plus_1[i-1] * Lower_Matrix_Diagonal_less_1;

        // Upper_Matrix_Diagonal_plus_1[i] = A_matrix_diagonal_plus_1[i];

		b[i] -= b[i-1] * Lower_Matrix_Diagonal_less_1;
    }
}

template<typename real_t>
void LUSolver<real_t>::getSolution(real_t *x)
{
	LU_Decomposition();

	x[N-1] = b[N-1] / A_matrix_diagonal[N-1];

    for (int i=N-2; i>=0; i--)

        x[i] = ( b[i] - A_matrix_diagonal_plus_1[i]*x[i+1] ) / A_matrix_diagonal[i];
}

template class LUSolver<long double>;
template class LUSolver<double>;
template class LUSolver<float>;