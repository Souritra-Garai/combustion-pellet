/**
 * @file Lower_Triangular_Matrix.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief 
 * @version 0.1
 * @date 2021-06-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "Lower_Triangular_Matrix.hpp"

#include "stdexcept"

template<typename real_t>
LowerTriangularMatrix<real_t>::LowerTriangularMatrix(
    unsigned int n
) : // Initialize constant N_ in the mem initialization list
    N(n)
{
    // For each row i, we store only the non sero elements
    // from column 0 to column i
    // Thus allocate \f$ \Sum_{i=0}^N i = \frac{N*(N+1)}{2} \f$
    // memory to the array
    array = new real_t[(N * (N+1)) / 2];
}

template<typename real_t>
LowerTriangularMatrix<real_t>::~LowerTriangularMatrix()
{
    // Deallocate the memory allocated for array
    delete [] array;
}

template<typename real_t>
unsigned int LowerTriangularMatrix<real_t>::getIndex(
    unsigned int i,
    unsigned int j
)
{
    // For each row i, we store only the non sero elements
    // from column 0 to column i
    return ((i * (i + 1)) / 2) + j;
}

template<typename real_t>
bool LowerTriangularMatrix<real_t>::indexOfZeroElement(
    unsigned int i,
    unsigned int j
)
{
    // columns with index j greater than row index i are zero
    return j > i;
}

template<typename real_t>
real_t LowerTriangularMatrix<real_t>::getElement(
    unsigned int i,
    unsigned int j
)
{
    // Check if indices are within range of square matrix
    if (i >= N || j >= N) throw std::out_of_range("Index outside range of square matrix");

    // Check if indices are those of zero elements
    if (indexOfZeroElement(i, j)) return 0;

    else
        // Identify the index of i,j th element in the flattened array
        // and return that value
        return array[getIndex(i, j)];
}

template<typename real_t>
void LowerTriangularMatrix<real_t>::setElement(
    unsigned int i,
    unsigned int j,
    real_t val
)
{
    // Check if indices are within range of square matrix
    if (i >= N || j >= N) throw std::out_of_range("Index outside range of square matrix");
    
    // Check if indices are those of zero elements
    if (indexOfZeroElement(i, j)) 
        throw std::out_of_range(
            "Referencing to a zero element in a Lower Triangular Matrix"
        );

    // Identify the index of i,j th element in the flattened array
    // and set the value of that element to val
    array[getIndex(i, j)] = val;
}

template<typename real_t>
void LowerTriangularMatrix<real_t>::printMatrix()
{
    // Iterate through all rows
    for (int i = 0; i < N; i++)
    {
        // Iterate through all columns
        for (int j = 0; j < N; j++)
        {
            // Print i,j th element of the tridiagonal matrix
            printf("%Lf\t", (long double) getElement(i, j));
            // getElement takes care of zero elements of the matrix
        }
        // Print a newline character at teh end of each row
        printf("\n");
    }
}

template class LowerTriangularMatrix<long double>;
template class LowerTriangularMatrix<double>;
template class LowerTriangularMatrix<float>;





