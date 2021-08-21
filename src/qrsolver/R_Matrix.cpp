/**
 * @file    R_Matrix.cpp
 * @author  Souritra Gari (souritra.garai@iitgn.ac.in)
 * @brief   Member function definitions for RMatrix class
 * @version 0.1
 * @date    2021-06-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "qrsolver/R_Matrix.hpp"

// Required for out of range exception
#include <stdexcept>

template<typename real_t>
RMatrix<real_t>::RMatrix(
    unsigned int n
) : // Initialize constant N in the mem initialization list
    N(n)
{
    // Allocate memory for the flattened array implemented for
    // representing 4 diagonals of the matrix
    array = new real_t[4*N];
}

template<typename real_t>
RMatrix<real_t>::~RMatrix()
{
    // Deallocate the array
    delete [] array;
}

template<typename real_t>
real_t RMatrix<real_t>::getElement(
    unsigned int i,
    unsigned int j
)
{
    // Check if indices are within range of square matrix
    if ((i >= N) || (j >= N)) throw std::out_of_range("Index outside range of square matrix");

    // Check if indices are those of zero elements
    if (indexOfZeroElement(i, j)) return 0;

    else
        // Identify the index of i,j th element in the flattened array
        // and return that value
        return array[getIndex(i, j)];
}

template<typename real_t>
void RMatrix<real_t>::setElement(
    unsigned int i,
    unsigned int j,
    real_t val
)
{
    // Check if indices are within range of square matrix
    if ((i >= N) || (j >= N)) throw std::out_of_range("Index outside range of square matrix");
    
    // Check if indices are those of zero elements
    if (indexOfZeroElement(i, j)) 
        throw std::out_of_range(
            "Referencing to a zero element in a R Matrix"
        );

    // Identify the index of i,j th element in the flattened array
    // and set the value of that element to val
    array[getIndex(i, j)] = val;
}

template<typename real_t>
void RMatrix<real_t>::printMatrix()
{
    // Iterate through all rows
    for (int i = 0; i < N; i++)
    {
        // Iterate through all columns
        for (int j = 0; j < N; j++)
        {
            // Print i,j th element of the matrix
            printf("%LE\t", (long double) getElement(i, j));
            // getElement takes care of zero elements of the matrix
        }
        // Print a newline character at teh end of each row
        printf("\n");
    }
}


template<typename real_t>
void RMatrix<real_t>::print()
{
    // Iterate through each of the four diagonals
    for (int k = 0; k < 4; k++)
    {
        // Iterate through the diagonal
        for (int i = 0; i < N; i++) 
            // Print each element of diagonal
            printf("%Lf\t", (long double) array[i*4 + k]);
        // Print a newline character at the end of each diagonal
        printf("\n");
    }
}

template<typename real_t>
unsigned int RMatrix<real_t>::getIndex(
    unsigned int i,
    unsigned int j
)
{ 
    // In the current implementation of 2D R Matrix as 
    // a flattened array, the elements of i th row of the 2D R matrix
    // is present in the 4 consecutive elements starting at index # 4*i
    // 4*i is the element in the sub diagonal in the i th row \f$\Rightarrow (i,i-1)\f$
    // 4*i + 1 is the element in the main diagonal in the i th row \f$\Rightarrow (i,i)\f$
    // 4*i + 2 is the element in the super diagonal in the i th row \f$\Rightarrow (i,i+1)\f$
    // 4*i + 3 is the element in the super super diagonal in the i th row \f$\Rightarrow (i,i+2)\f$
    return 3*i + (j + 1);
}

template<typename real_t>
bool RMatrix<real_t>::indexOfZeroElement(
    unsigned int i,
    unsigned int j
)
{
    // In our R matrix, only the elements at position
    // (i,i-1), (i,i), (i,i+1) and (i,i+2) can be non zero
    return j + 1 < i || j > i + 2;
}

template class RMatrix<long double>;
template class RMatrix<double>;
template class RMatrix<float>;