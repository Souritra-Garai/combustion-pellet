/**
 * @file Q_Matrix.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Member function definitions for QMatrix class
 * @version 0.1
 * @date 2021-06-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "qrsolver/Q_Matrix.hpp"

// Required for out of range exception 
#include <stdexcept>

// Required for memset function
#include <bits/stdc++.h>

template<typename real_t>
QMatrix<real_t>::QMatrix(
    unsigned int n
) : // Initialize constant N in the mem initialization list
    N(n)
{
    // For each row i, we store only the non zero elements
    // from column 1 to column i+1
    // Thus allocate \f$ \Sum_{i=1}^N i+1 = \frac{N(N+1)}{2} + N = \frac{N(N+3)}{2} \f$
    // memory to the array
    array = new real_t[(N * (N+3)) / 2];
}

template<typename real_t>
QMatrix<real_t>::~QMatrix()
{
    // Deallocate the memory allocated for array
    delete [] array;
}

template<typename real_t>
unsigned int QMatrix<real_t>::getIndex(
    unsigned int i,
    unsigned int j
)
{
    // For each row i, we store only the non zero elements
    // from column 0 to column i
    return ((i * (i + 3)) / 2) + j;
}

template<typename real_t>
bool QMatrix<real_t>::indexOfZeroElement(
    unsigned int i,
    unsigned int j
)
{
    // columns with index j greater than row index i+1 are zero
    return j > i + 1;
}

template<typename real_t>
real_t QMatrix<real_t>::getElement(
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
void QMatrix<real_t>::setElement(
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
            "Referencing to a zero element in a Q Matrix"
        );

    // Identify the index of i,j th element in the flattened array
    // and set the value of that element to val
    array[getIndex(i, j)] = val;
}

template<typename real_t>
void QMatrix<real_t>::printMatrix()
{
    // Iterate through all rows
    for (int i = 0; i < N; i++)
    {
        // Iterate through all columns
        for (int j = 0; j < N; j++)
        {
            // Print i,j th element of the Q matrix
            printf("%Lf\t", (long double) getElement(i, j));
            // getElement takes care of zero elements of the matrix
        }
        // Print a newline character at teh end of each row
        printf("\n");
    }
}

template<typename real_t>
void QMatrix<real_t>::identity()
{
    // Set the N * sizeof(real_t) bytes to zero
    // starting at array
    memset(array, 0, ((N * (N+3)) / 2) * sizeof(real_t));

    // Parallelize the loop
    #pragma omp parallel for
        // Set main diagonal to 1
        for (int i = 0; i < N; i++) array[getIndex(i, i)] = 1;
}

template<typename real_t>
void QMatrix<real_t>::multiply(
    real_t *b,
    real_t *x
) {
    // The calculation for x[i] is independent of other rows
    // Thus parallelization is possitble
    #pragma omp parallel for

        // For each row of the Q matrix except the last row
        for (int i = 0; i < N-1; i++)
        {
            // The calculation performed at each column of Q
            // is independent of other columns
            // But x[i] is the cumulative sum of the multiplications
            // across all the columns
            #pragma omp parallel for reduction(+:x[i])

                // For each non zero element of Q in the row i
                for (int j = 0; j <= i+1; j++)
                
                    // Cumulative sum of products of the elements of Q and b
                    x[i] += getElement(i, j) * b[j];
        }

    // Again, calculations performed at each column of Q
    // is independent of other columns
    // But x[N-1] is the cumulative sum of the multiplications
    // across all the columns
    #pragma omp parallel for reduction(+:x[N-1])

        // For the last row the iterations across
        // columns go till N-1
        for (int j = 0; j < N; j++)
        
            // Cumulative sum of products of the elements of Q and b
            x[N-1] += getElement(N-1, j) * b[j];
}

template class QMatrix<long double>;
template class QMatrix<double>;
template class QMatrix<float>;
