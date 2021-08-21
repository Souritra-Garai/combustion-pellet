/**
 * @file QR_Solver.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Member function definitions for QRSolver class
 * @version 0.1
 * @date 2021-07-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "qrsolver/QR_Solver.hpp"

// Required for the printf function
#include <stdio.h>
// Required for the memset function
#include <bits/stdc++.h>

template<typename real_t>
QRSolver<real_t>::QRSolver(
    unsigned int n
) : // Initialise the QMatrix, RMatrix and N in the mem initialisation list
    Q(n),
    R(n),
    N(n)
{
    // Allocate memory for the array represent equation constants
    b = new real_t[N];
}

template<typename real_t>
QRSolver<real_t>::~QRSolver()
{
    // Deallocate the memory allocated for b
    delete [] b;
}

template<typename real_t>
void QRSolver<real_t>::setEquation(
    unsigned int i,
    real_t e_val,
    real_t f_val,
    real_t g_val,
    real_t b_val
) {
    // \f$ e x_{i-1} + f x_i + g x_{i+1} = b \f$

    // Set \f$ a_{i,i-1} \f$ to e_val
    R.setElement(i, i-1, e_val);
    // Set \f$ a_{i,i} \f$ to f_val
    R.setElement(i, i,   f_val);
    // Set \f$ a_{i,i+1} \f$ to g_val
    R.setElement(i, i+1, g_val);

    // Set the constant
    b[i] = b_val;
}

template<typename real_t>
void QRSolver<real_t>::setEquationFirstRow(
    real_t f_val,
    real_t g_val,
    real_t b_val
) {
    // \f$ f x_i + g x_{i+1} = b \f$

    // Set \f$ a_{i,i} \f$ to f_val
    R.setElement(0, 0, f_val);
    // Set \f$ a_{i,i+1} \f$ to g_val
    R.setElement(0, 1, g_val);

    // Set the constant
    b[0] = b_val;
}

template<typename real_t>
void QRSolver<real_t>::setEquationLastRow(
    real_t e_val,
    real_t f_val,
    real_t b_val
) {
    // \f$ e x_{i-1} f x_i = b \f$

    // Set \f$ a_{i,i-1} \f$ to e_val
    R.setElement(N-1, N-2, e_val);
    // Set \f$ a_{i,i} \f$ to f_val
    R.setElement(N-1, N-1, f_val);

    // Set the constant
    b[N-1] = b_val;
}

template<typename real_t>
void QRSolver<real_t>::printMatrixEquation()
{
    printf("Matrix Eqn A.x = b\n");
    printf("\nA_{%d \\times %d}\n", N, N);
    // Print \f$ A \f$ matrix
    R.printMatrix();

    printf("\nb_{%d \\times 1}\n", N);
    // Print \f$ b \f$ vector
    for (k = 0; k < N; k++) printf("%LE\n", (long double) b[k]);
}

template<typename real_t>
void QRSolver<real_t>::QRFactorize()
{
    // Initialize the Q matrix to an identity matrix
    Q.identity();

    // Iteratively vanish \f$ r_{k+1,k} \f$
    // Needs to be performed serially
    for (k = 0; k < N-2; k++)
    {
        // Create a Givens Rotation Matrix to vanish
        // \f$ r_{k+1, k} \f$
        GMatrix<real_t> G(R, k);

        // Multiply Givens Rotation Matrix with \f$ R \f$
        G.multiply(R);
        // Multiply Givens Rotation Matrix with \f$ Q \f$
        G.multiply(Q);
    }

    // Last iteration needs to be handled differently as
    // N+1 column doesnot exist

    // Create a Givens Rotation Matrix to vanish
    // \f$ r_{k+1,k} = r_{N-1,N-2} \f$
    GMatrix<real_t> G(R, k);
    
    // Multiply Givens Rotation Matrix to \f$ R \f$
    G.multiplyLastRow(R);
    // Multiply Givens Rotation Matrix to \f$ Q \f$
    G.multiply(Q);
}

template<typename real_t>
void QRSolver<real_t>::printQRMatrices()
{
    printf("A = Q.R\n");

    printf("\nQ^T_{%d \\times %d}\n", N, N);
    // Print \f$ Q \f$ matrix
    Q.printMatrix();

    printf("\nR_{%d \\times %d}\n", N, N);
    // Print \f$ R \f$ matrix
    R.printMatrix();
}

template<typename real_t>
void QRSolver<real_t>::getSolution(
    real_t *x
) {
    // Initialise x to a zero vector
    memset(x, 0, N * sizeof(real_t));
    
    // Factorize \f$ A \f$ matrix to \f$ Q \cdot R \f$
    QRFactorize();

    // Multiply \f$ Q \cdot b \f$ and store the result in \f$ x \f$
    Q.multiply(b, x);

    // Solve the linear equations row wise starting from the bottom
    // serially

    // Last row
    x[N-1] = x[N-1] / R.getElement(N-1, N-1);

    // Second last row
    x[N-2] = (x[N-2] - R.getElement(N-2, N-1) * x[N-1]) / R.getElement(N-2, N-2);

    // For the next N-2 rows, solving the equations is similar
    for (k = N-3; k > 0; k--)
    {
        x[k] = (x[k] - R.getElement(k, k+1) * x[k+1] - R.getElement(k, k+2) * x[k+2]) / R.getElement(k, k);
        R.setElement(k, k+2, 0);
    }
    
    // k is an unsigned int, hence cannot be less than 0
    // decrementing k when its value is 0, leads to it attaining the maximum value
    // and it bypasses the k >= 0 criteria to stop iterations
    // Hence, the iteration corresponding to k = 0 is handled outside the loop
    x[0] = (x[0] - R.getElement(0, 1) * x[1] - R.getElement(0, 2) * x[2]) / R.getElement(0, 0);
    R.setElement(0, 2, 0);
}

template class QRSolver<long double>;
template class QRSolver<double>;
template class QRSolver<float>;
