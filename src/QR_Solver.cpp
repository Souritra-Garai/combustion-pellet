/**
 * @file QR_Solver.cpp
 * @author Souritra Garai (souritra,garai@iitgn.ac.in)
 * @brief Implementation of QR_Solver class
 * @version 0.1
 * @date 2021-07-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "QR_Solver.hpp"

#include <stdio.h>
#include <bits/stdc++.h>

template<typename real_t>
QRSolver<real_t>::QRSolver(
    unsigned int n
) : // Initialise the base classes and N in the mem initialisation list
    Q(n),
    R(n),
    N(n)
{
    b = new real_t[N];
}

template<typename real_t>
QRSolver<real_t>::~QRSolver()
{
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
    R.setElement(i, i-1, e_val);
    R.setElement(i, i,   f_val);
    R.setElement(i, i+1, g_val);

    b[i] = b_val;
}

template<typename real_t>
void QRSolver<real_t>::setEquationFirstRow(
    real_t f_val,
    real_t g_val,
    real_t b_val
) {
    R.setElement(0, 0, f_val);
    R.setElement(0, 1, g_val);

    b[0] = b_val;
}

template<typename real_t>
void QRSolver<real_t>::setEquationLastRow(
    real_t e_val,
    real_t f_val,
    real_t b_val
) {
    R.setElement(N-1, N-2, e_val);
    R.setElement(N-1, N-1, f_val);

    b[N-1] = b_val;
}

template<typename real_t>
void QRSolver<real_t>::printMatrixEquation()
{
    printf("Matrix Eqn A.x = b\n");
    printf("\nA_{%d \\times %d}\n", N, N);

    R.printMatrix();

    printf("\nb_{%d \\times 1}\n", N);

    for (k = 0; k < N; k++) printf("%Lf\n", (long double) b[k]);
}

template<typename real_t>
void QRSolver<real_t>::QRFactorize()
{
    Q.identity();

    for (k = 0; k < N-2; k++)
    {
        GMatrix<real_t> G(R, k);

        G.multiply(R);
        G.multiply(Q);
    }

    GMatrix<real_t> G(R, k);

    G.multiplyLastRow(R);
    G.multiply(Q);
}

template<typename real_t>
void QRSolver<real_t>::printQRMatrices()
{
    printf("A = Q.R\n");

    printf("\nQ^T_{%d \\times %d}\n", N, N);
    Q.printMatrix();

    printf("\nR_{%d \\times %d}\n", N, N);
    R.printMatrix();
}

template<typename real_t>
void QRSolver<real_t>::getSolution(
    real_t *x
) {
    memset(x, 0, N * sizeof(real_t));

    Q.multiply(b, x);

    x[N-1] = x[N-1] / R.getElement(N-1, N-1);

    x[N-2] = (x[N-2] - R.getElement(N-2, N-1) * x[N-1]) / R.getElement(N-2, N-2);

    for (k = N-3; k > 0; k--)
    {
        x[k] = (x[k] - R.getElement(k, k+1) * x[k+1] - R.getElement(k, k+2) * x[k+2]) / R.getElement(k, k);
    }

    x[0] = (x[0] - R.getElement(0, 1) * x[1] - R.getElement(0, 2) * x[2]) / R.getElement(0, 0);

}

template class QRSolver<long double>;
template class QRSolver<double>;
template class QRSolver<float>;
