/**
 * @file QR_Solver_Example.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Example cpp file to test out QRSolver class
 * @version 0.1
 * @date 2021-06-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "qrsolver/QR_Solver.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char const *argv[])
{
    const unsigned int N = 10000;

    QMatrix<long double> A(N);

    long double x[N], b[N];

    srand(time(0));

    A.setElement(0, 0, rand());
    A.setElement(0, 1, rand());

    x[0] = rand();
    b[0] = 0;

    for (int i = 1; i < N-1; i++)
    {
        A.setElement(i, i-1,    rand());
        A.setElement(i, i,      rand());
        A.setElement(i, i+1,    rand());

        x[i] = rand();
        b[i] = 0;
    }

    A.setElement(N-1, N-2, rand());
    A.setElement(N-1, N-1, rand());

    x[N-1] = rand();
    b[N-1] = 0;

    A.multiply(x, b);

    QRSolver<long double> my_solver(N);

    my_solver.setEquationFirstRow(A.getElement(0, 0), A.getElement(0, 1), b[0]);

    for (int i = 1; i < N-1; i++)
    {
        my_solver.setEquation(i, A.getElement(i, i-1), A.getElement(i, i), A.getElement(i, i+1), b[i]);
    }

    my_solver.setEquationLastRow(A.getElement(N-1, N-2), A.getElement(N-1, N-1), b[N-1]);

    // my_solver.printMatrixEquation();

    // printf("\n");

    // my_solver.printQRMatrices();

    long double x_soln[N];

    my_solver.getSolution(x_soln);

    long double MSE = 0;

    for (int i = 0; i < N; i++) MSE += (x[i] - x_soln[i]) * (x[i] - x_soln[i]);

    printf("\nMSE : %Lf\n", MSE);

    return 0;
}

