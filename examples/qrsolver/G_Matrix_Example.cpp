/**
 * @file G_Matrix_Example.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Example to test out GMatrix class
 * @version 0.1
 * @date 2021-07-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <stdio.h>

#include "qrsolver/G_Matrix.hpp"

int main(int argc, char const *argv[])
{
    RMatrix<long double> R(5);
    int k = 1;
    R.setElement(0, 0, k++);
    R.setElement(0, 1, k++);

    for (int i = 1; i < 4; i++)
    {
        R.setElement(i, i-1, k++);
        R.setElement(i, i, k++);
        R.setElement(i, i+1, k++);
    }

    R.setElement(4, 3, k++);
    R.setElement(4, 4, k++);

    printf("R Matrix\n");
    R.printMatrix();

    GMatrix<long double> G(R, 0);

    G.multiply(R);
    printf("\nR Matrix\n");
    R.printMatrix();

    QMatrix<long double> Q(5);
    Q.identity();

    G.multiply(Q);
    printf("\nQ Matrix\n");
    Q.printMatrix();

    return 0;
}
