/**
 * @file QR_Solver_Example.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Example cpp file to test out QR_Solver functions
 * @version 0.1
 * @date 2021-06-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "QR_Solver.hpp"

QR_Solver my_solver(5);

real_t A[3*5];

int main(int argc, char const *argv[])
{
    int k = 1;
    for (int i = 1; i < 5; i++)
    {
        A[getTridiagonalIndex(i, i-1)] = k;
        k++;

        A[getTridiagonalIndex(i, i)] = k;
    }
    return 0;
}
