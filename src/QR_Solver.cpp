#include "QR_Solver.hpp"

#include <bits/stdc++.h>

QRSolver::QRSolver(unsigned int n) :
    // Call base class TridiagonalMatrix constructor in
    // mem initialization list 
    TridiagonalMatrix(n)
{
    // Allocate memory for Q matrix
    Q = new real_t[N*N];
}

QRSolver::~QRSolver()
{
    // Deallocate memory assigned to Q matrix
    delete [] Q;
}

void QRSolver::loadR()
{
    // TridiagonalMatrix holds the R Matrix
    // Get the corresponding elements from the Tridiagonal Matrix

    // Row k
    // r_{k,k-1} is already zero
    R_k_k       = getElement(k, k);
    R_k_kp1     = getElement(k, k+1);
    // Row k+1
    R_kp1_k     = getElement(k+1, k);
    R_kp1_kp1   = getElement(k+1, k+1);
    R_kp1_kp2   = getElement(k+1, k+2);
}

void QRSolver::setupGivensRotationMatrix()
{
    // 
    real_t hypotenuse = hypotl(R_k_k, R_kp1_k);
        
    real_t cos_theta = R_k_k    / hypotenuse;
    real_t sin_theta = R_kp1_k  / hypotenuse;

    G_k_k    = cos_theta;
    G_k_kp1  = sin_theta;

    G_kp1_k      = - sin_theta;
    G_kp1_kp1    =   cos_theta;
}

void QRSolver::multiplyGivensMatrixWithR()
{
    // A = [a_i_j]_nxn
    // B = [b_i_j]_nxn
    // C = [c_i_j]_nxn = A . B
    // c_i_j = sum_over_l a_i_l * b_l_j

    // a_i_j =  1,          i = j != k, k+1
    //          a_k_k,      i = j = k
    //          a_k_kp1,    i = k, j = k + 1
    //          a_kp1_k,    i = k + 1, j = k
    //          a_kp1_kp1,  i = j = k + 1

    // c_i_j =  b_i_j,                                  i != k, k+1
    //          a_k_k * b_k_j   + a_k_kp1 * b_kp1_j     i = k
    //          a_kp1_k * b_k_j + a_kp1_kp1 * b_kp1_j   i = k + 1

    // b_i_j =  b_i_im1,    j = i - 1
    //          b_i_i,      j = i
    //          b_i_ip1,    j = 1 + 1

    // c_i_j =  b_i_im1,                                    j = i - 1, i != k, k+1
    //          b_i_i,                                      j = i, i != k, k+1
    //          b_i_ip1,                                    j = i + 1, i != k, k+1
    //          a_k_k * b_k_km1                             i = k, j = k - 1
    //          a_k_k * b_k_k + a_k_kp1 * b_kp1_k           i = j = k
    //          a_k_k * b_k_kp1 + a_k_kp1 * b_kp1_kp1       i = k, j = k + 1
    //          a_k_kp1 * b_kp1_kp2                         i = k, j = k + 2
    //          a_kp1_k * b_k_km1                           i = k + 1, j = k - 1
    //          a_kp1_k * b_k_k + a_kp1_kp1 * b_kp1_k       i = k + 1, j = k
    //          a_kp1_k * b_k_kp1 + a_kp1_kp1 * b_kp1_kp1   i = j = k + 1
    //          a_kp1_k * b_k_kp2 + a_kp1_kp1 * b_kp1_kp2   i = k + 1, j = k + 2

    setElement(k,   k,      G_k_k * R_k_k   + G_k_kp1 * R_kp1_k);
    setElement(k,   k+1,    G_k_k * R_k_kp1 + G_k_kp1 * R_kp1_kp1);
    setElement(k,   k+2,    G_k_kp1 * R_kp1_kp2);

    setElement(k+1, k+1,    G_kp1_k * R_k_kp1 + G_kp1_kp1 * R_kp1_kp1);
    setElement(k+1, k+2,    G_kp1_kp1 * R_kp1_kp2);
}

void QRSolver::multiplyGivensMatrixWithQ()
{
    #pragma omp parallel for

        for (int j = 0; j < k+1; j++)
        {
            real_t Q_k_j    = Q[getIndex(k, j)];
            real_t Q_kp1_j  = Q[getIndex(k+1, j)];
            
            Q[getIndex(k, j)]      = G_k_k * Q_k_j     + G_k_kp1 * Q_kp1_j;
            Q[getIndex(k+1, j)]    = G_kp1_k * Q_k_j   + G_kp1_kp1 * Q_kp1_j;
        }
}

void QRSolver::QRfactorize(
)
{
    initQ();

    for (int k = 0; k < N-2; k++)
    {
        ;
    }
}



void QRSolver::initQ()
{
    memset(Q, 0, sizeof(Q));

    #pragma omp parallel for

        for (int i = 0; i < N; i++)

            Q[getIndex(i, i)] = 0.0;
}