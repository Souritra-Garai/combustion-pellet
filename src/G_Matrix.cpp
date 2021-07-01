/**
 * @file GivensRotationMatrix.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Member function definitions for GivensRotationMatrix class
 * @version 0.1
 * @date 2021-07-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "G_Matrix.hpp"

// Required for the hypotl function
#include <math.h>

template<typename real_t>
GMatrix<real_t>::GMatrix(
    RMatrix<real_t> &R,
    unsigned int index
) {
    k = index;

    real_t r_k_k    = R.getElement(k, k);
    real_t r_kp1_k  = R.getElement(k+1, k);
    
    real_t hypotenuse = hypotl(r_k_k, r_kp1_k);

    real_t cos_theta = r_k_k / hypotenuse;
    real_t sin_theta = r_kp1_k / hypotenuse;

    g_k_k     = cos_theta;
    g_k_kp1   = sin_theta;

    g_kp1_k   = - sin_theta;
    g_kp1_kp1 =   cos_theta;
}

template<typename real_t>
void GMatrix<real_t>::multiply(
    QMatrix<real_t> &Q
) {
    #pragma omp parallel for

        for (int j = 0; j <= k+1; j++)
        {
            real_t q_k_j   = Q.getElement(k, j);
            real_t q_kp1_j = Q.getElement(k+1, j);
            
            Q.setElement(
                k, j,
                g_k_k * q_k_j + g_k_kp1 * q_kp1_j
            );

            Q.setElement(
                k+1, j,
                g_kp1_k * q_k_j + g_kp1_kp1 * q_kp1_j
            );
        }
}

template<typename real_t>
void GMatrix<real_t>::multiply(
    RMatrix<real_t> &R
) {
    real_t r_k_k        = R.getElement(k, k);
    real_t r_k_kp1      = R.getElement(k, k+1);
    real_t r_kp1_k      = R.getElement(k+1, k);
    real_t r_kp1_kp1    = R.getElement(k+1, k+1);
    real_t r_kp1_kp2    = R.getElement(k+1, k+2);

    R.setElement(
        k, k,
        g_k_k * r_k_k + g_k_kp1 * r_kp1_k
    );

    R.setElement(
        k, k+1,
        g_k_k * r_k_kp1 + g_k_kp1 * r_kp1_kp1
    );

    R.setElement(
        k, k+2,
        g_k_kp1 * r_kp1_kp2
    );

    R.setElement(
        k+1, k,
        g_kp1_k * r_k_k + g_kp1_kp1 * r_kp1_k
    );

    R.setElement(
        k+1, k+1,
        g_kp1_k * r_k_kp1 + g_kp1_kp1 * r_kp1_kp1
    );

    R.setElement(
        k+1, k+2,
        g_kp1_kp1 * r_kp1_kp2
    );
}

template class GMatrix<long double>;
template class GMatrix<double>;
template class GMatrix<float>;
