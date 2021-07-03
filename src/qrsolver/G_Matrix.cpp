/**
 * @file G_Matrix.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Member function definitions for GMatrix class
 * @version 0.1
 * @date 2021-07-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "qrsolver/G_Matrix.hpp"

// Required for the hypotl function
#include <math.h>

template<typename real_t>
GMatrix<real_t>::GMatrix(
    RMatrix<real_t> &R,
    unsigned int index
) {
    // To vanish element \f$ r_{k+1,k} \f$ of matrix \f$ R \f$
    k = index;

    // Retrieve the elements \f$ r_{k,k} \f$ ans \f$ r_{k+1,k} \f$
    real_t r_k_k    = R.getElement(k, k);
    real_t r_kp1_k  = R.getElement(k+1, k);
    
    // Calculate the hypotenuse \f$ \sqrt{r_{k,k}^2 + r_{k+1,k}^2} \f$
    real_t hypotenuse = hypotl(r_k_k, r_kp1_k);

    // Calculate \f$ \cos{\theta} = r_{k,k} / \sqrt{r_{k,k}^2 + r_{k+1,k}^2} \f$
    real_t cos_theta = r_k_k / hypotenuse;
    // Calculate \f$ \sin{\theta} = r_{k+1,k} / \sqrt{r_{k,k}^2 + r_{k+1,k}^2} \f$
    real_t sin_theta = r_kp1_k / hypotenuse;

    // \f$ g_{k,k} = \cos{\theta} \f$
    g_k_k     = cos_theta;
    // \f$ g_{k,k+1} = \sin({\theta}) \f$
    g_k_kp1   = sin_theta;

    // \f$ g_{k+1,k} = - \sin(\theta) \f$
    g_kp1_k   = - sin_theta;
    // \f$ g_{k+1,k+1} = \cos{\theta} \f$
    g_kp1_kp1 =   cos_theta;
}

template<typename real_t>
void GMatrix<real_t>::multiply(
    QMatrix<real_t> &Q
) {
    // Since the multiplication along each column of
    // \f$ Q \f$ is independent of the other, we can 
    // safely parallelize
    #pragma omp parallel for

        // Along each row j of the \f$ Q \f$
        for (int j = 0; j <= k+1; j++)
        {
            // Make copies of \f$ q_{k,j} \f$ and \f$ q_{k+1,j} \f$
            // as they will be changed
            real_t q_k_j   = Q.getElement(k, j);
            real_t q_kp1_j = Q.getElement(k+1, j);
            
            // Set \f$ q_{k,j} = g_{k,k} q_{k,j} + g_{k,k+1} q_{k+1,j} \f$
            Q.setElement(
                k, j,
                g_k_k * q_k_j + g_k_kp1 * q_kp1_j
            );

            // Set \f$ q_{k+1,j} = g_{k+1,k} q_{k,j} + g_{k+1,k+1} q_{k+1,j} \f$
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
    // Create copies of elements
    // \f$ r_{k,k} \f$, \f$ r_{k,k+1} \f$, \f$ r_{k+1,k} \f$
    // \f$ r_{k+1,k+1} \f$ and \f$ r_{k+1,k+2} \f$
    real_t r_k_k        = R.getElement(k, k);
    real_t r_k_kp1      = R.getElement(k, k+1);
    real_t r_kp1_k      = R.getElement(k+1, k);
    real_t r_kp1_kp1    = R.getElement(k+1, k+1);
    real_t r_kp1_kp2    = R.getElement(k+1, k+2);

    // \f$ r_{k,k-1} \f$ and \f$ r_{k+1,k-1} \f$ are not required
    // they are and will remain zero

    // Set \f$ r_{k,k} = g_{k,k} r_{k,k} + g_{k,k+1} r_{k+1,k} \f$
    R.setElement(
        k, k,
        g_k_k * r_k_k + g_k_kp1 * r_kp1_k
    );

    // Set \f$ r_{k,k+1} = g_{k,k} r_{k,k+1} + g_{k,k+1} r_{k+1,k+1} \f$
    R.setElement(
        k, k+1,
        g_k_k * r_k_kp1 + g_k_kp1 * r_kp1_kp1
    );

    // Set \f$ r_{k,k+2} = g_{k,k} r_{k,k+2} + g_{k,k+1} r_{k+1,k+2} \f$
    // But \f$ r_{k,k+2} \f$ is initially zero
    R.setElement(
        k, k+2,
        g_k_kp1 * r_kp1_kp2
    );

    // Set \f$ r_{k+1,k} = g_{k+1,k} r_{k,k} + g_{k+1,k+1} r_{k+1,k} = 0 \f$
    R.setElement(
        k+1, k,
        0   // g_kp1_k * r_k_k + g_kp1_kp1 * r_kp1_k
    );

    // Set \f$ r_{k+1,k+1} = g_{k+1,k} r_{k,k+1} + g_{k+1,k+1} r_{k+1,k+1} \f$
    R.setElement(
        k+1, k+1,
        g_kp1_k * r_k_kp1 + g_kp1_kp1 * r_kp1_kp1
    );

    // Set \f$ r_{k+1,k+2} = g_{k+1,k} r_{k,k+2} + g_{k+1,k+1} r_{k+1,k+2} \f$
    // But \f$ r_{k,k+2} \f$ is zero initially
    R.setElement(
        k+1, k+2,
        g_kp1_kp1 * r_kp1_kp2
    );
}

template<typename real_t>
void GMatrix<real_t>::multiplyLastRow(
    RMatrix<real_t> &R
) {
    // \f$ r_{k,k+2} \f$ and \f$ r_{k+1,k+2} \f$
    // do not exist for the last row

    // Create copies of elements
    // \f$ r_{k,k} \f$, \f$ r_{k,k+1} \f$, \f$ r_{k+1,k} \f$
    // \f$ r_{k+1,k+1} \f$ \f$
    real_t r_k_k        = R.getElement(k, k);
    real_t r_k_kp1      = R.getElement(k, k+1);
    real_t r_kp1_k      = R.getElement(k+1, k);
    real_t r_kp1_kp1    = R.getElement(k+1, k+1);

    // \f$ r_{k,k-1} \f$ and \f$ r_{k+1,k-1} \f$ are not required
    // they are and will remain zero

    // Set \f$ r_{k,k} = g_{k,k} r_{k,k} + g_{k,k+1} r_{k+1,k} \f$
    R.setElement(
        k, k,
        g_k_k * r_k_k + g_k_kp1 * r_kp1_k
    );

    // Set \f$ r_{k,k+1} = g_{k,k} r_{k,k+1} + g_{k,k+1} r_{k+1,k+1} \f$
    R.setElement(
        k, k+1,
        g_k_k * r_k_kp1 + g_k_kp1 * r_kp1_kp1
    );

    // Set \f$ r_{k+1,k} = g_{k+1,k} r_{k,k} + g_{k+1,k+1} r_{k+1,k} = 0 \f$
    R.setElement(
        k+1, k,
        0   // g_kp1_k * r_k_k + g_kp1_kp1 * r_kp1_k
    );

    // Set \f$ r_{k+1,k+1} = g_{k+1,k} r_{k,k+1} + g_{k+1,k+1} r_{k+1,k+1} \f$
    R.setElement(
        k+1, k+1,
        g_kp1_k * r_k_kp1 + g_kp1_kp1 * r_kp1_kp1
    );
}

template class GMatrix<long double>;
template class GMatrix<double>;
template class GMatrix<float>;
