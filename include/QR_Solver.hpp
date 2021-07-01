/**
 * @file QR_Solver.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header file defines a class for solving 2D matrix equations
 * of the form A.x = b (where A is an n x n matrix and x and b are n x 1 vectors)
 * using QR factorization technique.
 * Also the class is implemented in such a way that it may be parallelized easily
 * using openmp constructs.
 * @version 0.1
 * @date 2021-06-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __QR_SOLVER__
#define __QR_SOLVER__

#include "Lower_Triangular_Matrix.hpp"
#include "GivensRotationMatrix.hpp"
#include "Tridiagonal_Matrix.hpp"

/**
 * @brief Class to implement QR factorization algorithm
 * for solving matrix equations of the A.x = b
 * where A is a n x n tridiagonal matrix and
 * x and b are n x 1 vectors
 */
template <typename real_t>
class QRSolver :
{
    private :

        // Orthogonal 2D Matrix Q in the A = Q.R factorization
        // implemented as a pointer to linear array of size \f$ N^2 \f$
        real_t *Q;

        // Variables to temporary hold non-trivial elements of 
        // Givens' rotation matrix used to vanish element \f$ r_{k+1,k} \f$
        real_t G_k_k;       // Value of \f$ g_{k,k} \f$
        real_t G_k_kp1;     // Value of \f$ g_{k,k+1} \f$

        real_t G_kp1_k;     // Value of \f$ g_{k+1,k} \f$
        real_t G_kp1_kp1;   // Value of \f$ g_{k+1,k+1} \f$

        // Variables to temporary hold values of R matrix for multiplying
        // with Givens' rotation matrix used to vanish element \f$ r_{k+1,k} \f$
        real_t R_k_k;       // Value of \f$ r_{k,k} \f$
        real_t R_k_kp1;     // Value of \f$ r_{k,k+1} \f$

        real_t R_kp1_k;     // Value of \f$ r_{k+1,k} \f$
        real_t R_kp1_kp1;   // Value of \f$ r_{k+1,k+1} \f$
        real_t R_kp1_kp2;   // Value of \f$ r_{k+1,k+2} \f$

        // Iterator
        unsigned int k;

        /**
         * @brief Loads the values of R matrix that will change during multiplication
         * with Givens' rotation matrix to temporary variables
         */
        void loadR();

        /**
         * @brief Setup Givens' rotation matrix
         */
        void setupGivensRotationMatrix();

        /**
         * @brief Multiplies the Givens rotation matrix with
         * R matrix and updates its value
         */
        void multiplyGivensMatrixWithR();
        
        /**
         * @brief Multiplies the Givens rotation matrix with
         * R matrix and updates its value
         */
        void multiplyGivensMatrixWithQ();

        /**
         * @brief Get the index in a flattened linear array representation
         * of a 2D matrix
         * 
         * @param row_index Row index i of the desired element
         * @param column_index Column index j of the desired element
         * @return Index in a linear array of the i,j element of
         * a 2D matrix implemented using the linear array
         */
        inline unsigned int getIndex(
            unsigned int row_index,
            unsigned int column_index    
        ) { 
            // Row major implementation of Q matrix
            return row_index * N + column_index;
        }
       
    public :

        /**
         * @brief Construct a new QRSolver object
         * 
         * @param N 
         */
        QRSolver(unsigned int N);

        /**
         * @brief Destroy the QRSolver object
         * 
         */
        ~QRSolver();

        void QRfactorize();

        void initQ();
};

#endif
