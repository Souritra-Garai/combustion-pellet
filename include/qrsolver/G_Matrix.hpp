/**
 * @file G_Matrix.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header serves the definition of an implementation of 
 * Givens' Rotation matrix. The rotation matrix is used to solve
 * matrix equations through QR factorization method, particularly 
 * tridiagonal matrix equation.
 * @version 0.1
 * @date 2021-06-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __GIVENS_ROTATION_MATRIX__
#define __GIVENS_ROTATION_MATRIX__

// Contains class definition for R Matrix 
// required for performing multiplication
#include "qrsolver/R_Matrix.hpp"
// Contains class definition for Q Matrix
// required for performing multiplicaiton
#include "qrsolver/Q_Matrix.hpp"

/**
 * @brief Class to implement a memory efficient 
 * Givens' Rotation Matrix that can vanish the element 
 * \f$ a_{k+1, k} \f$ of a matrix \f$ A \f$
 * 
 * A Givens rotation matrix is an orthogonal such that
 * \f$ G_{N \times N} = [g_{ij}]_{N \times N} \f$
 * where  
 * \f{equation}{
 * g_{ij}   = \begin{array}{ll}
 *              1   &   i = j \ne k, k+1 \\
 *  \cos{\theta}    &   i = j = k, k+1 \\
 *  \sin{\theta}    &   i = k, j = k+1 \\
 *  -\sin{\theta}   &   i = k+1, j = k \\
 *              0   &   \text{otherwise}
 * \end{array}
 * \f}
 * 
 * When the Givens Rotation Matrix is multiplied to another
 * matrix \f$ B_{N \times N} = [b_{i,j}]_{N \times N} \f$
 * \f{equation}{ C_{N \times N} = G_{N \times N} \cdot B_{N \times N} \f}
 * \f{equation} \Rightarrow c_{ij} = \sum_{l=1}^N g_{i,l} b_{l,j} \f}
 * \f{equation}{
 * \Rightarrow c_{ij}   = \begin{array}{ll}
 *      g_{k,k} b_{k,j} + g_{k,k+1} b_{k+1,j}   &   i=k \\
 *  g_{k+1,k} b_{k,j} + g_{k+1,k+1} b_{k+1,j}   &   i=k+1 \\
 *                                  b_{i,j}     &   \text{otherwise}
 * \end{array}
 * \f}
 *  
 * @tparam real_t float, double or long double data types to represent real numbers
 */
template<typename real_t>
class GMatrix
{
    public :

        /**
         * @brief Construct a new Givens Rotation Matrix to
         * vanish the element at position \f$ (k+1, k) \f$
         * of R matrix.
         * 
         * @param R R matrix whose element at position 
         * \f$ (k+1, k) \f$ needs to be vanished
         * @param index k
         */
        GMatrix(
            RMatrix<real_t> &matrix,
            unsigned int index
        );

        /**
         * @brief Multiplies the Givens Rotation Matrix to the 
         * R Matrix and updates the R matrix in place
         * 
         * @param R R Matrix that will be converted into
         * upper tridiagonal matrix after multiplication
         * 
         * R matrix can be represented as \f$ R_{N \times N} = [r_{i,j}]_{N \times N} \f$
         * where
         * \f{equation}{
         *  r_{i,j}   = \begin{array}{lllll}
         *      0 & \text{if} & j < i-1 & \text{or} & j > i+1
         * \end{array}
         * \f}
         * 
         * Thus, upon multiplying the R matrix with the rotation matrix,
         * we get matrix \f$ C_{N \times N} \f$
         * 
         * \f{equation}{
         *  c_{ij}   = \begin{array}{ll}
         *      g_{k,k} r_{k,j} + g_{k,k+1} r_{k+1,j}       &   i=k \\
         *      g_{k+1,k} r_{k,j} + g_{k+1,k+1} r_{k+1,j}   &   i=k+1 \\
         *                                      r_{i,j}     &   \text{otherwise}
         * \end{array}
         * \f}
         * 
         * \f{equation}{
         *  \Rightarrow c_{ij}   = \begin{array}{lll}
         *                              g_{k,k} r_{k,k-1}   &   i=k &   j=k-1   \\
         *          g_{k,k} r_{k,k} + g_{k,k+1} r_{k+1,k}   &   i=k &   j=k     \\
         *      g_{k,k} r_{k,k+1} + g_{k,k+1} r_{k+1,k+1}   &   i=k &   j=k+1   \\
         *                          g_{k,k+1} r_{k+1,k+2}   &   i=k &   j=k+2   \\
         * 
         *                          g_{k+1,k} r_{k,k-1}     &   i=k+1   &   j=k-1   \\
         *      g_{k+1,k} r_{k,k} + g_{k+1,k+1} r_{k+1,k}   &   i=k+1   &   j=k     \\
         *  g_{k+1,k} r_{k,k+1} + g_{k+1,k+1} r_{k+1,k+1}   &   i=k+1   &   j=k+1   \\
         *                      g_{k+1,k+1} r_{k+1,k+2}     &   i=k+1   &   j=k+2   \\
         * 
         *                                      r_{i,j}     &   \text{otherwise}
         * \end{array}
         * \f}
         * 
         * But as we iterate over k and multiply the corresponding rotation marix,
         * the sub-diagonal (\f$ r_{i,i-1} \f$) becomes zero and the super-super-diagonal
         * (\f$ r_{i,i+2} \f$) fills up with non-zero values. 
         * 
         * \f{equation}{
         *  \Rightarrow c_{ij}   = \begin{array}{lll}
         *                                              0   &   i=k &   j=k-1   \\
         *          g_{k,k} r_{k,k} + g_{k,k+1} r_{k+1,k}   &   i=k &   j=k     \\
         *      g_{k,k} r_{k,k+1} + g_{k,k+1} r_{k+1,k+1}   &   i=k &   j=k+1   \\
         *                          g_{k,k+1} r_{k+1,k+2}   &   i=k &   j=k+2   \\
         * 
         *                                              0   &   i=k+1   &   j=k-1   \\
         *                                              0   &   i=k+1   &   j=k     \\
         *  g_{k+1,k} r_{k,k+1} + g_{k+1,k+1} r_{k+1,k+1}   &   i=k+1   &   j=k+1   \\
         *                      g_{k+1,k+1} r_{k+1,k+2}     &   i=k+1   &   j=k+2   \\
         * 
         *                                      r_{i,j}     &   \text{otherwise}
         * \end{array}
         * \f}
         */
        void multiply(RMatrix<real_t> &R);

        /**
         * @brief Multiplies the Givens rotation matrix to the R matrix
         * where the rotation matrix is formed to vanish the sub diagonal
         * element of the last row \f$ (N, N-1) \f$ of the \f$ N \times N \f$
         * R matrix 
         * 
         * \f{equation}{
         *  \Rightarrow c_{ij}   = \begin{array}{lll}
         *                                              0   &   i=k &   j=k-1   \\
         *  g_{N-1,N-1} r_{N-1,N-1} + g_{N-1,N} r_{N,N-1}   &   i=k &   j=k     \\
         *      g_{N-1,N-1} r_{N-1,N} + g_{N-1,N} r_{N,N}   &   i=k &   j=k+1   \\
         * 
         *                                              0   &   i=k+1   &   j=k-1   \\
         *                                              0   &   i=k+1   &   j=k     \\
         *          g_{N,N-1} r_{N-1,N} + g_{N,N} r_{N,N}   &   i=k+1   &   j=k+1   \\
         * 
         *                                      r_{i,j}     &   \text{otherwise}
         * \end{array}
         * \f}
         * 
         * @param R R Matrix that will be converted into
         * upper tridiagonal matrix after multiplication
         */
        void multiplyLastRow(RMatrix<real_t> &R);

        /**
         * @brief Multiplies the Givens Rotation Matrix to the 
         * Q Matrix and updates the Q matrix in place
         * 
         * @param Q Q Matrix passed as reference
         * 
         * Q matrix may be represented as
         * \f$ Q_{N \times N} = [q_{i,j}]_{N \times N} \f$ where
         * \f{eqnarray*}{
         *  q_{i,j} = 0, & j > i
         * \f}
         * 
         * Thus, when a rotation matrix is multiplied to Q matrix
         * \f{equation}{ C_{N \times N} = G_{N \times N} \cdot Q_{N \times N} \f}
         * \f{equation} \Rightarrow c_{ij} = \sum_{l=1}^N g_{i,l} q_{l,j} \f}
         * \f{equation}{
         * \Rightarrow c_{ij}   = \begin{array}{ll}
         *      g_{k,k} l_{k,j} + g_{k,k+1} q_{k+1,j}   &   i=k \\
         *  g_{k+1,k} l_{k,j} + g_{k+1,k+1} q_{k+1,j}   &   i=k+1 \\
         *                                  q_{i,j}     &   \text{otherwise}
         * \end{array}
         * \f}
         */
        void multiply(QMatrix<real_t> &Q);

    private :

        /**
         * @brief Element of Givens Rotation Matrix at position \f$ k, k \f$
         */
        real_t g_k_k;
        
        /**
         * @brief Element of Givens Rotation Matrix at position \f$ k, k+1 \f$
         */
        real_t g_k_kp1;

        /**
         * @brief Element of Givens Rotation Matrix at position \f$ k+1, k \f$
         */
        real_t g_kp1_k;
        
        /**
         * @brief Element of Givens Rotation Matrix at position \f$ k+1, k+1 \f$
         */
        real_t g_kp1_kp1;

        /**
         * @brief Index k, where we want to vanish the element 
         * at position \f$ k+1, k \f$
         */
        unsigned int k;
};

#endif