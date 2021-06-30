/**
 * @file GivensRotationMatrix.hpp
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

// Contains class definition for Tridiagonal Matrix
#include <Tridiagonal_Matrix.hpp>
// Contains class definition for Lower Triangular Matrix
#include <Lower_Triangular_Matrix.hpp>

/**
 * @brief Class to implement a memory efficient 
 * Givens' Rotation Matrix that can vanish the element 
 * \f$ a_{k+1, k} \f$ of a tridiagonal matrix \f$ A \f$
 * 
 * A Givens rotation matrix is an orthogonal such that
 * \f$ G_{N \times N} = [g_{ij}]_{N \times N} \f$
 * where  
 * \f{equation}{
 * g_{ij}   = \begin{array}{cc}
 *              1   &   i = j \ne k, k+1 \\
 *  \cos{\theta}    &   i = j = k, k+1 \\
 *  \sin{\theta}    &   i = k, j = k+1 \\
 *  -\sin{\theta}   &   i = k+1, j = k \\
 *              0   &   \textbf{otherwise}
 * \end{array}
 * \f}
 * 
 * When the Givens Rotation Matrix is multiplied to another
 * matrix \f$ B_{N \times N} = [b_{i,j}]_{N \times N} \f$
 * \f{equation}{ C_{N \times N} = G_{N \times N} \cdot B_{N \times N} \f}
 * \f{equation} \Rightarrow c_{ij} = \sum_{l=1}^N g_{i,l} b_{l,j} \f}
 * \f{equation}{
 * \Rightarrow c_{ij}   = \begin{array}{cc}
 *  g_{k,k} b_{k,j} + g_{k,k+1} b_{k+1,j}       &   i=k \\
 *  g_{k+1,k} b_{k,j} + g_{k+1,k+1} b_{k+1,j}   &   i=k+1 \\
 *                                      b_{i,j} &   \textbf{otherwise}
 * \end{array}
 * \f}
 *  
 * @tparam real_t float, double or long double data types to represent real numbers
 * 
 */
template<typename real_t>
class GivensRotationMatrix
{
    public :

        /**
         * @brief Construct a new Givens Rotation Matrix
         * 
         */
        GivensRotationMatrix();

        /**
         * @brief Sets up the Givens Rotation Matrix to
         * vanish the element at position \f$ (i+1, i) \f$
         * of tridiagonal matrix.
         * 
         * @param matrix Tridiagonal matrix whose element at position 
         * \f$ (i+1, i) \f$ needs to be vanished
         * @param index i
         */
        void setupRotationMatrix(
            TridiagonalMatrix<real_t> &matrix,
            unsigned int index
        );

        /**
         * @brief Multiplies the Givens Rotation Matrix to the 
         * Tridiagonal Matrix and updates the tridiagonal matrix in place
         * 
         * @param matrix Tridiagonal Matrix that will be converted into
         * upper tridiagonal matrix after multiplication
         * 
         * For multiplication with a Tridiagonal matrix, \f$ T_{N \times N} \f$
         * where
         * \f{eqnarray*}{
         *  t_{i,j}   = 0, & \text{if} j < i-1 \text{or} j > i+1
         * }
         * 
         * Thus, upon multiplying the tridiagonal matrix with the rotation matrix,
         * we get matrix \f$ C_{N \times N} \f$
         * 
         * \f{equation}{
         *  c_{ij}   = \begin{array}{cc}
         *      g_{k,k} t_{k,j} + g_{k,k+1} t_{k+1,j}       &   i=k \\
         *      g_{k+1,k} t_{k,j} + g_{k+1,k+1} t_{k+1,j}   &   i=k+1 \\
         *                                          b_{i,j} &   \textbf{otherwise}
         * \end{array}
         * 
         * \f{equation}{
         *  \Rightarrow c_{ij}   = \begin{array}{cc}
         *                              g_{k,k} t_{k,k-1}   &   i=k &   j=k-1   \\
         *          g_{k,k} t_{k,k} + g_{k,k+1} t_{k+1,k}   &   i=k &   j=k     \\
         *      g_{k,k} t_{k,k+1} + g_{k,k+1} t_{k+1,k+1}   &   i=k &   j=k+1   \\
         *                          g_{k,k+1} t_{k+1,k+2}   &   i=k &   j=k+2   \\
         * 
         *      g_{k+1,k} t_{k,j} + g_{k+1,k+1} t_{k+1,j}   &   i=k+1 \\
         *      g_{k+1,k} t_{k,j} + g_{k+1,k+1} t_{k+1,j}   &   i=k+1 \\
         *      g_{k+1,k} t_{k,j} + g_{k+1,k+1} t_{k+1,j}   &   i=k+1 \\
         *      g_{k+1,k} t_{k,j} + g_{k+1,k+1} t_{k+1,j}   &   i=k+1 \\
         *                                          b_{i,j} &   \textbf{otherwise}
         * \end{array}
         */
        void multiply(TridiagonalMatrix<real_t> &matrix);
        void multiply(LowerTriangularMatrix<real_t> &matrix);

    private :

        real_t element_i_i;
        real_t element_i_ip1;

        real_t element_ip1_i;
        real_t element_ip1_ip1;

        unsigned int i;
};

#endif