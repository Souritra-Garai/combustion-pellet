/**
 * @file QR_Solver.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header file defines a class for solving 2D matrix equations
 * of the form \f$ A \cdot x = b \f$ (where \f$ A \f$ is an \f$ N \times N \f$ matrix 
 * and \f$ x \f$ and \f$ b \f$ are \f$ N \times 1 \f$ vectors)
 * using QR factorization technique.
 * 
 * @version 0.1
 * @date 2021-06-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __QR_SOLVER__
#define __QR_SOLVER__

#include "qrsolver/Q_Matrix.hpp"
#include "qrsolver/R_Matrix.hpp"
#include "qrsolver/G_Matrix.hpp"

/**
 * @brief Class to implement QR factorization algorithm
 * for solving matrix equations of the \f$ A \cdot x = b \f$
 * where \f$ A \f$ is a \f$ N \times N \f$ tridiagonal matrix and
 * \f$ x \f$ and \f$ b \f$ are \f$ N \times 1 \f$ vectors
 * 
 * @tparam real_t 
 */
template <typename real_t>
class QRSolver
{
    private :

        /**
         * @brief \f$ Q \f$ matrix for QR Factorization
         */
        QMatrix<real_t> Q;

        /**
         * @brief \f$ R \f$ matrix for QR Factorization
         * 
         * The \f$ R \f$ matrix also serves as the initial tridiagonal
         * \f$ A \f$ matrix to save memory and more important reduce
         * redundant memory copy operations
         */
        RMatrix<real_t> R;

        /**
         * @brief Number of rows \f$ N \f$ of the matrix \f$ A \f$
         */
        const unsigned int N;

        /**
         * @brief One dimensional array to store \f$ N \times N \f$
         * vector \f$ b \f$
         */
        real_t *b;

        /**
         * @brief Iterator
         */
        unsigned int k;

        /**
         * @brief Factorizes the tridiagonal \f$ A \f$ matrix
         * stored in R to \f$ Q \cdot R \f$
         */
        void QRFactorize();
       
    public :

        /**
         * @brief Construct a new QRSolver object
         * 
         * @param N Size of the matrix equations
         */
        QRSolver(unsigned int N);

        /**
         * @brief Destroy the QRSolver object
         */
        ~QRSolver();

        /**
         * @brief Set up equation represented by ith
         * row of the matrix equation
         * \f$ e x_{i-1} + f x_i + g x_{i+1} = b \f$
         * 
         * @param index i
         * @param e Coefficient to \f$ x_{i-1} \f$
         * @param f Coefficient to \f$ x_{i} \f$
         * @param g Coefficient to \f$ x_{i+1} \f$
         * @param b Constant
         */
        void setEquation(
            unsigned int index,
            real_t e,
            real_t f,
            real_t g,
            real_t b
        );

        /**
         * @brief Set up equation represented by the first row
         * of the matrix equation
         * \f$ f x_i + g x_{i+1} = b \f$
         * 
         * @param f Coefficient to \f$ x_{i} \f$
         * @param g Coefficient to \f$ x_{i+1} \f$
         * @param b Constant
         */
        void setEquationFirstRow(
            real_t f,
            real_t g,
            real_t b
        );

        /**
         * @brief up equation represented by the last row
         * of the matrix equation
         * \f$ e x_{i-1} f x_i = b \f$
         * 
         * @param e Coefficient to \f$ x_{i-1} \f$
         * @param f Coefficient to \f$ x_{i} \f$
         * @param b Constant
         */
        void setEquationLastRow(
            real_t e,
            real_t f,
            real_t b
        );

        /**
         * @brief Prints the matrix \f$ A \f$ and vector \f$ b \f$
         * 
         * To be used before QR factorization is performed in getSolution
         */
        void printMatrixEquation();

        /**
         * @brief Prints the factorized matrices \f$ Q \cdot R \f$
         * 
         * To be used after QR factorization is performed in getSolution
         */
        void printQRMatrices();

        /**
         * @brief Finds the solution to matrix equation and saves
         * it to array x
         * 
         * @param x Array to store the solution of the matrix equation
         */
        void getSolution(real_t *x);
};

#endif
