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

#include "Q_Matrix.hpp"
#include "R_Matrix.hpp"
#include "G_Matrix.hpp"

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

        QMatrix<real_t> Q;

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

        void setEquation(
            unsigned int index,
            real_t e,
            real_t f,
            real_t g,
            real_t b
        );

        void setEquationFirstRow(
            real_t f,
            real_t g,
            real_t b
        );

        void setEquationLastRow(
            real_t e,
            real_t f,
            real_t b
        );

        void printMatrixEquation();

        void QRFactorize();

        void printQRMatrices();

        void getSolution(real_t *x);
};

#endif
