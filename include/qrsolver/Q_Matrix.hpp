/**
 * @file Q_Matrix.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header file defines a class for memory efficient 
 * implementation of \f$ Q \f$ matrix used for QR factorisation of tridiagonal matrix
 * using Givens rotation matrix
 * @version 0.1
 * @date 2021-06-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __LOWER_TRIANGULAR_MATRIX__
#define __LOWER_TRIANGULAR_MATRIX__

/**
 * @brief Class to implement a memory efficient model of
 * \f$ N \times N \f$ Q Matrix.
 * 
 * @tparam real_t float, double or long double data types
 * to represent real numbers
 */
template <typename real_t>
class QMatrix
{
    public:
        
        /**
         * @brief Construct a new Q Matrix
         * 
         * @param n Number of rows in the \f$ N \times N \f$ square matrix
         */
        QMatrix(unsigned int n);

        /**
         * @brief Destroy the Q Matrix object
         */
        ~QMatrix();

        /**
         * @brief Get the i,j th element of the Q Matrix
         * 
         * @param row_index Row index i
         * @param column_index Column index j
         * @return Value of the i,j th element of the Q Matrix
         */
        real_t getElement(
            unsigned int row_index,
            unsigned int column_index
        );

        /**
         * @brief Set the value of the i,j th element of the Q Matrix
         * 
         * @param row_index Row index i
         * @param column_index Column index j
         * @param value Value to be set at the i,j th position
         */
        void setElement(
            unsigned int row_index,
            unsigned int column_index,
            real_t value
        );

        /**
         * @brief Prints the Q Matrix in form of a 2D array
         */
        void printMatrix();

        /**
         * @brief Makes Q matrix an identity matrix
         */
        void identity();
        
        /**
         * @brief Multiplies the Q matrix with the column vector b
         * and stores the result in the column vector x
         * 
         * @param b 
         * @param x 
         */
        void multiply(real_t *b, real_t *x);
        
    private :

        /**
         * @brief One dimensional array to store only non zero
         * elements of the lower triangular matrix
         */
        real_t *array;

        /**
         * @brief Number of rows in a \f$ N \times N \f$
         * square matrix
         */
        const unsigned int N;

        /**
         * @brief Get the index of i,j th element of Q matrix in the
         * flattened array
         * 
         * @param row_index Row index i
         * @param column_index Column index j
         * @return Index of the i,j th element in the flattened array
         */
        unsigned int getIndex(
            unsigned int row_index,
            unsigned int column_index
        );
        
        /**
         * @brief Checks if the row index i and the column index j
         * belong to a zero element of the Q matrix
         * 
         * @param row_index Row index i
         * @param column_index Column index j
         * @return true if i,j are indices of zero elements in a Q matrix
         * @return false if i,j are indices of non zero elements in a Q matrix
         */
        bool indexOfZeroElement(
            unsigned int row_index,
            unsigned int column_index
        );
};

#endif