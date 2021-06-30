/**
 * @file Lower_Triangular_Matrix.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header file defines a class for memory efficient 
 * implementation of lower triangualer square matrices
 * @version 0.1
 * @date 2021-06-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __LOWER_TRIANGULAR_MATRIX__
#define __LOWER_TRIANGULAR_MATRIX__

/**
 * @brief Class to implement a memory efficient model of Lower Triangular
 * Square Matrix.
 * 
 * @tparam real_t float, double or long double data types to represent real numbers
 */
template <typename real_t>
class LowerTriangularMatrix
{
    public:
        
        /**
         * @brief Construct a new Lower Triangular Matrix object
         * 
         * @param n Number of rows in the \f$ N \times N \f$ square matrix
         */
        LowerTriangularMatrix(unsigned int n);

        /**
         * @brief Destroy the Lower Triangular Matrix object
         * 
         */
        ~LowerTriangularMatrix();

        /**
         * @brief Get the i,j th element of the Lower Triangular Matrix
         * 
         * @param row_index Row index i
         * @param column_index Column index j
         * @return Value of the i,j th element of the Lower Triangular Matrix
         */
        real_t getElement(
            unsigned int row_index,
            unsigned int column_index
        );

        /**
         * @brief Set the value of the i,j th element of the Lower Triangular Matrix
         * 
         * @param row_index Row index i
         * @param column_index Column index j
         * @param value Value to be set at the i,j th element
         */
        void setElement(
            unsigned int row_index,
            unsigned int column_index,
            real_t value
        );

        /**
         * @brief Prints the Lower Triangular Matrix in form of a 2D array
         * 
         */
        void printMatrix();
        
    private :

        /**
         * @brief One dimensional array to store only non zero
         * elements of the lower triangular matrix
         */
        real_t *array;

        /**
         * @brief Size of main diagonal of \f$ N \times N \f$
         * square matrix
         */
        const unsigned int N;

        /**
         * @brief Get the index of i,j th element of lower triangular
         * matrix
         * 
         * @param row_index Row index i
         * @param column_index Column index j
         * @return Index of the i,j th element in the flattened array representation
         */
        unsigned int getIndex(
            unsigned int row_index,
            unsigned int column_index
        );
        
        /**
         * @brief Checks if the row index i and the column index j
         * belong to a zero element of the lower triangular matrix
         * 
         * @param row_index Row index i
         * @param column_index Column index j
         * @return true if i,j are indices of zero elements in a lower triangular matrix
         * @return false if i,j are indices of non zero elements in a lower triangular matrix
         */
        bool indexOfZeroElement(
            unsigned int row_index,
            unsigned int column_index
        );
};

#endif