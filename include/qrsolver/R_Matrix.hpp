/**
 * @file R_Matrix.hpp
 * @author Souritra Gari (souritra.garai@iitgn.ac.in)
 * @brief This header file defines a class for memory efficient 
 * implementation of \f$ R \f$ matrix used for QR factorisation of tridiagonal matrix
 * @version 0.1
 * @date 2021-06-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __TRIDIAGONAL_MATRIX__
#define __TRIDIAGONAL_MATRIX__

/**
 * @brief Class to implement a memory efficient
 * \f$ N \times N \f$ R matrix
 * 
 * The class is specifically built for implementation in a
 * QR factorization algorithm for solving matrix equations of the form
 * \f$ A \cdot x = b \f$. The QR algo converts a normal tridiagonal matrix
 * (a matrix with non zero entries only at indices \f$ (i,i-1) \f$, \f$ (i,i) \f$
 * and \f$ (i,i+1) \f$) to an upper tridiagonal matrix (a matrix with 
 * non zero entries only at indices \f$ (i,i) \f$, \f$ (i,i+1) \f$ and \f$ (i,i+2) \f$).
 * Thus only indices \f$ (i,i-1) \f$, \f$ (i,i) \f$, \f$ (i,i+1) \f$ and
 * \f$ (i,i+2) \f$ are stored in memory for this matrix
 */
template<typename real_t>
class RMatrix
{
    public :

        /**
         * @brief Construct a new R Matrix
         * 
         * @param n Size of main diagonal of the 
         * \f$ N \times N \f$ R Matrix
         */
        RMatrix(unsigned int n);
        
        /**
         * @brief Destroy the R Matrix
         * 
         */
        ~RMatrix();

        /**
         * @brief Get the i,j th element of R Matrix
         * 
         * @param row_index Row index i
         * @param column_index Column index j
         * @return Value of the i,j th element of a R Matrix
         */
        real_t getElement(
            unsigned int row_index,
            unsigned int column_index
        );

        /**
         * @brief Set the value of the i,j th element of R Matrix
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
         * @brief Prints the R Matrix in form of a 2D array
         * 
         */
        void printMatrix();

        /**
         * @brief Prints the R Matrix in form of a flattened array
         * 
         */
        void print();

        
    private :
        
        /**
         * @brief Flattened array of size 4 * N to represent
         * R matrix of size \f$N \times N\f$
         * 
         */
        real_t *array;

        /**
         * @brief Size of main diagonal of the 
         * \f$N \times N\f$ R Matrix
         * 
         */
        const unsigned int N;

        /**
         * @brief Get the index of the i,j th element of R Matrix in
         * the flattened array
         * 
         * @param row_index Row index i
         * @param column_index Column index j
         * @return Returns the index in the the flattened array
         */
        unsigned int getIndex(
            unsigned int row_index,
            unsigned int column_index
        );

        /**
         * @brief Checks if the row index i and the column index j
         * belong to a zero element of the R matrix
         * 
         * @param row_index 
         * @param column_index 
         * @return true if i,j are indices of zero elements in a R matrix
         * @return false if i,j are indicess of non-zero elements in a R matrix
         */
        bool indexOfZeroElement(
            unsigned int row_index,
            unsigned int column_index
        );
};


#endif