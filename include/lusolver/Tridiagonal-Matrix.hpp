#ifndef __TRIDIAGONAL_MATRIX__
#define __TRIDIAGONAL_MATRIX__

#include "math/Data-Type.hpp"

class TridiagonalMatrix
{
	private:

		const unsigned int _n;

		real_t * _array;

		// i - row index
		// j - column index
		// i and j start from 0 to _n-1
		static inline unsigned int getIndex(
			unsigned int i,
			unsigned int j
		) {
			// if ( (j < (i-1)) || (j > (i+1)) ) throw std::out_of_range("Indices cannot be j < i-1 or j > i+1 for tridiagonal matrices.");
			return 2 * i + j;
		}
		
	public:

		// i - row index
		// j - column index
		// i and j start from 0 to _n-1
		inline real_t operator() (unsigned int i, unsigned int j) const
		{
			return _array[getIndex(i, j)];
		}

		// i - row index
		// j - column index
		// i and j start from 0 to _n-1
		inline real_t& operator() (unsigned int i, unsigned int j)
		{
			return _array[getIndex(i, j)];
		}
		
		// n - size of the square matrix
		TridiagonalMatrix(unsigned int n);
		
		~TridiagonalMatrix();

		// i - row index
		// j - column index
		// i and j start from 0 to _n-1
		inline real_t getElement(
			unsigned int i,
			unsigned int j
		) const 
		{
			return _array[getIndex(i, j)];
		}

		// i - row index
		// j - column index
		// i and j start from 0 to _n-1
		inline void setElement(
			unsigned int i,
			unsigned int j,
			real_t value
		) {
			_array[getIndex(i, j)] = value;
		}

		// Prints the Tridiagonal Matrix in the form of a 2D array
		void printMatrix();
		
		// Multiplies the Tridiagonal matrix with the column vector b
		// and stores the result in the column vector x
		// both arrays should be of size n
		void multiply(real_t *b, real_t *x) const;
};

#endif