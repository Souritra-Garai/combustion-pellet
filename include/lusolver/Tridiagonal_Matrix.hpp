#ifndef __TRIDIAGONAL_MATRIX__
#define __TRIDIAGONAL_MATRIX__

template<typename real_t>
class TridiagonalMatrix
{
	private:

		const unsigned int _n;

		real_t * _array;

		// row_index and column_index start from 0 to _n-1
		inline unsigned int getIndex(
			unsigned int row_index,
			unsigned int column_index
		) {
			return 2 * row_index + column_index;
		}
		
	public:
		
		// n - size of the square matrix
		TridiagonalMatrix(unsigned int n);
		
		~TridiagonalMatrix();

		// row_index and column_index start from 0 to _n-1
		real_t getElement(
			unsigned int row_index,
			unsigned int column_index
		) {
			return _array[getIndex(row_index, column_index)];
		}

		// row_index and column_index start from 0 to _n-1
		void setElement(
			unsigned int row_index,
			unsigned int column_index,
			real_t value
		) {
			_array[getIndex(row_index, column_index)] = value;
		}

		// Prints the Tridiagonal Matrix in form of a 2D array
		void printMatrix();
		
		// Multiplies the Tridiagonal matrix with the column vector b
        // and stores the result in the column vector x
		// both arrays should be of size n
        void multiply(real_t *b, real_t *x);
};

#endif