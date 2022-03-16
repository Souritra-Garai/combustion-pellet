#include "lusolver/Tridiagonal_Matrix.hpp"

#include <iostream>

template<typename real_t>
TridiagonalMatrix<real_t>::TridiagonalMatrix(unsigned int n)
{
	_n = n;
	
	_array = new real_t[3*n - 2];
}

template<typename real_t>
TridiagonalMatrix<real_t>::~TridiagonalMatrix()
{
	delete [] _array;
}

template<typename real_t>
void TridiagonalMatrix<real_t>::printMatrix()
{
	unsigned int i, j;

	
	std::cout << getElement(0,0) << '\t' << getElement(0,1);

	for (i = 2; i < _n; i++) std::cout << "\t0.0";

	std::cout << std::endl;

	
	for (i = 1; i < _n-1; i++)
	{
		for (j = 0; j < i-1; j++) std::cout << "0.0\t";

		std::cout << getElement(i, i-1) << '\t' << getElement(i, i) << '\t' << getElement(i, i+1);

		for (j = i+2; j < _n; j++) std::cout << "\t0.0";

		std::cout << std::endl;
	}


	for (i = 0; i < _n-2; i++) std::cout << "0.0\t";

	std::cout << getElement(_n-1, _n-2) << '\t' << getElement(_n-1, _n-1);

	std::cout << std::endl;
}

template<typename real_t>
void TridiagonalMatrix<real_t>::multiply(real_t *x, real_t *b)
{
	b[0] = getElement(0, 0) * x[0] + getElement(0, 1) * x[1];

	for (unsigned int i = 1; i < _n-1; i++)
	
		b[i] = getElement(i, i-1) * x[i-1] + getElement(i, i) * x[i] + getElement(i, i+1) * x[i+1];

	b[_n-1] = getElement(_n-1, _n-2) * x[_n-2] + getElement(_n-1, _n-1) * x[_n-1];
}

template class TridiagonalMatrix<long double>;
template class TridiagonalMatrix<double>;
template class TridiagonalMatrix<float>;