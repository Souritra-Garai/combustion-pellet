#include "lusolver/LU_Solver.hpp"

#include <iostream>

template<typename real_t>
LUSolver<real_t>::LUSolver(unsigned int n) : _n(n), _A(n)
{
	b = new real_t[_n];
}

template<typename real_t>
LUSolver<real_t>::~LUSolver()
{
	delete [] b;
}

template<typename real_t>
void LUSolver<real_t>::printMatrixEquation()
{
    std::cout << "Matrix Eqn A.x = b" << std::endl;
	
	std::cout << std::endl;
    
	std::cout << "A\t" << _n << " x " << _n << std::endl;    
	_A.printMatrix();

	std::cout << std::endl;

    std::cout << "b\t" << _n << " x 1" << std::endl;
    for (int k = 0; k < _n; k++) std::cout << b[k] << std::endl;
}

template class LUSolver<long double>;
template class LUSolver<double>;
template class LUSolver<float>;