#include "lusolver/LU-Solver.hpp"

#include <iostream>

LUSolver::LUSolver(unsigned int n) : _n(n), _A(n)
{
	_b = new real_t[_n];
}

LUSolver::~LUSolver()
{
	delete [] _b;
}

void LUSolver::printMatrixEquation()
{
	std::cout << "Matrix Eqn A.x = b" << std::endl;
	
	std::cout << std::endl;
	
	std::cout << "A\t" << _n << " x " << _n << std::endl;    
	_A.printMatrix();

	std::cout << std::endl;

	std::cout << "b\t" << _n << " x 1" << std::endl;
	for (int k = 0; k < _n; k++) std::cout << _b[k] << std::endl;
}
