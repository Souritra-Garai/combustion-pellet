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
void LUSolver<real_t>::setEquation(
	unsigned int i,
	real_t e_val,
    real_t f_val,
    real_t g_val,
    real_t b_val
) {
	// \f$ e x_{i-1} + f x_i + g x_{i+1} = b \f$

	// Set \f$ a_{i,i-1} \f$ to e_val
    _A.setElement(i, i-1, e_val);
    // Set \f$ a_{i,i} \f$ to f_val
    _A.setElement(i, i,   f_val);
    // Set \f$ a_{i,i+1} \f$ to g_val
    _A.setElement(i, i+1, g_val);

    // Set the constant
    b[i] = b_val;
}

template<typename real_t>
void LUSolver<real_t>::setEquationFirstRow(
    real_t f_val,
    real_t g_val,
    real_t b_val
) {
    // \f$ f x_i + g x_{i+1} = b \f$

    // Set \f$ a_{i,i} \f$ to f_val
    _A.setElement(0, 0, f_val);
    // Set \f$ a_{i,i+1} \f$ to g_val
    _A.setElement(0, 1, g_val);

    // Set the constant
    b[0] = b_val;
}

template<typename real_t>
void LUSolver<real_t>::setEquationLastRow(
    real_t e_val,
    real_t f_val,
    real_t b_val
) {
    // \f$ e x_{i-1} f x_i = b \f$

    // Set \f$ a_{i,i-1} \f$ to e_val
    _A.setElement(_n-1, _n-2, e_val);
    // Set \f$ a_{i,i} \f$ to f_val
    _A.setElement(_n-1, _n-1, f_val);

    // Set the constant
    b[_n-1] = b_val;
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

template<typename real_t>
void LUSolver<real_t>::LU_DecompositionAndForwardSubstitution()
{
	real_t Lower_Matrix_Diagonal_less_1;
    
    for (int i=1; i<_n; i++)
    {
		// Set l_{i,i-1} = a_{i,i-1} / u_{i-1,i-1}
        Lower_Matrix_Diagonal_less_1 = _A.getElement(i, i-1) / _A.getElement(i-1, i-1);

		// Set u_{i,i} = a_{i,i} - a_{i-1,i} * l_(i,i-1)
        _A.setElement(i, i,
			_A.getElement(i, i) - _A.getElement(i-1, i) * Lower_Matrix_Diagonal_less_1
		);
		
		// Set d_{i} = b_{i} - d_{i-1} * l_{i,i-1}
		b[i] -= b[i-1] * Lower_Matrix_Diagonal_less_1;
    }
}

template<typename real_t>
void LUSolver<real_t>::getSolution(real_t *x)
{
	LU_DecompositionAndForwardSubstitution();

	// Backward substitution
	// Last element is simply x_{n-1} = d_{n-1} / U_{n-1, n-1}
	x[_n-1] = b[_n-1] / _A.getElement(_n-1, _n-1);

    for (int i=_n-2; i>=0; i--)

		// Set x_{i} = ( d_{i} - U_{i,i+1} * x_{i+1} ) / U_{i,i}
        x[i] = ( b[i] - _A.getElement(i, i+1) * x[i+1] ) / _A.getElement(i, i);
}

template class LUSolver<long double>;
template class LUSolver<double>;
template class LUSolver<float>;