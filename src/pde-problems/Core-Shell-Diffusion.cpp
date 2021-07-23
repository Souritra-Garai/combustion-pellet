/**
 * @file Core-Shell-Diffusion.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Implementation of CoreShellDiffusion class
 * @version 0.1
 * @date 2021-07-16
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "pde-problems/Core-Shell-Diffusion.hpp"

#include <math.h>

template<typename real_t> real_t CoreShellDiffusion<real_t>::_delta_t = 0;

template<typename real_t> size_t CoreShellDiffusion<real_t>::_n = 1;
template<typename real_t> real_t CoreShellDiffusion<real_t>::_delta_r = 0;

template<typename real_t> real_t CoreShellDiffusion<real_t>::_pre_exponential_factor = 0;
template<typename real_t> real_t CoreShellDiffusion<real_t>::_activation_energy = 0;

template<typename real_t>
void CoreShellDiffusion<real_t>::setTimeStep(real_t Dt) {_delta_t = Dt;}

template<typename real_t>
void CoreShellDiffusion<real_t>::setDiffusivityParameters(
    real_t pre_exponential_factor,
    real_t activaiton_energy
) {
    _pre_exponential_factor = pre_exponential_factor;
    _activation_energy = activaiton_energy;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::setGridSize(size_t n)
{
    _n = n;
    _delta_r = CoreShellCombustionParticle<real_t>::_overall_radius / (_n - 1);
}

template<typename real_t>
real_t CoreShellDiffusion<real_t>::getRadialCoordinate(size_t i)
{
    return CoreShellCombustionParticle<real_t>::_overall_radius * ((real_t) i / (real_t) (_n-1));
}

template<typename real_t>
real_t CoreShellDiffusion<real_t>::getDiffusivity(real_t T)
{
    return _pre_exponential_factor * exp(- _activation_energy / (8.314 * T));
}

template<typename real_t>
CoreShellDiffusion<real_t>::CoreShellDiffusion() : 
    CoreShellCombustionParticle<real_t>(),
    _solver_A(_n),
    _solver_B(_n)
{
    _concentration_array_A = new real_t[_n];
    _concentration_array_B = new real_t[_n];

    #pragma omp parallel for

        for (size_t i = 0; i < _n; i++)
        {
            if (getRadialCoordinate(i) < this->_core_radius)
            {
                _concentration_array_A[i] = this->_core_material.getMolarDensity();
                _concentration_array_B[i] = 0;
            }
            else
            {
                _concentration_array_A[i] = 0;
                _concentration_array_B[i] = this->_shell_material.getMolarDensity();
            }
        }
}

template<typename real_t>
CoreShellDiffusion<real_t>::~CoreShellDiffusion()
{
    delete [] _concentration_array_A;
    delete [] _concentration_array_B;
}

template<typename real_t>
inline real_t CoreShellDiffusion<real_t>::getRxnConcA(size_t i)
{
    return std::max(_concentration_array_A[i] - _concentration_array_B[i], (real_t) 0.0);
}

template<typename real_t>
inline real_t CoreShellDiffusion<real_t>::getRxnConcB(size_t i)
{
    return std::max(_concentration_array_B[i] - _concentration_array_A[i], (real_t) 0.0);
}

template<typename real_t>
inline real_t CoreShellDiffusion<real_t>::getRxnConcAB(size_t i)
{
    return std::min(_concentration_array_A[i], _concentration_array_B[i]);
}

template<typename real_t>
void CoreShellDiffusion<real_t>::calcMassFractions()
{
    real_t Y_A  = 0.5 * getRxnConcA(_n-1)  * pow(this->_overall_radius, 2);
    real_t Y_B  = 0.5 * getRxnConcB(_n-1)  * pow(this->_overall_radius, 2);
    real_t Y_AB = 0.5 * getRxnConcAB(_n-1) * pow(this->_overall_radius, 2);
    
    #pragma omp parallel for reduction(+:Y_A, Y_B, Y_AB)
    
        for (size_t i = 1; i < _n-1; i++)
        {
            Y_A  += getRxnConcA(i)  * pow(getRadialCoordinate(i), 2);
            Y_B  += getRxnConcB(i)  * pow(getRadialCoordinate(i), 2);
            Y_AB += getRxnConcAB(i) * pow(getRadialCoordinate(i), 2);
        }
    
    this->_mass_fraction_core_material    = 4.0 * M_PI * _delta_r * Y_A  * this->_core_material.getMolecularWeight()    / this->_mass;
    this->_mass_fraction_shell_material   = 4.0 * M_PI * _delta_r * Y_B  * this->_shell_material.getMolecularWeight()   / this->_mass;
    this->_mass_fraction_product_material = 4.0 * M_PI * _delta_r * Y_AB * this->_product_material.getMolecularWeight() / this->_mass;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::setUpEquations(real_t T)
{
    real_t coefficient1 = getDiffusivity(T) / pow(_delta_r, 2);
    real_t coefficient2 = 1.0 / _delta_t;

    _solver_A.setEquationFirstRow(1, -1, 0);
    _solver_B.setEquationFirstRow(1, -1, 0);

    #pragma omp parallel for
    
        for (size_t i = 1; i < _n - 1; i++)
        {
            real_t coefficient3 = coefficient1 * pow(getRadialCoordinate(i+1) / getRadialCoordinate(i), 2);

            real_t e = - coefficient3;
            real_t f = coefficient1 + coefficient2 + coefficient3;
            real_t g = - coefficient1;

            _solver_A.setEquation(i, e, f, g, coefficient2 * _concentration_array_A[i]);
            _solver_B.setEquation(i, e, f, g, coefficient2 * _concentration_array_B[i]);
        }

    _solver_A.setEquationLastRow(-1, 1, 0);
    _solver_B.setEquationLastRow(-1, 1, 0);
}

template<typename real_t>
void CoreShellDiffusion<real_t>::solveEquations()
{
    _solver_A.getSolution(_concentration_array_A);
    _solver_B.getSolution(_concentration_array_B);   
}

template class CoreShellDiffusion<float>;
template class CoreShellDiffusion<double>;
template class CoreShellDiffusion<long double>;