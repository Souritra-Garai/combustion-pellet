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

template<typename real_t>
CoreShellDiffusion<real_t>::CoreShellDiffusion() : CoreShellCombustionParticle<real_t>()
{
    _concentration_array_A = new real_t[_n];
    _concentration_array_B = new real_t[_n];

    for (size_t i = 0; i < N; i++)
    {
        if (getRadialCoordinate(i) < this->_core_radius)
        {
            _concentration_array_A[i] = this->_core_material.getDensity() / this->_core_material.getMolecularWeight();
            _concentration_array_B[i] = 0;
        }
        else
        {
            _concentration_array_A[i] = 0;
            _concentration_array_B[i] = this->_shell_material.getDensity() / this->_shell_material.getMolecularWeight();
        }
    }
    
    // calcMassFractions();
}

template<typename real_t>
CoreShellDiffusion<real_t>::~CoreShellDiffusion()
{
    delete [] concentration_array_A;
    delete [] concentration_array_B;
}

template<typename real_t>
real_t CoreShellDiffusion<real_t>::getDiffusivity(real_t T)
{
    return pre_exponential_factor * exp(- activation_energy / (8.314 * T));
}

template<typename real_t>
real_t CoreShellDiffusion<real_t>::getRadialCoordinate(size_t i)
{
    return radius * ((real_t) i / (real_t) N);
}

template<typename real_t>
void CoreShellDiffusion<real_t>::calcMassFractions()
{
    real_t Y_A  = 0;
    real_t Y_B  = 0;
    real_t Y_AB = 0;
    
    #pragma omp parallel for reduction(+:Y_A, Y_B, Y_AB)
    
        for (size_t i = 0; i < N-1; i++)
        {
            Y_A += 
                std::max(_concentration_array_A[i]   - _concentration_array_B[i],   (real_t) 0.0) * pow(getRadialCoordinate(i), 2) +
                std::max(_concentration_array_A[i+1] - _concentration_array_B[i+1], (real_t) 0.0) * pow(getRadialCoordinate(i+1), 2);

            Y_B +=
                std::max(_concentration_array_B[i]   - _concentration_array_A[i],   (real_t) 0.0) *  pow(getRadialCoordinate(i), 2) +
                std::max(_concentration_array_B[i+1] - _concentration_array_A[i+1], (real_t) 0.0l) * pow(getRadialCoordinate(i+1), 2);

            Y_AB += 
                std::min(_concentration_array_A[i],   _concentration_array_B[i])   * pow(getRadialCoordinate(i), 2) +
                std::min(_concentration_array_A[i+1], _concentration_array_B[i+1]) * pow(getRadialCoordinate(i+1), 2);
        }
    
    this->_mass_fraction_core_material    = 0.5 * 4.0 * M_PI * _delta_r * Y_A  * this->_core_material.getMolecularWeight()    / this->_mass;
    this->_mass_fraction_shell_material   = 0.5 * 4.0 * M_PI * _delta_r * Y_B  * this->_shell_material.getMolecularWeight()   / this->_mass;
    this->_mass_fraction_product_material = 0.5 * 4.0 * M_PI * _delta_r * Y_AB * this->_product_material.getMolecularWeight() / this->_mass;
}

template class CoreShellDiffusion<float>;
template class CoreShellDiffusion<double>;
template class CoreShellDiffusion<long double>;