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

#include <iostream>

template<typename real_t>
CoreShellDiffusion<real_t>::CoreShellDiffusion(
    CoreShellCombustionParticle<real_t> combustion_particle,
    size_t number_of_grid_points
) : CoreShellCombustionParticle<real_t>(combustion_particle),
    N(number_of_grid_points)
{
    concentration_array_A = new real_t[N];
    concentration_array_B = new real_t[N];

    Delta_r = (*(this->overall_radius)) / (N-1);

    for (size_t i = 0; i < N; i++)
    {
        if (getRadialCoordinate(i) <= (*(this->core_radius)))
        {
            concentration_array_A[i] = this->reactant_A->getDensity() / this->reactant_A->getMolecularWeight();
            concentration_array_B[i] = 0;
        }
        else
        {
            concentration_array_A[i] = 0;
            concentration_array_B[i] = this->reactant_B->getDensity() / this->reactant_B->getMolecularWeight();
        }
    }
    
    calcMassFractions();
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
                std::max(concentration_array_A[i] - concentration_array_B[i], (real_t) 0.0) * pow(getRadialCoordinate(i), 2) +
                std::max(concentration_array_A[i+1] - concentration_array_B[i+1], (real_t) 0.0) * pow(getRadialCoordinate(i+1), 2);

            Y_B +=
                std::max(concentration_array_B[i] - concentration_array_A[i], (real_t) 0.0) * pow(getRadialCoordinate(i), 2) +
                std::max(concentration_array_B[i+1] - concentration_array_A[i+1], (real_t) 0.0l) * pow(getRadialCoordinate(i+1), 2);

            Y_AB += 
                std::min(concentration_array_A[i], concentration_array_B[i]) * pow(getRadialCoordinate(i), 2) +
                std::min(concentration_array_A[i+1], concentration_array_B[i+1]) * pow(getRadialCoordinate(i+1), 2);
        }

    std::cout << this->mass;
    
    this->mass_fraction_reactant_A = 0.5 * 4.0 * M_PI * Delta_r * Y_A  * this->reactant_A->getMolecularWeight() / this->mass;
    this->mass_fraction_reactant_B = 0.5 * 4.0 * M_PI * Delta_r * Y_B  * this->reactant_B->getMolecularWeight() / this->mass;
    this->mass_fraction_product_AB = 0.5 * 4.0 * M_PI * Delta_r * Y_AB * this->product_AB->getMolecularWeight() / this->mass;
}

template class CoreShellDiffusion<float>;
template class CoreShellDiffusion<double>;
template class CoreShellDiffusion<long double>;