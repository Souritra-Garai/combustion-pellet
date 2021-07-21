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
CoreShellDiffusion<real_t>::CoreShellDiffusion(
    CoreShellCombustionParticle<real_t> combustion_particle,
    size_t number_of_grid_points,
    real_t particle_radius,
    real_t core_radius
) : CoreShellCombustionParticle<real_t>(combustion_particle),
    N(number_of_grid_points)
{
    concentration_array_A = new real_t[N];
    concentration_array_B = new real_t[N];

    radius = particle_radius;

    Delta_r = radius / (N-1);
    
    mass = (4.0 * M_PI / 3.0) * (
        this->_core_material.getDensity() * pow(core_radius, 3) +
        this->_shell_material.getDensity() * (pow(particle_radius, 3) - pow(core_radius, 3))
    );

    for (size_t i = 0; i < N; i++)
    {
        if (getRadialCoordinate(i) <= core_radius)
        {
            concentration_array_A[i] = this->_core_material.getDensity() / this->_core_material.getMolecularWeight();
            concentration_array_B[i] = 0;
        }
        else
        {
            concentration_array_A[i] = 0;
            concentration_array_B[i] = this->_shell_material.getDensity() / this->_shell_material.getMolecularWeight();
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
    
    this->_mass_fraction_core_material =    0.5 * 4.0 * M_PI * Delta_r * Y_A  * this->_core_material.getMolecularWeight()    / mass;
    this->_mass_fraction_shell_material =   0.5 * 4.0 * M_PI * Delta_r * Y_B  * this->_shell_material.getMolecularWeight()   / mass;
    this->_mass_fraction_product_material = 0.5 * 4.0 * M_PI * Delta_r * Y_AB * this->_product_material.getMolecularWeight() / mass;
}

template class CoreShellDiffusion<float>;
template class CoreShellDiffusion<double>;
template class CoreShellDiffusion<long double>;