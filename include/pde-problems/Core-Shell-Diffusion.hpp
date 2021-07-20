/**
 * @file Core-Shell-Diffusion.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Definition for Core-Shell Diffusion problem class
 * @version 0.1
 * @date 2021-07-16
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __CORE_SHELL_DIFFUSION__
#define __CORE_SHELL_DIFFUSION__

#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"
#include "qrsolver/QR_Solver.hpp"

template<typename real_t>
class CoreShellDiffusion : public CoreShellCombustionParticle<real_t>
{
    private :

        const size_t  N;

        real_t * concentration_array_A;
        real_t * concentration_array_B;

        real_t radius;
        real_t Delta_r;
        real_t mass;

        /**
         * @brief Pre exponential factor for 
         * Arrhenius Diffusivity model
         */
        real_t pre_exponential_factor;
        /**
         * @brief Activation energy for 
         * Arrhenius Diffusivity model
         */
        real_t activation_energy;

        real_t getRadialCoordinate(size_t index);

        void calcMassFractions();

    public :

        CoreShellDiffusion(
            CoreShellCombustionParticle<real_t> combustion_particle,
            size_t number_of_grid_points
        );

        ~CoreShellDiffusion();

        /**
         * @brief Get the Diffusivity for the interdiffusion of the core and shell
         * material
         * 
         * @param temperature Overall temperature of the particle
         * @return real_t Diffusivity for interdiffusion model at the specified temperature
         */
        real_t getDiffusivity(real_t temperature);
};

#endif