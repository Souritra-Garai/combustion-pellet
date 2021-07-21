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

        static size_t _n;
        static real_t _delta_r;

        /**
         * @brief Pre exponential factor for 
         * Arrhenius Diffusivity model
         */
        static real_t _pre_exponential_factor;
        /**
         * @brief Activation energy for 
         * Arrhenius Diffusivity model
         */
        static real_t _activation_energy;

        real_t * _concentration_array_A;
        real_t * _concentration_array_B;

        QRSolver<real_t> _solver_A;
        QRSolver<real_t> _solver_B;

    public :

        CoreShellDiffusion();
        ~CoreShellDiffusion();

        static void setGridSize(size_t n);
        static void setDiffusivityParameters(
            real_t pre_exponential_constant,
            real_t activation_energy
        );

        /**
         * @brief Get the Diffusivity for the interdiffusion of the core and shell
         * material
         * 
         * @param temperature Overall temperature of the particle
         * @return real_t Diffusivity for interdiffusion model at the specified temperature
         */
        static real_t getDiffusivity(real_t temperature);

        static real_t getRadialCoordinate(size_t index);
        
        void calcMassFractions();
};

#endif