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

        static real_t _delta_t;

        static size_t _n;
        static real_t _delta_r;

        static real_t _pre_exponential_factor;
        static real_t _activation_energy;

        real_t * _concentration_array_A;
        real_t * _concentration_array_B;

        QRSolver<real_t> _solver_A;
        QRSolver<real_t> _solver_B;

        static real_t getDiffusivity(real_t temperature);

        static real_t getRadialCoordinate(size_t index);
        
        void calcMassFractions();

        real_t getRxnConcA  (size_t index);
        real_t getRxnConcB  (size_t index);
        real_t getRxnConcAB (size_t index);

    public :

        CoreShellDiffusion();
        ~CoreShellDiffusion();

        static void setGridSize(size_t n);
        static void setTimeStep(real_t Delta_t);

        static void setDiffusivityParameters(
            real_t pre_exponential_factor,
            real_t activation_energy
        );

        void setUpEquations(real_t temperature);
        void solveEquations();
};

#endif