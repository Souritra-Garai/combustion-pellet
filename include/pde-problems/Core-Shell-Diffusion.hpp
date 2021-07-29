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

// Required for definition of CoreShellCombustionParticle class
#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"
// Required for definition of QRSolver class
#include "qrsolver/QR_Solver.hpp"

/**
 * @brief Class with functionalities to perform diffusion in a
 * spherical core shell particle
 * 
 * The class is derived from CoreShellCombustionParticle class
 * to retain the thrmo-physical properties of the core-shell particle
 * @tparam real_t float, double or long double data types
 * to represent real numbers
 */
template<typename real_t>
class CoreShellDiffusion : public CoreShellCombustionParticle<real_t>
{
    private :

        /**
         * @brief Duration of interval between consecutive time steps
         */
        static real_t _delta_t;

        /**
         * @brief Number of grid points along the radial axis
         * includeing the first point at \f$ r = 0 \f$ and last point
         * at \f$ r = r_P \f$
         */
        static size_t _n;
        /**
         * @brief Distance between consecutive grid points
         */
        static real_t _delta_r;

        /**
         * @brief Array to store molar concentration of substance A in \f$ mol/m^3 \f$,
         * initially present in the core, at each grid point
         */
        real_t * _concentration_array_A;
        /**
         * @brief Array to store molar concentration of substance B in \f$ mol/m^3 \f$,
         * initially present in the shell, at each grid point
         */
        real_t * _concentration_array_B;

        /**
         * @brief Tridiagonal matrix equation solver for concentration
         * profile of substance A
         */
        QRSolver<real_t> _solver_A;
        /**
         * @brief Tridiagonal matrix equation solver for concentration
         * profile of substance B
         */
        QRSolver<real_t> _solver_B;

        /**
         * @brief Get the Radial Coordinate of the grid point 
         * at the specified index
         * @param index Index of the grid point
         * @return real_t Radial distance of the grid point from origin in \f$ m \f$
         */
        static real_t getRadialCoordinate(size_t index);
        
        /**
         * @brief Calculate mass fractions of the substances A, B and AB
         * in the particle and set the mass fraction varibles to the calculated values
         */
        void calcMassFractions();

        /**
         * @brief Get the reaction concentration of substance A
         * 
         * Reaction concentration assumes all the amount of the limiting agent
         * has been used up to form the product AB
         * @param index Index of grid point where the concentration of A is desired
         * @return real_t Molar concentration of substance A in \f$ mol/m^3 \f$
         */
        real_t getRxnConcA  (size_t index);
        /**
         * @brief Get the reaction concentration of substance B
         * 
         * Reaction concentration assumes all the amount of the limiting agent
         * has been used up to form the product AB
         * @param index Index of grid point where the concentration of B is desired
         * @return real_t Molar concentration of substance B in \f$ mol/m^3 \f$
         */
        real_t getRxnConcB  (size_t index);
        /**
         * @brief Get the reaction concentration of substance AB
         * 
         * Reaction concentration assumes all the amount of the limiting agent
         * has been used up to form the product AB
         * @param index Index of grid point where the concentration of AB is desired
         * @return real_t Molar concentration of substance AB in \f$ mol/m^3 \f$
         */
        real_t getRxnConcAB (size_t index);

    public :

        /**
         * @brief Construct a new Core Shell Diffusion object
         */
        CoreShellDiffusion();
        /**
         * @brief Destroy the Core Shell Diffusion object
         */
        ~CoreShellDiffusion();

        /**
         * @brief Set the Grid Size
         * @param n Number of grid points to consider,
         * including the points \f$ r = 0 \f$ and \f$ r = r_P \f$
         */
        static void setGridSize(size_t n);
        /**
         * @brief Set the duration of time step interval
         * @param Delta_t Duration of time step intervals in \f$ s \f$
         */
        static void setTimeStep(real_t Delta_t);

        /**
         * @brief Set up the matrix equation for diffusion over one time step
         * @param diffusivity Value of diffusivity in \f$ m^2 / s \f$
         */
        void setUpEquations(real_t diffusivity);
        /**
         * @brief Solve the matrix equations and upadte the particle state internally
         */
        void solveEquations();

        /**
         * @brief Print the diffusion concentration profile of substance A to the output stream
         * @param output_stream ostream object
         * @param delimiter Delimiter to separate consecutive values
         */
        void printConcentrationProfileA(std::ostream &output_stream, char delimiter = '\t');
        /**
         * @brief Print the diffusion concentration profile of substance B to the output stream
         * @param output_stream ostream object
         * @param delimiter Delimiter to separate consecutive values
         */
        void printConcentrationProfileB(std::ostream &output_stream, char delimiter = '\t');
        
        /**
         * @brief Numerically evaluate the diffusion mass of the substance A
         * initially at the core of the particle
         * @return real_t Mass of the substance initially at the core, at the
         * current state in \f$ kg \f$
         */
        real_t numCalcCoreMass();
        /**
         * @brief Numerically evaluate the diffusion mass of the substance B
         * initially in the shell
         * @return real_t Mass of the substance initially in the shell, at the
         * current state in \f$ kg \f$
         */
        real_t numCalcShellMass();
};

#endif