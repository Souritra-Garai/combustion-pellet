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
#include "lusolver/LU_Solver.hpp"

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
        LUSolver<real_t> _solver_A;
        /**
         * @brief Tridiagonal matrix equation solver for concentration
         * profile of substance B
         */
        LUSolver<real_t> _solver_B;

        /**
         * @brief Get the Radial Coordinate of the grid point 
         * at the specified index
         * @param index Index of the grid point
         * @return real_t Radial distance of the grid point from origin in \f$ m \f$
         */
        static real_t getRadialCoordinate(size_t index);
        
        /**
         * @brief Calculate reaction mass fractions of the substances A, B and AB
         * in the particle and set the mass fraction variables to the calculated values
         */
        void calcRxnMassFractions();

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
         * @brief Print the PDE solver configuration to the output stream
         * 
         * @param output_stream Stream to which configuration is printed
         */
        static void printConfiguration(std::ostream &output_stream);

        /**
         * @brief Initialize the concentration of the species in the core and shell
         * regions of the particle
         */
        void initializeParticle();

        /**
         * @brief Set up the matrix equation for diffusion over one time step
         * @param diffusivity Value of diffusivity in \f$ m^2 / s \f$
         */
        void setUpEquations(real_t diffusivity);
		/**
         * @brief Set up the matrix equation for diffusion over one time step
         * @param diffusivity Value of diffusivity in \f$ m^2 / s \f$
         */
        void setUpEquations(real_t diffusivity, CoreShellDiffusion<real_t> &diffusion_problem, real_t Delta_t = _delta_t);

        /**
         * @brief Solve the matrix equations and update the particle state internally
         */
        void solveEquations();

        /**
         * @brief Get the Diffusion Mass of substance A (initially present in core)
         * 
         * @return real_t Numerically evaluate and return the diffusion mass 
         * of the substance B in \f$ kg \f$
         */
        real_t getDiffusionMassA();

        /**
         * @brief Get the Diffusion Mass of substatance B (initially present in shell)
         * 
         * @return real_t Numerically evaluate and return the diffusion mass 
         * of the substance B in \f$ kg \f$
         */
        real_t getDiffusionMassB();

        /**
         * @brief Copy the state of diffusion and reaction from the argument
         * 
         * @param diffusion_problem Copy from the diffusion_problem object to the calling object
         */
        void copyFrom(CoreShellDiffusion<real_t> & diffusion_problem);

        /**
         * @brief Copy the state of diffusion and reaction to the argument
         * 
         * @param diffusion_problem Copy to the diffusion_problem object from the calling object
         */
        void copyTo(CoreShellDiffusion<real_t> & diffusion_problem);

        /**
         * @brief Print the diffusion concentration profile of substance A to the output stream
         * @param output_stream ostream object
         * @param delimiter Delimiter to separate consecutive values
         */
        void printConcentrationProfileA(std::ostream &output_stream, char delimiter = '\t', real_t curr_time = 0);
        /**
         * @brief Print the diffusion concentration profile of substance B to the output stream
         * @param output_stream ostream object
         * @param delimiter Delimiter to separate consecutive values
         */
        void printConcentrationProfileB(std::ostream &output_stream, char delimiter = '\t', real_t curr_time = 0);

        /**
         * @brief Prints the r-coordinates of the grid points in the discretized spherical particle
         * 
         * @param output_stream Output stream where the coordinates will be printed
         * @param delimiter Character to separate two consecutive values of temperature
         */
        void printGridPoints(std::ostream &output_stream, char delimiter = '\t');
};

#endif