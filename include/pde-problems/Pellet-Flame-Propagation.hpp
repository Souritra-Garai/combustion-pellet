/**
 * @file Pellet-Flame-Propagation.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Definition of class representing energy transport problem in
 * a packed pellet of core shell particles
 * @version 0.1
 * @date 2021-08-02
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __PELLET_FLAME_PROPAGATION__
#define __PELLET_FLAME_PROPAGATION__

// Required for output stream class 
// (to print temperature / other properties to screen / file)
#include <ostream>

// Required for PackedPellet class representing thermo-physical properties of 
// a pellet packed with energetic particles and degassed with an inert
#include "thermo-physical-properties/Packed-Pellet.hpp"
// Required for CoreShellDiffusion class representing the 
// diffusion in a core-shell particle
#include "pde-problems/Core-Shell-Diffusion.hpp"
// Required for QRSolver class used to set up and solve matrix equations
#include "lusolver/LU_Solver.hpp"

/**
 * @brief Class to represent a linear expression of the form
 * \f$ a_0 + a_1 \cdot x \f$
 * @tparam real_t float, double or long double data types
 * to represent real numbers
 */
template<typename real_t>
class LinearExpression {
    
    public :
        /**
         * @brief Constant \f$ a_0 \f$ in the linear expression
         */
        real_t constant;
        /**
         * @brief Coefficient of \f$ x \f$, \f$ a_1 \f$ in
         * the linear expression
         */
        real_t coefficient;
};

/**
 * @brief Class to represent combustion front propagation problem in
 * a pellet packed with energetic particles and degassed with an inert
 * @tparam real_t float, double or long double data types
 * to represent real numbers
 */
template<typename real_t>
class PelletFlamePropagation : public PackedPellet<real_t>
{
    private:

        /**
         * @brief Stores the time evolved since the zeroth iteration
         */
        real_t _time;

        /**
         * @brief Array to store temperature at each grid point of the pellet
         */
        real_t * _temperature_array;

		real_t * _thermal_conductivity;

		real_t * _prev_enthalpy_particle;
        
        /**
         * @brief Array of Objects to solve the core-shell diffusion problem at
         * each grid point in the pellet. It stores the state of the particles
         * at the current time step and is updated once the temperature is updated
         */
        CoreShellDiffusion<real_t> * _particles_array;

        /**
         * @brief Array of Objects to solve the core-shell diffusion problem and
         * store the state of the pellet at the previous time step, for each 
         * grid point in the pellet
         */
        CoreShellDiffusion<real_t> * _particles_array_const_temperature_evolution;
        /**
         * @brief Objects to solve the core-shell diffusion problem and
         * store the state of the particles at the current time step, for each
         * grid point in the pellet
         */
        CoreShellDiffusion<real_t> * _particles_array_raised_temperature_evolution;

        /**
         * @brief Matrix equation solver for solving the implicit temperature
         * update equation
         */
        LUSolver<real_t> _solver;

        /**
         * @brief Get the x-coordinate of the grid point
         * 
         * @param index Index of the grid point whose x-coordinate is desired
         * @return real_t x-coordiante of grid point # index
         */
        static real_t getXCoordinate(size_t index);

        /**
         * @brief Determine if the grid point is inside reaction zone or not
         * 
         * @param index Index of the grid point
         * @return true If the grid point is within the reaction zone
         * @return false If the grid point is out of the reaction zone
         */
        bool inReactionZone(size_t index);

        /**
         * @brief Evolve the constant temperature evolution particle and the
         * raised temperature evolution particle for determining the time and
         * temperature derivatives of enthalpy
         * 
         * @param index Index of the grid point where particle needs to be evolved
         */
        void evolveParticleForEnthalpyDerivative(size_t index);

        /**
         * @brief Calculates the linear expression in terms of Temperature
         * for the transient term in the energy transport equation for the
         * pellet
         * 
         * @param index Index of grid point where the linearized transient term is desired
         * @return LinearExpression<real_t> Linear expression in terms of temperature 
         * to describe the discretized transient term at the present grid point
         */
        LinearExpression<real_t> calcTransientTerm(size_t index);
        /**
         * @brief Calculates the linear expression in terms of Temperature
         * for the heat loss term in the energy transport equation for the
         * pellet
         * 
         * @param index Index of grid point where the linearized transient term is desired
         * @return LinearExpression<real_t> Linear expression in terms of temperature 
         * to describe the discretized heat loss term at the present grid point
         */
        LinearExpression<real_t> calcHeatLossTerm(size_t index);

        /**
         * @brief Set up the boundary condition at grid point # 0
         */
        void setUpBoundaryConditionX0();
        /**
         * @brief Set up the boundary condition at grid point # M
         */
        void setUpBoundaryConditionXN();

        /**
         * @brief Update the state of the particles using
         * the present temperature profile
         */
        void updateParticlesState();

		real_t getInterstitialGasTransientTermCoefficient(size_t index);

    public:
	
		static const real_t kappa;

		static const real_t gamma;
    
        /**
         * @brief Number of grid points in the pellet
         * including both points on the boundaries
         */
        static const size_t m;
        /**
         * @brief Distance between consecutive grid points
         * on the x axis, \f$ \Delta x \f$
         */
        static const real_t delta_x;

        /**
         * @brief Duration of time step, \f$ \Delta t \f$
         */
        static const real_t delta_t;

        /**
         * @brief Infinitesimal change in temperature used to
         * find derivative of a function with respect to temperature
         * using first principle
         */
        static const real_t delta_T;

        /**
         * @brief Construct a new Pellet Flame Propagation object
         * 
         * @param particle_volume_fraction Volume fraction of the pellet
         * occupied by energetic particles
         */
        PelletFlamePropagation(real_t particle_volume_fraction);
        /**
         * @brief Destroy the Pellet Flame Propagation object
         */
        ~PelletFlamePropagation();

        /**
         * @brief Initialize the pellet by setting the temperature profile and energetic
         * particles to the initial conditions
         */
        void initializePellet(
			real_t initial_ignition_temperature = 1500.,
			real_t initial_ignition_length_fraction = 0.1
		);

        /**
         * @brief Set up the discretized energy transport equations
         * for solving the temperature profile in the pellet
         */
        void setUpEquations();
        /**
         * @brief Solve the discretized energy transport equation
         * to find the temperature profile in the pellet. Also, updates
         * the reaction state of the energetic particles in the reaction zone
         */
        void solveEquations();

        /**
         * @brief Determine whether combustion throughout the pellet is complete
         * 
         * @return true If combustion throughout the pellet is complete
         * @return false If any grid is still undergoing combustion
         */
        bool isCombustionComplete();

        /**
         * @brief Print to file / output stream the temperature profile of the pellet
         * 
         * @param output_stream ostream object where the temperature profile is output
         * @param delimiter Character to separate two consecutive values of temperature
         */
        void printTemperatureProfile(std::ostream &output_stream, char delimiter = '\t');
        
        /**
         * @brief Prints the x-coordinates of the grid points in the discretized pellet
         * 
         * @param output_stream Output stream where the coordinates will be printed
         * @param delimiter Character to separate two consecutive values of temperature
         */
        void printGridPoints(std::ostream &output_stream, char delimiter = '\t');

		void printDiffusionParticleGridPoints(std::ostream &output_stream, unsigned int particle_index, char delimiter = '\t');

		void printDiffusionParticleConcentationProfiles(
			std::ostream &output_stream_A,
			std::ostream &output_stream_B, 
			unsigned int particle_index,
			char delimiter = ','
		);
};

#endif