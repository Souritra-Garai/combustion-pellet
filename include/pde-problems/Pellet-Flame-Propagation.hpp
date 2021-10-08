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

// Required for ArrheniusDiffusivityModel class that returns value of diffusivity
// for various temperatures
#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.hpp"
// Required for PackedPellet class representing thermo-physical properties of 
// a pellet packed with energetic particles and degassed with an inert
#include "thermo-physical-properties/Packed-Pellet.hpp"
// Required for CoreShellDiffusion class representing the 
// diffusion in a core-shell particle
#include "pde-problems/Core-Shell-Diffusion.hpp"
// Required for QRSolver class used to set up and solve matrix equations
#include "qrsolver/QR_Solver.hpp"

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
         * @brief Number of grid points in the pellet
         * including both points on the boundaries
         */
        static size_t _m;
        /**
         * @brief Distance between consecutive grid points
         * on the x axis, \f$ \Delta x \f$
         */
        static real_t _delta_x;

        /**
         * @brief Duration of time step, \f$ \Delta t \f$
         */
        static real_t _delta_t;

        /**
         * @brief Infinitesimal change in temperature used to
         * find derivative of a function with respect to temperature
         * using first principle
         */
        static real_t _delta_T;

        /**
         * @brief Length \f$ l \f$ of the pellet that is at temperature
         * \f$ T_{ign,0} \f$ at time \f$ t = 0 \t$
         * \f$ T \left( t = 0, x \right) = \begin{cases}
         *      T_{ign,0}   &   0 \leq x \leq l\\
         *      T_a         &   l < x \leq L
         * \f$
         */
        static real_t _initial_ignition_length;
        /**
         * @brief Temperature \f$ T_{ign,0} \f$ of the initial small segment 
         * of the pellet at time \f$ t = 0 \f$, to start the combustion front propagation
         * \f$ T \left( t = 0, x \right) = \begin{cases}
         *      T_{ign,0}   &   0 \leq x \leq l\\
         *      T_a         &   l < x \leq L
         * \f$
         */
        static real_t _initial_ignition_temperature;

        /**
         * @brief Stores the time evolved since the zeroth iteration
         */
        real_t _time;

        /**
         * @brief Array to store temperature at each grid point of the pellet
         */
        real_t * _temperature_array;

		real_t * _thermal_conductivity;
        
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
        QRSolver<real_t> _solver;

        /**
         * @brief Arrhenius Diffusivity model to get value of diffusivity in core-shell
         * particle at different diffusivities
         */
        ArrheniusDiffusivityModel<real_t> _diffusivity_model;

        /**
         * @brief Get the x-coordinate of the grid point
         * 
         * @param index Index of the grid point whose x-coordinate is desired
         * @return real_t x-coordiante of grid point # index
         */
        real_t getXCoordinate(size_t index);

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
         * @brief Get the partial derivative of enthalpy of a particle with respect to temperature
         * \f$ \left\delimiter0\frac{\partial Y_k}{\partial T}\right|_{x,t} \left(t_n, x_j\right)
         * \approx \left\delimiter0\frac{\Delta Y_k}{\Delta T}\right|_j^n \f$
         * 
         * @param index Index of the grid point where the particle whose enthalpy derivative is desired, is
         * @return real_t The partial derivative of enthalpy of a particle at the specified grid point 
         * with respect to temperature
         */
        real_t getParticleEnthalpyTemperatureDerivative(size_t index);
        /**
         * @brief Get the partial derivative of enthalpy of a particle with respect to time
         * \left\delimiter0\frac{\partial Y_k}{\partial t}\right|_{x,T} \left(t_n, x_j\right)
         * \approx \left\delimiter0\frac{\Delta Y_k}{\Delta t}\right|_j^n
         * 
         * @param index Index of the grid point where the particle whose enthalpy derivative is desired, is
         * @return real_t The partial derivative of enthalpy of a particle at the specified grid point 
         * with respect to time 
         */
        real_t getParticleEnthalpyTimeDerivative(size_t index);

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

    public:

        /**
         * @brief Set the Grid Size to M including boundary points
         * 
         * @param M Number of grid points including the boundary points
         */
        static void setGridSize(size_t M);
        /**
         * @brief Set the time step duration
         * 
         * @param delta_t 
         */
        static void setTimeStep(real_t delta_t);
        /**
         * @brief Set the infinitesimal change in temperature used to
         * find derivative of a function with respect to temperature
         * using first principle
         * 
         * @param delta_T Infinitesimal change in temperature
         */
        static void setInfinitesimalChangeTemperature(real_t delta_T);

        /**
         * @brief Set the Initial Ignition Parameters
         * 
         * @param initial_ignition_temperature 
         * @param initial_ignition_length 
         */
        static void setInitialIgnitionParameters(
            real_t initial_ignition_temperature,
            real_t initial_ignition_length
        );

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
         * @brief Set the Diffusivity Model for core-shell diffusion
         * of the energetic particle
         * 
         * @param diffusivity_model ArrheniusDiffusivityModel type object
         * that returns diffusivity at a given temperature
         */
        void setDiffusivityModel(ArrheniusDiffusivityModel<real_t> diffusivity_model);

        /**
         * @brief Initialize the pellet by setting the temperature profile and energetic
         * particles to the initial conditions
         */
        void initializePellet();

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
         * @brief Print the configuration of the PDE solver to the output stream
         * 
         * @param output_stream Output stream where the configuration details will be printed
         */
        void printConfiguration(std::ostream &output_stream);

        /**
         * @brief Prints the x-coordinates of the grid points in the discretized pellet
         * 
         * @param output_stream Output stream where the coordinates will be printed
         * @param delimiter Character to separate two consecutive values of temperature
         */
        void printGridPoints(std::ostream &output_stream, char delimiter = '\t');
};

#endif