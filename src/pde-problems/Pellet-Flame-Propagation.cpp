/**
 * @file Pellet-Flame-Propagation.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Implementation details for PelletFlamePropagation class
 * @version 0.1
 * @date 2021-08-02
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "pde-problems/Pellet-Flame-Propagation.hpp"

// Required for pow function
#include <math.h>

#define STEFAN_BOLTZMANN_CONSTANT 5.670374419E-8 // W / m2 - K4

// #include <iostream>

/******************************************************************************************************************/
// Instantiating static member variables of PelletFlamePropagation class

template<typename real_t> real_t PelletFlamePropagation<real_t>::_theta = 1.0;
template<typename real_t> real_t PelletFlamePropagation<real_t>::_gamma = 1.0;

// Number of grid points in the pellet
template<typename real_t> size_t PelletFlamePropagation<real_t>::_m = 2;
// Distance between consecutive grid points
template<typename real_t> real_t PelletFlamePropagation<real_t>::_delta_x = 0;

// Duration of time steps
template<typename real_t> real_t PelletFlamePropagation<real_t>::_delta_t = 0.001;

// Infinitesimal change in temperature
template<typename real_t> real_t PelletFlamePropagation<real_t>::_delta_T = 0.001;

// Initial temperature in the ignited part of the pellet
template<typename real_t> real_t PelletFlamePropagation<real_t>::_initial_ignition_temperature = 273.15;
// Length of initially ignited part of the pellet
template<typename real_t> real_t PelletFlamePropagation<real_t>::_initial_ignition_length = 1;
/******************************************************************************************************************/

/******************************************************************************************************************/
// Static member function definitions

template<typename real_t>
void PelletFlamePropagation<real_t>::setImplicitnessSourceTerm(real_t gamma)
{
	if (0.0 <= gamma && gamma <= 1.0) _gamma = gamma;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setImplicitnessDiffusionTerm(real_t theta)
{
	if (0.0 <= theta && theta <= 1.0) _theta = theta;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setGridSize(size_t M)
{
    // Set the number of grid points, including both boundary points, to M
    _m = M;
    // Set the distance between two consecutive grid points
    _delta_x = PackedPellet<real_t>::_length / (real_t) (_m - 1);
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setTimeStep(real_t delta_t)
{
    // Set the duration of time step
    _delta_t = delta_t;
    // Set the duration of time step for Core-Shell Diffusion problem
    CoreShellDiffusion<real_t>::setTimeStep(_delta_t);
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setInfinitesimalChangeTemperature(real_t delta_T)
{
    // Set the infinitesimal change in temperature used to
    // find derivative of a function with respect to temperature
    // using first principle
    _delta_T = delta_T;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setInitialIgnitionParameters(
    real_t initial_ignition_temperature,
    real_t initial_ignition_length
) {
    // Set length of the pellet that is initial ignited to high temperatures
    _initial_ignition_length = initial_ignition_length;
    // Set the temperature to which a small region of the pellet is initially ignited
    _initial_ignition_temperature = initial_ignition_temperature;
}
/******************************************************************************************************************/

/******************************************************************************************************************/
// Constructors and destructors

template<typename real_t>
PelletFlamePropagation<real_t>::PelletFlamePropagation(
    real_t particle_volume_fraction
) : PackedPellet<real_t>(particle_volume_fraction),
    _solver(_m)
{
    // Allocate memory for holding the temperature profile
    _temperature_array = new real_t[_m];

	_thermal_conductivity = new real_t[_m];

	_prev_internal_energy = new real_t[_m];

    // Allocate memory for core shell diffusion problems at each
    // grid point
    _particles_array = new CoreShellDiffusion<real_t>[_m];

    // Allocate memory for storing core shell diffusion problems
    // in a state one time step back
    _particles_array_const_temperature_evolution = new CoreShellDiffusion<real_t>[_m];
    // Allocate memory for storing a copy of core shell diffusion problems
    // at each grid point
    _particles_array_raised_temperature_evolution = new CoreShellDiffusion<real_t>[_m];

    // Initialize the temperature profile and the core-shell diffusion problems
    // to the initial conditions
    initializePellet();
}

template<typename real_t>
PelletFlamePropagation<real_t>::~PelletFlamePropagation()
{
    // Deallocate memory for temperature array
    delete [] _temperature_array;
	delete [] _thermal_conductivity;
	delete [] _prev_internal_energy;

    // Deallocate memory for core shell diffusion problems
    delete [] _particles_array;
    delete [] _particles_array_const_temperature_evolution;
    delete [] _particles_array_raised_temperature_evolution;
}
/******************************************************************************************************************/

/******************************************************************************************************************/
// Defining private member functions

template<typename real_t>
inline real_t PelletFlamePropagation<real_t>::getXCoordinate(size_t index) {
    // Linearly map M grid points to the length of the pellet
    return (real_t) index * this->_length / (real_t) (_m - 1);
}

template<typename real_t>
inline bool PelletFlamePropagation<real_t>::inReactionZone(size_t index)
{
    // A grid point is inside the reaction zone when
    return
        // the temperature at that grid point is equal to or more than 
        // the ignition temperature of energetic particle
        _temperature_array[index] >= this->_ignition_temperature &&
        // and the energetic particle has not undergone complete combustion
        !_particles_array[index].isCombustionComplete();
}

template<typename real_t>
real_t PelletFlamePropagation<real_t>::getParticleEnthalpyTemperatureDerivative(size_t i)
{
    // Calculate \f$ \sum_{k \in \text{Particle}} \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right) \right\}
    // \cdot \left\delimiter0\frac{\Delta Y_k}{\Delta T}\right|_j^n \f$
    
    // Here, \f$ \left\delimiter0\frac{\Delta Y_k}{\Delta T}\right|_j^n \f$ is calculated by
    // \f$ \left\delimiter0\frac{\Delta Y_k}{\Delta T}\right|_j^n =
    // \frac{   Y_{k,j}^{n}\left(T_j^{n-1} + \Delta T, \left\{C_i^{n-1}\right\}_j\right) - 
    //          Y_{k,j}^{n}\left(T_j^{n-1}, \left\{C_i^{n-1}\right\}_j\right)   }
    //      {   \Delta T    }
    // \f$

    // Return the difference between enthalpy of particle after evolving
    // the particle at raised temperature and that of the particle after evolving
    // at the same temperature, divided the by the incremental change in temperature
    return
    (   // \f$  \sum_{k \in \text{Particle}} \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right) \right\}
        //      \cdot Y_{k,j}^{n}\left(T_j^{n-1} + \Delta T, \left\{C_i^{n-1}\right\}_j\right)  \f$
        _particles_array_raised_temperature_evolution[i].getEnthalpy(_temperature_array[i]) - 
        // \f$  \sum_{k \in \text{Particle}} \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right) \right\}
        //      \cdot Y_{k,j}^{n}\left(T_j^{n-1}, \left\{C_i^{n-1}\right\}_j\right) \f$
        _particles_array_const_temperature_evolution[i].getEnthalpy(_temperature_array[i])
    ) / _delta_T;
}

template<typename real_t>
real_t PelletFlamePropagation<real_t>::getParticleEnthalpyTimeDerivative(size_t i)
{
    // Calculate \f$ \sum_{k \in \text{Particle}} \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right) \right\}
    // \cdot \left\delimiter0\frac{\Delta Y_k}{\Delta t}\right|_j^n \f$

    // Here, \f$ \left\delimiter0\frac{\Delta Y_k}{\Delta t}\right|_j^n \f$ is calculated by
    // \f$ \left\delimiter0\frac{\Delta Y_k}{\Delta t}\right|_j^n =
    // \frac{   Y_{k,j}^{n}\left(T_j^{n-1}, \left\{C_i^{n-1}\right\}_j\right) -
    //          Y_{k,j}^{n-1}   }
    //      {   \Delta t    }
    // \f$

    // Return the difference between enthalpy of particle after evolving
    // the particle at same temperature and that of the particle at the initial state,
    // divided by the time step
    return
    (   // \f$  \sum_{k \in \text{Particle}} \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right) \right\}
        //      \cdot Y_{k,j}^{n}\left(T_j^{n-1}, \left\{C_i^{n-1}\right\}_j\right) \f$
        _particles_array_const_temperature_evolution[i].getEnthalpy(_temperature_array[i]) -
        // \f$  \sum_{k \in \text{Particle}} \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right) \right\}
        //      \cdot Y_{k,j}^{n-1} \f$
        _particles_array[i].getEnthalpy(_temperature_array[i])
    ) / _delta_t;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::evolveParticleForEnthalpyDerivative(size_t i)
{
    // Since the evolution of the particles with constant temperature and raised temperature
    // are independent of each other, they can be parallelized
    #pragma omp parallel sections default(shared)
    {
        // Parallel section for constant temperature evolution
        #pragma omp section
        {
            // // Create copies of the current particle to evolve
            // _particles_array_const_temperature_evolution[i].copyFrom(_particles_array[i]);
            // Set up equations to evolve particle from \f$ \left( t_{n-1}, T_j^{n-1}, \left\{ C_j^{n-1} \right\} \right) \f$
            // to \f$ \left( t_n, T_j^n, \left\{ C_j^n \right\} \right) \f$
            _particles_array_const_temperature_evolution[i].setUpEquations(
                _diffusivity_model.getDiffusivity(_temperature_array[i]),
				_particles_array[i]
            );
            // Solve the equations to update the state of the particle
            _particles_array_const_temperature_evolution[i].solveEquations();
        }

        // Parallel section for evolution at an infinitesimally higher temperature
        #pragma omp section
        {
            // // Create copies of the current particle to evolve
            // _particles_array_raised_temperature_evolution[i].copyFrom(_particles_array[i]);
            // Set up equations to evolve particle from \f$ \left( t_{n-1}, T_j^{n-1} + \Delta T, \left\{ C_j^{n-1} \right\} \right) \f$
            // to \f$ \left( t_n, T_j^n, \left\{ C_j^n \right\} \right) \f$
            _particles_array_raised_temperature_evolution[i].setUpEquations(
                _diffusivity_model.getDiffusivity(_temperature_array[i] + _delta_T),
				_particles_array[i]
            );
            // Solve the equations to update the state of the particle
            _particles_array_raised_temperature_evolution[i].solveEquations();
        }
    }
}

template<typename real_t>
LinearExpression<real_t> PelletFlamePropagation<real_t>::calcTransientTerm(size_t i)
{
    // Develop a linear expression for the trasient term in terms of temperature
    // at the present time step, \f$ \alpha_{0,j}^n + \alpha_{1,j}^n \cdot T_j^n \f$
    LinearExpression<real_t> expression;

    // Here,
    // \f$ \alpha_{0,j}^n = - \sum_{k \in \text{Pellet}} Y_{k,Pellet} \cdot c_k \cdot \dfrac{T^{n-1}_j}{\Delta t} + 
    // \sum_{k \in \text{Particle}} \left( 1 - Y_{Ar,Pellet} \right) \cdot \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right)\right\}
    // \cdot \left( \left\delimiter0\frac{\Delta Y_k}{\Delta t}\right|_j^n - 
    // \left\delimiter0\frac{\Delta Y_k}{\Delta T}\right|_j^n \cdot \dfrac{T^{n-1}_j}{\Delta t} \right) \f$
    // and,
    // \f$ \alpha_{1,j}^n = \sum_{k \in \text{Pellet}} Y_{k,Pellet} \cdot c_k \cdot \dfrac{1}{\Delta t} +
    // \sum_{k \in \text{Particle}} \left( 1 - Y_{Ar,Pellet} \right) \cdot \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right)\right\}
    // \cdot \left( \left\delimiter0\frac{\Delta Y_k}{\Delta T}\right|_j^n \cdot \dfrac{1}{\Delta t} \right) \f$

    // Note, \f$ \alpha_{0,j}^n = - \alpha_{1,j}^n \cdot T_j^{n-1} +
    // \sum_{k \in \text{Particle}} \left( 1 - Y_{Ar,Pellet} \right) \cdot \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right)\right\}
    // \cdot \left\delimiter0\frac{\Delta Y_k}{\Delta t}\right|_j^n \f$
    
    // \f$ \alpha_{1,j}^n = \sum_{k \in \text{Pellet}} Y_{k,Pellet} \cdot c_k \cdot \dfrac{1}{\Delta t} \f$
    expression.coefficient = this->getHeatCapacity(_particles_array + i, _temperature_array[i]) / _delta_t;
    // Instantiate \f$ \alpha_{0,j}^n \f$ and initialize with 0
    expression.constant = 0.0;

    // If the grid point is within the reaction zone, only then activate the reaction term
    if (inReactionZone(i))
    {
        // Evolve the particles used for determining change in enthalpy at constant and
        // raised temperatures. This must be done before calling 
        evolveParticleForEnthalpyDerivative(i);
        
        // std::cout << "\tParticle # " << i << std::endl;
        // std::cout << "\t\tdudT : " << getParticleEnthalpyTemperatureDerivative(i) / _delta_t << std::endl;
        // std::cout << "\t\tdudt : " << getParticleEnthalpyTimeDerivative(i) << std::endl;

        // \f$ \alpha_{1,j}^n += \sum_{k \in \text{Particle}} \left( 1 - Y_{Ar,Pellet} \right) \cdot \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right)\right\}
        // \cdot \left( \left\delimiter0\frac{\Delta Y_k}{\Delta T}\right|_j^n \cdot \dfrac{1}{\Delta t} \right) \f$
        expression.coefficient += _gamma * this->_particle_mass_fractions * getParticleEnthalpyTemperatureDerivative(i) / _delta_t;

        // \f$ \alpha_{0,j}^n += \sum_{k \in \text{Particle}} \left( 1 - Y_{Ar,Pellet} \right) \cdot \left\{ \Delta H_{f,k}^0 + c_k \left(T^{n-1}_j - T_{ref}\right)\right\}
        // \cdot \left\delimiter0\frac{\Delta Y_k}{\Delta t}\right|_j^n \f$
        expression.constant += this->_particle_mass_fractions * (
			(1 - _gamma) * (_particles_array[i].getEnthalpy(_temperature_array[i]) - _prev_internal_energy[i]) / _delta_t +
			_gamma  * getParticleEnthalpyTimeDerivative(i));
    }
    
    // Return the linearized expression for transient term
    return expression;
}

template<typename real_t>
LinearExpression<real_t> PelletFlamePropagation<real_t>::calcHeatLossTerm(size_t i)
{
    // Develop a linear expression for the heat loss term in terms of temperature
    // at the present time step, \f$ \beta_{0,j}^n + \beta_{1,j}^n \cdot T_j^n \f$
    LinearExpression<real_t> expression;

    // Here,
    // \f$ \beta_{0,j}^n = - \frac{4}{D} \left[ h \cdot T_a + 
    // \varepsilon \sigma \left\{ 3 \left(T_j^{n-1}\right)^4 + T_a^4 \right\} \right] \f$
    expression.constant =
		this->_convective_heat_transfer_coefficient_curved_surface * (_temperature_array[i] - this->_ambient_temperature) +
		this->_radiative_emissivity * STEFAN_BOLTZMANN_CONSTANT * (pow(_temperature_array[i], 4) - pow(this->_ambient_temperature, 4))
	;

    // and,
    // \f$ \beta_{1,j}^n = \frac{4}{D} \left\{ h + 4 \varepsilon \sigma \left(T_j^{n-1}\right)^3 \right\} \f$
    expression.coefficient = _gamma * (
		this->_convective_heat_transfer_coefficient_curved_surface +
		4.0 * this->_radiative_emissivity * STEFAN_BOLTZMANN_CONSTANT * pow(_temperature_array[i], 3)
	);

    // Return the linearized expression for heat loss term
    return expression;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setUpBoundaryConditionX0()
{
    // Set up the discretized boundary equation
    // \f$ - \frac{\lambda}{\Delta x} \cdot T_1^n + 
    // \left( \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,0}^n \right) T_0^n 
    // = - \frac{D}{4} \cdot \beta_{0,0}^n \f$

    // Get the linear expression for heat loss term at boundary \f$ x = 0 \f$, \f$ \beta_{0,0}^n + \beta_{1,0}^n \cdot T_0^n \f$ 
    LinearExpression<real_t> heat_loss_term = calcHeatLossTerm(0);
    
    // Get the effective heat conductivity of particle - gas mixture, divided by the grid size
    real_t lambda_by_delta_x = _thermal_conductivity[0] / _delta_x;
    // Since, the particle at \f$ x = 0 \f$ grid point is not evolved, the effective heat conductivity
    // the next grid point is used

    // Set up the first row for the matrix equation
    _solver.setEquationFirstRow(
        // Coefficient of \f$ T_0^n \f$, \f$ \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,0}^n \f$
        (this->_diameter / 4.0) * heat_loss_term.coefficient + lambda_by_delta_x,
        // Coefficient of \f$ T_1^n \f$, \f$ - \frac{\lambda}{\Delta x} \f$
        - lambda_by_delta_x,
        // Constant term, \f$ - \frac{D}{4} \cdot \beta_{0,0}^n \f$
        - (this->_diameter / 4.0) * heat_loss_term.constant
    );
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setUpBoundaryConditionXN()
{
    // Set up the discretized boundary equation
    // \f$ \left( \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,M}^n \right) T_M^n 
    // - \frac{\lambda}{\Delta x} \cdot T_{M-1}^n
    // = - \frac{D}{4} \cdot \beta_{0,M}^n

    // Get the linear expression for heat loss term at boundary \f$ x = 1 \f$, \f$ \beta_{0,M}^n + \beta_{1,M}^n \cdot T_0^n \f$ 
    LinearExpression<real_t> heat_loss_term = calcHeatLossTerm(_m-1);

    // Get the effective heat conductivity of particle - gas mixture, divided by the grid size
    real_t lambda_by_delta_x = _thermal_conductivity[_m - 1] / _delta_x;
    // Since, the particle at \f$ x = 0 \f$ grid point is not evolved, the effective heat conductivity
    // the next grid point is used

    // Set up the first row for the matrix equation
    _solver.setEquationLastRow(
        // Coefficient of \f$ T_{M-1}^n \f$, \f$ - \frac{\lambda}{\Delta x} \f$
        - lambda_by_delta_x,
        // Coefficient of \f$ T_M^n \f$, \f$ \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,M}^n \f$
        (this->_diameter / 4.0) * heat_loss_term.coefficient + lambda_by_delta_x,
        // Constant term, \f$ - \frac{D}{4} \cdot \beta_{0,M}^n \f$
        - (this->_diameter / 4.0) * heat_loss_term.constant
    );
}

template<typename real_t>
void PelletFlamePropagation<real_t>::updateParticlesState()
{
    // Updating the state of the energetic particle at each grid point
    // parallely as they are independent of each other
    #pragma omp parallel for default(shared) schedule(dynamic)
        // For each grid point update the state of the energetic particle
        for (size_t i = 1; i < _m-1; i++)
        {
			_prev_internal_energy[i] = _particles_array[i].getEnthalpy(_temperature_array[i]);

            // Update the state of energetic particles only in the reaction zone
            if (inReactionZone(i))
            {
                // Set up the matrix equation representing diffusion at the
                // updated temperature of that grid point
                _particles_array[i].setUpEquations(_diffusivity_model.getDiffusivity(_temperature_array[i]));
                // Solve the equations to update the state of the particle
                _particles_array[i].solveEquations();
            }

			_thermal_conductivity[i] = this->getThermalConductivity(_particles_array + i, _temperature_array[i]);
        }

	_thermal_conductivity[0]	= this->getThermalConductivity(_particles_array + 1, 		_temperature_array[0]);
	_thermal_conductivity[_m-1] = this->getThermalConductivity(_particles_array + (_m - 2), _temperature_array[_m-1]);
}

/******************************************************************************************************************/

/******************************************************************************************************************/
// Defining public member functions

template<typename real_t>
void PelletFlamePropagation<real_t>::initializePellet()
{
    // Initialize the temperature and energetic particle at each grid point
    // with the specified initial conditions

    // Set time to 0
    _time = 0;

    // Set temperature at grid point \f$ x = 0 \f$ to the initial ignition temperature
    _temperature_array[0] = _initial_ignition_temperature;

    // Parallelly initialize the grid points
    #pragma omp parallel for
    
        // For each grid point
        for (size_t i = 1; i < _m-1; i++)
        {
            // If grid point is within the initial ignition length
            if (getXCoordinate(i) < _initial_ignition_length) _temperature_array[i] = _initial_ignition_temperature;
            // Set the temperature to the initial ignition temperature

            // Else set the temperature to ambient temperature
            else _temperature_array[i] = this->_ambient_temperature;

            // Initialize the main particle at each grid point
            _particles_array[i].initializeParticle();
            
            // Initialize the particles used for 
            _particles_array_raised_temperature_evolution[i].initializeParticle();
            _particles_array_const_temperature_evolution[i].initializeParticle();

			_thermal_conductivity[i] = this->getThermalConductivity(_particles_array + i, _temperature_array[i]);
			_prev_internal_energy[i] = _particles_array[i].getEnthalpy(_temperature_array[i]);
        }
    // Set temperature at grid point at \f$ x = L \f$ to ambient temperature
    _temperature_array[_m-1] = this->_ambient_temperature;

	_thermal_conductivity[0]	= this->getThermalConductivity(_particles_array + 1,		_temperature_array[0]);
	_thermal_conductivity[_m-1] = this->getThermalConductivity(_particles_array + (_m - 2), _temperature_array[_m-1]);
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setDiffusivityModel(
    ArrheniusDiffusivityModel<real_t> diffusivity_model
) {
    // Set the diffusivity model to the input model
    _diffusivity_model = diffusivity_model;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setUpEquations()
{
    // Set up the matrix equation representing energy transport in the pellet

    // For the first row of the matrix equation, set up the boundary condition at \f$ x = 0 \f$
    // \f$ - \frac{\lambda}{\Delta x} \cdot T_1^n + \left( \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,0}^n \right) T_0^n 
    // = - \frac{D}{4} \cdot \beta_{0,0}^n \f$
    setUpBoundaryConditionX0();

    // As the 
    #pragma omp parallel for default(shared) schedule(dynamic)

        for (size_t i = 1; i < _m-1; i++)
        {
            // For matrix equations corresponding to grid points inside the pellet,
            // set up the discretized energy transport equation
            // \f$ - \frac{\lambda}{\left( \Delta x \right)^2} \cdot T_{j+1}^n + 
            // \left\{ \rho_0 \alpha_{1,j}^n + \frac{2 \lambda}{\left( \Delta x \right)^2} + \beta_{1,j}^n \right\} \cdot T_j^n -
            // \frac{\lambda}{\left( \Delta x \right)^2} \cdot T_{j-1}^n
            // = - \rho_0 \alpha_{0,j}^n - \beta_{0,j}^n \f$

            // Get the linearized expression for the transient term, \f$ \alpha_{0,j}^n + \alpha_{1,j}^n \cdot T_j^n \f$
            LinearExpression<real_t> transient_term = calcTransientTerm(i);

            // Get the linearized expression for the heat loss term, \f$ \beta_{0,j}^n + \beta_{1,j}^n \cdot T_j^n \f$
            LinearExpression<real_t> heat_loss_term = calcHeatLossTerm(i);

			real_t lambda_forward_by_delta_x_sqr  = 0.5 * (_thermal_conductivity[i+1] + _thermal_conductivity[i]) / pow(_delta_x, 2);
			real_t lambda_backward_by_delta_x_sqr = 0.5 * (_thermal_conductivity[i] + _thermal_conductivity[i-1]) / pow(_delta_x, 2);

            // Set up the matrix equation
            _solver.setEquation(
                // Matrix row number
                i,
                // Coefficient of \f$ T_{j-1}^n \f$, \f$ - \frac{\lambda}{\left( \Delta x \right)^2} \f$
                - _theta * lambda_backward_by_delta_x_sqr,
                // Coefficient of \f$ T_j^n \f$, \f$ \left\{ \rho_0 \alpha_{1,j}^n + \frac{2 \lambda}{\left( \Delta x \right)^2} + \beta_{1,j}^n \right\} \f$
                this->_density * transient_term.coefficient + (4 * heat_loss_term.coefficient / this->_diameter) + _theta * (lambda_forward_by_delta_x_sqr + lambda_backward_by_delta_x_sqr),
                // Coefficient of \f$ T_{j+1}^n \f$, \f$ - \frac{\lambda}{\left( \Delta x \right)^2} \f$
                - _theta * lambda_forward_by_delta_x_sqr,
                (1 - _theta) * lambda_forward_by_delta_x_sqr * _temperature_array[i+1] +
				(this->_density * transient_term.coefficient - (1 - _theta) * (lambda_forward_by_delta_x_sqr + lambda_backward_by_delta_x_sqr) + (4 * heat_loss_term.coefficient / this->_diameter)) * _temperature_array[i] +
				(1 - _theta) * lambda_backward_by_delta_x_sqr * _temperature_array[i-1] -
				(this->_density * transient_term.constant + (4 * heat_loss_term.constant / this->_diameter))                
            );
        }

    // For the last row of the matrix equation, set up the boundary condition at \f$ x = M \f$
    // \f$ \left( \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,M}^n \right) T_M^n 
    // - \frac{\lambda}{\Delta x} \cdot T_{M-1}^n
    // = - \frac{D}{4} \cdot \beta_{0,M}^n
    setUpBoundaryConditionXN();

    // _solver.printMatrixEquation();
}

template<typename real_t>
void PelletFlamePropagation<real_t>::solveEquations()
{
    _time += _delta_t;
    // Solve the matrix equation
    _solver.getSolution(_temperature_array);
    // Update the state of the energetic particles
    // using the temperature at the present step
    updateParticlesState();
}

template<typename real_t>
bool PelletFlamePropagation<real_t>::isCombustionComplete()
{
    // Assuming combustion is complete by default
    bool flag = true;
    // If a point is in the reaction zone, then combustion is not complete

    // For all the points inside the pellet where
    // energetic particle undergoes combustion
    for (size_t i=1; i < _m-1 && flag; i++)
    {
        // If the particle is in the reaction zone then 
        // combustion in the pellet is incomplete
        flag = !inReactionZone(i);
    }
    // Return whether combustion is complete or not
    return flag;

	// return _temperature_array[_m - 1] > 1000;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::printTemperatureProfile(
    std::ostream &output_stream,
    char delimiter
) {
    // First print the current time, followed by the temperature profile
    output_stream << _time << delimiter;
    
    // Print the temperature followed by the delimiter for all except
    // the last grid point
    for (size_t i = 0; i < _m -1; i++) output_stream << _temperature_array[i] << delimiter;
    // For the last grid point print the temperature followed by endline
    output_stream << _temperature_array[_m-1] << std::endl;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::printGridPoints(
    std::ostream &output_stream,
    char delimiter
) {
    // First value is NaN since first column stores the time
    output_stream << NAN << delimiter;

    // Print the x-coordinate of the grid point followed by the delimiter for all
    // except the last grid point
    for (size_t i = 0; i < _m - 1; i++) output_stream << getXCoordinate(i) << delimiter;
    // For the last grid point print the x-coordinate followed by endline
    output_stream << getXCoordinate(_m-1) << std::endl;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::printConfiguration(
    std::ostream &output_stream
) {
	output_stream << "Pellet Flame Propagation PDE Solver Configuration\n\n";
    output_stream << "Time Step\t:\t" << _delta_t << "\ts" << std::endl;

    output_stream << "Number of Grid Points\t:\t" << _m << std::endl;

    output_stream << "Grid Size\t:\t" << _delta_x << "\tm" << std::endl;

    output_stream << "Delta T\t:\t" << _delta_T << "\tK" << std::endl;

	output_stream << "Diffusion Term Implicitness\t:\t" << _theta << std::endl;

	output_stream << "Source Term Implicitness\t:\t" << _gamma << std::endl;

    output_stream << "\nInitial Ignitiion Parameters" << std::endl;
    output_stream << "Temperature\t:\t" << _initial_ignition_temperature << "\tK" << std::endl;
    output_stream << "Length\t:\t" << _initial_ignition_length << "\tm" << std::endl;

    output_stream << "\nCore-Shell Particle Diffusivity Parameter" << std::endl;
    output_stream << "Pre Exponential Factor\t:\t" << _diffusivity_model.getPreExponentialFactor() << "\tm2 / s" << std::endl;
    output_stream << "Activation Energy\t:\t" << _diffusivity_model.getActivationEnergy() << "\tJ / mol." << std::endl;

    output_stream << "\nCore-Shell Particle Diffusion Solver Parameters" << std::endl;
    _particles_array->printConfiguration(output_stream);

	output_stream << std::endl;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::printDiffusionParticleGridPoints(
	std::ostream &output_stream,
	unsigned int particle_index,
	char delimiter
) {
	_particles_array[particle_index].printGridPoints(output_stream, delimiter);
}

template<typename real_t>
void PelletFlamePropagation<real_t>::printDiffusionParticleConcentationProfiles(
	std::ostream &output_stream_A,
	std::ostream &output_stream_B, 
	unsigned int particle_index,
	char delimiter
) {
	_particles_array[particle_index].printConcentrationProfileA(output_stream_A, delimiter, _time);
	_particles_array[particle_index].printConcentrationProfileB(output_stream_B, delimiter, _time);
}

/******************************************************************************************************************/

/******************************************************************************************************************/
// Create classes for float, double and long double data types
template class PelletFlamePropagation<float>;
template class PelletFlamePropagation<double>;
template class PelletFlamePropagation<long double>;