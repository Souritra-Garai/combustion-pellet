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

#include "utilities/Read-Data.hpp"

// Required for pow function
#include <math.h>

#define STEFAN_BOLTZMANN_CONSTANT 5.670374419E-8 // W / m2 - K4

/******************************************************************************************************************/
// Instantiating static member variables of PelletFlamePropagation class

template<typename real_t> const real_t PelletFlamePropagation<real_t>::kappa = readScalarData<real_t>("data/PDE-solver-config", "kappa.txt");
template<typename real_t> const real_t PelletFlamePropagation<real_t>::gamma = readScalarData<real_t>("data/PDE-solver-config", "gamma.txt");

// Number of grid points in the pellet
template<typename real_t> const size_t PelletFlamePropagation<real_t>::m = readScalarData<real_t>("data/PDE-solver-config", "number-of-grid-points-pellet.txt");
// Distance between consecutive grid points
template<typename real_t> const real_t PelletFlamePropagation<real_t>::delta_x = readScalarData<real_t>("data/pellet", "length.txt") / ((real_t) PelletFlamePropagation<real_t>::m - 1.);

// Duration of time steps
template<typename real_t> const real_t PelletFlamePropagation<real_t>::delta_t = readScalarData<real_t>("data/PDE-solver-config", "delta_t.txt");

// Infinitesimal change in temperature
template<typename real_t> const real_t PelletFlamePropagation<real_t>::delta_T = readScalarData<real_t>("data/PDE-solver-config", "delta_T.txt");

/******************************************************************************************************************/

/******************************************************************************************************************/
// Static member function definitions

template<typename real_t>
inline real_t PelletFlamePropagation<real_t>::getXCoordinate(size_t index) {
    // Linearly map M grid points to the length of the pellet
    return (real_t) index * PackedPellet<real_t>::length / (real_t) (m - 1);
}

/******************************************************************************************************************/

/******************************************************************************************************************/
// Constructors and destructors

template<typename real_t>
PelletFlamePropagation<real_t>::PelletFlamePropagation(
    real_t particle_volume_fraction
) : PackedPellet<real_t>(particle_volume_fraction),
    _solver(m)
{
    // Allocate memory for holding the temperature profile
    _temperature_array = new real_t[m];

	_thermal_conductivity = new real_t[m];

	_prev_enthalpy_particle = new real_t[m];

    // Allocate memory for core shell diffusion problems at each
    // grid point
    _particles_array = new CoreShellDiffusion<real_t>[m];

    // Allocate memory for storing core shell diffusion problems
    // in a state one time step back
    _particles_array_const_temperature_evolution = new CoreShellDiffusion<real_t>[m];
    // Allocate memory for storing a copy of core shell diffusion problems
    // at each grid point
    _particles_array_raised_temperature_evolution = new CoreShellDiffusion<real_t>[m];
	
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
	delete [] _prev_enthalpy_particle;

    // Deallocate memory for core shell diffusion problems
    delete [] _particles_array;
    delete [] _particles_array_const_temperature_evolution;
    delete [] _particles_array_raised_temperature_evolution;
}

/******************************************************************************************************************/

/******************************************************************************************************************/
// Defining private member functions

template<typename real_t>
inline bool PelletFlamePropagation<real_t>::inReactionZone(size_t index)
{
    // A grid point is inside the reaction zone when
    return
		// the energetic particle has not undergone complete combustion
        !_particles_array[index].isCombustionComplete();
}

template<typename real_t>
void PelletFlamePropagation<real_t>::evolveParticleForEnthalpyDerivative(size_t i)
{
	// Set up equations to evolve particle from \f$ \left( t_{n-1}, T_j^{n-1}, \left\{ C_j^{n-1} \right\} \right) \f$
	// to \f$ \left( t_n, T_j^n, \left\{ C_j^n \right\} \right) \f$
	_particles_array_const_temperature_evolution[i].setUpEquations(
		_temperature_array[i],
		_particles_array[i]
	);

	// Solve the equations to update the state of the particle
	_particles_array_const_temperature_evolution[i].solveEquations();

	// Set up equations to evolve particle from \f$ \left( t_{n-1}, T_j^{n-1} + \Delta T, \left\{ C_j^{n-1} \right\} \right) \f$
	// to \f$ \left( t_n, T_j^n, \left\{ C_j^n \right\} \right) \f$
	_particles_array_raised_temperature_evolution[i].setUpEquations(
		_temperature_array[i] + delta_T,
		_particles_array[i]
	);

	// Solve the equations to update the state of the particle
	_particles_array_raised_temperature_evolution[i].solveEquations();
}

template<typename real_t>
inline LinearExpression<real_t> PelletFlamePropagation<real_t>::calcTransientTerm(size_t i)
{
	static const real_t delta_t_inverse				= 1. / delta_t;
	static const real_t gamma_by_delta_t			= gamma / delta_t;
	static const real_t one_minus_gamma_by_delta_t	= (1. - gamma) / delta_t;
	static const real_t gamma_by_delta_T_delta_t	= gamma / (delta_t * delta_T);
    
	LinearExpression<real_t> expression;
	
	expression.coefficient = _particles_array[i].getHeatCapacity(_temperature_array[i]) * delta_t_inverse;
    
    expression.constant = 0.0;

    // If the grid point is within the reaction zone, only then activate the reaction term
    if (inReactionZone(i))
    {
		evolveParticleForEnthalpyDerivative(i);

		real_t enthalpy						= _particles_array[i].getEnthalpy(_temperature_array[i]);
		real_t enthalpy_const_T_evolution	= _particles_array_const_temperature_evolution[i].getEnthalpy(_temperature_array[i]);
		real_t enthalpy_raised_T_evolution	= _particles_array_raised_temperature_evolution[i].getEnthalpy(_temperature_array[i]);
        
        expression.coefficient	+=	gamma_by_delta_T_delta_t * (enthalpy_raised_T_evolution - enthalpy_const_T_evolution);

        expression.constant		+=	one_minus_gamma_by_delta_t * (enthalpy - _prev_enthalpy_particle[i])
								+ 	gamma_by_delta_t * (enthalpy_const_T_evolution - enthalpy);

		_prev_enthalpy_particle[i] = enthalpy;
    }

	expression.coefficient	*= PackedPellet<real_t>::overall_particle_density;
	expression.constant		*= PackedPellet<real_t>::overall_particle_density;
    
    // Return the linearized expression for transient term
    return expression;
}

template<typename real_t> 
inline real_t PelletFlamePropagation<real_t>::getInterstitialGasTransientTermCoefficient(size_t i)
{
	static const real_t constant = PackedPellet<real_t>::interstitial_volume_fractions / delta_t;

	return
		constant * 
		PackedPellet<real_t>::interstitial_gas.getDensity(_temperature_array[i]) *
		PackedPellet<real_t>::interstitial_gas.getCp(_temperature_array[i]);
}

template<typename real_t>
inline LinearExpression<real_t> PelletFlamePropagation<real_t>::calcHeatLossTerm(size_t i)
{
	static const real_t constant1 = PackedPellet<real_t>::radiative_emissivity * STEFAN_BOLTZMANN_CONSTANT;
	static const real_t constant2 = 4. * PackedPellet<real_t>::radiative_emissivity * STEFAN_BOLTZMANN_CONSTANT;

    // Develop a linear expression for the heat loss term in terms of temperature
    // at the present time step, \f$ \beta_{0,j}^n + \beta_{1,j}^n \cdot T_j^n \f$
    LinearExpression<real_t> expression;

    expression.constant =
		PackedPellet<real_t>::convective_heat_transfer_coefficient_curved_surface * (_temperature_array[i] - PackedPellet<real_t>::ambient_temperature) +
		constant1 * (pow(_temperature_array[i], 4) - pow(PackedPellet<real_t>::ambient_temperature, 4))
	;

    expression.coefficient = gamma * (
		PackedPellet<real_t>::convective_heat_transfer_coefficient_curved_surface +
		constant2 * pow(_temperature_array[i], 3)
	);

    // Return the linearized expression for heat loss term
    return expression;
}

template<typename real_t>
inline void PelletFlamePropagation<real_t>::setUpBoundaryConditionX0()
{
	static const real_t constant = 1. / delta_x;
    // Set up the discretized boundary equation
    // \f$ - \frac{\lambda}{\Delta x} \cdot T_1^n + 
    // \left( \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,0}^n \right) T_0^n 
    // = - \frac{D}{4} \cdot \beta_{0,0}^n \f$

    // Get the linear expression for heat loss term at boundary \f$ x = 0 \f$, \f$ \beta_{0,0}^n + \beta_{1,0}^n \cdot T_0^n \f$ 
    LinearExpression<real_t> beta = calcHeatLossTerm(0);
    
    // Get the effective heat conductivity of particle - gas mixture, divided by the grid size
    real_t lambda_by_delta_x = _thermal_conductivity[0] * constant;
    // Since, the particle at \f$ x = 0 \f$ grid point is not evolved, the effective heat conductivity
    // the next grid point is used

    // Set up the first row for the matrix equation
    _solver.setEquationFirstRow(
        // Coefficient of \f$ T_0^n \f$, \f$ \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,0}^n \f$
        beta.coefficient + lambda_by_delta_x,
        // Coefficient of \f$ T_1^n \f$, \f$ - \frac{\lambda}{\Delta x} \f$
        - lambda_by_delta_x,
        // Constant term, \f$ - \frac{D}{4} \cdot \beta_{0,0}^n \f$

        beta.coefficient * _temperature_array[0] - beta.constant
    );
}

template<typename real_t>
inline void PelletFlamePropagation<real_t>::setUpBoundaryConditionXN()
{
	static const real_t constant = 1. / delta_x;
    // Set up the discretized boundary equation
    // \f$ \left( \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,M}^n \right) T_M^n 
    // - \frac{\lambda}{\Delta x} \cdot T_{M-1}^n
    // = - \frac{D}{4} \cdot \beta_{0,M}^n

    // Get the linear expression for heat loss term at boundary \f$ x = 1 \f$, \f$ \beta_{0,M}^n + \beta_{1,M}^n \cdot T_0^n \f$ 
    LinearExpression<real_t> beta = calcHeatLossTerm(m-1);

    // Get the effective heat conductivity of particle - gas mixture, divided by the grid size
    real_t lambda_by_delta_x = _thermal_conductivity[m - 1] * constant;
    // Since, the particle at \f$ x = 0 \f$ grid point is not evolved, the effective heat conductivity
    // the next grid point is used

    // Set up the first row for the matrix equation
    _solver.setEquationLastRow(
        // Coefficient of \f$ T_{M-1}^n \f$, \f$ - \frac{\lambda}{\Delta x} \f$
        - lambda_by_delta_x,
        // Coefficient of \f$ T_M^n \f$, \f$ \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,M}^n \f$
        beta.coefficient + lambda_by_delta_x,
        // Constant term, \f$ - \frac{D}{4} \cdot \beta_{0,M}^n \f$
        beta.coefficient * _temperature_array[m - 1] - beta.constant
    );
}

template<typename real_t>
void PelletFlamePropagation<real_t>::updateParticlesState()
{
    // Updating the state of the energetic particle at each grid point
    // parallely as they are independent of each other
    #pragma omp parallel for default(shared) schedule(static, 1)
        // For each grid point update the state of the energetic particle
        for (size_t i = 1; i < m-1; i++)
        {
            // Update the state of energetic particles only in the reaction zone
            if (inReactionZone(i))
            {
                // Set up the matrix equation representing diffusion at the
                // updated temperature of that grid point
                _particles_array[i].setUpEquations(_temperature_array[i]);
                // Solve the equations to update the state of the particle
                _particles_array[i].solveEquations();
            }
			
			_thermal_conductivity[i] = PackedPellet<real_t>::getThermalConductivity(_particles_array + i, _temperature_array[i]);
        }

	_thermal_conductivity[0]	= PackedPellet<real_t>::getThermalConductivity(_particles_array + 1, 		_temperature_array[0]);
	_thermal_conductivity[m-1]	= PackedPellet<real_t>::getThermalConductivity(_particles_array + (m - 2),_temperature_array[m-1]);
}

/******************************************************************************************************************/

/******************************************************************************************************************/
// Defining public member functions

template<typename real_t>
void PelletFlamePropagation<real_t>::initializePellet(
	real_t initial_ignition_temperature,
	real_t initial_ignition_length_fraction
)
{
    // Initialize the temperature and energetic particle at each grid point
    // with the specified initial conditions

    // Set time to 0
    _time = 0;

    // Set temperature at grid point \f$ x = 0 \f$ to the initial ignition temperature
    _temperature_array[0] = initial_ignition_temperature;

	real_t initial_ignition_length = PackedPellet<real_t>::length * initial_ignition_length_fraction;

    // Parallelly initialize the grid points
    #pragma omp parallel for
    
        // For each grid point
        for (size_t i = 1; i < m-1; i++)
        {
            // If grid point is within the initial ignition length
            if (getXCoordinate(i) < initial_ignition_length) _temperature_array[i] = initial_ignition_temperature;
            // Set the temperature to the initial ignition temperature

            // Else set the temperature to ambient temperature
            else _temperature_array[i] = PackedPellet<real_t>::ambient_temperature;

            // Initialize the main particle at each grid point
            _particles_array[i].initializeParticle();

            // Store the present enthalpy
			_prev_enthalpy_particle[i] = _particles_array[i].getEnthalpy(_temperature_array[i]);
            
            // Initialize the particles used for
            _particles_array_raised_temperature_evolution[i].initializeParticle();
            _particles_array_const_temperature_evolution[i].initializeParticle();
        }
    	
	// Set temperature at grid point at \f$ x = L \f$ to ambient temperature
    _temperature_array[m-1] = PackedPellet<real_t>::ambient_temperature;

	updateParticlesState();
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setUpEquations()
{
    // Set up the matrix equation representing energy transport in the pellet

    // For the first row of the matrix equation, set up the boundary condition at \f$ x = 0 \f$
    // \f$ - \frac{\lambda}{\Delta x} \cdot T_1^n + \left( \frac{\lambda}{\Delta x} + \frac{D}{4} \cdot \beta_{1,0}^n \right) T_0^n 
    // = - \frac{D}{4} \cdot \beta_{0,0}^n \f$
    setUpBoundaryConditionX0();

	static const real_t one_minus_kappa			= 1. - kappa;
	static const real_t half_by_delta_x_sqr		= 0.5 / pow(delta_x, 2);
	static const real_t four_by_pellet_diameter = 4. / PackedPellet<real_t>::diameter;

    // As the 
    #pragma omp parallel for default(shared) schedule(static, 1)

        for (size_t i = 1; i < m-1; i++)
        {
            // For matrix equations corresponding to grid points inside the pellet,
            // set up the discretized energy transport equation
            // \f$ - \frac{\lambda}{\left( \Delta x \right)^2} \cdot T_{j+1}^n + 
            // \left\{ \rho_0 \alpha_{1,j}^n + \frac{2 \lambda}{\left( \Delta x \right)^2} + \beta_{1,j}^n \right\} \cdot T_j^n -
            // \frac{\lambda}{\left( \Delta x \right)^2} \cdot T_{j-1}^n
            // = - \rho_0 \alpha_{0,j}^n - \beta_{0,j}^n \f$

            // Get the linearized expression for the transient term, \f$ \alpha_{0,j}^n + \alpha_{1,j}^n \cdot T_j^n \f$
            LinearExpression<real_t> alpha = calcTransientTerm(i);

			real_t coeff_fluid = getInterstitialGasTransientTermCoefficient(i);

            // Get the linearized expression for the heat loss term, \f$ \beta_{0,j}^n + \beta_{1,j}^n \cdot T_j^n \f$
            LinearExpression<real_t> beta = calcHeatLossTerm(i);

			real_t lambda_forward_by_delta_x_sqr  = half_by_delta_x_sqr * (_thermal_conductivity[i+1] + _thermal_conductivity[i]);
			real_t lambda_backward_by_delta_x_sqr = half_by_delta_x_sqr * (_thermal_conductivity[i] + _thermal_conductivity[i-1]);

			real_t kappa_lambda_forward_by_delta_x_sqr  = kappa * lambda_forward_by_delta_x_sqr;
			real_t kappa_lambda_backward_by_delta_x_sqr = kappa * lambda_backward_by_delta_x_sqr;

            // Set up the matrix equation
            _solver.setEquation(
                // Matrix row number
                i,
                // Coefficient of \f$ T_{j-1}^n \f$, \f$ - \frac{\lambda}{\left( \Delta x \right)^2} \f$
                - kappa_lambda_backward_by_delta_x_sqr,
                // Coefficient of \f$ T_j^n \f$
                alpha.coefficient + coeff_fluid - four_by_pellet_diameter * beta.coefficient + kappa_lambda_backward_by_delta_x_sqr + kappa_lambda_forward_by_delta_x_sqr,
                // Coefficient of \f$ T_{j+1}^n \f$, \f$ - \frac{\lambda}{\left( \Delta x \right)^2} \f$
                - kappa_lambda_forward_by_delta_x_sqr,
				// Constant
				alpha.coefficient * _temperature_array[i] - alpha.constant + coeff_fluid * _temperature_array[i]
				+ four_by_pellet_diameter * (beta.constant - beta.coefficient * _temperature_array[i])
				+ one_minus_kappa * (
					lambda_forward_by_delta_x_sqr  * (_temperature_array[i+1] - _temperature_array[i]) -
					lambda_backward_by_delta_x_sqr * (_temperature_array[i] - _temperature_array[i-1])
				)
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
    _time += delta_t;
	// _solver.printMatrixEquation();
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
    for (size_t i=1; i < m-1 && flag; i++)
    {
        // If the particle is in the reaction zone then 
        // combustion in the pellet is incomplete
        flag = !inReactionZone(i);
    }
    // Return whether combustion is complete or not
    return flag;
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
    for (size_t i = 0; i < m -1; i++) output_stream << _temperature_array[i] << delimiter;
    // For the last grid point print the temperature followed by endline
    output_stream << _temperature_array[m-1] << '\n';
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
    for (size_t i = 0; i < m - 1; i++) output_stream << getXCoordinate(i) << delimiter;
    // For the last grid point print the x-coordinate followed by endline
    output_stream << getXCoordinate(m-1) << '\n';
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