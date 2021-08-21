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

#include <math.h>

#include <iostream>

#include "pde-problems/Pellet-Flame-Propagation.hpp"

#define STEFAN_BOLTZMANN_CONSTANT 5.670374419E-8 // W / m2 - K4

/**
 * @brief Get an infinitesimal increment change in the value of variable
 * 
 * @tparam real_t float, double or long double data types
 * to represent real numbers
 * @param variable Value in which incremental change is required
 * @param tol Tolerance for incremental change
 * @return real_t Value of infinitesimal change in the variable
 */
template<typename real_t>
real_t getInfinitesimalIncrement(
    real_t variable,
    real_t tol = 1E-5
) {
    // If variable is greater than 1 return the tolerance,
    // else multiply the value with tolerance and return the value
    return variable > 1.0 ? tol : variable * tol;
}

template<typename real_t> size_t PelletFlamePropagation<real_t>::_n = 0;
template<typename real_t> real_t PelletFlamePropagation<real_t>::_delta_x = 0;
template<typename real_t> real_t PelletFlamePropagation<real_t>::_delta_t = 0;
template<typename real_t> real_t PelletFlamePropagation<real_t>::_ignition_temperature = 273.15;
template<typename real_t> real_t PelletFlamePropagation<real_t>::_ignition_length = 0;

template<typename real_t>
void PelletFlamePropagation<real_t>::setGridSize(size_t N)
{
    _n = N;
    _delta_x = PackedPellet<real_t>::_length / (real_t) (_n - 1);
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setTimeStep(real_t delta_t)
{
    _delta_t = delta_t;
    CoreShellDiffusion<real_t>::setTimeStep(_delta_t);
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setIgnitionParameters(
    real_t ignition_temperature,
    real_t ignition_length
) {
    _ignition_temperature = ignition_temperature;
    _ignition_length = ignition_length;
}

template<typename real_t>
PelletFlamePropagation<real_t>::PelletFlamePropagation(
    real_t particle_volume_fraction
) : PackedPellet<real_t>(particle_volume_fraction),
    _solver(_n)
{
    _temperature_array = new real_t[_n];
    _particles_array   = new CoreShellDiffusion<real_t>[_n];

    _prev_particles_array = new CoreShellDiffusion<real_t>[_n];
    _curr_particles_array = new CoreShellDiffusion<real_t>[_n];

    initializePellet();
}

template<typename real_t>
PelletFlamePropagation<real_t>::~PelletFlamePropagation()
{
    delete [] _temperature_array;
    delete [] _particles_array;

    delete [] _prev_particles_array;
    delete [] _curr_particles_array;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setDiffusivityModel(
    ArrheniusDiffusivityModel<real_t> diffusivity_model
) {
    _diffusivity_model = diffusivity_model;
}

template<typename real_t>
real_t PelletFlamePropagation<real_t>::getParticleEnthalpyTemperatureDerivative(size_t i)
{
    real_t delta_T = getInfinitesimalIncrement(_temperature_array[i], (real_t) 1);

    real_t D = _diffusivity_model.getDiffusivity(_temperature_array[i] + delta_T);

    // #pragma omp critical    
    // {
    //     std::cout << "\t\tParticle # " << i << std::endl;
    //     std::cout << "\t\t\tDiffusivity : " << D << "\t";
    //     std::cout << "\t\t\tNormal : " << _diffusivity_model.getDiffusivity(_temperature_array[i]) << std::endl;
    // }

    _prev_particles_array[i].setUpEquations(D);
    _prev_particles_array[i].solveEquations();

    return
    (   _prev_particles_array[i].getEnthalpy(_temperature_array[i]) - 
        _particles_array[i].getEnthalpy(_temperature_array[i])
    ) / delta_T;
}

template<typename real_t>
real_t PelletFlamePropagation<real_t>::getParticleEnthalpyTimeDerivative(size_t i)
{
    real_t D = _diffusivity_model.getDiffusivity(_temperature_array[i]);

    _curr_particles_array[i].setUpEquations(D);
    _curr_particles_array[i].solveEquations();

    return
    (   _curr_particles_array[i].getEnthalpy(_temperature_array[i]) -
        _particles_array[i].getEnthalpy(_temperature_array[i])
    ) / _delta_t;
}

template<typename real_t>
LinearExpression<real_t> PelletFlamePropagation<real_t>::calcTransientTerm(size_t i)
{
    LinearExpression<real_t> expression;

    expression.coefficient = this->getHeatCapacity(_particles_array + i) / _delta_t;
    expression.constant = 0.0;

    if (inReactionZone(i))
    {
        real_t particle_enthalpy_temperature_derivative = getParticleEnthalpyTemperatureDerivative(i);
        real_t particle_enthalpy_time_derivative = getParticleEnthalpyTimeDerivative(i);
        
        // std::cout << "\tParticle # " << i << std::endl;
        // std::cout << "\t\tdudT : " << particle_enthalpy_temperature_derivative / _delta_t << std::endl;
        // std::cout << "\t\tdudt : " << particle_enthalpy_time_derivative << std::endl;

        expression.coefficient += this->_particle_mass_fractions * particle_enthalpy_temperature_derivative / _delta_t;
        expression.constant += - this->_particle_mass_fractions * particle_enthalpy_time_derivative;

        // expression.constant += - this->_particle_mass_fractions * (_curr_particles_array[i].getEnthalpy(_temperature_array[i]) - _prev_particles_array[i].getEnthalpy(_temperature_array[i])) / _delta_t;
    }

    expression.constant += expression.coefficient * _temperature_array[i];
    
    return expression;
}

template<typename real_t>
LinearExpression<real_t> PelletFlamePropagation<real_t>::calcHeatLossTerm(size_t i)
{
    LinearExpression<real_t> expression;

    expression.coefficient =
        this->_convective_heat_transfer_coefficient + 
        4.0 * STEFAN_BOLTZMANN_CONSTANT * this->_radiative_emissivity * pow(_temperature_array[i], 3);

    expression.constant = 
        this->_convective_heat_transfer_coefficient * this->_ambient_temperature +
        STEFAN_BOLTZMANN_CONSTANT * this->_radiative_emissivity * (
            3.0 * pow(_temperature_array[i], 4) - 
            pow(this->_ambient_temperature, 4)
        );

    return expression;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setUpBoundaryConditionX0()
{
    LinearExpression<real_t> heat_loss_term = calcHeatLossTerm(0);
    
    real_t lambda_by_delta_x = this->getHeatConductivity(_particles_array + 1) / _delta_x;

    _solver.setEquationFirstRow(
        heat_loss_term.coefficient + lambda_by_delta_x,
        - lambda_by_delta_x,
        heat_loss_term.constant
    );
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setUpBoundaryConditionXN()
{
    LinearExpression<real_t> heat_loss_term = calcHeatLossTerm(_n-1);

    real_t lambda_by_delta_x = this->getHeatConductivity(_particles_array + _n - 2) / _delta_x;

    _solver.setEquationLastRow(
        - lambda_by_delta_x,
        heat_loss_term.coefficient + lambda_by_delta_x,
        heat_loss_term.constant
    );
}

template<typename real_t>
void PelletFlamePropagation<real_t>::setUpEquations()
{
    setUpBoundaryConditionX0();

    #pragma omp parallel for

        for (size_t i = 1; i < _n-1; i++)
        {
            LinearExpression<real_t> transient_term = calcTransientTerm(i);

            LinearExpression<real_t> heat_loss_term = calcHeatLossTerm(i);

            real_t minus_lambda_by_delta_x_sqr = - this->getHeatConductivity(_particles_array + i) / pow(_delta_x, 2);

            _solver.setEquation( i,
                minus_lambda_by_delta_x_sqr,
                this->_density * transient_term.coefficient + 4.0 * heat_loss_term.coefficient / this->_diameter - 2.0 * minus_lambda_by_delta_x_sqr,
                minus_lambda_by_delta_x_sqr,
                this->_density * transient_term.constant + 4.0 * heat_loss_term.constant / this->_diameter
            );
        }

    setUpBoundaryConditionXN();

    // _solver.printMatrixEquation();
}

template<typename real_t>
void PelletFlamePropagation<real_t>::updateParticlesState()
{
    #pragma omp parallel for

        for (size_t i = 1; i < _n-1; i++)
        {
            _particles_array[i].copyTo(_prev_particles_array[i]);

            _particles_array[i].setUpEquations(_diffusivity_model.getDiffusivity(_temperature_array[i]));
            _particles_array[i].solveEquations();

            _particles_array[i].copyTo(_curr_particles_array[i]);
        }
}

template<typename real_t>
void PelletFlamePropagation<real_t>::solveEquations()
{
    _solver.getSolution(_temperature_array);

    updateParticlesState();
}

template<typename real_t>
void PelletFlamePropagation<real_t>::printTemperatureProfile(
    std::ostream &output_stream,
    char delimiter
) {
    for (size_t i = 0; i < _n -1; i++) output_stream << _temperature_array[i] << delimiter;

    output_stream << _temperature_array[_n-1] << std::endl;
}

template<typename real_t>
bool PelletFlamePropagation<real_t>::isCombustionComplete()
{
    bool flag = true;

    for (size_t i=1; i < _n-1 && flag; i++)
    {
        if (!_particles_array[i].isCombustionComplete()) 
        {
            flag = false;
        }
    }

    return flag;
}

template<typename real_t>
void PelletFlamePropagation<real_t>::initializePellet()
{
    _temperature_array[0] = 1200;

    #pragma omp parallel for
    
        for (size_t i = 1; i < _n-1; i++)
        {
            if (getxCoordinate(i) < _ignition_length) _temperature_array[i] = 1200;

            else _temperature_array[i] = this->_ambient_temperature;

            _particles_array[i].initializeParticle();
            
            _curr_particles_array[i].initializeParticle();
            _prev_particles_array[i].initializeParticle();
        }

    _temperature_array[_n-1] = this->_ambient_temperature;
}

template class PelletFlamePropagation<float>;
template class PelletFlamePropagation<double>;
template class PelletFlamePropagation<long double>;