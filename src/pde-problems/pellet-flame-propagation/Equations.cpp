#include "pde-problems/Pellet-Flame-Propagation.hpp"

#define STEFAN_BOLTZMANN_CONSTANT 5.670374419E-8 // W / m2 - K4

inline
bool PelletFlamePropagation::inReactionZone(size_t index)
{
	return !_particles_array[index].isCombustionComplete();
}

bool PelletFlamePropagation::isCombustionComplete()
{
	bool flag = true;
	
	for (size_t i=1; i < m-1 && flag; i++)
	{
		flag = !inReactionZone(i);
	}
    
	return flag;
}

void PelletFlamePropagation::evolveParticleForEnthalpyDerivative(size_t i)
{
	_particles_array_const_temperature_evolution[i].setUpEquations(
		_temperature_array[i],
		_particles_array[i]
	);

	_particles_array_const_temperature_evolution[i].solveEquations();

	_particles_array_raised_temperature_evolution[i].setUpEquations(
		_temperature_array[i] + delta_T,
		_particles_array[i]
	);

	_particles_array_raised_temperature_evolution[i].solveEquations();
}

inline LinearExpression PelletFlamePropagation::calcTransientTerm(size_t i)
{
	static const real_t delta_t_inverse				= 1. / delta_t;
	static const real_t gamma_by_delta_t			= gamma / delta_t;
	static const real_t one_minus_gamma_by_delta_t	= (1. - gamma) / delta_t;
	static const real_t gamma_by_delta_T_delta_t	= gamma / (delta_t * delta_T);
    
	LinearExpression expression;
	
	expression.a_1 = _particles_array[i].getHeatCapacity(_temperature_array[i]) * delta_t_inverse;
    
    expression.a_0 = 0.0;

    if (inReactionZone(i))
    {
		evolveParticleForEnthalpyDerivative(i);

		real_t enthalpy						= _particles_array[i].getEnthalpy(_temperature_array[i]);
		real_t enthalpy_const_T_evolution	= _particles_array_const_temperature_evolution[i].getEnthalpy(_temperature_array[i]);
		real_t enthalpy_raised_T_evolution	= _particles_array_raised_temperature_evolution[i].getEnthalpy(_temperature_array[i]);
        
        expression.a_1	+=	gamma_by_delta_T_delta_t * (enthalpy_raised_T_evolution - enthalpy_const_T_evolution);

        expression.a_0	+=	one_minus_gamma_by_delta_t * (enthalpy - _prev_enthalpy_particle[i])
						+ 	gamma_by_delta_t * (enthalpy_const_T_evolution - enthalpy);

		_prev_enthalpy_particle[i] = enthalpy;
    }

	expression *= PackedPellet::overall_particle_density;
    
    return expression;
}

inline LinearExpression PelletFlamePropagation::calcHeatLossTerm(size_t i)
{
	static const real_t constant1 = PackedPellet::radiative_emissivity * STEFAN_BOLTZMANN_CONSTANT;
	static const real_t constant2 = 4. * PackedPellet::radiative_emissivity * STEFAN_BOLTZMANN_CONSTANT;

    LinearExpression expression;

    expression.a_0 =
		PackedPellet::convective_heat_transfer_coefficient_curved_surface * (_temperature_array[i] - PackedPellet::ambient_temperature) +
		constant1 * (pow(_temperature_array[i], 4) - pow(PackedPellet::ambient_temperature, 4))
	;

    expression.a_1 = gamma * (
		PackedPellet::convective_heat_transfer_coefficient_curved_surface +
		constant2 * pow(_temperature_array[i], 3)
	);

	return expression;
}

inline real_t PelletFlamePropagation::getInterstitialGasTransientTermCoefficient(size_t i)
{
	static const real_t constant = PackedPellet::interstitial_volume_fractions / delta_t;

	return
		constant * 
		PackedPellet::interstitial_gas.getDensity(_temperature_array[i]) *
		PackedPellet::interstitial_gas.getCp(_temperature_array[i]);
}

inline void PelletFlamePropagation::setUpBoundaryConditionX0()
{
	static const real_t constant = 1. / delta_x;
    
    LinearExpression beta = calcHeatLossTerm(0);
    
    real_t lambda_by_delta_x = _thermal_conductivity[0] * constant;
    
    _solver.setEquationFirstRow(
          beta.a_1 + lambda_by_delta_x,
        - lambda_by_delta_x,
		  beta.evaluateExpression(_temperature_array[0])
    );
}

inline void PelletFlamePropagation::setUpBoundaryConditionXN()
{
	static const real_t constant = 1. / delta_x;
    
	LinearExpression beta = calcHeatLossTerm(m-1);

    real_t lambda_by_delta_x = _thermal_conductivity[m - 1] * constant;
    
	_solver.setEquationLastRow(
    	- lambda_by_delta_x,
		beta.a_1 + lambda_by_delta_x,
		beta.evaluateExpression(_temperature_array[m - 1])
    );
}

void PelletFlamePropagation::setUpEquations()
{
	setUpBoundaryConditionX0();

	static const real_t one_minus_kappa			= 1.0 - kappa;
	static const real_t half_by_delta_x_sqr		= 0.5 / (delta_x * delta_x);
	static const real_t four_by_pellet_diameter = 4.0 / PackedPellet::diameter;
	
	#pragma omp parallel for default(shared) schedule(static, 1)

		for (size_t i = 1; i < m-1; i++)
		{
			LinearExpression alpha = calcTransientTerm(i);

			LinearExpression beta = calcHeatLossTerm(i);

			real_t coeff_fluid = getInterstitialGasTransientTermCoefficient(i);

			real_t lambda_forward_by_delta_x_sqr  = half_by_delta_x_sqr * (_thermal_conductivity[i+1] + _thermal_conductivity[i]);
			real_t lambda_backward_by_delta_x_sqr = half_by_delta_x_sqr * (_thermal_conductivity[i] + _thermal_conductivity[i-1]);

			real_t kappa_lambda_forward_by_delta_x_sqr  = kappa * lambda_forward_by_delta_x_sqr;
			real_t kappa_lambda_backward_by_delta_x_sqr = kappa * lambda_backward_by_delta_x_sqr;

            _solver.setEquation(
				i,
				- kappa_lambda_backward_by_delta_x_sqr,
				
				  alpha.a_1 + coeff_fluid - four_by_pellet_diameter * beta.a_1
				+ kappa_lambda_backward_by_delta_x_sqr + kappa_lambda_forward_by_delta_x_sqr,
				
				- kappa_lambda_forward_by_delta_x_sqr,
				
				  alpha.a_1 * _temperature_array[i] - alpha.a_0 + coeff_fluid * _temperature_array[i]
				+ four_by_pellet_diameter * (beta.a_0 - beta.a_1 * _temperature_array[i])
				+ one_minus_kappa * (
					lambda_forward_by_delta_x_sqr  * (_temperature_array[i+1] - _temperature_array[i]) -
					lambda_backward_by_delta_x_sqr * (_temperature_array[i] - _temperature_array[i-1])
				)
            );
        }
	
	setUpBoundaryConditionXN();
}

void PelletFlamePropagation::solveEquations()
{
	_time += delta_t;
	
	_solver.getSolution(_temperature_array);
	
	updateParticles();
}

void PelletFlamePropagation::updateParticles()
{
	#pragma omp parallel for default(shared) schedule(static, 1)
        
		for (size_t i = 1; i < m-1; i++)
		{
			if (inReactionZone(i))
			{
				_particles_array[i].setUpEquations(_temperature_array[i]);
				_particles_array[i].solveEquations();
			}
			
			_thermal_conductivity[i] = PackedPellet::getThermalConductivity(_particles_array + i, _temperature_array[i]);
		}

	_thermal_conductivity[0]	= PackedPellet::getThermalConductivity(_particles_array + 1, 		_temperature_array[0]);
	_thermal_conductivity[m-1]	= PackedPellet::getThermalConductivity(_particles_array + (m - 2),_temperature_array[m-1]);
}
