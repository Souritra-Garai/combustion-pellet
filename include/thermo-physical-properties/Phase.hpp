/**
 * @file Phase.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief 
 * @version 0.1
 * @date 2021-09-17
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __PHASE__
#define __PHASE__

#include "math/Data-Type.hpp"
#include "math/Sigmoid-Function.hpp"

#include "math/Shomate-Expression.hpp"
#include "math/Quadratic-Expression.hpp"

class Phase
{
    private:

        ShomateExpression	_enthalpy;
        QuadraticExpression _thermal_conductivity;
    
        const real_t _density;
        
        const real_t _temperature_lower_bound;
        const real_t _temperature_upper_bound;

    public:

		static const real_t sharpness_coefficient;

        Phase(
            real_t density,
            ShomateExpression	&enthalpy,
            QuadraticExpression &thermal_conductivity,
            real_t temperature_lower_bound = 0,
            real_t temeprature_upper_bound = INFINITY
        ) : _enthalpy(enthalpy),
            _thermal_conductivity(thermal_conductivity),
			_density(density),
			_temperature_lower_bound(temperature_lower_bound),
			_temperature_upper_bound(temeprature_upper_bound)
        {
			;
        }

        inline real_t getDensity(real_t temperature) const 
        {
            return _density * (
                getSigmoid(temperature, _temperature_lower_bound, sharpness_coefficient) -
                getSigmoid(temperature, _temperature_upper_bound, sharpness_coefficient)
            );
        }

        inline real_t getStandardEnthalpy(real_t temperature) const
        {
            return
			ShomateExpression::
			normalizeIntegralOutput(
				_enthalpy.evaluateExpressionIntegral(
					ShomateExpression::normalizeInput(temperature)))
			* ( getSigmoid(temperature, _temperature_lower_bound, sharpness_coefficient) -
                getSigmoid(temperature, _temperature_upper_bound, sharpness_coefficient) );
        }

        inline real_t getHeatCapacity(real_t temperature) const
        {
            return
			_enthalpy.evaluateExpression(ShomateExpression::normalizeInput(temperature))
			* ( getSigmoid(temperature, _temperature_lower_bound, sharpness_coefficient) -
                getSigmoid(temperature, _temperature_upper_bound, sharpness_coefficient) ) +			
			ShomateExpression::
			normalizeIntegralOutput(
				_enthalpy.evaluateExpressionIntegral(
					ShomateExpression::normalizeInput(temperature)))
			* ( getSigmoidDerivative(temperature, _temperature_lower_bound, sharpness_coefficient) -
                getSigmoidDerivative(temperature, _temperature_upper_bound, sharpness_coefficient) );
        }

        inline real_t getThermalConductivity(real_t temperature) const
        {
            return _thermal_conductivity.evaluateExpression(temperature) * (
                getSigmoid(temperature, _temperature_lower_bound, sharpness_coefficient) -
                getSigmoid(temperature, _temperature_upper_bound, sharpness_coefficient)
            );
        }
};

#endif