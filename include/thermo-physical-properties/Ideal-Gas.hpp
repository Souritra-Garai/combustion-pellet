/**
 * @file Ideal-Gas.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief 
 * @version 2.0
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __IDEAL_GAS__
#define __IDEAL_GAS__

#include "math/Data-Type.hpp"

#include "math/Shomate-Expression.hpp"
#include "math/Quadratic-Expression.hpp"

class IdealGas
{
	private :

		ShomateExpression	_enthalpy;
		QuadraticExpression	_thermal_conductivity;

		const real_t molar_mass_inverse;

	public :

		const real_t molar_mass;
		const real_t gamma;

		IdealGas(
			real_t molar_mass,
			real_t gamma,
			ShomateExpression	&enthalpy,
			QuadraticExpression	&thermal_conductivity
		) :	_enthalpy(enthalpy),
			_thermal_conductivity(thermal_conductivity),
			molar_mass(molar_mass),
			gamma(gamma),
			molar_mass_inverse(1./molar_mass)
		{
			;
		}

		inline real_t getMolarDensity(real_t temperature, real_t pressure = 1.01325E5) const
        { 
            return (1./8.314) * pressure / temperature;
        }

		inline real_t getDensity(real_t temperature, real_t pressure = 1.01325E5) const
		{
			return molar_mass * getMolarDensity(temperature, pressure);
		}

		inline real_t getEnthalpy(real_t temperature) const
		{
			return 
			ShomateExpression::
			normalizeIntegralOutput(
				_enthalpy.evaluateExpressionIntegral(
					ShomateExpression::normalizeInput(temperature))) * molar_mass_inverse;
		}

		inline real_t getCp(real_t temperature) const
		{
			return _enthalpy.evaluateExpression(ShomateExpression::normalizeInput(temperature)) * molar_mass_inverse;
		}

		inline real_t getThermalConductivity(real_t temperature) const
		{
			return _thermal_conductivity.evaluateExpression(temperature);
		}
};

#endif