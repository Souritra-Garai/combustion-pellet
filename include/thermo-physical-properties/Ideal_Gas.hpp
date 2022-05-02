/**
 * @file IdealGas.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef __IDEAL_GAS__
#define __IDEAL_GAS__

// Required for printing to file / screen
#include <ostream>

#include "thermo-physical-properties/Enthalpy.hpp"
#include "thermo-physical-properties/Thermal_Conductivity.hpp"

template<typename real_t>
class IdealGas
{
	private :

		const real_t _molar_mass;
		const real_t _gamma;

		Enthalpy<real_t> &_enthalpy;
		ThermalConductivityQuadraticPolynomial<real_t> &_thermal_conductivity;

	public :

		IdealGas(
			real_t molar_mass,
			real_t gamma,
			Enthalpy<real_t> &enthalpy,
			ThermalConductivityQuadraticPolynomial<real_t> &thermal_conductivity
		) :	_enthalpy(enthalpy),
			_thermal_conductivity(thermal_conductivity),
			_molar_mass(molar_mass),
			_gamma(gamma)
		{
			;
		}

		IdealGas() : _enthalpy(), _thermal_conductivity(), _molar_mass(1.0), _gamma(1.4)
		{;}

		inline real_t getDensity(real_t temperature, real_t pressure = 1.01325E5)
		{
			return _molar_mass * pressure / (8.314 * temperature);
		}

		inline real_t getEnthalpy(real_t temperature)
		{
			return _enthalpy.getStandardEnthalpy(temperature) / _molar_mass;
		}

		inline real_t getCp(real_t temperature)
		{
			return _enthalpy.getHeatCapacity(temperature) / _molar_mass;
		}

		inline real_t getThermalConductivity(real_t temperature)
		{
			return _thermal_conductivity.getThermalConductivity(temperature);
		}

		inline real_t getMolarMass()  { return _molar_mass; }

		inline real_t getMolarDensity(real_t temperature, real_t pressure = 1.01325E5)
        { 
            return pressure / (8.314 * temperature);
        }
};

#endif