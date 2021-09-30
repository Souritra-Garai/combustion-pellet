/**
 * @file Phase.hpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-09-17
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __PHASE__
#define __PHASE__

#include <math.h>

#include "thermo-physical-properties/Enthalpy.hpp"
#include "thermo-physical-properties/Thermal_Conductivity.hpp"

template<typename real_t>
class Phase
{
    private:

        Enthalpy<real_t> &_enthalpy;

        ThermalConductivityQuadraticPolynomial<real_t> &_thermal_conductivity;
    
        real_t _density;
        
        real_t _temperature_lower_bound;
        real_t _temperature_upper_bound;

        real_t _sharpness_coefficient;

        static real_t _getSigmoid(real_t x, real_t origin = 0, real_t scale = 0)
        {
            return 1 / (1 + exp( - scale * (x - origin)));
        }

        static real_t _getSigmoidDerivative(real_t x, real_t origin = 0, real_t scale = 0)
        {
            real_t sigma = _getSigmoid(x, origin, scale);

            return scale * sigma * (1 - sigma);
        }

    public:

        Phase(
            real_t density,
            Enthalpy<real_t> &enthalpy,
            ThermalConductivityQuadraticPolynomial<real_t> &thermal_conductivity,
            real_t temperature_lower_bound = 0,
            real_t temeprature_upper_bound = INFINITY,
            real_t sharpness_coefficient = 1
        ) : _enthalpy(enthalpy),
            _thermal_conductivity(thermal_conductivity)
        {
            _density = density;

            _temperature_lower_bound = temperature_lower_bound;
            _temperature_upper_bound = temeprature_upper_bound;

            _sharpness_coefficient = sharpness_coefficient;
        }

        real_t getDensity(real_t temperature) 
        {
            return _density * (
                _getSigmoid(temperature, _temperature_lower_bound, _sharpness_coefficient) -
                _getSigmoid(temperature, _temperature_upper_bound, _sharpness_coefficient)
            );
        }

        real_t getStandardEnthalpy(real_t temperature)
        {
            return _enthalpy.getStandardEnthalpy(temperature) * (
                _getSigmoid(temperature, _temperature_lower_bound, _sharpness_coefficient) -
                _getSigmoid(temperature, _temperature_upper_bound, _sharpness_coefficient)
            );
        }

        real_t getHeatCapacity(real_t temperature)
        {
            return _enthalpy.getHeatCapacity(temperature) * (
                _getSigmoid(temperature, _temperature_lower_bound, _sharpness_coefficient) -
                _getSigmoid(temperature, _temperature_upper_bound, _sharpness_coefficient)
            ) + _enthalpy.getStandardEnthalpy(temperature) * (
                _getSigmoidDerivative(temperature, _temperature_lower_bound, _sharpness_coefficient) -
                _getSigmoidDerivative(temperature, _temperature_upper_bound, _sharpness_coefficient)
            );
        }

        real_t getThermalConductivity(real_t temperature)
        {
            return _thermal_conductivity.getThermalConductivity(temperature) * (
                _getSigmoid(temperature, _temperature_lower_bound, _sharpness_coefficient) -
                _getSigmoid(temperature, _temperature_upper_bound, _sharpness_coefficient)
            );
        }
};

#endif