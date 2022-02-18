/**
 * @file Substance.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header defines a class to represent pure substances
 * with their thermodynamic properties like density, heat capacity etc.
 * @version 0.1
 * @date 2021-07-06
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __SUBSTANCE__
#define __SUBSTANCE__

// Required for printing to file / screen
#include <ostream>

#include <thermo-physical-properties/Phase.hpp>

template <typename real_t>
class Substance
{
    private :

        unsigned int _num_phases;

        Phase<real_t> *_phases;

        real_t _molar_mass;

    public :

        Substance(
            unsigned int number_of_phases,
            Phase<real_t> array_of_phases[],
            real_t molar_mass
        ) {
            _num_phases = number_of_phases;
            _phases = array_of_phases;
            
            _molar_mass = molar_mass;
        }

        Substance() {;}

        real_t getMolarMass()  { return _molar_mass; }
        
        real_t getDensity(real_t temperature)
        {
            real_t density = 0;
            
            for (unsigned int i = 0; i < _num_phases; i++)
            {
                density += _phases[i].getDensity(temperature);
            }

            return density;
        }
        
        real_t getMolarDensity(real_t temperature)
        { 
            return getDensity(temperature) / _molar_mass;
        }

        real_t getHeatCapacity(real_t temperature)
        {
            real_t heat_capacity = 0;

            for (unsigned int i = 0; i < _num_phases; i++)
            {
                heat_capacity += _phases[i].getHeatCapacity(temperature);
            }

            return heat_capacity / _molar_mass;
        }

        real_t getEnthalpy(real_t temperature)
        {
            real_t enthalpy = 0;

            for (unsigned int i = 0; i < _num_phases; i++)
            {
                enthalpy += _phases[i].getStandardEnthalpy(temperature);
            }

            return enthalpy / _molar_mass;
        }

        real_t getThermalConductivity(real_t temperature)
        { 
            real_t heat_conductivity = 0;

            for (unsigned int i = 0; i < _num_phases; i++)
            {
                heat_conductivity += _phases[i].getThermalConductivity(temperature);
            }

            return heat_conductivity;
        }
};

#endif