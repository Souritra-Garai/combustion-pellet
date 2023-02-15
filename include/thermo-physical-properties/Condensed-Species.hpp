/**
 * @file Condensed-Species.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header defines a class to represent pure substances
 * with their thermodynamic properties like density, heat capacity etc.
 * @version 0.1
 * @date 2021-07-06
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __CONDENSED_SPECIES__
#define __CONDENSED_SPECIES__

#include <cstdlib>
#include <cstring>

#include "math/Data-Type.hpp"
#include "thermo-physical-properties/Phase.hpp"

class CondensedSpecies
{
    private :

        Phase *_phases;

    public :

        const unsigned int num_phases;
        const real_t molar_mass;

        CondensedSpecies(
            unsigned int number_of_phases,
            Phase array_of_phases[],
            real_t molar_mass
        ) : num_phases(number_of_phases),
			molar_mass(molar_mass) 
		{
			_phases = (Phase*) std::malloc(number_of_phases * sizeof(Phase));
			std::memcpy(_phases, array_of_phases, number_of_phases * sizeof(Phase));
        }

        ~CondensedSpecies()
		{
			std::free(_phases);
		}

        inline real_t getDensity(real_t temperature) const
        {
            real_t density = 0;
            
			#pragma unroll

				for (unsigned int i = 0; i < num_phases; i++)
				{
					density += _phases[i].getDensity(temperature);
				}

            return density;
        }
        
        inline real_t getMolarDensity(real_t temperature) const
        { 
            return getDensity(temperature) / molar_mass;
        }

        inline real_t getHeatCapacity(real_t temperature) const
        {
            real_t heat_capacity = 0;

			#pragma unroll

				for (unsigned int i = 0; i < num_phases; i++)
				{
					heat_capacity += _phases[i].getHeatCapacity(temperature);
				}

            return heat_capacity / molar_mass;
        }

        inline real_t getEnthalpy(real_t temperature) const
        {
            real_t enthalpy = 0;

			#pragma unroll

				for (unsigned int i = 0; i < num_phases; i++)
				{
					enthalpy += _phases[i].getStandardEnthalpy(temperature);
				}

            return enthalpy / molar_mass;
        }

        inline real_t getThermalConductivity(real_t temperature) const
        { 
            real_t heat_conductivity = 0;

			#pragma unroll

				for (unsigned int i = 0; i < num_phases; i++)
				{
					heat_conductivity += _phases[i].getThermalConductivity(temperature);
				}

            return heat_conductivity;
        }
};

#endif