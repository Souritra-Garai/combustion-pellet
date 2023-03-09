#ifndef __CONDENSED_SPECIES__
#define __CONDENSED_SPECIES__

#include <cstdlib>	// malloc & free
#include <cstring>	// memcpy

#include "math/Data-Type.hpp"
#include "thermo-physical-properties/Phase.hpp"

class CondensedSpecies
{
	private :

		Phase *_phases;

		const real_t _molar_mass_inverse;

	public :

		const unsigned int num_phases;
		const real_t molar_mass;

		CondensedSpecies(
			unsigned int number_of_phases,
			Phase array_of_phases[],
			real_t molar_mass
		) : num_phases(number_of_phases),
			molar_mass(molar_mass),
			_molar_mass_inverse(1./molar_mass)
		{
			_phases = (Phase*) std::malloc(number_of_phases * sizeof(Phase));
			std::memcpy(_phases, array_of_phases, number_of_phases * sizeof(Phase));
		}

		~CondensedSpecies()
		{
			std::free(_phases);
		}

		// Input temperature T in K
		// Returns density at temperature T and standard pressure in kg/m^3
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
		
		// Input temperature T in K
		// Returns molar density at temperature T and standard pressure in mol./m^3
		inline real_t getMolarDensity(real_t temperature) const
		{ 
			return getDensity(temperature) * _molar_mass_inverse;
		}

		// Input temperature T in K
		// Returns specific heat capacity at temperature T and standard pressure in J/kg-K
		inline real_t getHeatCapacity(real_t temperature) const
		{
			real_t heat_capacity = 0;

			#pragma unroll

				for (unsigned int i = 0; i < num_phases; i++)
				{
					heat_capacity += _phases[i].getHeatCapacity(temperature);
				}

			return heat_capacity * _molar_mass_inverse;
		}

		// Input temperature T in K
		// Returns specific enthalpy at temperature T and standard pressure in J/kg
		inline real_t getEnthalpy(real_t temperature) const
		{
			real_t enthalpy = 0;

			#pragma unroll

				for (unsigned int i = 0; i < num_phases; i++)
				{
					enthalpy += _phases[i].getStandardEnthalpy(temperature);
				}

			return enthalpy * _molar_mass_inverse;
		}

		// Input temperature T in K
		// Returns thermal conductivity at temperature T and standard pressure in W/m-K
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