#ifndef __CORE_SHELL_PARTICLE__
#define __CORE_SHELL_PARTICLE__

#include "math/Data-Type.hpp"
#include "thermo-physical-properties/Condensed-Species.hpp"
#include "thermo-physical-properties/Arrhenius-Diffusivity-Model.hpp"

class CoreShellParticle
{
	protected :
	
		real_t _mass_fraction_core_material;
		real_t _mass_fraction_shell_material;
		real_t _mass_fraction_product_material;

	public :

		static CondensedSpecies core_species;
		static CondensedSpecies shell_species;
		static CondensedSpecies product_species;

		static ArrheniusDiffusivityModel diffusivity_model;

		static const real_t overall_radius;
		static const real_t core_radius;

		static const real_t mass;

		CoreShellParticle();

		// Input temperature T in K
		// Returns density at temperature T and standard pressure in kg/m^3
		inline real_t getDensity(real_t temperature) const
		{
			real_t volume_core_material		= _mass_fraction_core_material	  / core_species.getDensity(temperature);
			real_t volume_shell_material	= _mass_fraction_shell_material	  / shell_species.getDensity(temperature);
			real_t volume_product_material	= _mass_fraction_product_material / product_species.getDensity(temperature);

			return 1. / (
				volume_core_material +
				volume_shell_material +
				volume_product_material
			);
		}

		// Input temperature T in K
		// Returns specific enthalpy at temperature T and standard pressure in J/kg
		inline real_t getEnthalpy(real_t temperature) const
		{
			return
				_mass_fraction_core_material	* core_species.getEnthalpy(temperature) +
				_mass_fraction_shell_material	* shell_species.getEnthalpy(temperature) +
				_mass_fraction_product_material	* product_species.getEnthalpy(temperature);
		}

		// Input temperature T in K
		// Returns specific heat capacity at temperature T and standard pressure in J/kg-K
		inline real_t getHeatCapacity(real_t temperature) const
		{
			return 
				_mass_fraction_core_material	* core_species.getHeatCapacity(temperature) +
				_mass_fraction_shell_material	* shell_species.getHeatCapacity(temperature) +
				_mass_fraction_product_material	* product_species.getHeatCapacity(temperature); 
		}

		// Input temperature T in K
		// Returns thermal conductivity at temperature T and standard pressure in W/m-K
		inline real_t getThermalConductivity(real_t temperature) const
		{
			real_t volume_fraction_core_material	= _mass_fraction_core_material	  / core_species.getDensity(temperature);
			real_t volume_fraction_shell_material	= _mass_fraction_shell_material	  / shell_species.getDensity(temperature);
			real_t volume_fraction_product_material	= _mass_fraction_product_material / product_species.getDensity(temperature);

			real_t sum = volume_fraction_core_material + volume_fraction_shell_material + volume_fraction_product_material;

			volume_fraction_core_material	 /= sum;
			volume_fraction_shell_material	 /= sum;
			volume_fraction_product_material /= sum;

			return
				volume_fraction_core_material	 * core_species.getThermalConductivity(temperature)     +
				volume_fraction_shell_material	 * shell_species.getThermalConductivity(temperature)    +
				volume_fraction_product_material * product_species.getThermalConductivity(temperature)  ;
		}

		inline real_t getMassFractionsCoreMaterial() const
		{
			return _mass_fraction_core_material;
		}

		inline real_t getMassFractionsShellMaterial() const
		{
			return _mass_fraction_shell_material;
		}

		inline real_t getMassFractionsProductMaterial() const
		{
			return _mass_fraction_product_material;
		}

		inline bool isCombustionComplete(real_t tolerance = 1E-3) const
		{
			return _mass_fraction_core_material < tolerance || _mass_fraction_shell_material < tolerance;
		}

};

#endif