/**
 * @file Core-Shell-Particle.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header defines a class for core shell particle
 * @version 0.1
 * @date 2021-07-08
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __CORE_SHELL_PARTICLE__
#define __CORE_SHELL_PARTICLE__

// Required for Substance class
#include "thermo-physical-properties/Condensed-Species.hpp"

#include "thermo-physical-properties/Arrhenius-Diffusivity-Model.hpp"

/**
 * @brief Class to represent Core-Shell Particle with functions to
 * estimate thermodynamic properties for varying composition 
 * and temperature
 * @tparam real_t float, double or long double data types
 * to represent real numbers
 */
template<typename real_t>
class CoreShellParticle
{
    protected :

        /**
         * @brief Mass fraction of substance present at
         * the core of the core-shell particle
         */
        real_t _mass_fraction_core_material;
        /**
         * @brief Mass fraction of substance present in
         * the shell of the core-shell particle
         */
        real_t _mass_fraction_shell_material;
        /**
         * @brief Mass fraction of substance formed as
         * a product of the reaction between core and
         * shell substances
         */
        real_t _mass_fraction_product_material;

        /**
         * @brief Substance present in the core of 
         * the core-shell particle
         */
        static CondensedSpecies<real_t> _core_species;
        /**
         * @brief Substance present in the shell of
         * the core-shell particle
         */
        static CondensedSpecies<real_t> _shell_species;
        /**
         * @brief Substance formed upon reaction of
         * the core and shell substances
         */
        static CondensedSpecies<real_t> _product_species;

        /**
         * @brief Overall radius of the Core-Shell Particle
         */
        static const real_t _overall_radius;
        /**
         * @brief Core radius of the Core-Shell Particle
         */
        static const real_t _core_radius;

        /**
         * @brief Mass of the Core-Shell Particle
         */
        static const real_t _mass;

		static ArrheniusDiffusivityModel<real_t> _diffusivity_model;

    public :

        /**
         * @brief Function to calculate volume of the core
         * of the core-shell particle
         * @return real_t Volume of the core in \f$ m^3 \f$
         */
        static real_t calcCoreVolume();
        /**
         * @brief Function to calculate volume of the shell
         * of the core-shell particle
         * @return real_t Volume of the shell in \f$ m^3 \f$
         */
        static real_t calcShellVolume();

        /**
         * @brief Function to calculate the mass of the core
         * of the core-shell particle
         * @return real_t Mass of the core in kg
         */
        static real_t calcCoreMass();
        /**
         * @brief Function to calculate the mass of the shell
         * of the core-shell particle
         * @return real_t Mass of the shell in kg
         */
        static real_t calcShellMass();
        /**
         * @brief Function to calculate the mass of the
         * core-shell particle
         * @return real_t Mass of the core-shell particle in kg
         */
        static real_t calcParticleMass();

        /**
         * @brief Construct a new Core Shell Combustion Particle
         */
        CoreShellParticle();

        /**
         * @brief Get the Density of the Particle
         * 
         * @return real_t Density of the particle
         */
        inline real_t getDensity(real_t temperature)
		{
			real_t volume_fraction_core_material    = _mass_fraction_core_material      / _core_species.getDensity(temperature);
			real_t volume_fraction_shell_material   = _mass_fraction_shell_material     / _shell_species.getDensity(temperature);
			real_t volume_fraction_product_material = _mass_fraction_product_material   / _product_species.getDensity(temperature);

			return 1 / (
				volume_fraction_core_material +
				volume_fraction_shell_material +
				volume_fraction_product_material
			);
		}

        /**
         * @brief Get the Heat Capacity of the particle
         * 
         * @return real_t Heat capacity of the particle
         */
        inline real_t getHeatCapacity(real_t temperature)
		{
			// \f$ c = \sum_{k \in \left\{ Core, Shell, Product \right}} Y_k c_k \f$
			return 
				_mass_fraction_core_material    * _core_species.getHeatCapacity(temperature) +
				_mass_fraction_shell_material   * _shell_species.getHeatCapacity(temperature) +
				_mass_fraction_product_material * _product_species.getHeatCapacity(temperature); 
		}

        /**
         * @brief Get the Heat Conductivity of the particle
         * 
         * @return real_t Heat capacity of the particle
         */
        inline real_t getThermalConductivity(real_t temperature)
		{
			// \f$ \lambda = \sum_{k \in \left\{ Core, Shell, Product \right}} Y_k \lambda_k \f$

			real_t volume_fraction_core_material    = _mass_fraction_core_material      / _core_species.getDensity(temperature);
			real_t volume_fraction_shell_material   = _mass_fraction_shell_material     / _shell_species.getDensity(temperature);
			real_t volume_fraction_product_material = _mass_fraction_product_material   / _product_species.getDensity(temperature);

			real_t sum = volume_fraction_core_material + volume_fraction_shell_material + volume_fraction_product_material;

			volume_fraction_core_material       /= sum;
			volume_fraction_shell_material      /= sum;
			volume_fraction_product_material    /= sum;

			return
				volume_fraction_core_material       * _core_species.getThermalConductivity(temperature)     +
				volume_fraction_shell_material      * _shell_species.getThermalConductivity(temperature)    +
				volume_fraction_product_material    * _product_species.getThermalConductivity(temperature)  ;
		}

        /**
         * @brief Get the Enthalpy of the particle
         * 
         * @param temperature Overall temperature of the particle
         * @return real_t Enthalpy of the particle at the specified temperature
         */
        inline real_t getEnthalpy(real_t temperature)
		{
			// \f$ h \left( T \right) 
			// = \sum_{k \in \left\{ Core, Shell, Product \right}} Y_k h_k \left( T \right)\f$
			return
				_mass_fraction_core_material    * _core_species.getEnthalpy(temperature) +
				_mass_fraction_shell_material   * _shell_species.getEnthalpy(temperature) +
				_mass_fraction_product_material * _product_species.getEnthalpy(temperature);
		}

        /**
         * @brief Determine whether combustion is complete based on
         * remaining reaction mass fractions of reactant substances
         * 
         * @param tolerance Mass fraction of reactant is assumed zero 
         * if less than this value
         * @return true Combustion is complete (one or more reactants have depleted)
         * @return false Combustion is incomplete
         */
        inline bool isCombustionComplete(real_t tolerance = 1E-3)
		{
    		return _mass_fraction_core_material < tolerance || _mass_fraction_shell_material < tolerance;
		}

        inline real_t getMassFractionsCoreMaterial()
		{
			return _mass_fraction_core_material;
		}

        inline real_t getMassFractionsShellMaterial()
		{
			return _mass_fraction_shell_material;
		}

        inline real_t getMassFractionsProductMaterial()
		{
			return _mass_fraction_product_material;
		}
};

#endif