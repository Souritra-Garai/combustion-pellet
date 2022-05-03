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
#include "thermo-physical-properties/Condensed_Species.hpp"

// Required for << operator for printing to file / screen
#include <ostream>

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
        static CondensedSpecies<real_t> *_core_material;
        /**
         * @brief Substance present in the shell of
         * the core-shell particle
         */
        static CondensedSpecies<real_t> *_shell_material;
        /**
         * @brief Substance formed upon reaction of
         * the core and shell substances
         */
        static CondensedSpecies<real_t> *_product_material;

        /**
         * @brief Overall radius of the Core-Shell Particle
         */
        static real_t _overall_radius;
        /**
         * @brief Core radius of the Core-Shell Particle
         */
        static real_t _core_radius;

        /**
         * @brief Mass of the Core-Shell Particle
         */
        static real_t _mass;

    public :

        /**
         * @brief Set up Core Shell Combustion Particle 
         * 
         * @param core_material Substance that forms the core of the core-shell particle
         * @param shell_material Substance that forms the shell of the core-shell particle
         * @param product_material Substance that is produced upon reaction of the core and shell materials
         * @param overall_radius Overall radius of the core-shell particle in m
         * @param core_radius Radius of the core of the particle in m
         */
        static void setUpCoreShellParticle(
            CondensedSpecies<real_t> &core_material,
            CondensedSpecies<real_t> &shell_material,
            CondensedSpecies<real_t> &product_material,
            real_t overall_radius,
            real_t core_radius
        );

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
			real_t volume_fraction_core_material    = _mass_fraction_core_material      / _core_material->getDensity(temperature);
			real_t volume_fraction_shell_material   = _mass_fraction_shell_material     / _shell_material->getDensity(temperature);
			real_t volume_fraction_product_material = _mass_fraction_product_material   / _product_material->getDensity(temperature);

			real_t sum = volume_fraction_core_material + volume_fraction_shell_material + volume_fraction_product_material;

			volume_fraction_core_material       /= sum;
			volume_fraction_shell_material      /= sum;
			volume_fraction_product_material    /= sum;

			// \f$ \rho = \sum_{k \in \left\{ Core, Shell, Product \right}} Y_k \rho_k \f$
			return 1.0 / (
				volume_fraction_core_material    / _core_material->getDensity(temperature) +
				volume_fraction_shell_material   / _shell_material->getDensity(temperature) +
				volume_fraction_product_material / _product_material->getDensity(temperature)
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
				_mass_fraction_core_material    * _core_material->getHeatCapacity(temperature) +
				_mass_fraction_shell_material   * _shell_material->getHeatCapacity(temperature) +
				_mass_fraction_product_material * _product_material->getHeatCapacity(temperature); 
		}

        /**
         * @brief Get the Heat Conductivity of the particle
         * 
         * @return real_t Heat capacity of the particle
         */
        inline real_t getThermalConductivity(real_t temperature)
		{
			// \f$ \lambda = \sum_{k \in \left\{ Core, Shell, Product \right}} Y_k \lambda_k \f$

			real_t volume_fraction_core_material    = _mass_fraction_core_material      / _core_material->getDensity(temperature);
			real_t volume_fraction_shell_material   = _mass_fraction_shell_material     / _shell_material->getDensity(temperature);
			real_t volume_fraction_product_material = _mass_fraction_product_material   / _product_material->getDensity(temperature);

			real_t sum = volume_fraction_core_material + volume_fraction_shell_material + volume_fraction_product_material;

			volume_fraction_core_material       /= sum;
			volume_fraction_shell_material      /= sum;
			volume_fraction_product_material    /= sum;

			return
				volume_fraction_core_material       * _core_material->getThermalConductivity(temperature)     +
				volume_fraction_shell_material      * _shell_material->getThermalConductivity(temperature)    +
				volume_fraction_product_material    * _product_material->getThermalConductivity(temperature)  ;
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
				_mass_fraction_core_material    * _core_material->getEnthalpy(temperature) +
				_mass_fraction_shell_material   * _shell_material->getEnthalpy(temperature) +
				_mass_fraction_product_material * _product_material->getEnthalpy(temperature);
		}

        /**
         * @brief Print the properties of the substance to the given stream
         * 
         * @param output_stream Stream to which the properties are printed
         */
        void printProperties(std::ostream &output_stream);

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
			// printf("%lf, %lf", _mass_fraction_shell_material, _mass_fraction_core_material);
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


/**
 * @brief Function to calculate mass of core shell type particle
 * 
 * @tparam real_t 
 * @param core_material Substance forming the core of the core shell particle
 * @param shell_material Substance forming the shell of the core shell particle
 * @param overall_radius Overall radius of the core shell particle
 * @param core_radius Radius of the core of the core sheel particle
 * @return real_t Mass of the core sheel particle
 */
template<typename real_t>
real_t calcMassCoreShellParticle(
    CondensedSpecies<real_t> core_material,
    CondensedSpecies<real_t> shell_material,
    real_t overall_radius,
    real_t core_radius
);

#endif