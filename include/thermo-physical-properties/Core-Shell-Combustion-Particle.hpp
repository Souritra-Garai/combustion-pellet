/**
 * @file Core-Shell-Combustion-Particle.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header defines a class for core shell particle
 * @version 0.1
 * @date 2021-07-08
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __CORE_SHELL_COMBUSTION_PARTICLE__
#define __CORE_SHELL_COMBUSTION_PARTICLE__

// Required for Substance class
#include "thermo-physical-properties/Substance.hpp"

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
class CoreShellCombustionParticle
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
        static Substance<real_t> *_core_material;
        /**
         * @brief Substance present in the shell of
         * the core-shell particle
         */
        static Substance<real_t> *_shell_material;
        /**
         * @brief Substance formed upon reaction of
         * the core and shell substances
         */
        static Substance<real_t> *_product_material;

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
        static void setUpCoreShellCombustionParticle(
            Substance<real_t> &core_material,
            Substance<real_t> &shell_material,
            Substance<real_t> &product_material,
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
        CoreShellCombustionParticle();

        /**
         * @brief Get the Density of the Particle
         * 
         * @return real_t Density of the particle
         */
        real_t getDensity(real_t temperature);

        /**
         * @brief Get the Heat Capacity of the particle
         * 
         * @return real_t Heat capacity of the particle
         */
        real_t getHeatCapacity(real_t temperature);

        /**
         * @brief Get the Heat Conductivity of the particle
         * 
         * @return real_t Heat capacity of the particle
         */
        real_t getThermalConductivity(real_t temperature);

        /**
         * @brief Get the Enthalpy of the particle
         * 
         * @param temperature Overall temperature of the particle
         * @return real_t Enthalpy of the particle at the specified temperature
         */
        real_t getEnthalpy(real_t temperature);

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
        bool isCombustionComplete(real_t tolerance = 1E-3);

        real_t getMassFractionsCoreMaterial();
        real_t getMassFractionsShellMaterial();
        real_t getMassFractionsProductMaterial();
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
    Substance<real_t> core_material,
    Substance<real_t> shell_material,
    real_t overall_radius,
    real_t core_radius
);

#endif