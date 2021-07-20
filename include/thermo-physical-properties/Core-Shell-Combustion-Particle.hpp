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
 * estimate thermodynamic properties for varying composition and temperature
 * 
 * @tparam real_t 
 */
template<typename real_t>
class CoreShellCombustionParticle
{
    protected :

        /**
         * @brief Mass fraction of reactant A
         */
        real_t mass_fraction_core_material;
        /**
         * @brief Mass fraction of reactant B
         */
        real_t mass_fraction_shell_material;
        /**
         * @brief Mass fraction of product AB
         */
        real_t mass_fraction_product_material;

    public :

        /**
         * @brief Reactant A substance
         * initially present in the core
         */
        static Substance<real_t> core_material;
        /**
         * @brief Reactant B substance
         * initially present in the shell
         */
        static Substance<real_t> shell_material;
        /**
         * @brief Product AB substance 
         */
        static Substance<real_t> product_material;

        /**
         * @brief Overall radius of the Core-Shell Particle
         */
        static real_t overall_radius;
        /**
         * @brief Core radius of the Core-Shell Particle
         * 
         */
        static real_t core_radius;

        /**
         * @brief Mass of the Core-Shell Particle
         */
        static real_t mass;

        CoreShellCombustionParticle();

        /**
         * @brief Get the Density of the Particle
         * 
         * @return real_t Density of the particle
         */
        real_t getDensity();

        /**
         * @brief Get the Heat Capacity of the particle
         * 
         * @return real_t Heat capacity of the particle
         */
        real_t getHeatCapacity();

        /**
         * @brief Get the Heat Conductivity of the particle
         * 
         * @return real_t Heat capacity of the particle
         */
        real_t getHeatConductivity();

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

        static void setUpCoreShellParticle(
            Substance<real_t> core_material,
            Substance<real_t> shell_material,
            real_t overall_radius,
            real_t core_radius
        );
};

template<typename real_t>
real_t calcMassCoreShellParticle(
    Substance<real_t> core_material,
    Substance<real_t> shell_material,
    real_t overall_radius,
    real_t core_radius
);

template<typename real_t> Substance<real_t> CoreShellCombustionParticle<real_t>::core_material = Substance<real_t>(0, 0, 0, 0, 0);
template<typename real_t> Substance<real_t> CoreShellCombustionParticle<real_t>::shell_material = Substance<real_t>(0, 0, 0, 0, 0);
template<typename real_t> Substance<real_t> CoreShellCombustionParticle<real_t>::product_material = Substance<real_t>(0, 0, 0, 0, 0);

template<typename real_t> real_t CoreShellCombustionParticle<real_t>::overall_radius = 0;
template<typename real_t> real_t CoreShellCombustionParticle<real_t>::core_radius = 0;
template<typename real_t> real_t CoreShellCombustionParticle<real_t>::mass = 0;

#endif