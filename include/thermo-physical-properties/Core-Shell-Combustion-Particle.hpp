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
        real_t mass_fraction_reactant_A;
        /**
         * @brief Mass fraction of reactant B
         */
        real_t mass_fraction_reactant_B;
        /**
         * @brief Mass fraction of product AB
         */
        real_t mass_fraction_product_AB;

        /**
         * @brief Reactant A substance
         * initially present in the core
         */
        Substance<real_t> *const reactant_A;
        /**
         * @brief Reactant B substance
         * initially present in the shell
         */
        Substance<real_t> *const reactant_B;
        /**
         * @brief Product AB substance 
         */
        Substance<real_t> *const product_AB;

        /**
         * @brief Overall radius of the core-sheel particle
         */
        real_t *const overall_radius;
        /**
         * @brief Radius of the core of the core sheel particle
         */
        real_t *const core_radius;

        /**
         * @brief Mass of the core-shell particle
         */
        real_t *const mass;

    public :

        /**
         * @brief Construct a new Core Shell Combustion Particle object
         * 
         * @param core_material Substance that forms the core material
         * @param shell_material Substance that forms the shell material
         * @param product_material Substance that is produced as a result of the reaction of the core and shell material
         * @param particle_radius Overall radius of the particle
         * @param core_radius Radius of the core of the Core Shell Particle
         * @param mass Mass of the particle
         */
        CoreShellCombustionParticle(
            Substance<real_t> &core_material,
            Substance<real_t> &shell_material,
            Substance<real_t> &product_material,
            real_t &particle_radius,
            real_t &core_radius,
            real_t &mass
        );

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