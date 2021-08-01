/**
 * @file Packed-Pellet.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Definition of class to represent thermo-physical properties
 * of a pellet packed with core-shell type combustion particles
 * @version 0.1
 * @date 2021-07-30
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __PACKED_PELLET__
#define __PACKED_PELLET__

// Required for << operator for printing to file / screen
#include <ostream>

// Rquired for Substance class
#include "thermo-physical-properties/Substance.hpp"
// Required for CoreShellCombustionParticle class
#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"

/**
 * @brief Class to represent thermo - physical properties of a cylindrical
 *  pellet packed with core - shell type combustion particles
 * 
 * @tparam real_t float, double or long double data types
 * to represent real numbers
 */
template<typename real_t>
class PackedPellet
{
    protected :

        /**
         * @brief Length of cylindrical pellet in \f$ m \f$
         */
        static real_t _length;
        /**
         * @brief Diameter of cylindrical pellet in \f$ m \f$
         */
        static real_t _diameter;

        /**
         * @brief Convective heat transfer coefficient
         * for heat loss to surrounding atmosphere via 
         * convection. Units - \f$ W / m^2 - K \f$
         */
        static real_t _convective_heat_transfer_coefficient;
        /**
         * @brief Emissivity for the surface of the cylinder
         * through which radiative heat loss to the surrounding
         * atmosphere occurs
         */
        static real_t _radiative_emissivity;
        /**
         * @brief Temperature of surrounding atmosphere in \f$ K \f$
         */
        static real_t _ambient_temperature;

        /**
         * @brief Fluid (gas) filling the voids of the packed pellet
         */
        static Substance<real_t> _degassing_fluid;

        /**
         * @brief Fraction of the pellet volume occupied by
         * core-shell type combustion particles
         */
        const real_t _particle_volume_fractions;
        /**
         * @brief Fraction of the pellet mass contributed by
         * core-shell type combustion particles
         */
        const real_t _particle_mass_fractions;
        /**
         * @brief Density of the pellet in \f$ kg / m^3 \f$
         */
        const real_t _density;

        /**
         * @brief Calculates the particle mass fraction given the 
         * volume fraction assuming a fresh instance of CoreShellCombustionParticle
         * for the particle composition
         * @param particle_volume_fractions Fraction of the pellet volume occupied by
         * core-shell type combustion particles
         * @return real_t Fraction of the pellet mass contributed by
         * core-shell type combustion particles
         */
        static real_t calcParticleMassFractions(real_t particle_volume_fractions);
        /**
         * @brief Calculates the mean density of the particle - degassing fluid
         * mixture, assuming a fresh instance of CoreShellCombustionParticle
         * for the particle composition
         * @param particle_volume_fractions Fraction of the pellet volume occupied by
         * core-shell type combustion particles
         * @return real_t Mean density of the pellet in \f$ kg / m^3 \f$
         */
        static real_t calcDensity(real_t particle_volume_fractions);

    public :

        /**
         * @brief Construct a new Packed Pellet object
         * 
         * @param particle_volume_fractions Fraction of the pellet volume occupied by
         * core-shell type combustion particles
         */
        PackedPellet(real_t particle_volume_fractions);

        /**
         * @brief Set the dimensions of the cylindrical pellet
         * 
         * @param length Length of the cylindrical pellet
         * @param diameter Diameter of the cylindrical pellet
         */
        static void setPelletDimensions(
            real_t length,
            real_t diameter
        );

        /**
         * @brief Set the parameters for Heat Loss to the ambient
         * 
         * @param convective_heat_transfer_coefficient Convective heat transfer 
         * coefficient for heat loss to surrounding atmosphere via convection
         * @param radiative_emissivity Emissivity for the surface of the cylinder
         * through which radiative heat loss to the surrounding atmosphere occurs
         * @param ambient_temperature Temperature of surrounding atmosphere
         */
        static void setAmbientHeatLossParameters(
            real_t convective_heat_transfer_coefficient,
            real_t radiative_emissivity,
            real_t ambient_temperature
        );

        /**
         * @brief Set the fluid substance that fills the voids in the pellet
         * 
         * @param degassing_fluid Fluid (gas) substance degassing the pellet
         */
        static void setDegassingFluid(Substance<real_t> degassing_fluid);

        /**
         * @brief Get the mean density of the pellet (evaluated at the initial condition)
         * 
         * @return real_t Mean Density of the pellet in \f$ kg /m^3 \f$
         */
        real_t getDensity();

        /**
         * @brief Get the mean heat capacity at a point in the pellet
         * where state of the core-shell combustion particle is represented by
         * the particle in the argument
         * 
         * @param ptr_2_particle Pointer to a core-shell combustion particle
         * that represents the state of all particles in the pellet
         * @return real_t Mean heat capacity of the pellet in \f$ J / kg - K \f$
         */
        real_t getHeatCapacity(CoreShellCombustionParticle<real_t> *ptr_2_particle);

        /**
         * @brief Get he mean heat conductivity at a point in the pellet
         * where state of the core-shell combustion particle is represented by
         * the particle in the argument
         * 
         * @param ptr_2_particle Pointer to a core-shell combustion particle
         * that represents the state of all particles in the pellet
         * @return real_t Mean heat conductivity of the pellet in \f$ W / m - K \f$
         */
        real_t getHeatConductivity(CoreShellCombustionParticle<real_t> *ptr_2_particle);

        /**
         * @brief Get he mean enthalpy at a point in the pellet, at the specified temperature,
         * where state of the core-shell combustion particle is represented by
         * the particle in the argument
         * 
         * @param ptr_2_particle Pointer to a core-shell combustion particle
         * that represents the state of all particles in the pellet
         * @param temperature Temperature at which enthalpy is evaluated \f$ K \f$
         * @return real_t Mean enthalpy of the pellet in \f$ J / kg \f$
         */
        real_t getEnthalpy(CoreShellCombustionParticle<real_t> *ptr_2_particle, real_t temperature);

        /**
         * @brief Print the properties of the pellet (evaluated at the initial condition)
         * to the output stream
         * @param output_stream Stream to which the properties are to be printed 
         */
        void printProperties(std::ostream &output_stream);
};

#endif