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

// Required for IdealGas class
#include "thermo-physical-properties/IdealGas.hpp"
// Required for Species class
#include "thermo-physical-properties/Substance.hpp"
// Required for CoreShellCombustionParticle class
#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"

/**
 * @brief Class to represent thermo - physical properties of a cylindrical
 * pellet packed with core - shell type combustion particles
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
         * convection at the curved surface of the pellet.
		 * Units - \f$ W / m^2 - K \f$
         */
        static real_t _convective_heat_transfer_coefficient_curved_surface;
		/**
         * @brief Convective heat transfer coefficient
         * for heat loss to surrounding atmosphere via 
         * convection at the flat surface of the pellet.
		 * Units - \f$ W / m^2 - K \f$
         */
        static real_t _convective_heat_transfer_coefficient_flat_surface;
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
         * @brief Temperature at which reactions begin in the pellet
         */
        static real_t _ignition_temperature;

        /**
         * @brief Fluid (gas) filling the voids of the packed pellet
         */
        static IdealGas<real_t> *_degassing_fluid;

        /**
         * @brief Fraction of the pellet volume occupied by
         * the degassing fluid
         */
        const real_t _degassing_fluid_volume_fractions;
        
        /**
         * @brief Density of particle in the pellet, i.e.,
         * \f$ \rho_P \phi_P \f$ in \f$ kg / m^3 \f$
         */
		const real_t _overall_particle_density;

        /**
         * @brief Calculates the overall density of particles in the pellet,
         * assuming a fresh instance of CoreShellCombustionParticle
         * for the particle composition
         * @param particle_volume_fractions Fraction of the pellet volume occupied by
         * core-shell type combustion particles
         * @return real_t Overall density of particles in the pellet in \f$ kg / m^3 \f$
         */
        static real_t calcOverallParticleDensity(real_t particle_volume_fractions);

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
         */

		/**
		 * @brief Set the parameters for Heat Loss to the ambient
		 * 
		 * @param convective_heat_transfer_coefficient_curved_surface Convective heat transfer 
         * coefficient for heat loss to surrounding atmosphere via convection at the curved surface
		 * @param convective_heat_transfer_coefficient_flat_surface Convective heat transfer 
         * coefficient for heat loss to surrounding atmosphere via convection at the flat surface
		 * @param radiative_emissivity Emissivity for the surface of the cylinder
         * through which radiative heat loss to the surrounding atmosphere occurs
		 */
        static void setAmbientHeatLossParameters(
            real_t convective_heat_transfer_coefficient_curved_surface,
			real_t convective_heat_transfer_coefficient_flat_surface,
            real_t radiative_emissivity
        );

        /**
         * @brief Set the Temperature Parameters for the pellet
         * 
         * @param ignition_temperature Temperature at which reactions of energetic 
         * particles begin inside the pellet
         * @param ambient_temperature Temperature of surrounding atmosphere
         */
        static void setTemperatureParameters(
            real_t ignition_temperature,
            real_t ambient_temperature
        );

        /**
         * @brief Set the fluid that fills the voids in the pellet
         * 
         * @param degassing_fluid Fluid degassing the pellet
         */
        static void setDegassingFluid(IdealGas<real_t> &degassing_fluid);

        /**
         * @brief Get he mean heat conductivity at a point in the pellet
         * where state of the core-shell combustion particle is represented by
         * the particle in the argument
         * 
         * @param ptr_2_particle Pointer to a core-shell combustion particle
         * that represents the state of all particles in the pellet
         * @return real_t Mean heat conductivity of the pellet in \f$ W / m - K \f$
         */
        real_t getThermalConductivity(CoreShellCombustionParticle<real_t> *ptr_2_particle, real_t temperature);

        /**
         * @brief Print the properties of the pellet (evaluated at the initial condition)
         * to the output stream
         * @param output_stream Stream to which the properties are to be printed 
         */
        void printProperties(std::ostream &output_stream);
};

#endif