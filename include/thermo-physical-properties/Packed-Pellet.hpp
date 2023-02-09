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
#include "thermo-physical-properties/Ideal-Gas.hpp"
// Required for CoreShellParticle class
#include "thermo-physical-properties/Core-Shell-Particle.hpp"
// Required for functions to calculate heat conductivity
#include "thermo-physical-properties/Thermal-Conductivity-Pellet.hpp"

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
    public :

        /**
         * @brief Length of cylindrical pellet in \f$ m \f$
         */
        static const real_t length;
        /**
         * @brief Diameter of cylindrical pellet in \f$ m \f$
         */
        static const real_t diameter;

        /**
         * @brief Convective heat transfer coefficient
         * for heat loss to surrounding atmosphere via 
         * convection at the curved surface of the pellet.
		 * Units - \f$ W / m^2 - K \f$
         */
        static const real_t convective_heat_transfer_coefficient_curved_surface;
		/**
         * @brief Convective heat transfer coefficient
         * for heat loss to surrounding atmosphere via 
         * convection at the flat surface of the pellet.
		 * Units - \f$ W / m^2 - K \f$
         */
        static const real_t convective_heat_transfer_coefficient_flat_surface;
        /**
         * @brief Emissivity for the surface of the cylinder
         * through which radiative heat loss to the surrounding
         * atmosphere occurs
         */
        static const real_t radiative_emissivity;
        
        /**
         * @brief Temperature of surrounding atmosphere in \f$ K \f$
         */
        static const real_t ambient_temperature;

        /**
         * @brief Fluid (gas) filling the voids of the packed pellet
         */
        static IdealGas<real_t> interstitial_gas;

        /**
         * @brief Fraction of the pellet volume occupied by
         * the degassing fluid
         */
        const real_t interstitial_volume_fractions;
        
        /**
         * @brief Density of particle in the pellet, i.e.,
         * \f$ \rho_P \phi_P \f$ in \f$ kg / m^3 \f$
         */
		const real_t overall_particle_density;

        /**
         * @brief Calculates the overall density of particles in the pellet,
         * assuming a fresh instance of CoreShellParticle
         * for the particle composition
         * @param particle_volume_fractions Fraction of the pellet volume occupied by
         * core-shell type combustion particles
         * @return real_t Overall density of particles in the pellet in \f$ kg / m^3 \f$
         */
        static real_t calcOverallParticleDensity(real_t particle_volume_fractions);

        /**
         * @brief Construct a new Packed Pellet object
         * 
         * @param particle_volume_fractions Fraction of the pellet volume occupied by
         * core-shell type combustion particles
         */
        PackedPellet(real_t particle_volume_fractions);

        /**
         * @brief Get he mean heat conductivity at a point in the pellet
         * where state of the core-shell combustion particle is represented by
         * the particle in the argument
         * 
         * @param ptr_2_particle Pointer to a core-shell combustion particle
         * that represents the state of all particles in the pellet
         * @return real_t Mean heat conductivity of the pellet in \f$ W / m - K \f$
         */
        inline real_t getThermalConductivity(CoreShellParticle<real_t> *ptr_2_particle, real_t temperature)
		{    
			// Return the heat conductivity determined using Bruggeman model
			return getThermalConductivityMEB(
				interstitial_volume_fractions,
				interstitial_gas.getThermalConductivity(temperature),
				ptr_2_particle->getThermalConductivity(temperature)
			);
		}
};

#endif