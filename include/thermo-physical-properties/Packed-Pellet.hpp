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

#include "math/Data-Type.hpp"

#include "thermo-physical-properties/Ideal-Gas.hpp"
#include "thermo-physical-properties/Core-Shell-Particle.hpp"
#include "thermo-physical-properties/Thermal-Conductivity-Pellet.hpp"

class PackedPellet
{
    public :

		static const real_t length;
		static const real_t diameter;

		static const real_t convective_heat_transfer_coefficient_curved_surface;
		static const real_t convective_heat_transfer_coefficient_flat_surface;
		static const real_t radiative_emissivity;
		
		static const real_t ambient_pressure;
    	static const real_t ambient_temperature;

		static IdealGas interstitial_gas;

		const real_t interstitial_volume_fractions;
        
        const real_t overall_particle_density;

		PackedPellet(real_t particle_volume_fractions);

		inline real_t getThermalConductivity(CoreShellParticle *ptr_2_particle, real_t temperature) const
		{
			return getThermalConductivityMEB(
				interstitial_volume_fractions,
				interstitial_gas.getThermalConductivity(temperature),
				ptr_2_particle->getThermalConductivity(temperature)
			);
		}
};

#endif