#ifndef __PACKED_PELLET__
#define __PACKED_PELLET__

#include <ostream>

#include "thermo-physical-properties/Ideal_Gas.cuh"
#include "thermo-physical-properties/Species.cuh"
#include "thermo-physical-properties/Core-Shell-Particle.cuh"

namespace PackedPellet
{
	__device__ double length;
	__device__ double diameter;

	__device__ double convective_heat_transfer_coefficient_curved_surface;
	__device__ double convective_heat_transfer_coefficient_flat_surface;
	__device__ double radiative_emissivity;

	__device__ double ambient_temperature;
	__device__ double ignition_temperature;

	__device__ IdealGas *degassing_fluid;

	__device__ void initializePelletDimensions(double length, double diameter)
	{
		PackedPellet::length = length;
		PackedPellet::diameter = diameter;
	}

	__device__ void setHeatLossParameters(
		double convective_heat_transfer_coefficient_curved_surface,
		double convective_heat_transfer_coefficient_flat_surface,
		double radiative_emissivity
	) {
		PackedPellet::convective_heat_transfer_coefficient_curved_surface = convective_heat_transfer_coefficient_curved_surface;
		PackedPellet::convective_heat_transfer_coefficient_flat_surface = convective_heat_transfer_coefficient_flat_surface;
		PackedPellet::radiative_emissivity = radiative_emissivity;
	}

	__device__ void setTemperatureParameters(
		double ambient_temperature,
		double ignition_temperature
	) {
		PackedPellet::ambient_temperature = ambient_temperature;
		PackedPellet::ignition_temperature = ignition_temperature;
	}

	__device__ void setDegassingFluid(
		IdealGas *degassing_fluid
	) {
		PackedPellet::degassing_fluid = degassing_fluid;
	}

	class Pellet
	{
		private:
	
			double _degassing_fluid_volume_fractions;
			double _overall_particle_density;

		public:

			Pellet(double particle_volume_fractions)
			{
				_degassing_fluid_volume_fractions = 1.0 - particle_volume_fractions;
				
				_overall_particle_density = particle_volume_fractions * CoreShellParticle::Particle().getDensity(ambient_temperature);
			}

			__device__ __forceinline__ double getPelletThermalConductivity(CoreShellParticle::Particle *particle, double temperature)
			{
				return 1.0;
			}
	};
} // namespace PackedPellet


class PackedPellet
{
    protected :

        static double _length;
        /**
         * @brief Diameter of cylindrical pellet in \f$ m \f$
         */
        static double _diameter;

        /**
         * @brief Convective heat transfer coefficient
         * for heat loss to surrounding atmosphere via 
         * convection at the curved surface of the pellet.
		 * Units - \f$ W / m^2 - K \f$
         */
        static double _convective_heat_transfer_coefficient_curved_surface;
		/**
         * @brief Convective heat transfer coefficient
         * for heat loss to surrounding atmosphere via 
         * convection at the flat surface of the pellet.
		 * Units - \f$ W / m^2 - K \f$
         */
        static double _convective_heat_transfer_coefficient_flat_surface;
        /**
         * @brief Emissivity for the surface of the cylinder
         * through which radiative heat loss to the surrounding
         * atmosphere occurs
         */
        static double _radiative_emissivity;
        
        /**
         * @brief Temperature of surrounding atmosphere in \f$ K \f$
         */
        static double _ambient_temperature;

        /**
         * @brief Temperature at which reactions begin in the pellet
         */
        static double _ignition_temperature;

        /**
         * @brief Fluid (gas) filling the voids of the packed pellet
         */
        static IdealGas *_degassing_fluid;

        /**
         * @brief Fraction of the pellet volume occupied by
         * the degassing fluid
         */
        const double _degassing_fluid_volume_fractions;
        
        /**
         * @brief Density of particle in the pellet, i.e.,
         * \f$ \rho_P \phi_P \f$ in \f$ kg / m^3 \f$
         */
		const double _overall_particle_density;

        /**
         * @brief Calculates the overall density of particles in the pellet,
         * assuming a fresh instance of CoreShellCombustionParticle
         * for the particle composition
         * @param particle_volume_fractions Fraction of the pellet volume occupied by
         * core-shell type combustion particles
         * @return double Overall density of particles in the pellet in \f$ kg / m^3 \f$
         */
        static double calcOverallParticleDensity(double particle_volume_fractions);

    public :

        /**
         * @brief Construct a new Packed Pellet object
         * 
         * @param particle_volume_fractions Fraction of the pellet volume occupied by
         * core-shell type combustion particles
         */
        PackedPellet(double particle_volume_fractions);

        /**
         * @brief Set the dimensions of the cylindrical pellet
         * 
         * @param length Length of the cylindrical pellet
         * @param diameter Diameter of the cylindrical pellet
         */
        static void setPelletDimensions(
            double length,
            double diameter
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
            double convective_heat_transfer_coefficient_curved_surface,
			double convective_heat_transfer_coefficient_flat_surface,
            double radiative_emissivity
        );

        /**
         * @brief Set the Temperature Parameters for the pellet
         * 
         * @param ignition_temperature Temperature at which reactions of energetic 
         * particles begin inside the pellet
         * @param ambient_temperature Temperature of surrounding atmosphere
         */
        static void setTemperatureParameters(
            double ignition_temperature,
            double ambient_temperature
        );

        /**
         * @brief Set the fluid that fills the voids in the pellet
         * 
         * @param degassing_fluid Fluid degassing the pellet
         */
        static void setDegassingFluid(IdealGas &degassing_fluid);

        /**
         * @brief Get he mean heat conductivity at a point in the pellet
         * where state of the core-shell combustion particle is represented by
         * the particle in the argument
         * 
         * @param ptr_2_particle Pointer to a core-shell combustion particle
         * that represents the state of all particles in the pellet
         * @return double Mean heat conductivity of the pellet in \f$ W / m - K \f$
         */
        double getThermalConductivity(CoreShellParticle::Particle *ptr_2_particle, double temperature);

        /**
         * @brief Print the properties of the pellet (evaluated at the initial condition)
         * to the output stream
         * @param output_stream Stream to which the properties are to be printed 
         */
        void printProperties(std::ostream &output_stream);
};

#endif