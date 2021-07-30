/**
 * @file Packed-Pellet.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief 
 * @version 0.1
 * @date 2021-07-30
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __PACKED_PELLET__
#define __PACKED_PELLET__

#include <ostream>

#include "thermo-physical-properties/Substance.hpp"
#include "pde-problems/Core-Shell-Diffusion.hpp"

template<typename real_t>
class PackedPellet
{
    protected :

        static real_t _length;
        static real_t _diameter;

        static real_t _convective_heat_transfer_coefficient;
        static real_t _radiative_emissivity;

        static real_t _ambient_temperature;

        static Substance<real_t> _degassing_fluid;

        const real_t _particle_volume_fractions;
        
        const real_t _particle_mass_fractions;
        
        const real_t _density;

        static real_t calcParticleMassFractions(real_t particle_volume_fractions);
        static real_t calcDensity(real_t particle_volume_fractions);

    public :

        PackedPellet(real_t particle_volume_fractions);

        static void setPelletDimensions(
            real_t length,
            real_t diameter
        );

        static void setAmbientHeatLossProperties(
            real_t convective_heat_transfer_coefficient,
            real_t radiative_emissivity,
            real_t ambient_temperature
        );

        static void setDegassingFluid(Substance<real_t> degassing_fluid);

        real_t getDensity();

        real_t getHeatCapacity(CoreShellDiffusion<real_t> &particle);

        real_t getHeatConductivity(CoreShellDiffusion<real_t> &particle);

        real_t getEnthalpy(CoreShellDiffusion<real_t> &particle, real_t temperature);

        void printProperties(std::ostream &output_stream);
};

#endif