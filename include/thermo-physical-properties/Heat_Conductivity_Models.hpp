/**
 * @file Heat_Conductivity_Models.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Declares functions for determining heat conductivity for a mixture
 * of solid particles and gases
 * @version 0.1
 * @date 2021-08-01
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __HEAT_CONDUCTIVITY__
#define __HEAT_CONDUCTIVITY__

// Required for sqrt and pow functions
#include <math.h>

/**
 * @brief Get the heat conductivity for particle - gas mixture using 
 * Bruggeman model
 * 
 * @tparam real_t 
 * @param particle_volume_fraction 
 * @param heat_conductivity_particle 
 * @param heat_conductivity_fluid 
 * @return real_t 
 */
template<typename real_t>
real_t getBruggemanHeatConductivity(
    real_t particle_volume_fraction,
    real_t heat_conductivity_particle,
    real_t heat_conductivity_fluid
) {
    real_t discriminant = 
        pow((3.0 * particle_volume_fraction - 1.0) * heat_conductivity_particle / heat_conductivity_fluid, 2) +
        pow(2.0 - 3.0 * particle_volume_fraction, 2) +
        (2.0 + 9.0 * particle_volume_fraction - 9.0 * pow(particle_volume_fraction, 2));

    return (
        (3.0 * particle_volume_fraction - 1.0) * heat_conductivity_particle / heat_conductivity_fluid +
        2.0 - 3.0 * particle_volume_fraction +
        sqrt(discriminant)
    ) * heat_conductivity_fluid / 4.0;
}

template<typename real_t>
real_t getCoContinuousHeatConductivity(
    real_t particle_volume_fraction,
    real_t heat_conductivity_particle,
    real_t heat_conductivity_fluid
) {
    real_t heat_conductivity_series = 1.0 / (
        particle_volume_fraction / heat_conductivity_particle +
        (1.0 - particle_volume_fraction) / heat_conductivity_fluid
    );

    real_t heat_conductivity_parallel = 
        particle_volume_fraction * heat_conductivity_particle + 
        (1.0 - particle_volume_fraction) * heat_conductivity_fluid
    ;

    return (heat_conductivity_series / 2.0) * (
        sqrt(1 + 8.0 * heat_conductivity_parallel / heat_conductivity_series) -
        1.0
    );
}

template <typename real_t>
real_t getMaxwellEucken1HeatConudctivity(
    real_t particle_volume_fraction,
    real_t heat_conductivity_particle,
    real_t heat_conductivity_fluid
) {
    return (
        heat_conductivity_particle * particle_volume_fraction * (2 * heat_conductivity_particle + heat_conductivity_fluid) +
        heat_conductivity_fluid * (1.0 - particle_volume_fraction) * 3.0 * heat_conductivity_particle
    ) / (
        particle_volume_fraction * (2 * heat_conductivity_particle + heat_conductivity_fluid) +
        (1.0 - particle_volume_fraction) * 3.0 * heat_conductivity_particle
    );
}

template <typename real_t>
real_t getMaxwellEucken2HeatConudctivity(
    real_t particle_volume_fraction,
    real_t heat_conductivity_particle,
    real_t heat_conductivity_fluid
) {
    return getMaxwellEucken1HeatConudctivity(
        1.0 - particle_volume_fraction,
        heat_conductivity_fluid,
        heat_conductivity_particle
    );
}

#endif