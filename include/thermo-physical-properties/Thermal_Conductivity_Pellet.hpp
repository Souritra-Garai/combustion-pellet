/**
 * @file Thermal_Conductivity_Pellet.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief 
 * @version 1.0
 * @date 2021-10-08
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __THERMAL_CONDUCTIVITY_PELLET__
#define __THERMAL_CONDUCTIVITY_PELLET__

template<typename real_t>
real_t getThermalConductivityEMT(
    real_t porosity,
    real_t thermal_conductivity_fluid,
    real_t thermal_conductivity_particle
);

template<typename real_t>
real_t getThermalConductivityCC(
    real_t porosity,
    real_t thermal_conductivity_fluid,
    real_t thermal_conductivity_particle
);

template<typename real_t>
real_t getThermalConductivityME1(
    real_t porosity,
    real_t thermal_conductivity_fluid,
    real_t thermal_conductivity_particle
);

template<typename real_t>
real_t getThermalConductivityME2(
    real_t porosity,
    real_t thermal_conductivity_fluid,
    real_t thermal_conductivity_particle
);

template<typename real_t>
real_t getThermalConductivityMEB(
    real_t porosity,
    real_t thermal_conductivity_fluid,
    real_t thermal_conductivity_particle
);

template<typename real_t>
real_t getThermalConductivityCCB(
    real_t porosity,
    real_t thermal_conductivity_fluid,
    real_t thermal_conductivity_particle
);

#endif