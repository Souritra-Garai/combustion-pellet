#ifndef __ALUMINIUM__
#define __ALUMINIUM__

#include "thermo-physical-properties/Thermal_Conductivity.hpp"
#include "thermo-physical-properties/Internal_Energy.hpp"
#include "thermo-physical-properties/Phase.hpp"

#include "thermo-physical-properties/Substance.hpp"

#define SHARPNESS_COEFFICIENT_AL 100.0

InternalEnergy<double> internal_energy_solid_Al(
    28.08920,
   -5.414849,
    8.560423,
    3.427370,
   -0.277375,
   -9.147187
);

InternalEnergy<double> internal_energy_liquid_Al(
    31.75104,
    3.935826E-8,
   -1.786515E-8,
    2.694171E-9,
    5.480037E-9,
   -0.945684
);

ThermalConductivityQuadraticPolynomial<double> thermal_conductivity_solid_Al(248.0, -0.067, 0.0);

ThermalConductivityQuadraticPolynomial<double> thermal_conductivity_liquid_Al(33.9, 7.892E-2, -2.099E-5);

Phase<double> solid_Al(
    2700.0,
    internal_energy_solid_Al,
    thermal_conductivity_solid_Al,
    0,
    933.47,
    SHARPNESS_COEFFICIENT_AL
);

Phase<double> liquid_Al(
    2375.0,
    internal_energy_liquid_Al,
    thermal_conductivity_liquid_Al,
    933.47,
    INFINITY,
    SHARPNESS_COEFFICIENT_AL
);

Phase<double> phases_Al[] = {solid_Al, liquid_Al};

Substance<double> Aluminium(2, phases_Al, 26.9815386);

#endif