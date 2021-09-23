#ifndef __NICKEL_ALUMINIDE__
#define __NICKEL_ALUMINIDE__

#include "thermo-physical-properties/Thermal_Conductivity.hpp"
#include "thermo-physical-properties/Internal_Energy.hpp"
#include "thermo-physical-properties/Phase.hpp"

#include "thermo-physical-properties/Substance.hpp"

#define SHARPNESS_COEFFICIENT_NIAL 100.0

InternalEnergy<double> internal_energy_solid_NiAl(
    45.44,
    10.4,
    0.0,
    0.0,
    0.0,
   -147.710
);

InternalEnergy<double> internal_energy_liquid_NiAl(
    71.128,
    0.0,
    0.0,
    0.0,
    0.0,
   -84.95
);

ThermalConductivityQuadraticPolynomial<double> thermal_conductivity_solid_NiAl(72.8, 0.0156, -1.35E-5);

ThermalConductivityQuadraticPolynomial<double> thermal_conductivity_liquid_NiAl(23.9, 0, 0);

Phase<double> solid_NiAl(
    6050.0,
    internal_energy_solid_NiAl,
    thermal_conductivity_solid_NiAl,
    273.0,
    1912.0,
    SHARPNESS_COEFFICIENT_NIAL
);

Phase<double> liquid_NiAl(
    4798.3,
    internal_energy_liquid_NiAl,
    thermal_conductivity_liquid_NiAl,
    1912.0,
    INFINITY,
    SHARPNESS_COEFFICIENT_NIAL
);

Phase<double> phases_NiAl[] = {solid_NiAl, liquid_NiAl};

Substance<double> NickelAluminide(2, phases_NiAl, 85.6748386E-3);

#endif