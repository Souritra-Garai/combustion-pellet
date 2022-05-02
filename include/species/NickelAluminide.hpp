#ifndef __NICKEL_ALUMINIDE__
#define __NICKEL_ALUMINIDE__

#include "thermo-physical-properties/Thermal_Conductivity.hpp"
#include "thermo-physical-properties/Enthalpy.hpp"
#include "thermo-physical-properties/Phase.hpp"

#include "thermo-physical-properties/Condensed_Species.hpp"

extern long double sharpness_coefficient;

Enthalpy<long double> enthalpy_solid_NiAl(
    41.86,
    13.81,
    0.0,
    0.0,
    0.0,
   -131.5
);

Enthalpy<long double> enthalpy_liquid_NiAl(
    71.16,
    0.0,
    0.0,
    0.0,
    0.0,
   -99.54
);

ThermalConductivityQuadraticPolynomial<long double> thermal_conductivity_solid_NiAl(72.8, 0.0156, -1.35E-5);

ThermalConductivityQuadraticPolynomial<long double> thermal_conductivity_liquid_NiAl(23.9, 0, 0);

Phase<long double> solid_NiAl(
    5650.0,
    enthalpy_solid_NiAl,
    thermal_conductivity_solid_NiAl,
    -INFINITY,
    1912.0
);

Phase<long double> liquid_NiAl(
    4798.3,
    enthalpy_liquid_NiAl,
    thermal_conductivity_liquid_NiAl,
    1912.0,
    INFINITY
);

Phase<long double> phases_NiAl[] = {solid_NiAl, liquid_NiAl};

CondensedSpecies<long double> NickelAluminide(2, phases_NiAl, 85.6748386E-3);

#endif