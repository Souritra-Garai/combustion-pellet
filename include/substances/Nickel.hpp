#ifndef __NICKEL__
#define __NICKEL__

#include "thermo-physical-properties/Thermal_Conductivity.hpp"
#include "thermo-physical-properties/Enthalpy.hpp"
#include "thermo-physical-properties/Phase.hpp"

#include "thermo-physical-properties/Substance.hpp"

#define SHARPNESS_COEFFICIENT_NI 1.0

Enthalpy<long double> enthalpy_solid_Ni_1(
    13.69160,
    82.46509,
   -174.9548,
    161.6011,
   -0.092417,
   -6.833644
);

Enthalpy<long double> enthalpy_solid_Ni_2(
    1248.045,
   -1257.510,
    0,
    0,
   -165.1266,
   -788.8263
);

Enthalpy<long double> enthalpy_solid_Ni_3(
    16.49839,
    18.74913,
   -6.639841,
    1.717278,
    1.872051,
   -0.467675
);

Enthalpy<long double> enthalpy_liquid_Ni(
    38.91103,
    0.0,
    0.0,
    0.0,
    0.0,
   -2.722630
);

ThermalConductivityQuadraticPolynomial<long double> thermal_conductivity_solid_Ni_1(107, -0.096, 3.2E-5);

ThermalConductivityQuadraticPolynomial<long double> thermal_conductivity_solid_Ni_2_3(59.5, -7.67E-3, 1.7E-5);

ThermalConductivityQuadraticPolynomial<long double> thermal_conductivity_liquid_Ni(17.95, 2.097E-2, 0.0);

Phase<long double> solid_Ni_1(
    8902.0,
    enthalpy_solid_Ni_1,
    thermal_conductivity_solid_Ni_1,
    273.0,
    600.0,
    SHARPNESS_COEFFICIENT_NI
);

Phase<long double> solid_Ni_2(
    8902.0,
    enthalpy_solid_Ni_2,
    thermal_conductivity_solid_Ni_2_3,
    600.0,
    700.0,
    SHARPNESS_COEFFICIENT_NI
);

Phase<long double> solid_Ni_3(
    8902.0,
    enthalpy_solid_Ni_3,
    thermal_conductivity_solid_Ni_2_3,
    700.0,
    1728.0,
    SHARPNESS_COEFFICIENT_NI
);

Phase<long double> liquid_Ni(
    7810.0,
    enthalpy_liquid_Ni,
    thermal_conductivity_liquid_Ni,
    1728.0,
    INFINITY,
    SHARPNESS_COEFFICIENT_NI
);

Phase<long double> phases_Ni[] = {solid_Ni_1, solid_Ni_2, solid_Ni_3, liquid_Ni};

Substance<long double> Nickel(4, phases_Ni, 58.6934E-3);

#endif