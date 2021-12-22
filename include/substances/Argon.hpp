#ifndef __ARGON__
#define __ARGON__

#include "thermo-physical-properties/Thermal_Conductivity.hpp"
#include "thermo-physical-properties/Enthalpy.hpp"
#include "thermo-physical-properties/Phase.hpp"

#include "thermo-physical-properties/Substance.hpp"

#define SHARPNESS_COEFFICIENT_AR 100.0

Enthalpy<long double> internal_energy_gaseous_Ar(
    20.78600,
    2.825911E-7,
   -1.464191E-7,
    1.092131E-8,
   -3.661371E-8,
   -6.197350
);

ThermalConductivityQuadraticPolynomial<long double> thermal_conductivity_gaseous_Ar(1.49E-3, 5.98E-5, -1.92E-8);

Phase<long double> gaseous_Ar(
    1.6025,
    internal_energy_gaseous_Ar,
    thermal_conductivity_gaseous_Ar,
    273.0,
    INFINITY,
    SHARPNESS_COEFFICIENT_AR
);

Phase<long double> phases_Ar[] = {gaseous_Ar};

Substance<long double> Argon(1, phases_Ar, 39.948E-3, 5.0/3.0);

#endif