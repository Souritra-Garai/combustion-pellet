#ifndef __ALUMINIUM__
#define __ALUMINIUM__

#include "thermo-physical-properties/Thermal_Conductivity.hpp"
#include "thermo-physical-properties/Enthalpy.hpp"
#include "thermo-physical-properties/Phase.hpp"

#include "thermo-physical-properties/Condensed_Species.hpp"

extern long doub;

Enthalpy<long double> enthalpy_solid_Al(
    28.08920,
   -5.414849,
    8.560423,
    3.427370,
   -0.277375,
   -9.147187
);

Enthalpy<long double> enthalpy_liquid_Al(
    31.75104,
    3.935826E-8,
   -1.786515E-8,
    2.694171E-9,
    5.480037E-9,
   -0.945684
);

ThermalConductivityQuadraticPolynomial<long double> thermal_conductivity_solid_Al(248.0, -0.067, 0.0);

ThermalConductivityQuadraticPolynomial<long double> thermal_conductivity_liquid_Al(33.9, 7.892E-2, -2.099E-5);

Phase<long double> solid_Al(
    2700.0,
    enthalpy_solid_Al,
    thermal_conductivity_solid_Al,
    -INFINITY,
    933.47
);

Phase<long double> liquid_Al(
    2375.0,
    enthalpy_liquid_Al,
    thermal_conductivity_liquid_Al,
    933.47,
    INFINITY
);

Phase<long double> phases_Al[] = {solid_Al, liquid_Al};

CondensedSpecies<long double> Aluminium(2, phases_Al, 26.9815386E-3);

#endif