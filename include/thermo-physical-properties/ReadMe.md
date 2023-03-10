# Thermo-physical Properties

This group of libraries define classes to model thermo-physical properties of species and systems. All properties are returned in their SI units. Furthermore, the intensice properties are calculated on a mass basis (unless explicitly specified).

## Ideal Gas
The header `Ideal-Gas.hpp` defines the class `IdealGas` to model the density (ideal gas equation of state), specific heat capacity (Shomate equation), specific enthalpy (integral of Shomate equation) and thermal conductivity (quadratic polynomial fit) of an ideal gas. Specific heat capacity, specific enthalpy and thermal conductivity are assumed to be independent of pressure. Hence, the corresponding values at standard pressure conditions, i.e., 101325 Pa, are calculated.

## Phase
The header `Phase.hpp` defines the class `Phase` to model the density (constant), specific heat capacity (Shomate equation), specific enthalpy (integral of Shomate equation) and thermal conductivity (quadratic polynomial fit) of a condensed phase of a species. All properties are assumed to be independent of pressure. Hence, the corresponding values at standard pressure conditions, i.e., 101325 Pa, are calculated.

Additionally, `Phase` objects have upper and lower temperature bounds. Enthalpy and thermal conductivity smoothly transition to zero above and below the upper and lower temperature bounds respectively (using the sigmoid function). Specific heat rises in magnitude near the bounds representing the latent heat during phase transition.

## Condensed Species
The header `Condensed-Species.hpp` defines the class `CondensedSpecies` to model the density, specific heat capacity, specific enthalpy and thermal conductivity of a species. It strings together a number of phases of the species. The overall thermo-physical property of the species at a temperature is determined from the sum of that thermo-physical property for all the phases at that temperature. Using the upper and lower bound for the phases guarantees only the thermo-physical property of the predominant phases of the species at that temperature is returned.

## Arrhenius Diffusivity Model
The header `Arrhenius-Diffusivity-Model.hpp` defines the class `ArrheniusDiffusivityModel` to model intermetallic diffusivity as an exponential function of temperature
```
D(T) = D_0 * exp(- E_a / R_u * T)
```
where `D_0` is the pre-exponential factor in m^2/s, `E_a` is the activation energy in J/mol. and `R_u` is the universal gas constant in J/mol.-K.

The class has a `critical_temperature` above which one set of Arrhenius model parameters are used (`pre_exponential_factor_high` and `activation_energy_high`); and another set of parameters (`pre_exponential_factor_low` and `activation_energy_low`) are used for temperatures below critical temperature.

## Core-Shell Particle
The header `Core-Shell-Particle.hpp` defines the class `CoreShellParticle` to represent the macroscopic composition of a core-shell structured particle. All `CoreShellParticle` objects are composed of three `Species` objects, namely, `core_species`, `shell_species` and `product_species`. The thermo-physical properties are determined from the composition fractions of the particle - specific heat capacity and specific enthalpy are mass fraction weighted; and density and thermal conductivity are volume fraction weighted.

The combustion reaction is completed when the mass of one of the core or shell species is depleted (mass fractions is below an user-defined limit).

## Thermal Conductivity of Pellet
The header `Thermal-Conductivity-Pellet.hpp` declares a set of functions to calculate effective thermal conductivity of a heterogenous mixture. Provided models are - 
- Bruggeman (EMT)
- Co-continuous (CC)
- Maxwell-Eucken 1 (ME1)
- Maxwell-Eucken 2 (ME2)
- Combined EMT + ME2 (MEB)
The thermal conductivities input to the functions in this library should have consistent units. The output will then be returned in the same unit.

## Packed Pellet
The header `Packed-Pellet.hpp` defines the class `PackedPellet` to represent a pellet packed with particles. It mostly contains constants and functions essential for calculating thermo-physical properties at pellet level.