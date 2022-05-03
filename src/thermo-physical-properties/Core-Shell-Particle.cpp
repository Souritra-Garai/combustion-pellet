/**
 * @file Core-Shell-Combustion-Particle.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief 
 * @version 0.1
 * @date 2021-07-08
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "thermo-physical-properties/Core-Shell-Particle.hpp"

// Required pow() function and value of \f$ \pi \f$
#include <math.h>

// Instantiate the static private members of the CoreShellParticle class
// Creating the Substance members
template<typename real_t> CondensedSpecies<real_t> * CoreShellParticle<real_t>::_core_material;
template<typename real_t> CondensedSpecies<real_t> * CoreShellParticle<real_t>::_shell_material;
template<typename real_t> CondensedSpecies<real_t> * CoreShellParticle<real_t>::_product_material;
// Creating the real_t members
template<typename real_t> real_t CoreShellParticle<real_t>::_overall_radius = 0;
template<typename real_t> real_t CoreShellParticle<real_t>::_core_radius    = 0;
template<typename real_t> real_t CoreShellParticle<real_t>::_mass = 0;

template<typename real_t> real_t Phase<real_t>::sharpness_coefficient = 0.1;

// Set up the static members of the CoreShellParticle class
template<typename real_t>
void CoreShellParticle<real_t>::setUpCoreShellParticle(
    CondensedSpecies<real_t> &core_material,
    CondensedSpecies<real_t> &shell_material,
    CondensedSpecies<real_t> &product_material,
    real_t overall_radius,
    real_t core_radius
) {
    // Copy the substances to the private static members
    // Copy constructor for CondensedSpecies<real_t> is called
    _core_material    = &core_material;
    _shell_material   = &shell_material;
    _product_material = &product_material;

    // Set the overall and core radii of the core-shell particle
    _overall_radius = overall_radius;
    _core_radius    = core_radius;

    // Calculate the total mass of the core-shell particle
    // and save it in memory in the private static member _mass
    _mass = calcParticleMass();
}

// Calculate and return volume of sphere 
// with static member _core_radius as the radius
template<typename real_t>
real_t CoreShellParticle<real_t>::calcCoreVolume()
{
    // Volume of sphere = \f$ \frac{4}{3} \pi r^3
    return 4.0 * M_PI * pow(_core_radius, 3) / 3.0;
}

// Calculate and return volume fo a shell
// with static members _overall_radius and _core_radius
// as the outer and inner radii respectively
template<typename real_t>
real_t CoreShellParticle<real_t>::calcShellVolume()
{
    // Volume of spherical shell = \f$ \frac{4}{3} \pi \left( r_{outer}^3 - r_{inner}^3 \right) \f$
    return 4.0 * M_PI * (pow(_overall_radius, 3) - pow(_core_radius, 3)) / 3.0;
}

// Calculate and return mass of uniform sphere
// with density of core material
template<typename real_t>
real_t CoreShellParticle<real_t>::calcCoreMass()
{
    // Mass of core = Density of core material \f$ \times \f$ Volume of core
    return _core_material->getDensity(298.15) * calcCoreVolume();
}

// Calculate and return mass of spherical shell
// with uniform density equal to that of shell material
template<typename real_t>
real_t CoreShellParticle<real_t>::calcShellMass()
{
    // Mass of shell = Density of shell material \f$ \times \f$ Volume of shell
    return _shell_material->getDensity(298.15) * calcShellVolume();
}

// Sum and return the masses of core part and the shell part
template<typename real_t>
real_t CoreShellParticle<real_t>::calcParticleMass()
{
    // Mass of Core-Shell Particle = Mass of core + Mass of shell
    return calcCoreMass() + calcShellMass();
}

// Create an instance of CoreShellParticle
// and initialize the mass fraction variables
template<typename real_t>
CoreShellParticle<real_t>::CoreShellParticle()
{
    // Calculate and set mass fractions of the 
    // core, shell and product material
    _mass_fraction_core_material    = calcCoreMass()  / _mass;
    _mass_fraction_shell_material   = calcShellMass() / _mass;
    _mass_fraction_product_material = 0.0;
}

template<typename real_t>
void CoreShellParticle<real_t>::printProperties(
    std::ostream &output_stream
) {
    output_stream << "Core Diameter\t\t:\t" << 2.0 * _core_radius << "\tm" << std::endl;
    output_stream << "Particle Diameter\t:\t" << 2.0 * _overall_radius << "\tm" << std::endl;

    output_stream << "\nMass\t\t\t:\t" << _mass << "\tkg" << std::endl;

    output_stream << "Density\t\t\t:\t" << getDensity(298.15) << "\tkg/m3" << std::endl;
    output_stream << "Heat Capacity\t\t:\t" << getHeatCapacity(298.15) << "\tJ/kg-K" << std::endl;
    output_stream << "Enthalpy\t\t:\t" << getEnthalpy(298.15) << "\tJ/kg" << std::endl;
    output_stream << "Heat Conductivity\t:\t" << getThermalConductivity(298.15) << "\tW/m" << std::endl;
}

template class CoreShellParticle<float>;
template class CoreShellParticle<double>;
template class CoreShellParticle<long double>;

template<typename real_t>
real_t calcMassCoreShellParticle(
    CondensedSpecies<real_t> core_material,
    CondensedSpecies<real_t> shell_material,
    real_t overall_radius,
    real_t core_radius
) {
    return (4.0 * M_PI / 3.0) * (
        core_material.getDensity(298.15) * pow(core_radius, 3) +
        shell_material.getDensity(298.15) * (pow(overall_radius, 3) - pow(core_radius, 3))
    );
}

template
float calcMassCoreShellParticle(
    CondensedSpecies<float> core_material,
    CondensedSpecies<float> shell_material,
    float overall_radius,
    float core_radius
);

template
double calcMassCoreShellParticle(
    CondensedSpecies<double> core_material,
    CondensedSpecies<double> shell_material,
    double overall_radius,
    double core_radius
);

template
long double calcMassCoreShellParticle(
    CondensedSpecies<long double> core_material,
    CondensedSpecies<long double> shell_material,
    long double overall_radius,
    long double core_radius
);