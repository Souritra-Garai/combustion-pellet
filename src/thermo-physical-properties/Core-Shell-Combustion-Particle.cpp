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

#include "thermo-physical-properties/Core-Shell-Combustion-Particle.hpp"

// Required pow() function and value of \f$ \pi \f$
#include <math.h>

// Instantiate the static private members of the CoreShellCombustionParticle class
// Creating the Substance members
template<typename real_t> Substance<real_t> CoreShellCombustionParticle<real_t>::_core_material    = Substance<real_t>();
template<typename real_t> Substance<real_t> CoreShellCombustionParticle<real_t>::_shell_material   = Substance<real_t>();
template<typename real_t> Substance<real_t> CoreShellCombustionParticle<real_t>::_product_material = Substance<real_t>();
// Creating the real_t members
template<typename real_t> real_t CoreShellCombustionParticle<real_t>::_overall_radius = 0;
template<typename real_t> real_t CoreShellCombustionParticle<real_t>::_core_radius    = 0;
template<typename real_t> real_t CoreShellCombustionParticle<real_t>::_mass = 0;

// Set up the static members of the CoreShellCombustionParticle class
template<typename real_t>
void CoreShellCombustionParticle<real_t>::setUpCoreShellCombustionParticle(
    Substance<real_t> core_material,
    Substance<real_t> shell_material,
    Substance<real_t> product_material,
    real_t overall_radius,
    real_t core_radius
) {
    // Copy the substances to the private static members
    // Copy constructor for Substance<real_t> is called
    _core_material    = core_material;
    _shell_material   = shell_material;
    _product_material = product_material;

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
real_t CoreShellCombustionParticle<real_t>::calcCoreVolume()
{
    // Volume of sphere = \f$ \frac{4}{3} \pi r^3
    return 4.0 * M_PI * pow(_core_radius, 3) / 3.0;
}

// Calculate and return volume fo a shell
// with static members _overall_radius and _core_radius
// as the outer and inner radii respectively
template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::calcShellVolume()
{
    // Volume of spherical shell = \f$ \frac{4}{3} \pi \left( r_{outer}^3 - r_{inner}^3 \right) \f$
    return 4.0 * M_PI * (pow(_overall_radius, 3) - pow(_core_radius, 3)) / 3.0;
}

// Calculate and return mass of uniform sphere
// with density of core material
template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::calcCoreMass()
{
    // Mass of core = Density of core material \f$ \times \f$ Volume of core
    return _core_material.getDensity() * calcCoreVolume();
}

// Calculate and return mass of spherical shell
// with uniform density equal to that of shell material
template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::calcShellMass()
{
    // Mass of shell = Density of shell material \f$ \times \f$ Volume of shell
    return _shell_material.getDensity() * calcShellVolume();
}

// Sum and return the masses of core part and the shell part
template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::calcParticleMass()
{
    // Mass of Core-Shell Particle = Mass of core + Mass of shell
    return calcCoreMass() + calcShellMass();
}

// Create an instance of CoreShellCombustionParticle
// and initialize the mass fraction variables
template<typename real_t>
CoreShellCombustionParticle<real_t>::CoreShellCombustionParticle()
{
    // Calculate and set mass fractions of the 
    // core, shell and product material
    _mass_fraction_core_material    = calcCoreMass()  / _mass;
    _mass_fraction_shell_material   = calcShellMass() / _mass;
    _mass_fraction_product_material = 0.0;
}

// Take weighted sum of densities of core, shell and product materials
// with mass fraction as the weights
template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getDensity()
{
    // \f$ \rho = \sum_{k \in \left\{ Core, Shell, Product \right}} Y_k \rho_k \f$
    return 1.0 / (
        _mass_fraction_core_material    / _core_material.getDensity() +
        _mass_fraction_shell_material   / _shell_material.getDensity() +
        _mass_fraction_product_material / _product_material.getDensity()
    );
}

// Take weighted sum of heat capacities of core, shell and product materials
// with mass fraction as the weights
template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getHeatCapacity()
{
    // \f$ c = \sum_{k \in \left\{ Core, Shell, Product \right}} Y_k c_k \f$
    return 
        _mass_fraction_core_material    * _core_material.getHeatCapacity() +
        _mass_fraction_shell_material   * _shell_material.getHeatCapacity() +
        _mass_fraction_product_material * _product_material.getHeatCapacity(); 
}

// Take weighted sum of heat conductivities of core, shell and product materials
// with mass fraction as the weights
template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getHeatConductivity()
{
    // \f$ \lambda = \sum_{k \in \left\{ Core, Shell, Product \right}} Y_k \lambda_k \f$
    return
        _mass_fraction_core_material    * _core_material.getHeatConductivity() +
        _mass_fraction_shell_material   * _shell_material.getHeatConductivity() +
        _mass_fraction_product_material * _product_material.getHeatConductivity();
}

// Take weighted sum of enthalpies of core, shell and product materials
// at the specified temperature with mass fraction as the weights
template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getEnthalpy(real_t T)
{
    // \f$ h \left( T \right) 
    // = \sum_{k \in \left\{ Core, Shell, Product \right}} Y_k h_k \left( T \right)\f$
    return
        _mass_fraction_core_material    * _core_material.getEnthalpy(T) +
        _mass_fraction_shell_material   * _shell_material.getEnthalpy(T) +
        _mass_fraction_product_material * _product_material.getEnthalpy(T);
}

template<typename real_t>
void CoreShellCombustionParticle<real_t>::printProperties(
    std::ostream &output_stream
) {
    output_stream << "Core Diameter\t\t:\t" << 2.0 * _core_radius << "\tm" << std::endl;
    output_stream << "Particle Diameter\t:\t" << 2.0 * _overall_radius << "\tm" << std::endl;

    output_stream << "\nMass\t\t\t:\t" << _mass << "\tkg" << std::endl;

    output_stream << "Density\t\t\t:\t" << getDensity() << "\tkg/m3" << std::endl;
    output_stream << "Heat Capacity\t\t:\t" << getHeatCapacity() << "\tJ/kg-K" << std::endl;
    output_stream << "Heat Conductivity\t:\t" << getHeatConductivity() << "\tW/m" << std::endl;

    output_stream << "\nCore Material\nMass Fraction\t:\t" << _mass_fraction_core_material << std::endl;
    _core_material.printProperties(output_stream);

    output_stream << "\nShell Material\nMass Fraction\t:\t" << _mass_fraction_shell_material << std::endl;
    _shell_material.printProperties(output_stream);

    output_stream << "\nProduct Material\nMass Fraction\t:\t" << _mass_fraction_product_material << std::endl;
    _product_material.printProperties(output_stream);
}

template class CoreShellCombustionParticle<float>;
template class CoreShellCombustionParticle<double>;
template class CoreShellCombustionParticle<long double>;