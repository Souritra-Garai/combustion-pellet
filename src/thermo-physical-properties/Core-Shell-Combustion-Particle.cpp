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

#include <math.h>

template<typename real_t>
CoreShellCombustionParticle<real_t>::CoreShellCombustionParticle() {
    mass_fraction_core_material = 0.5;
    mass_fraction_shell_material = 0.5;
    mass_fraction_product_material = 0;
}

template<typename real_t>
void CoreShellCombustionParticle<real_t>::printProperties(
    std::ostream &output_stream
) {
    output_stream << "Density\t\t\t:\t" << getDensity() << "\tkg/m3" << std::endl;
    
    output_stream << "Heat Capacity\t\t:\t" << getHeatCapacity() << "\tJ/kg-K" << std::endl;
    
    output_stream << "Heat Conductivity\t:\t" << getHeatConductivity() << "\tW/m" << std::endl;
    
    output_stream << "Core Material\nMass Fraction\t:\t" << mass_fraction_core_material << std::endl;

    core_material.printProperties(output_stream);

    output_stream << "Shell Material\nMass Fraction\t:\t" << mass_fraction_shell_material << std::endl;

    shell_material.printProperties(output_stream);

    output_stream << "Product Material\nMass Fraction\t:\t" << mass_fraction_product_material << std::endl;

    product_material.printProperties(output_stream);
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getDensity()
{
    return 
        mass_fraction_core_material * core_material.getDensity() +
        mass_fraction_shell_material * shell_material.getDensity() +
        mass_fraction_product_material * product_material.getDensity();
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getHeatCapacity()
{
    return 
        mass_fraction_core_material * core_material.getHeatCapacity() +
        mass_fraction_shell_material * shell_material.getHeatCapacity() +
        mass_fraction_product_material * product_material.getHeatCapacity(); 
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getHeatConductivity()
{
    return
        mass_fraction_core_material * core_material.getHeatConductivity() +
        mass_fraction_shell_material * shell_material.getHeatConductivity() +
        mass_fraction_product_material * product_material.getHeatConductivity();
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getEnthalpy(real_t T)
{
    return
        mass_fraction_core_material * core_material.getEnthalpy(T) +
        mass_fraction_shell_material * shell_material.getEnthalpy(T) +
        mass_fraction_product_material * product_material.getEnthalpy(T);
}

template<typename real_t>
void CoreShellCombustionParticle<real_t>::setUpCoreShellParticle(
    Substance<real_t> core,
    Substance<real_t> shell,
    real_t r_P,
    real_t r_C
) {
    Substance<real_t> core_material(core);
    Substance<real_t> shell_material(shell);

    core_radius = r_C;
    overall_radius = r_P;

    mass = calcMassCoreShellParticle<real_t>(
        core_material,
        shell_material,
        overall_radius,
        core_radius
    );
}

// template<typename real_t> Substance<real_t> CoreShellCombustionParticle<real_t>::core_material = Substance<real_t>(0, 0, 0, 0, 0);
// template<typename real_t> Substance<real_t> CoreShellCombustionParticle<real_t>::shell_material = Substance<real_t>(0, 0, 0, 0, 0);
// template<typename real_t> Substance<real_t> CoreShellCombustionParticle<real_t>::product_material = Substance<real_t>(0, 0, 0, 0, 0);

// template<typename real_t> real_t CoreShellCombustionParticle<real_t>::overall_radius = 0;
// template<typename real_t> real_t CoreShellCombustionParticle<real_t>::core_radius = 0;
// template<typename real_t> real_t CoreShellCombustionParticle<real_t>::mass = 0;

template class CoreShellCombustionParticle<float>;
template class CoreShellCombustionParticle<double>;
template class CoreShellCombustionParticle<long double>;

template<typename real_t>
real_t calcMassCoreShellParticle(
    Substance<real_t> core_material,
    Substance<real_t> shell_material,
    real_t overall_radius,
    real_t core_radius
) {
    return (4.0 * M_PI / 3.0) * (
        core_material.getDensity() * pow(core_radius, 3) +
        shell_material.getDensity() * (pow(overall_radius, 3) - pow(core_radius, 3))
    );
}