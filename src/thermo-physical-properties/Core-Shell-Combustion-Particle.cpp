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
CoreShellCombustionParticle<real_t>::CoreShellCombustionParticle(
    Substance<real_t> &core_material,
    Substance<real_t> &shell_material,
    Substance<real_t> &product_material,
    real_t &r_P,
    real_t &r_C,
    real_t &m
) : mass(&m),
    core_radius(&r_C),
    overall_radius(&r_P),
    reactant_A(&core_material),
    reactant_B(&shell_material),
    product_AB(&product_material)
{
    mass_fraction_reactant_A = reactant_A->getDensity() * 4 * M_PI * pow(*core_radius, 3) / (3 * (*mass));
    mass_fraction_reactant_B = reactant_B->getDensity() * 4 * M_PI * (pow(*overall_radius, 3) - pow(*core_radius, 3)) / (3 * (*mass));
    mass_fraction_product_AB = 0;
}

template<typename real_t>
void CoreShellCombustionParticle<real_t>::printProperties(
    std::ostream &output_stream
) {
    output_stream << "Density\t\t\t:\t" << getDensity() << "\tkg/m3" << std::endl;
    
    output_stream << "Heat Capacity\t\t:\t" << getHeatCapacity() << "\tJ/kg-K" << std::endl;
    
    output_stream << "Heat Conductivity\t:\t" << getHeatConductivity() << "\tW/m" << std::endl;
    
    output_stream << "Core Material\nMass Fraction\t:\t" << mass_fraction_reactant_A << std::endl;

    reactant_A->printProperties(output_stream);

    output_stream << "Shell Material\nMass Fraction\t:\t" << mass_fraction_reactant_B << std::endl;

    reactant_B->printProperties(output_stream);

    output_stream << "Product Material\nMass Fraction\t:\t" << mass_fraction_product_AB << std::endl;

    product_AB->printProperties(output_stream);
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getDensity()
{
    return 
        mass_fraction_reactant_A * reactant_A->getDensity() +
        mass_fraction_reactant_B * reactant_B->getDensity() +
        mass_fraction_product_AB * product_AB->getDensity();
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getHeatCapacity()
{
    return 
        mass_fraction_reactant_A * reactant_A->getHeatCapacity() +
        mass_fraction_reactant_B * reactant_B->getHeatCapacity() +
        mass_fraction_product_AB * product_AB->getHeatCapacity(); 
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getHeatConductivity()
{
    return
        mass_fraction_reactant_A * reactant_A->getHeatConductivity() +
        mass_fraction_reactant_B * reactant_B->getHeatConductivity() +
        mass_fraction_product_AB * product_AB->getHeatConductivity();
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getEnthalpy(real_t T)
{
    return
        mass_fraction_reactant_A * reactant_A->getEnthalpy(T) +
        mass_fraction_reactant_B * reactant_B->getEnthalpy(T) +
        mass_fraction_product_AB * product_AB->getEnthalpy(T);
}

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

template
float calcMassCoreShellParticle(
    Substance<float> core_material,
    Substance<float> shell_material,
    float overall_radius,
    float core_radius
);

template
double calcMassCoreShellParticle(
    Substance<double> core_material,
    Substance<double> shell_material,
    double overall_radius,
    double core_radius
);

template
long double calcMassCoreShellParticle(
    Substance<long double> core_material,
    Substance<long double> shell_material,
    long double overall_radius,
    long double core_radius
);