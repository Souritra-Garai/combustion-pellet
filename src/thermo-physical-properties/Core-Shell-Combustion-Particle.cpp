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
    Substance<real_t> core_material,
    Substance<real_t> shell_material,
    Substance<real_t> product_material
) : reactant_A(core_material),
    reactant_B(shell_material),
    product(product_material)
{
    setMassFractionCoreMaterial(0.5);
    setMassFractionShellMaterial(0.5);
    setMassFractionProductMaterial(0);
}

template<typename real_t>
void CoreShellCombustionParticle<real_t>::setMassFractionCoreMaterial(
    real_t mass_fraction
) {
    mass_fraction_reactant_A = mass_fraction;
}

template<typename real_t>
void CoreShellCombustionParticle<real_t>::setMassFractionShellMaterial(
    real_t mass_fraction
) {
    mass_fraction_reactant_B = mass_fraction;
}

template<typename real_t>
void CoreShellCombustionParticle<real_t>::setMassFractionProductMaterial(
    real_t mass_fraction
) {
    mass_fraction_product = mass_fraction;
}

template<typename real_t>
void CoreShellCombustionParticle<real_t>::printProperties(
    std::ostream &output_stream
) {
    output_stream << "Density\t\t\t:\t" << getDensity() << "\tkg/m3" << std::endl;
    
    output_stream << "Heat Capacity\t\t:\t" << getHeatCapacity() << "\tJ/kg-K" << std::endl;
    
    output_stream << "Heat Conductivity\t:\t" << getHeatConductivity() << "\tW/m" << std::endl;
    
    output_stream << "Core Material\nMass Fraction\t:\t" << getMassFractionCoreMaterial() << std::endl;

    reactant_A.printProperties(output_stream);

    output_stream << "Shell Material\nMass Fraction\t:\t" << getMassFractionShellMaterial() << std::endl;

    reactant_B.printProperties(output_stream);

    output_stream << "Product Material\nMass Fraction\t:\t" << getMassFractionProductMaterial() << std::endl;

    product.printProperties(output_stream);
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getDensity()
{
    return 
        getMassFractionCoreMaterial()    * reactant_A.getDensity() +
        getMassFractionShellMaterial()   * reactant_B.getDensity() +
        getMassFractionProductMaterial() * product.getDensity();
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getDiffusionConstant(real_t T)
{
    return pre_exponential_factor * exp(- activation_energy / (8.314 * T));
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getHeatCapacity()
{
    return 
        getMassFractionCoreMaterial()    * reactant_A.getHeatCapacity() +
        getMassFractionShellMaterial()   * reactant_B.getHeatCapacity() +
        getMassFractionProductMaterial() * product.getHeatCapacity(); 
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getHeatConductivity()
{
    return
        getMassFractionCoreMaterial()    * reactant_A.getHeatConductivity() +
        getMassFractionShellMaterial()   * reactant_B.getHeatConductivity() +
        getMassFractionProductMaterial() * product.getHeatConductivity();
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getMassFractionCoreMaterial()
{
    return mass_fraction_reactant_A;
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getMassFractionProductMaterial()
{
    return mass_fraction_product;
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getMassFractionShellMaterial()
{
    return mass_fraction_reactant_B;
}

template<typename real_t>
real_t CoreShellCombustionParticle<real_t>::getEnthalpy(real_t T)
{
    return
        getMassFractionCoreMaterial()    * reactant_A.getEnthalpy(T) +
        getMassFractionShellMaterial()   * reactant_B.getEnthalpy(T) +
        getMassFractionProductMaterial() * product.getEnthalpy(T);
}

template class CoreShellCombustionParticle<float>;
template class CoreShellCombustionParticle<double>;
template class CoreShellCombustionParticle<long double>;