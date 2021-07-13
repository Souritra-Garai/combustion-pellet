/**
 * @file Core-Shell-Combustion-Particle.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header defines a class for core shell particle
 * @version 0.1
 * @date 2021-07-08
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __CORE_SHELL_COMBUSTION_PARTICLE__
#define __CORE_SHELL_COMBUSTION_PARTICLE__

#include "thermo-physical-properties/Substance.hpp"

#include <ostream>

template<typename real_t>
class CoreShellCombustionParticle
{
    private :

        real_t mass_fraction_reactant_A;
        real_t mass_fraction_reactant_B;
        real_t mass_fraction_product;

        real_t pre_exponential_factor;
        real_t activation_energy;

    protected :

        Substance<real_t> reactant_A;
        Substance<real_t> reactant_B;
        Substance<real_t> product;

        void setMassFractionCoreMaterial(real_t mass_fraction);
        void setMassFractionShellMaterial(real_t mass_fraction);
        void setMassFractionProductMaterial(real_t mass_fraction);

    public :

        CoreShellCombustionParticle(
            Substance<real_t> core_material,
            Substance<real_t> shell_material,
            Substance<real_t> product_material
        );

        real_t getDensity();
        real_t getHeatCapacity();
        real_t getHeatConductivity();

        real_t getEnthalpy(real_t temperature);
        real_t getDiffusionConstant(real_t temperature);

        real_t getMassFractionCoreMaterial();
        real_t getMassFractionShellMaterial();
        real_t getMassFractionProductMaterial();

        void printProperties(std::ostream &output_stream);
};

#endif