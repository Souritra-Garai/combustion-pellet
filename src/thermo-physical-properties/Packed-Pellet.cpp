#include "thermo-physical-properties/Packed-Pellet.hpp"

#include <math.h>

template<typename real_t> real_t PackedPellet<real_t>::_length = 0;
template<typename real_t> real_t PackedPellet<real_t>::_diameter = 0;

template<typename real_t> real_t PackedPellet<real_t>::_convective_heat_transfer_coefficient = 0;
template<typename real_t> real_t PackedPellet<real_t>::_radiative_emissivity = 0;
template<typename real_t> real_t PackedPellet<real_t>::_ambient_temperature = 0;

template<typename real_t> Substance<real_t> PackedPellet<real_t>::_degassing_fluid;

template<typename real_t>
void PackedPellet<real_t>::setPelletDimensions(
    real_t length,
    real_t diameter
) {
    _length = length;
    _diameter = diameter;
}

template<typename real_t>
void PackedPellet<real_t>::setAmbientHeatLossProperties(
    real_t convective_heat_transfer_coefficient,
    real_t radiative_emissivity,
    real_t ambient_temperature
) {
    _convective_heat_transfer_coefficient = convective_heat_transfer_coefficient;
    _radiative_emissivity = radiative_emissivity;
    _ambient_temperature = ambient_temperature;
}

template<typename real_t>
void PackedPellet<real_t>::setDegassingFluid(Substance<real_t> degassing_fluid)
{
    _degassing_fluid = degassing_fluid;
}

template<typename real_t>
real_t PackedPellet<real_t>::calcParticleMassFractions(
    real_t particle_volume_fractions
) {
    real_t initial_particle_density = CoreShellCombustionParticle<real_t>().getDensity();
    
    return 
    
        particle_volume_fractions * initial_particle_density / (
            
            particle_volume_fractions           * initial_particle_density + 
            
            (1.0 - particle_volume_fractions)   * _degassing_fluid.getDensity()
        );
}

template<typename real_t>
real_t PackedPellet<real_t>::calcDensity(real_t particle_volume_fractions)
{
    real_t initial_particle_density = CoreShellCombustionParticle<real_t>().getDensity();

    return particle_volume_fractions * initial_particle_density +
    (1.0 - particle_volume_fractions) * _degassing_fluid.getDensity();
}

template<typename real_t>
PackedPellet<real_t>::PackedPellet(
    real_t particle_volume_fractions
) : _density(calcDensity(particle_volume_fractions)),
    _particle_volume_fractions(particle_volume_fractions),
    _particle_mass_fractions(calcParticleMassFractions(particle_volume_fractions))
{ ; }

template<typename real_t>
real_t PackedPellet<real_t>::getDensity()
{
    return _density;
}

template<typename real_t>
real_t PackedPellet<real_t>::getHeatCapacity(CoreShellDiffusion<real_t> &particle)
{
    return _particle_mass_fractions * particle.getHeatCapacity() +
    (1.0 - _particle_mass_fractions) * _degassing_fluid.getHeatCapacity();
}

template<typename real_t>
real_t PackedPellet<real_t>::getHeatConductivity(CoreShellDiffusion<real_t> &particle)
{    
    real_t lambda_p = particle.getHeatConductivity();
    real_t lambda_f = _degassing_fluid.getHeatConductivity();

    real_t D =
        pow((lambda_p / lambda_f) * (3*_particle_volume_fractions - 1), 2) + 
        pow(2-3*_particle_volume_fractions, 2) + 
        2 * (lambda_p / lambda_f) * (
            2 + 
            9*_particle_volume_fractions - 
            9*_particle_volume_fractions*_particle_volume_fractions
        );

    return (lambda_f/4) * ((3*_particle_volume_fractions - 1)*(lambda_p / lambda_f) + 2 - 3*_particle_volume_fractions + sqrt(D));
}

template<typename real_t>
real_t PackedPellet<real_t>::getEnthalpy(CoreShellDiffusion<real_t> &particle, real_t T)
{
    return _particle_mass_fractions * particle.getEnthalpy(T) +
    (1.0 - _particle_mass_fractions) * _degassing_fluid.getEnthalpy(T);
}

template<typename real_t>
void PackedPellet<real_t>::printProperties(std::ostream &output_stream)
{
    output_stream << "Length\t:\t" << _length << "\tm" << std::endl;
    output_stream << "Diameter\t:\t" << _diameter << "\tm" << std::endl;

    output_stream << "\nAmbient Heat Loss Parameters" << std::endl;
    output_stream << "Convective Heat Transfer Coefficient\t:\t" << _convective_heat_transfer_coefficient << "\tW/m2-K" << std::endl;
    output_stream << "Radiative emissivity\t:\t" << _radiative_emissivity << std::endl;
    output_stream << "Ambient Temperature\t:\t" << _ambient_temperature << "\tK" << std::endl;

    CoreShellDiffusion<real_t> particle;

    output_stream << "Density\t\t\t\t:\t" << getDensity() << "\tkg/m3" << std::endl;
    output_stream << "Heat Capacity\t\t\t:\t" << getHeatCapacity(particle) << "\tJ/kg-K" << std::endl;
    output_stream << "Heat Conductivity\t\t:\t" << getHeatConductivity(particle) << "\tW/m-K" << std::endl;

    output_stream << "\nParticle Properties" << std::endl;
    output_stream << "Particle Volume Fraction\t:\t" << _particle_volume_fractions << std::endl;
    output_stream << "Particle Mass Fraction\t:\t" << _particle_mass_fractions << std::endl;
    particle.printProperties(output_stream);

    output_stream << "\nDegassing Fluid" << std::endl;
    output_stream << "Degassing Fluid Volume Fraction\t:\t" << 1.0 - _particle_volume_fractions << std::endl;
    output_stream << "Degassing Fluid Mass Fraction\t:\t" << 1.0 - _particle_mass_fractions << std::endl;
    _degassing_fluid.printProperties(output_stream);
}

template class PackedPellet<float>;
template class PackedPellet<double>;
template class PackedPellet<long double>;