#include <iostream>

#include "thermo-physical-properties/Packed-Pellet.hpp"

int main(int argc, char const *argv[])
{
	PackedPellet pellet(0.5);
    
    std::cout << "Overall Particle Density :\t" << pellet.overall_particle_density << "\tkg/m3" << std::endl;

    return 0;
}