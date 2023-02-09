#include <iostream>

#include "thermo-physical-properties/Packed-Pellet.hpp"

int main(int argc, char const *argv[])
{
	PackedPellet<long double> pellet(0.5);
    
    std::cout << pellet.length << std::endl;

    return 0;
}