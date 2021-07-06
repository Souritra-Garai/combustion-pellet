/**
 * @file Substance_Example.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Example to test functions of Substance class
 * @version 0.1
 * @date 2021-07-06
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include <iostream>

#include "thermo-physical-properties/Substance.hpp"

Substance<float> Water(1000, 4180, 18E-3, 10);

int main(int argc, char const *argv[])
{
    Water.printProperties(std::cout);
    
    return 0;
}
