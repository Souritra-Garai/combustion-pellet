/**
 * @file Core-Shell-Diffusion.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Implementation of CoreShellDiffusion class
 * @version 0.1
 * @date 2021-07-16
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "pde-problems/Core-Shell-Diffusion.hpp"

#include <math.h>

// Instatiating static member variables of CoreShellDiffusion Class
// Time step interval duration
template<typename real_t> real_t CoreShellDiffusion<real_t>::_delta_t = 0;
// Grid size
template<typename real_t> size_t CoreShellDiffusion<real_t>::_n = 1;
template<typename real_t> real_t CoreShellDiffusion<real_t>::_delta_r = 0;

// Defining static member functions of CoreShellDiffusion Class
// Setting time step interval duration
template<typename real_t>
void CoreShellDiffusion<real_t>::setTimeStep(real_t Dt) {_delta_t = Dt;}
// Setting grid size
template<typename real_t>
void CoreShellDiffusion<real_t>::setGridSize(size_t n)
{
    // Set number of grid points
    // Include the boundary \f$ r = 0 \f$ and \f$ r = r_p \f$ points
    _n = n;
    // Total n grid points (inclusive of boundary points) implies n-1 divisions
    // Thus \f$ \Delta r = r_P / (n-1) \f$
    _delta_r = CoreShellCombustionParticle<real_t>::_overall_radius / (real_t) (_n - 1);
}

template<typename real_t>
real_t CoreShellDiffusion<real_t>::getRadialCoordinate(size_t i)
{
    // Radial distance from origin of the grid point \f$ i \f$
    // \f$ r_i = r_p \cdot \frac{i}{n-1}
    return CoreShellCombustionParticle<real_t>::_overall_radius * ((real_t) i) / ((real_t) (_n-1));
    // For loops will iterate from \f$ i = 0 \f$ to \f$ i = n-1 \f$
}

template<typename real_t>
CoreShellDiffusion<real_t>::CoreShellDiffusion() : // Mem initialization list
    CoreShellCombustionParticle<real_t>(),  // Call the constructor of the base class
    // Initialize solvers
    _solver_A(_n),
    _solver_B(_n)
{
    // Allocate memory for concentration profiles
    _concentration_array_A = new real_t[_n];
    _concentration_array_B = new real_t[_n];

    // Parallelize for loops
    #pragma omp parallel for

        // Loop over the grid points to initialize the concentration of the
        // substances A and B
        for (size_t i = 0; i < _n; i++)
        {
            // Only substance A is present in the core (\f$ r < r_C \f$)
            if (getRadialCoordinate(i) < this->_core_radius)
            {
                // Concentration of pure substance A is its molar density
                _concentration_array_A[i] = this->_core_material.getMolarDensity();
                _concentration_array_B[i] = 0;
            }
            // Only substance B is present in the shell (\f$ r > r_V \f$)
            else
            {
                // Concentration of pure substance B is its molar density
                _concentration_array_A[i] = 0;
                _concentration_array_B[i] = this->_shell_material.getMolarDensity();
            }
        }
}

template<typename real_t>
CoreShellDiffusion<real_t>::~CoreShellDiffusion()
{
    // Deallocate the memory allocated for concentration profiles
    delete [] _concentration_array_A;
    delete [] _concentration_array_B;
}

template<typename real_t>
inline real_t CoreShellDiffusion<real_t>::getRxnConcA(size_t i)
{
    // If molar concentration of A < molar concentration of B
    // All of A exists as AB => 0 concentration of pure A
    // Else all of B exists as AB
    // => molar concentration of AB = molar concentration of B (w.r.t. diffusion)
    // => molar concentration of A = molar concentration of A (w.r.t. diffusion) - molar concentration of B (w.r.t. diffusion)
    return std::max(_concentration_array_A[i]  - _concentration_array_B[i], (real_t) 0.0);
}

template<typename real_t>
inline real_t CoreShellDiffusion<real_t>::getRxnConcB(size_t i)
{
    // If molar concentration of B < molar concentration of A
    // All of B exists as AB => 0 concentration of pure B
    // Else all of A exists as AB
    // => molar concentration of AB = molar concentration of A (w.r.t. diffusion)
    // => molar concentration of B = molar concentration of B (w.r.t. diffusion) - molar concentration of A (w.r.t. diffusion)
    return std::max(_concentration_array_B[i] - _concentration_array_A[i], (real_t) 0.0);
}

template<typename real_t>
inline real_t CoreShellDiffusion<real_t>::getRxnConcAB(size_t i)
{
    // If molar concentration of A < molar concentration of B
    // All of A exists as AB 
    // => molar concentration of AB = molar concentration of A (w.r.t. diffusion)
    // Else all of B exists as AB
    // => molar concentration of AB = molar concentration of B (w.r.t. diffusion)
    return std::min(_concentration_array_A[i], _concentration_array_B[i]);
}

template<typename real_t>
real_t CoreShellDiffusion<real_t>::numCalcCoreMass()
{
    /**
     * Using the trapezoidal rule
     * \f{equation}{ m_A = \left( MW \right)_A \sum_{i=0}^{n-1} 4 \pi \cdot \frac{1}{2} \left( C_{A,i} r_i^2 + C_{A,i+1} r_{i+1}^2 \right) \Delta r \f}
     * \f{equation}{     = 4 \pi \left( MW \right)_A \Delta r \left\{ \frac{1}{2} C_{A,n}r_p^2 + \sum_{i=1}^{n-1} C_{A,i} r_i^2 \right\} \f}
     */
    
    // Initialise variable to store the integral
    // with the value \f$ \frac{1}{2} C_{A,n}r_n^2 \f$
    real_t sum = 0.5 * _concentration_array_A[_n-1] * pow(this->_overall_radius, 2);

    // Parallelize the loop
    // create private copies for the variable sum
    // and add all of them together at the end (reduction clause)
    #pragma omp parallel for reduction(+:sum)
        // Iterate over all the grid points but the last
        for (size_t i = 0; i < _n-1; i++)
        {
            // Cumulatively add the term \f$ C_{A,i} r_i^2 \f$
            sum += _concentration_array_A[i] * pow(getRadialCoordinate(i), 2);
        }
    // Multiply the constant factors in the summation
    return 4.0 * M_PI * this->_core_material.getMolecularWeight() * _delta_r * sum;
}

template<typename real_t>
real_t CoreShellDiffusion<real_t>::numCalcShellMass()
{
    /**
     * Using the trapezoidal rule
     * \f{equation}{ m_A = \left( MW \right)_A \sum_{i=0}^{n-1} 4 \pi \cdot \frac{1}{2} \left( C_{A,i} r_i^2 + C_{A,i+1} r_{i+1}^2 \right) \Delta r \f}
     * \f{equation}{     = 4 \pi \left( MW \right)_A \Delta r \left\{ \frac{1}{2} C_{A,n}r_p^2 + \sum_{i=1}^{n-1} C_{A,i} r_i^2 \right\} \f}
     */
    
    // Initialise variable to store the integral
    // with the value \f$ \frac{1}{2} C_{A,n}r_n^2 \f$
    real_t sum = 0.5 * _concentration_array_B[_n-1] * pow(this->_overall_radius, 2);

    // Parallelize the loop
    // create private copies for the variable sum
    // and add all of them together at the end (reduction clause)
    #pragma omp parallel for reduction(+:sum)
        // Iterate over all the grid points but the last
        for (size_t i = 0; i < _n-1; i++)
        {
            // Cumulatively add the term \f$ C_{A,i} r_i^2 \f$
            sum += _concentration_array_B[i] * pow(getRadialCoordinate(i), 2);
        }
    // Multiply the constant factors in the summation
    return 4.0 * M_PI * this->_shell_material.getMolecularWeight() * _delta_r * sum;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::calcMassFractions()
{
    real_t Y_A  = 0.5 * getRxnConcA(_n-1)  * pow(this->_overall_radius, 2);
    real_t Y_B  = 0.5 * getRxnConcB(_n-1)  * pow(this->_overall_radius, 2);
    real_t Y_AB = 0.5 * getRxnConcAB(_n-1) * pow(this->_overall_radius, 2);
    
    #pragma omp parallel for reduction(+:Y_A, Y_B, Y_AB)
    
        for (size_t i = 1; i < _n-1; i++)
        {
            real_t r2 = pow(getRadialCoordinate(i), 2);
            
            Y_A  += getRxnConcA(i)  * r2;
            Y_B  += getRxnConcB(i)  * r2;
            Y_AB += getRxnConcAB(i) * r2;
        }
    
    this->_mass_fraction_core_material    = 4.0 * M_PI * _delta_r * Y_A  * this->_core_material.getMolecularWeight()    / this->_mass;
    this->_mass_fraction_shell_material   = 4.0 * M_PI * _delta_r * Y_B  * this->_shell_material.getMolecularWeight()   / this->_mass;
    this->_mass_fraction_product_material = 4.0 * M_PI * _delta_r * Y_AB * this->_product_material.getMolecularWeight() / this->_mass;

    // Y_A  *= 4.0 * M_PI * _delta_r * this->_core_material.getMolecularWeight()    / this->_mass;
    // Y_B  *= 4.0 * M_PI * _delta_r * this->_shell_material.getMolecularWeight()   / this->_mass;
    // Y_AB *= 4.0 * M_PI * _delta_r * this->_product_material.getMolecularWeight() / this->_mass;

    // float sum = Y_A + Y_B + Y_AB;

    // this->_mass_fraction_core_material    = Y_A  / sum;
    // this->_mass_fraction_shell_material   = Y_B  / sum;
    // this->_mass_fraction_product_material = Y_AB / sum;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::setUpEquations(real_t D)
{
    real_t coefficient1 = D / pow(_delta_r, 2);
    real_t coefficient2 = 1.0 / _delta_t;

    _solver_A.setEquationFirstRow(1, -1, 0);
    _solver_B.setEquationFirstRow(1, -1, 0);

    #pragma omp parallel for
    
        for (size_t i = 1; i < _n - 1; i++)
        {
            real_t coefficient3 = coefficient1 * pow(getRadialCoordinate(i+1) / getRadialCoordinate(i), 2);

            real_t e = - coefficient3;
            real_t f = coefficient1 + coefficient2 + coefficient3;
            real_t g = - coefficient1;

            _solver_A.setEquation(i, g, f, e, coefficient2 * _concentration_array_A[i]);
            _solver_B.setEquation(i, g, f, e, coefficient2 * _concentration_array_B[i]);
        }

    _solver_A.setEquationLastRow(-1, 1, 0);
    _solver_B.setEquationLastRow(-1, 1, 0);
}

template<typename real_t>
void CoreShellDiffusion<real_t>::solveEquations()
{
    _solver_A.getSolution(_concentration_array_A);
    _solver_B.getSolution(_concentration_array_B);

    calcMassFractions();
}

template<typename real_t>
void CoreShellDiffusion<real_t>::printConcentrationProfileA(std::ostream &output_stream, char delimiter)
{
    for (size_t i = 0; i < _n-1; i++) output_stream << _concentration_array_A[i] << delimiter;

    output_stream << _concentration_array_A[_n-1] << std::endl;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::printConcentrationProfileB(std::ostream &output_stream, char delimiter)
{
    for (size_t i = 0; i < _n-1; i++) output_stream << _concentration_array_B[i] << delimiter;

    output_stream << _concentration_array_B[_n-1] << std::endl;
}

template class CoreShellDiffusion<float>;
template class CoreShellDiffusion<double>;
template class CoreShellDiffusion<long double>;