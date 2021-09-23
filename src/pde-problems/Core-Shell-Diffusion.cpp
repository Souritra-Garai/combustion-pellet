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

// Required for pow function
#include <math.h>

/******************************************************************************************************************/
// Instatiating static member variables of CoreShellDiffusion Class

// Time step interval duration
template<typename real_t> real_t CoreShellDiffusion<real_t>::_delta_t = 0;
// Grid size
template<typename real_t> size_t CoreShellDiffusion<real_t>::_n = 2;
template<typename real_t> real_t CoreShellDiffusion<real_t>::_delta_r = 0;
/******************************************************************************************************************/

/******************************************************************************************************************/
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
void CoreShellDiffusion<real_t>::printConfiguration(std::ostream &output_stream)
{
    output_stream << "Time Step\t:\t" << _delta_t << "\ts" << std::endl;

    output_stream << "Number of Grid Points\t:\t" << _n << std::endl;

    output_stream << "Grid Size\t:\t" << _delta_r << "\tm" << std::endl;
}
/******************************************************************************************************************/

/******************************************************************************************************************/
// Constructors and destructors

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

    initializeParticle();
}

template<typename real_t>
CoreShellDiffusion<real_t>::~CoreShellDiffusion()
{
    // Deallocate the memory allocated for concentration profiles
    delete [] _concentration_array_A;
    delete [] _concentration_array_B;
}
/******************************************************************************************************************/

/******************************************************************************************************************/
// Defining private member functions

template<typename real_t>
void CoreShellDiffusion<real_t>::calcRxnMassFractions()
{
    // \f$ Y_A  = 4 \pi \Delta r \cdot \frac{\left(MW\right)_{A}}{m_P}  \left\{ \frac{1}{2} \max{\left(C_{A,N} - C_{B,N}, 0\right)} \cdot r_P^2 + \sum_{i=1}^{N-1} \max{\left(C_{A,i} - C_{B,i}, 0\right)} \cdot r_i^2 \right\} \f$
    // \f$ Y_B  = 4 \pi \Delta r \cdot \frac{\left(MW\right)_{B}}{m_P}  \left\{ \frac{1}{2} \max{\left(C_{B,N} - C_{A,N}, 0\right)} \cdot r_P^2 + \sum_{i=1}^{N-1} \max{\left(C_{B,i} - C_{A,i}, 0\right)} \cdot r_i^2 \right\} \f$
    // \f$ Y_AB = 4 \pi \Delta r \cdot \frac{\left(MW\right)_{AB}}{m_P} \left\{ \frac{1}{2} \min{\left(C_{A,N}, C_{B,N}\right)} \cdot r_P^2 + \sum_{i=1}^{N-1} \min{\left(C_{A,i}, C_{B,i}\right)} \cdot r_i^2 \right\} \f$
    
    // Performing the above summation

    // Initialising the sums with the term for grid point # N    
    
    // Calculate \f$ r^2 \f$
    real_t r2 = pow(this->_overall_radius, 2);
    // Intialise the sums without the common multiplicative factors
    real_t Y_A  = 0.5 * getRxnConcA(_n-1)  * pow(this->_overall_radius, 2);
    real_t Y_B  = 0.5 * getRxnConcB(_n-1)  * pow(this->_overall_radius, 2);
    real_t Y_AB = 0.5 * getRxnConcAB(_n-1) * pow(this->_overall_radius, 2);
    
    // Parallelizing the summation loop
    // Create private copies for Y_A, Y_B and Y_AB; and
    // sum the copies at the end
    #pragma omp parallel for reduction(+:Y_A, Y_B, Y_AB) private(r2) num_threads(4) schedule(static, 250)
        // Summing over all grid points except the first and the last
        for (size_t i = 1; i < _n-1; i++)
        {
            // Calculate \f$ r^2 \f$
            r2 = pow(getRadialCoordinate(i), 2);
            // Add the specific terms for the summation
            Y_A  += getRxnConcA(i)  * r2;
            Y_B  += getRxnConcB(i)  * r2;
            Y_AB += getRxnConcAB(i) * r2;
        }
    
    // Update the reaction mass fractions of the core-shell combustion particle

    // Multiply the summations with the required factors and divide
    // by the total mass of the particle calculated during initialisation
    this->_mass_fraction_core_material    = 4.0 * M_PI * _delta_r * Y_A  * this->_core_material->getMolarMass()    / this->_mass;
    this->_mass_fraction_shell_material   = 4.0 * M_PI * _delta_r * Y_B  * this->_shell_material->getMolarMass()   / this->_mass;
    this->_mass_fraction_product_material = 4.0 * M_PI * _delta_r * Y_AB * this->_product_material->getMolarMass() / this->_mass;

    // // Multiply the summations with the required factors
    // Y_A  *= 4.0 * M_PI * _delta_r * this->_core_material->getMolarMass();
    // Y_B  *= 4.0 * M_PI * _delta_r * this->_shell_material->getMolarMass();
    // Y_AB *= 4.0 * M_PI * _delta_r * this->_product_material->getMolarMass();
    
    // // Add the summations to get the total mass
    // float sum = Y_A + Y_B + Y_AB;
    // // Divide by the sum of the summations to get the mass fractions
    // this->_mass_fraction_core_material    = Y_A  / sum;
    // this->_mass_fraction_shell_material   = Y_B  / sum;
    // this->_mass_fraction_product_material = Y_AB / sum;
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
/******************************************************************************************************************/

/******************************************************************************************************************/
// Defining public member functions

template<typename real_t>
void CoreShellDiffusion<real_t>::initializeParticle()
{
    // Parallelize for loops
    #pragma omp parallel for num_threads(4) schedule(static, 250)

        // Loop over the grid points to initialize the concentration of the
        // substances A and B
        for (size_t i = 0; i < _n; i++)
        {
            // Only substance A is present in the core (\f$ r < r_C \f$)
            if (getRadialCoordinate(i) < this->_core_radius)
            {
                // Concentration of pure substance A is its molar density
                _concentration_array_A[i] = this->_core_material->getMolarDensity(298.15);
                _concentration_array_B[i] = 0;
            }
            // Only substance B is present in the shell (\f$ r > r_V \f$)
            else
            {
                // Concentration of pure substance B is its molar density
                _concentration_array_A[i] = 0;
                _concentration_array_B[i] = this->_shell_material->getMolarDensity(298.15);
            }
        }
}

template<typename real_t>
void CoreShellDiffusion<real_t>::setUpEquations(real_t D)
{
    // \f$ \Rightarrow -
    // \mathcal{D} \cdot \left( \frac{r_{i+1}}{\Delta r \cdot r_i} \right)^2 \cdot C_{k,i+1}^n +
    // \left\{ \frac{1}{\Delta t} + \mathcal{D} \cdot \left( \frac{r_{i+1}}{\Delta r \cdot r_i} \right)^2 + \frac{\mathcal{D}}{\left(\Delta r\right)^2} \right\} \cdot C_{k,i}^n \\ -
    // \frac{\mathcal{D}}{\left(\Delta r\right)^2} \cdot C_{k,i-1}^n
    // = \frac{C_{k,i}^{n-1}}{\Delta t} \f$

    // We can calculate the coefficients in three parts to avoid repetitive calculations

    // Coefficient 1 = \f$ \frac{\mathcal{D\left(T\right)}}{\left(\Delta r\right)^2} \f$
    real_t coefficient1 = D / pow(_delta_r, 2);
    // Coefficient 2 = \f$ \frac{1}{\Delta t} \f$
    real_t coefficient2 = 1.0 / _delta_t;
    // Coefficient 3 = \f$ \mathcal{D} \cdot \left( \frac{r_{i+1}}{\Delta r \cdot r_i} \right)^2 \f$
    real_t coefficient3;

    // Set zero flux boundary condition at grid point # 0
    // \f$ C_{k,1}^n - C_{k,0}^n = 0 \f$
    _solver_A.setEquationFirstRow(1, -1, 0);
    _solver_B.setEquationFirstRow(1, -1, 0);

    // Parallely setup the matrix for the linear algebra solver
    // As coefficient 3 is different for every grid point,
    // it needs to be calculated for each row 
    // and hence every thread needs a private copy
    #pragma omp parallel for private(coefficient3) num_threads(4) schedule(static, 250)
        // Iteratively set up the discretized governing diffusion equation
        // for all grid points except the grid points # 0 and N, where
        // zero flux boundary conditions apply
        for (size_t i = 1; i < _n - 1; i++)
        {
            // Calculating coefficient 3 for the current grid point
            coefficient3 = coefficient1 * pow(getRadialCoordinate(i+1) / getRadialCoordinate(i), 2);

            // Set up row for matrix equation representing 
            // diffusion of substance A
            _solver_A.setEquation(
                i, 
                - coefficient1, 
                coefficient1 + coefficient2 + coefficient3, 
                - coefficient3, 
                coefficient2 * _concentration_array_A[i]
            );

            // Set up row for matrix equation representing 
            // diffusion of substance B
            _solver_B.setEquation(
                i, 
                - coefficient1, 
                coefficient1 + coefficient2 + coefficient3, 
                - coefficient3, 
                coefficient2 * _concentration_array_B[i]
            );
        }

    // Set zero flux boundary condition at grid point # N
    // C_{k,N}^n - C_{k,N-1}^n = 0
    _solver_A.setEquationLastRow(-1, 1, 0);
    _solver_B.setEquationLastRow(-1, 1, 0);
}

template<typename real_t>
void CoreShellDiffusion<real_t>::solveEquations()
{
    // Parallelly solve the two linear algebra problems
    #pragma omp parallel sections
    {
        // Parallel section for substance A
        #pragma omp section
            // Solve the matrix equation representing diffusion of substance A
            // and save the solution in array A
            _solver_A.getSolution(_concentration_array_A);

        // Parallel section for substance B
        #pragma omp section
            // Solve the matrix equation representing diffusion of substance B
            // and save the solution in array B
            _solver_B.getSolution(_concentration_array_B);
    }

    // Calculate the reaction mass fractions
    // and update the same for the core-shell combustion particle
    calcRxnMassFractions();
}

template<typename real_t>
real_t CoreShellDiffusion<real_t>::getDiffusionMassA()
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
    return 4.0 * M_PI * this->_core_material->getMolarMass() * _delta_r * sum;
}

template<typename real_t>
real_t CoreShellDiffusion<real_t>::getDiffusionMassB()
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
    return 4.0 * M_PI * this->_shell_material->getMolarMass() * _delta_r * sum;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::copyFrom(CoreShellDiffusion<real_t> &diffusion_problem)
{
    // Parallelize the copying operation
    #pragma omp parallel for num_threads(4) schedule(static, 250)
        // For each element in the concentration array
        // copy the corresponding values
        for (size_t i = 0; i < _n; i++)
        {
            _concentration_array_A[i] = diffusion_problem._concentration_array_A[i];
            _concentration_array_B[i] = diffusion_problem._concentration_array_B[i];
        }

    // Update reaction mass fractions of the particle
    this->_mass_fraction_core_material    = diffusion_problem._mass_fraction_core_material;
    this->_mass_fraction_shell_material   = diffusion_problem._mass_fraction_shell_material;
    this->_mass_fraction_product_material = diffusion_problem._mass_fraction_product_material;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::copyTo(CoreShellDiffusion<real_t> &diffusion_problem)
{
    // Parallelize the copying operation
    #pragma omp parallel for num_threads(4) schedule(static, 250)
        // For each element in the concentration array
        // copy the corresponding values
        for (size_t i = 0; i < _n; i++)
        {
            diffusion_problem._concentration_array_A[i] = _concentration_array_A[i];
            diffusion_problem._concentration_array_B[i] = _concentration_array_B[i];
        }

    // Update reaction mass fractions of the particle
    diffusion_problem._mass_fraction_core_material    = this->_mass_fraction_core_material;
    diffusion_problem._mass_fraction_shell_material   = this->_mass_fraction_shell_material;
    diffusion_problem._mass_fraction_product_material = this->_mass_fraction_product_material;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::printConcentrationProfileA(std::ostream &output_stream, char delimiter, real_t curr_time)
{
    // First print the current time, followed by the temperature profile
    output_stream << curr_time << delimiter;

    // Shift the array elements into the output stream followed by a delimter
    for (size_t i = 0; i < _n-1; i++) output_stream << _concentration_array_A[i] << delimiter;
    // For the last element shift an endline instead of the delimiter
    output_stream << _concentration_array_A[_n-1] << std::endl;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::printConcentrationProfileB(std::ostream &output_stream, char delimiter, real_t curr_time)
{
    // First print the current time, followed by the temperature profile
    output_stream << curr_time << delimiter;
    
    // Shift the array elements into the output stream followed by a delimter
    for (size_t i = 0; i < _n-1; i++) output_stream << _concentration_array_B[i] << delimiter;
    // For the last element shift an endline instead of the delimiter
    output_stream << _concentration_array_B[_n-1] << std::endl;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::printGridPoints(
    std::ostream &output_stream,
    char delimiter
) {
    // First value is NaN since first column stores the time
    output_stream << NAN << delimiter;

    // Print the r-coordinate of the grid point followed by the delimiter for all
    // except the last grid point
    for (size_t i = 0; i < _n - 1; i++) output_stream << getRadialCoordinate(i) << delimiter;
    // For the last grid point print the x-coordinate followed by endline
    output_stream << getRadialCoordinate(_n-1) << std::endl;
}
/******************************************************************************************************************/

/******************************************************************************************************************/
// Create classes for float, double and long double data types
template class CoreShellDiffusion<float>;
template class CoreShellDiffusion<double>;
template class CoreShellDiffusion<long double>;