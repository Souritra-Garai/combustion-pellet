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

#include "utilities/Read-Data.hpp"

#include <iostream>

// Required for pow function
#include <math.h>

#include <string.h>

/******************************************************************************************************************/
// Instatiating static member variables of CoreShellDiffusion Class

// Time step interval duration
template<typename real_t> const real_t CoreShellDiffusion<real_t>::delta_t	= readScalarData<real_t>("data/PDE-solver-config", "delta_t.txt");
// Grid size
template<typename real_t> const size_t CoreShellDiffusion<real_t>::n		= readScalarData<size_t>("data/PDE-solver-config", "number-of-grid-points-particle.txt");
template<typename real_t> const real_t CoreShellDiffusion<real_t>::delta_r	= readScalarData<real_t>("data/core-shell-particle", "overall-radius.txt") / ((real_t) CoreShellDiffusion<real_t>::n - 1);

template<typename real_t> real_t * CoreShellDiffusion<real_t>::radial_coordinate_sqr;
template<typename real_t> real_t * CoreShellDiffusion<real_t>::radial_ratio;
/******************************************************************************************************************/

/******************************************************************************************************************/
// Defining static member functions of CoreShellDiffusion Class

// Setting grid size
template<typename real_t>
void CoreShellDiffusion<real_t>::setUpRadiusArray()
{
	// std::cout << delta_t << "\t" << delta_r << "\t" << n << std::endl;
	radial_coordinate_sqr = new real_t[n];
	radial_ratio = new real_t[n];

	for (size_t i = 0; i < n; i++)		radial_coordinate_sqr[i] = pow(getRadialCoordinate(i), 2);

	for (size_t i = 1; i < n-1; i++)	radial_ratio[i] = radial_coordinate_sqr[i+1] / radial_coordinate_sqr[i];

	// for (size_t i = 0; i < n; i++)	std::cout << i << "\t" << radial_ratio[i] << "\t" << radial_coordinate_sqr[i] << std::endl;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::deallocateRadiusArray()
{
	delete [] radial_coordinate_sqr;
	delete [] radial_ratio;
}

template<typename real_t>
inline real_t CoreShellDiffusion<real_t>::getRadialCoordinate(size_t i)
{
	static const real_t constant = CoreShellParticle<real_t>::_overall_radius / ((real_t) (n-1));
    // Radial distance from origin of the grid point \f$ i \f$
    // \f$ r_i = r_p \cdot \frac{i}{n-1}
    return constant * (real_t) i;
    // For loops will iterate from \f$ i = 0 \f$ to \f$ i = n-1 \f$
}

/******************************************************************************************************************/

/******************************************************************************************************************/
// Constructors and destructors

template<typename real_t>
CoreShellDiffusion<real_t>::CoreShellDiffusion() : // Mem initialization list
    CoreShellParticle<real_t>(),  // Call the constructor of the base class
    // Initialize solvers
    _solver_A(n),
    _solver_B(n)
{
    // Allocate memory for concentration profiles
    _concentration_array_A = new real_t[n];
    _concentration_array_B = new real_t[n];

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
    
    // Intialise the sums without the common multiplicative factors
    real_t Y_A  = 0.5 * getRxnConcA(n-1)  * radial_coordinate_sqr[n - 1];
    real_t Y_B  = 0.5 * getRxnConcB(n-1)  * radial_coordinate_sqr[n - 1];
    real_t Y_AB = 0.5 * getRxnConcAB(n-1) * radial_coordinate_sqr[n - 1];

	// Summing over all grid points except the first and the last
	for (size_t i = 1; i < n-1; i++)
	{
		// Add the specific terms for the summation
		Y_A  += getRxnConcA(i)  * radial_coordinate_sqr[i];
		Y_B  += getRxnConcB(i)  * radial_coordinate_sqr[i];
		Y_AB += getRxnConcAB(i) * radial_coordinate_sqr[i];
	}
    
    // Update the reaction mass fractions of the core-shell combustion particle

    // // Multiply the summations with the required factors and divide
    // // by the total mass of the particle calculated during initialisation
    // this->_mass_fraction_core_material    = 4.0 * M_PI * delta_r * Y_A  * this->_core_species->getMolarMass()    / this->_mass;
    // this->_mass_fraction_shell_material   = 4.0 * M_PI * delta_r * Y_B  * this->_shell_species->getMolarMass()   / this->_mass;
    // this->_mass_fraction_product_material = 4.0 * M_PI * delta_r * Y_AB * this->_product_species->getMolarMass() / this->_mass;

    // Multiply the summations with the required factors
    Y_A  *= CoreShellParticle<real_t>::_core_species.getMolarMass();
    Y_B  *= CoreShellParticle<real_t>::_shell_species.getMolarMass();
    Y_AB *= CoreShellParticle<real_t>::_product_species.getMolarMass();
    
    // Add the summations to get the total mass
    real_t sum = Y_A + Y_B + Y_AB;
    // Divide by the sum of the summations to get the mass fractions
    CoreShellParticle<real_t>::_mass_fraction_core_material    = Y_A  / sum;
    CoreShellParticle<real_t>::_mass_fraction_shell_material   = Y_B  / sum;
    CoreShellParticle<real_t>::_mass_fraction_product_material = Y_AB / sum;
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
	// Loop over the grid points to initialize the concentration of the
	// substances A and B
	for (size_t i = 0; i < n; i++)
	{
		// Only substance A is present in the core (\f$ r < r_C \f$)
		if (getRadialCoordinate(i) < CoreShellParticle<real_t>::_core_radius)
		{
			// Concentration of pure substance A is its molar density
			_concentration_array_A[i] = CoreShellParticle<real_t>::_core_species.getMolarDensity(298.15);
			_concentration_array_B[i] = 0;
		}
		// Only substance B is present in the shell (\f$ r > r_V \f$)
		else
		{
			// Concentration of pure substance B is its molar density
			_concentration_array_A[i] = 0;
			_concentration_array_B[i] = CoreShellParticle<real_t>::_shell_species.getMolarDensity(298.15);
		}
	}
}

template<typename real_t>
void CoreShellDiffusion<real_t>::setUpEquations(real_t T, CoreShellDiffusion<real_t> &diffusion_problem)
{
    // \f$ \Rightarrow -
    // \mathcal{D} \cdot \left( \frac{r_{i+1}}{\Delta r \cdot r_i} \right)^2 \cdot C_{k,i+1}^n +
    // \left\{ \frac{1}{\Delta t} + \mathcal{D} \cdot \left( \frac{r_{i+1}}{\Delta r \cdot r_i} \right)^2 + \frac{\mathcal{D}}{\left(\Delta r\right)^2} \right\} \cdot C_{k,i}^n \\ -
    // \frac{\mathcal{D}}{\left(\Delta r\right)^2} \cdot C_{k,i-1}^n
    // = \frac{C_{k,i}^{n-1}}{\Delta t} \f$

	static const real_t constant_r = 0.5 / pow(delta_r, 2);

    // We can calculate the coefficients in three parts to avoid repetitive calculations

    // Coefficient 1 = \f$ \frac{\mathcal{D\left(T\right)}}{\left(\Delta r\right)^2} \f$
    real_t coefficient1 = CoreShellParticle<real_t>::_diffusivity_model.getDiffusivity(T) * constant_r;
    // Coefficient 2 = \f$ \frac{1}{\Delta t} \f$
    static const real_t coefficient2 = 1 / delta_t;
    // Coefficient 3 = \f$ \mathcal{D} \cdot \left( \frac{r_{i+1}}{\Delta r \cdot r_i} \right)^2 \f$
    real_t coefficient3;

	real_t coeff_sum  = coefficient1 + coefficient2;
	real_t coeff_diff = coefficient1 - coefficient2;

    // Set zero flux boundary condition at grid point # 0
    // \f$ C_{k,1}^n - C_{k,0}^n = 0 \f$
    _solver_A.setEquationFirstRow(1, -1, 0);
    _solver_B.setEquationFirstRow(1, -1, 0);

	// Iteratively set up the discretized governing diffusion equation
	// for all grid points except the grid points # 0 and N, where
	// zero flux boundary conditions apply
	for (size_t i = 1; i < n - 1; i++)
	{
		// Calculating coefficient 3 for the current grid point
		coefficient3 = coefficient1 * radial_ratio[i];

		// Set up row for matrix equation representing 
		// diffusion of substance A
		_solver_A.setEquation(
			i, 
			- coefficient1, 
			coeff_sum + coefficient3, 
			- coefficient3, 
			coefficient1 * diffusion_problem._concentration_array_A[i-1]
			- (coeff_diff + coefficient3) * diffusion_problem._concentration_array_A[i]
			+ coefficient3 * diffusion_problem._concentration_array_A[i+1]
		);

		// Set up row for matrix equation representing 
		// diffusion of substance B
		_solver_B.setEquation(
			i, 
			- coefficient1, 
			coeff_sum + coefficient3, 
			- coefficient3,
			coefficient1 * diffusion_problem._concentration_array_B[i-1]
			- (coeff_diff + coefficient3) * diffusion_problem._concentration_array_B[i]
			+ coefficient3 * diffusion_problem._concentration_array_B[i+1]
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
    // #pragma omp parallel sections
    // {
        // Parallel section for substance A
        // #pragma omp section
            // Solve the matrix equation representing diffusion of substance A
            // and save the solution in array A
            _solver_A.getSolution(_concentration_array_A);

        // Parallel section for substance B
        // #pragma omp section
            // Solve the matrix equation representing diffusion of substance B
            // and save the solution in array B
            _solver_B.getSolution(_concentration_array_B);
    // }

	// for (size_t i = 0; i < n; i++) std::cout << i << "\t" << _concentration_array_A[i] << "\t" << _concentration_array_B[i] << "\n";

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
    real_t sum = 0.5 * _concentration_array_A[n-1] * radial_coordinate_sqr[n - 1];

	// Iterate over all the grid points but the last
	for (size_t i = 0; i < n-1; i++)
	{
		// Cumulatively add the term \f$ C_{A,i} r_i^2 \f$
		sum += _concentration_array_A[i] * radial_coordinate_sqr[i];
	}

    // Multiply the constant factors in the summation
    return 4.0 * M_PI * CoreShellParticle<real_t>::_core_species.getMolarMass() * delta_r * sum;
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
    real_t sum = 0.5 * _concentration_array_B[n-1] * radial_coordinate_sqr[n - 1];

	// Iterate over all the grid points but the last
	for (size_t i = 0; i < n-1; i++)
	{
		// Cumulatively add the term \f$ C_{A,i} r_i^2 \f$
		sum += _concentration_array_B[i] * radial_coordinate_sqr[i];
	}

    // Multiply the constant factors in the summation
    return 4.0 * M_PI * CoreShellParticle<real_t>::_shell_species.getMolarMass() * delta_r * sum;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::copyFrom(CoreShellDiffusion<real_t> &diffusion_problem)
{
	memcpy(_concentration_array_A, diffusion_problem._concentration_array_A, n * sizeof(real_t));
	memcpy(_concentration_array_B, diffusion_problem._concentration_array_B, n * sizeof(real_t));

    // Update reaction mass fractions of the particle
    this->_mass_fraction_core_material    = diffusion_problem._mass_fraction_core_material;
    this->_mass_fraction_shell_material   = diffusion_problem._mass_fraction_shell_material;
    this->_mass_fraction_product_material = diffusion_problem._mass_fraction_product_material;
}

template<typename real_t>
void CoreShellDiffusion<real_t>::copyTo(CoreShellDiffusion<real_t> &diffusion_problem)
{
	memcpy(diffusion_problem._concentration_array_A, _concentration_array_A, n * sizeof(real_t));
	memcpy(diffusion_problem._concentration_array_B, _concentration_array_B, n * sizeof(real_t));

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
    for (size_t i = 0; i < n-1; i++) output_stream << _concentration_array_A[i] << delimiter;
    // For the last element shift an endline instead of the delimiter
    output_stream << _concentration_array_A[n-1] << '\n';
}

template<typename real_t>
void CoreShellDiffusion<real_t>::printConcentrationProfileB(std::ostream &output_stream, char delimiter, real_t curr_time)
{
    // First print the current time, followed by the temperature profile
    output_stream << curr_time << delimiter;
    
    // Shift the array elements into the output stream followed by a delimter
    for (size_t i = 0; i < n-1; i++) output_stream << _concentration_array_B[i] << delimiter;
    // For the last element shift an endline instead of the delimiter
    output_stream << _concentration_array_B[n-1] << '\n';
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
    for (size_t i = 0; i < n - 1; i++) output_stream << getRadialCoordinate(i) << delimiter;
    // For the last grid point print the x-coordinate followed by endline
    output_stream << getRadialCoordinate(n-1) << '\n';
}
/******************************************************************************************************************/

/******************************************************************************************************************/
// Create classes for float, double and long double data types
template class CoreShellDiffusion<float>;
template class CoreShellDiffusion<double>;
template class CoreShellDiffusion<long double>;