#include "pde-problems/Core-Shell-Diffusion.hpp"

// Required for pow function
#include <math.h>

#include <string.h>

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
inline real_t CoreShellDiffusion<real_t>::getRadialCoordinate(size_t i)
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
	cudaMalloc(&_concentration_array_A, _n * sizeof(real_t));
	cudaMalloc(&_concentration_array_B, _n * sizeof(real_t));

    initializeParticle();
}

template<typename real_t>
CoreShellDiffusion<real_t>::~CoreShellDiffusion()
{
    // Deallocate the memory allocated for concentration profiles
    cudaFree(_concentration_array_A);
	cudaFree(_concentration_array_B);
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
