#include "pde-problems/Core-Shell-Diffusion.hpp"

#include <cmath>

inline real_t CoreShellDiffusion::getRxnConcA(size_t i) const
{
	return std::max(_concentration_array_A[i]  - _concentration_array_B[i], (real_t) 0.0);
}

inline real_t CoreShellDiffusion::getRxnConcB(size_t i) const
{
	return std::max(_concentration_array_B[i] - _concentration_array_A[i], (real_t) 0.0);
}

inline real_t CoreShellDiffusion::getRxnConcAB(size_t i) const
{
	return std::min(_concentration_array_A[i], _concentration_array_B[i]);
}

void CoreShellDiffusion::updateMassFractions()
{
	real_t m_A  = 0.5 * getRxnConcA(n-1)  * radial_coordinate_sqr[n - 1];
	real_t m_B  = 0.5 * getRxnConcB(n-1)  * radial_coordinate_sqr[n - 1];
	real_t m_AB = 0.5 * getRxnConcAB(n-1) * radial_coordinate_sqr[n - 1];

	for (size_t i = 1; i < n-1; i++)
	{
		m_A  += getRxnConcA(i)  * radial_coordinate_sqr[i];
		m_B  += getRxnConcB(i)  * radial_coordinate_sqr[i];
		m_AB += getRxnConcAB(i) * radial_coordinate_sqr[i];
	}
    
	m_A  *= CoreShellParticle::core_species.molar_mass;
	m_B  *= CoreShellParticle::shell_species.molar_mass;
	m_AB *= CoreShellParticle::product_species.molar_mass;
    
    real_t sum = m_A + m_B + m_AB;
	
    CoreShellParticle::_mass_fraction_core_material    = m_A  / sum;
    CoreShellParticle::_mass_fraction_shell_material   = m_B  / sum;
    CoreShellParticle::_mass_fraction_product_material = m_AB / sum;
}

void CoreShellDiffusion::setUpEquations(real_t T, CoreShellDiffusion &diffusion_problem)
{
	static const real_t one_over_delta_t = 1 / delta_t;
	static const real_t half_over_delta_r_sqr = 0.5 / std::pow(delta_r, 2);

	real_t coefficient1 = CoreShellParticle::diffusivity_model.getDiffusivity(T) * half_over_delta_r_sqr;

	real_t coefficient2 = coefficient1 + one_over_delta_t;
	real_t coefficient3 = coefficient1 - one_over_delta_t;

	real_t coefficient4;

    // Set zero flux boundary condition at grid point # 0
    _solver_A.setEquationFirstRow(1, -1, 0);
    _solver_B.setEquationFirstRow(1, -1, 0);

	for (size_t i = 1; i < n - 1; i++)
	{
		coefficient4 = coefficient1 * radial_ratio[i];

		_solver_A.setEquationSerially(
			  i, 
			- coefficient1, 
			  coefficient2 + coefficient4, 
			- coefficient4, 
			  coefficient1 * diffusion_problem._concentration_array_A[i-1]
			- (coefficient3 + coefficient4) * diffusion_problem._concentration_array_A[i]
			+ coefficient4 * diffusion_problem._concentration_array_A[i+1]
		);

		_solver_B.setEquationSerially(
			i, 
			- coefficient1,
			  coefficient2 + coefficient4,
			- coefficient4,
			  coefficient1 * diffusion_problem._concentration_array_B[i-1]
			- (coefficient3 + coefficient4) * diffusion_problem._concentration_array_B[i]
			+ coefficient4 * diffusion_problem._concentration_array_B[i+1]
		);
	}

    // Set zero flux boundary condition at grid point # N
    _solver_A.setEquationLastRowSerially(-1, 1, 0);
    _solver_B.setEquationLastRowSerially(-1, 1, 0);
}

void CoreShellDiffusion::solveEquations()
{
	_solver_A.getSolutionSerially(_concentration_array_A);
	_solver_B.getSolutionSerially(_concentration_array_B);

	updateMassFractions();
}

real_t CoreShellDiffusion::getAtomMassA() const
{
	real_t sum = 0.5 * _concentration_array_A[n-1] * radial_coordinate_sqr[n - 1];

	for (size_t i = 0; i < n-1; i++)
	{
		sum += _concentration_array_A[i] * radial_coordinate_sqr[i];
	}

    return 4.0 * M_PI * CoreShellParticle::core_species.molar_mass * delta_r * sum;
}

real_t CoreShellDiffusion::getAtomMassB() const
{
	real_t sum = 0.5 * _concentration_array_B[n-1] * radial_coordinate_sqr[n - 1];
	
	for (size_t i = 0; i < n-1; i++)
	{
		sum += _concentration_array_B[i] * radial_coordinate_sqr[i];
	}
	
	return 4.0 * M_PI * CoreShellParticle::shell_species.molar_mass * delta_r * sum;
}
