#include "pde-problems/Core-Shell-Diffusion.hpp"

#include "utilities/Read-Data.hpp"

#include <cmath>	// pow
#include <cstring>	// memcpy

const real_t CoreShellDiffusion::delta_t = readScalarData<real_t>("data/PDE-solver-config", "delta_t.txt");

const size_t CoreShellDiffusion::n		 = readScalarData<size_t>("data/PDE-solver-config", "number-of-grid-points-particle.txt");
const real_t CoreShellDiffusion::delta_r = readScalarData<real_t>("data/core-shell-particle", "overall-radius.txt") / ((real_t) CoreShellDiffusion::n - 1);

real_t * CoreShellDiffusion::radial_coordinate_sqr;
real_t * CoreShellDiffusion::radial_coordinate_sqr_ratio;

inline real_t CoreShellDiffusion::getRadialCoordinate(size_t i)
{
	return delta_r * (real_t) i;
}

void CoreShellDiffusion::setUpRadiusArray()
{
	radial_coordinate_sqr		= new real_t[n];
	radial_coordinate_sqr_ratio	= new real_t[n];

	for (size_t i = 0; i < n; i++)		radial_coordinate_sqr[i] = std::pow(getRadialCoordinate(i), 2);

	for (size_t i = 1; i < n-1; i++)	radial_coordinate_sqr_ratio[i] = radial_coordinate_sqr[i+1] / radial_coordinate_sqr[i];
}

void CoreShellDiffusion::deallocateRadiusArray()
{
	delete [] radial_coordinate_sqr;
	delete [] radial_coordinate_sqr_ratio;
}

CoreShellDiffusion::CoreShellDiffusion() : 
	CoreShellParticle(),
	_solver_A(n),
	_solver_B(n)
{
	// Allocate memory for concentration profiles
	_concentration_array_A = new real_t[n];
	_concentration_array_B = new real_t[n];

	initializeParticle();
}

CoreShellDiffusion::~CoreShellDiffusion()
{
	// Deallocate the memory allocated for concentration profiles
	delete [] _concentration_array_A;
	delete [] _concentration_array_B;
}

void CoreShellDiffusion::initializeParticle()
{
	for (size_t i = 0; i < n; i++)
	{
		if (getRadialCoordinate(i) < CoreShellParticle::core_radius)
		{
			_concentration_array_A[i] = CoreShellParticle::core_species.getMolarDensity(298.15);
			_concentration_array_B[i] = 0;
		}
		
		else
		{
			_concentration_array_A[i] = 0;
			_concentration_array_B[i] = CoreShellParticle::shell_species.getMolarDensity(298.15);
		}
	}

	updateMassFractions();
}

void CoreShellDiffusion::copyFrom(CoreShellDiffusion &diffusion_problem)
{
	std::memcpy(_concentration_array_A, diffusion_problem._concentration_array_A, n * sizeof(real_t));
	std::memcpy(_concentration_array_B, diffusion_problem._concentration_array_B, n * sizeof(real_t));

	CoreShellParticle::_mass_fraction_core_material    = diffusion_problem._mass_fraction_core_material;
	CoreShellParticle::_mass_fraction_shell_material   = diffusion_problem._mass_fraction_shell_material;
	CoreShellParticle::_mass_fraction_product_material = diffusion_problem._mass_fraction_product_material;
}

void CoreShellDiffusion::copyTo(CoreShellDiffusion &diffusion_problem)
{
	std::memcpy(diffusion_problem._concentration_array_A, _concentration_array_A, n * sizeof(real_t));
	std::memcpy(diffusion_problem._concentration_array_B, _concentration_array_B, n * sizeof(real_t));

	diffusion_problem._mass_fraction_core_material    = CoreShellParticle::_mass_fraction_core_material;
	diffusion_problem._mass_fraction_shell_material   = CoreShellParticle::_mass_fraction_shell_material;
	diffusion_problem._mass_fraction_product_material = CoreShellParticle::_mass_fraction_product_material;
}

void CoreShellDiffusion::printConcentrationProfileA(std::ostream &output_stream, char delimiter, real_t curr_time) const
{
	output_stream << curr_time << delimiter;

	for (size_t i = 0; i < n-1; i++) output_stream << _concentration_array_A[i] << delimiter;
	
	output_stream << _concentration_array_A[n-1] << '\n';
}

void CoreShellDiffusion::printConcentrationProfileB(std::ostream &output_stream, char delimiter, real_t curr_time) const
{
	output_stream << curr_time << delimiter;
	
	for (size_t i = 0; i < n-1; i++) output_stream << _concentration_array_B[i] << delimiter;
	
	output_stream << _concentration_array_B[n-1] << '\n';
}

void CoreShellDiffusion::printGridPoints(
	std::ostream &output_stream,
	char delimiter
) {
	output_stream << NAN << delimiter;
	
	for (size_t i = 0; i < n - 1; i++) output_stream << getRadialCoordinate(i) << delimiter;
	
	output_stream << getRadialCoordinate(n-1) << '\n';
}
