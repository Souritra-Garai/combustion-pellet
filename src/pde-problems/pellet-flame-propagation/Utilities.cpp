#include "pde-problems/Pellet-Flame-Propagation.hpp"

#include "utilities/Read-Data.hpp"

const real_t PelletFlamePropagation::kappa = readScalarData<real_t>("data/PDE-solver-config", "kappa.txt");
const real_t PelletFlamePropagation::gamma = readScalarData<real_t>("data/PDE-solver-config", "gamma.txt");

const size_t PelletFlamePropagation::m		 = readScalarData<real_t>("data/PDE-solver-config", "number-of-grid-points-pellet.txt");
const real_t PelletFlamePropagation::delta_x = readScalarData<real_t>("data/pellet", "length.txt") / ((real_t) PelletFlamePropagation::m - 1.);

const real_t PelletFlamePropagation::delta_t = readScalarData<real_t>("data/PDE-solver-config", "delta_t.txt");

const real_t PelletFlamePropagation::delta_T = readScalarData<real_t>("data/PDE-solver-config", "delta_T.txt");

PelletFlamePropagation::PelletFlamePropagation(
	real_t particle_volume_fraction
) : PackedPellet(particle_volume_fraction),
	_solver(m)
{
	_temperature_array = new real_t[m];

	_thermal_conductivity = new real_t[m];

	_prev_enthalpy_particle = new real_t[m];

	_particles_array = new CoreShellDiffusion[m];

	_particles_array_const_temperature_evolution = new CoreShellDiffusion[m];
	_particles_array_raised_temperature_evolution = new CoreShellDiffusion[m];

	initializePellet();
}

PelletFlamePropagation::~PelletFlamePropagation()
{
	// Deallocate memory for temperature array
	delete [] _temperature_array;
	delete [] _thermal_conductivity;
	delete [] _prev_enthalpy_particle;

	// Deallocate memory for core shell diffusion problems
	delete [] _particles_array;
	delete [] _particles_array_const_temperature_evolution;
	delete [] _particles_array_raised_temperature_evolution;
}


inline real_t PelletFlamePropagation::getXCoordinate(size_t index)
{
	return (real_t) index * delta_x;
}

void PelletFlamePropagation::initializePellet(
	real_t initial_ignition_temperature,
	real_t initial_ignition_length_fraction
) {
	_time = 0;
	_temperature_array[0] = initial_ignition_temperature;

	real_t initial_ignition_length = PackedPellet::length * initial_ignition_length_fraction;

	#pragma omp parallel for

		for (size_t i = 1; i < m-1; i++)
		{
			if (getXCoordinate(i) < initial_ignition_length) _temperature_array[i] = initial_ignition_temperature;
			
			else _temperature_array[i] = PackedPellet::ambient_temperature;
			
			_particles_array[i].initializeParticle();
			
			_particles_array_raised_temperature_evolution[i].initializeParticle();
			_particles_array_const_temperature_evolution[i].initializeParticle();
		}
		
	_temperature_array[m-1] = PackedPellet::ambient_temperature;

	updateParticles();
}

void PelletFlamePropagation::printTemperatureProfile(
	std::ostream &output_stream,
	char delimiter
) {
	output_stream << _time << delimiter;

	for (size_t i = 0; i < m -1; i++) output_stream << _temperature_array[i] << delimiter;

	output_stream << _temperature_array[m-1] << '\n';
}

void PelletFlamePropagation::printGridPoints(
	std::ostream &output_stream,
	char delimiter
) {
	output_stream << NAN << delimiter;

	for (size_t i = 0; i < m - 1; i++) output_stream << getXCoordinate(i) << delimiter;

	output_stream << getXCoordinate(m-1) << '\n';
}

void PelletFlamePropagation::printDiffusionParticleGridPoints(
	std::ostream &output_stream,
	unsigned int particle_index,
	char delimiter
) {
	_particles_array[particle_index].printGridPoints(output_stream, delimiter);
}

void PelletFlamePropagation::printDiffusionParticleConcentationProfiles(
	std::ostream &output_stream_A,
	std::ostream &output_stream_B, 
	unsigned int particle_index,
	char delimiter
) {
	_particles_array[particle_index].printConcentrationProfileA(output_stream_A, delimiter, _time);
	_particles_array[particle_index].printConcentrationProfileB(output_stream_B, delimiter, _time);
}