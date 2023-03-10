#ifndef __PELLET_FLAME_PROPAGATION__
#define __PELLET_FLAME_PROPAGATION__

#include <ostream>

#include "math/Data-Type.hpp"
#include "math/Linear-Expression.hpp"

#include "thermo-physical-properties/Packed-Pellet.hpp"
#include "pde-problems/Core-Shell-Diffusion.hpp"
#include "lusolver/LU-Solver.hpp"

class PelletFlamePropagation : public PackedPellet
{
	private:

		real_t _time;
		
		real_t * _temperature_array;

		real_t * _thermal_conductivity;

		real_t * _prev_enthalpy_particle;
		
		CoreShellDiffusion * _particles_array;

		CoreShellDiffusion * _particles_array_const_temperature_evolution;
		CoreShellDiffusion * _particles_array_raised_temperature_evolution;

		LUSolver _solver;

		static real_t getXCoordinate(size_t index);

		bool inReactionZone(size_t index);

		void evolveParticleForEnthalpyDerivative(size_t index);

		LinearExpression calcTransientTerm(size_t index);
		LinearExpression calcHeatLossTerm(size_t index);

		void setUpBoundaryConditionX0();
		void setUpBoundaryConditionXN();

		void updateParticles();

		real_t getInterstitialGasTransientTermCoefficient(size_t index);

	public:
	
		static const real_t kappa;

		static const real_t gamma;
	
		static const size_t m;
		static const real_t delta_x;

		static const real_t delta_t;

		static const real_t delta_T;

		PelletFlamePropagation(real_t particle_volume_fraction);
		~PelletFlamePropagation();

		void initializePellet(
			real_t initial_ignition_temperature = 1500.,
			real_t initial_ignition_length_fraction = 0.1
		);

		void setUpEquations();
		void solveEquations();

		bool isCombustionComplete();

		void printTemperatureProfile(std::ostream &output_stream, char delimiter = '\t');
		
		void printGridPoints(std::ostream &output_stream, char delimiter = '\t');

		void printDiffusionParticleGridPoints(std::ostream &output_stream, unsigned int particle_index, char delimiter = '\t');

		void printDiffusionParticleConcentationProfiles(
			std::ostream &output_stream_A,
			std::ostream &output_stream_B, 
			unsigned int particle_index,
			char delimiter = ','
		);
};

#endif