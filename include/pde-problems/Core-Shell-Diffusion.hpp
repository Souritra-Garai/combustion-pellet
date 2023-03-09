#ifndef __CORE_SHELL_DIFFUSION__
#define __CORE_SHELL_DIFFUSION__

#include <ostream>

#include "math/Data-Type.hpp"

#include "thermo-physical-properties/Core-Shell-Particle.hpp"
#include "lusolver/LU-Solver.hpp"

class CoreShellDiffusion : public CoreShellParticle
{
    private :

		real_t * _concentration_array_A;
		real_t * _concentration_array_B;
		
		LUSolver _solver_A;
		LUSolver _solver_B;

		static real_t *radial_coordinate_sqr;
		static real_t *radial_coordinate_sqr_ratio;

		static real_t getRadialCoordinate(size_t index);

		void updateMassFractions();
		
		real_t getRxnConcA  (size_t index) const;
		real_t getRxnConcB  (size_t index) const;
		real_t getRxnConcAB (size_t index) const;

    public :

		static const real_t delta_t;
		
		static const size_t n;
		static const real_t delta_r;

		CoreShellDiffusion();
		~CoreShellDiffusion();

		static void setUpRadiusArray();
		static void deallocateRadiusArray();
		
		void initializeParticle();
		
		inline
		void setUpEquations(real_t temperature) { setUpEquations(temperature, *this);}
		void setUpEquations(real_t temperature, CoreShellDiffusion &diffusion_problem);
		void solveEquations();

		real_t getAtomMassA() const;
		real_t getAtomMassB() const;

		void copyFrom(CoreShellDiffusion & diffusion_problem);
		void copyTo(CoreShellDiffusion & diffusion_problem);

		void printConcentrationProfileA(std::ostream &output_stream, char delimiter = '\t', real_t curr_time = 0) const;
		void printConcentrationProfileB(std::ostream &output_stream, char delimiter = '\t', real_t curr_time = 0) const;

		static void printGridPoints(std::ostream &output_stream, char delimiter = '\t');
};

#endif