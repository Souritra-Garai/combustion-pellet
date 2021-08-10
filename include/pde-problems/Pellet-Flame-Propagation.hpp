/**
 * @file Pellet-Flame-Propagation.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Definition of class representing energy transport problem in
 * a packed pellet of core shell particles
 * @version 0.1
 * @date 2021-08-02
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __PELLET_FLAME_PROPAGATION__
#define __PELLET_FLAME_PROPAGATION__

#include <ostream>

#include "thermo-physical-properties/Arrhenius_Diffusivity_Model.hpp"
#include "thermo-physical-properties/Packed-Pellet.hpp"
#include "pde-problems/Core-Shell-Diffusion.hpp"

#include "qrsolver/QR_Solver.hpp"

template<typename real_t>
class LinearExpression {
    
    public :

        real_t constant;
        real_t coefficient;
};

template<typename real_t>
class PelletFlamePropagation : public PackedPellet<real_t>
{
    private:
    
        static size_t _n;
        static real_t _delta_x;

        static real_t _delta_t;
        
        static real_t _ignition_temperature;
        static real_t _ignition_length;

        real_t * _temperature_array;
        
        CoreShellDiffusion<real_t> * _particles_array;

        CoreShellDiffusion<real_t> * _prev_particles_array;
        CoreShellDiffusion<real_t> * _curr_particles_array;

        QRSolver<real_t> _solver;

        ArrheniusDiffusivityModel<real_t> _diffusivity_model;

        real_t getParticleEnthalpyTemperatureDerivative(size_t index);
        real_t getParticleEnthalpyTimeDerivative(size_t index);

        real_t getxCoordinate(size_t index) { return (real_t) index * this->_length / (real_t) (_n - 1); }

        bool inReactionZone(size_t index)
        {
            return
                _temperature_array[index] >= _ignition_temperature &&
                !_particles_array[index].isCombustionComplete();
        }

        LinearExpression<real_t> calcTransientTerm(size_t index);
        LinearExpression<real_t> calcHeatLossTerm(size_t index);

        void setUpBoundaryConditionX0();
        void setUpBoundaryConditionXN();

        void updateParticlesState();

    public:

        static void setGridSize(size_t N);
        static void setTimeStep(real_t delta_t);

        static void setIgnitionParameters(
            real_t ignition_temperature,
            real_t ignition_length
        );

        void setDiffusivityModel(ArrheniusDiffusivityModel<real_t> diffusivity_model);

        PelletFlamePropagation(real_t particle_volume_fraction);

        ~PelletFlamePropagation();

        void setUpEquations();
        void solveEquations();

        bool isCombustionComplete();

        void printTemperatureProfile(std::ostream &output_stream, char delimiter = '\t');
};

#endif