/**
 * @file Arrhenius_Diffusivity_Model.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header defines a class for Arrhenius type model
 * for diffusivity
 * @version 0.1
 * @date 2021-07-26
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __ARRHENIUS_DIFFUSIVITY_MODEL__
#define __ARRHENIUS_DIFFUSIVITY_MODEL__

// Required for the exponential function
#include <math.h>

/**
 * @brief Class representing Arrhenius rate law type model for
 * diffusivity coefficient in Fick's Diffusion Law
 * 
 * @tparam real_t real_t float, double or long double data types
 * to represent real numbers
 */
template<typename real_t>
class ArrheniusDiffusivityModel
{
    private:
        
        /**
         * @brief Pre-exponential factor in Arrhenius model
         * for diffusivity dependence on temperature
         */
        real_t _pre_exponential_factor;
        /**
         * @brief Activation energy in Arrhenius model
         * for diffusivity dependence on temperature
         */
        real_t _activation_energy;

    public:
        
        /**
         * @brief Construct a new Arrhenius Diffusivity Model object
         * 
         * @param pre_exponential_factor Pre-exponential factor in Arrhenius model
         * @param activation_energy Activation energy in Arrhenius model
         */
        ArrheniusDiffusivityModel(
            real_t pre_exponential_factor,
            real_t activation_energy
        ) {
            // Set the Arrhenius model paramemters
            setParameters(pre_exponential_factor, activation_energy);
        }

        /**
         * @brief Construct a new Arrhenius Diffusivity Model object
         */
        ArrheniusDiffusivityModel() {
            // Set the Arrhenius model parameters to zero
            setParameters(0, 0);
        }

        /**
         * @brief Set the parameters in Arrhenius model for diffusivity
         * @param pre_exponential_factor Pre-exponential factor in \f$ m^2 / s \f$
         * @param activation_energy Activation energy in \f$ J / mol \f$
         */
        void setParameters(
            real_t pre_exponential_factor,
            real_t activation_energy
        ) {
            _pre_exponential_factor = pre_exponential_factor;
            _activation_energy = activation_energy;
        }

        /**
         * @brief Get the Diffusivity at the specified temperature
         * 
         * @param temperature Temperature at which diffusivity is desired
         * @return real_t Diffusivity at the specified temperature in \f$ m^2 / s \f$
         */
        real_t getDiffusivity(real_t temperature) {
            return _pre_exponential_factor * exp(- _activation_energy / (8.314 * temperature));
        }
};

#endif