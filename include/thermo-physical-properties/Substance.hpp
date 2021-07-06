/**
 * @file Substance.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This header defines a class to represent pure substances
 * with their thermodynamic properties like density, heat capacity etc.
 * @version 0.1
 * @date 2021-07-06
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __SUBSTANCE__
#define __SUBSTANCE__

// Required for << operator for printing to file / screen
#include <ostream>

/**
 * @brief Class to represent a pure substance
 * 
 * @tparam real_t 
 */
template <typename real_t>
class Substance
{
    public :

        /**
         * @brief Construct a new Substance
         * 
         * @param density Density of the substance in \f$ kg / m^3 \f$
         * @param heat_capacity Heat capacity of the substance in \f$ J / kg-K \f$
         * @param molecular_weight Molecular weight of the substance in \f$ kg / mol \f$
         * @param heat_conductivity Heat conductivity of the substance in \f$ W / m \f$
         */
        Substance(
            real_t density,
            real_t heat_capacity,
            real_t molecular_weight,
            real_t heat_conductivity
        ) : _density(density),
            _heat_capacity(heat_capacity),
            _molecular_weight(molecular_weight),
            _heat_conductivity(heat_conductivity)
        { ;}
        
        /**
         * @brief Get the Density
         * 
         * @return real_t Density of the substance
         */
        inline real_t getDensity()          { return _density; }

        /**
         * @brief Get the Heat Capacity
         * 
         * @return real_t Heat Capacity of the substance
         */
        inline real_t getHeatCapacity()     { return _heat_capacity; }

        /**
         * @brief Get the Molecular Weight
         * 
         * @return real_t Molecular Weight of the substance
         */
        inline real_t getMolecularWeight()  { return _molecular_weight; }

        /**
         * @brief Get the Heat Conductivity
         * 
         * @return real_t Heat Conductivity of the substance
         */
        inline real_t getHeatConductivity() { return _heat_conductivity; }

        /**
         * @brief 
         * 
         * @param output_stream 
         */
        void printProperties(std::ostream &output_stream)
        {
            output_stream << "Density\t\t\t:\t" << getDensity() << "\tkg/m3" << std::endl;
            output_stream << "Heat Capacity\t\t:\t" << getHeatCapacity() << "\tJ/kg-K" << std::endl;
            output_stream << "Molecular Weight\t:\t" << getMolecularWeight() << "\tkg/mol" << std::endl;
            output_stream << "Heat Conductivity\t:\t" << getHeatConductivity() << "\tW/m" << std::endl;
        }

    private :

        /**
         * @brief Density of the substance
         */
        const real_t _density;

        /**
         * @brief Heat capacity of the substance
         */
        const real_t _heat_capacity;

        /**
         * @brief Molecular Weight of the substance
         */
        const real_t _molecular_weight;

        /**
         * @brief Heat Conductivity of the substance
         */
        const real_t _heat_conductivity;
};

#endif