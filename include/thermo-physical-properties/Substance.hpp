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

#define T_REF 298.15 // Reference temperature for enthalpy in K

// Required for << operator for printing to file / screen
#include <ostream>

/**
 * @brief Class to represent a pure substance
 * 
 * @tparam real_t float, double or long double data types
 * to represent real numbers
 */
template <typename real_t>
class Substance
{
    public :

        /**
         * @brief uct a new Substance
         * 
         * @param density Density of the substance in \f$ kg / m^3 \f$
         * @param heat_capacity Heat capacity of the substance in \f$ J / kg-K \f$
         * @param molecular_weight Molecular weight of the substance in \f$ kg / mol \f$
         * @param heat_conductivity Heat conductivity of the substance in \f$ W / m-K \f$
         * @param standard_enthalpy_of_formation Standard enthalpy of formation of the substance at 298 K in \f$ J / kg \f$
         */
        Substance(
            real_t density = 0,
            real_t heat_capacity = 0,
            real_t molecular_weight = 0,
            real_t heat_conductivity = 0,
            real_t standard_enthalpy_of_formation = 0
        ) : _density(density),
            _heat_capacity(heat_capacity),
            _molecular_weight(molecular_weight),
            _heat_conductivity(heat_conductivity),
            _std_enthalpy_formation(standard_enthalpy_of_formation)
        { ;}
        
        /**
         * @brief Get the Density of the substance
         * 
         * @return real_t Density in \f$ kg / m^3 \f$
         */
        inline real_t getDensity()  { return _density; }

        /**
         * @brief Get the Heat Capacity of the substance
         * 
         * @return real_t Heat Capacity in \f$ J / kg-K \f$
         */
        inline real_t getHeatCapacity() { return _heat_capacity; }

        /**
         * @brief Get the Molecular Weight of the weight
         * 
         * @return real_t Molecular Weight in \f$ kg / mol. \f$
         */
        inline real_t getMolecularWeight()  { return _molecular_weight; }

        /**
         * @brief Get the Molar Density of the substance
         * 
         * @return real_t Density in \f$ mol. / m^3 \f$
         */
        inline real_t getMolarDensity() { return _density / _molecular_weight; }

        /**
         * @brief Get the Heat Conductivity of the substance
         * 
         * @return real_t Heat Conductivity in \f$ W / m-K \f$
         */
        inline real_t getHeatConductivity() { return _heat_conductivity; }

        /**
         * @brief Get the Standard Enthalpy Of Formation of the substance
         * 
         * @return real_t Standard Enthalpy of Formation of in \f$ J / kg \f$
         */
        inline real_t getStandardEnthalpyOfFormation() { return _std_enthalpy_formation; }

        /**
         * @brief Get the Enthalpy of the substance with respect to the T_REF
         * 
         * @param Temperature Temperature at which enthalpy of the object needs to be evaluated
         * @return real_t Enthalpy of the substance in J / kg
         */
        inline real_t getEnthalpy(real_t Temperature)
        {
            return getStandardEnthalpyOfFormation() + getHeatCapacity() * (Temperature - T_REF);
        }

        /**
         * @brief Print the properties of the substance to the given stream
         * 
         * @param output_stream Stream to which the properties are printed
         */
        void printProperties(std::ostream &output_stream)
        {
            output_stream << "Density\t\t\t\t:\t" << getDensity() << "\tkg/m3" << std::endl;
            output_stream << "Heat Capacity\t\t\t:\t" << getHeatCapacity() << "\tJ/kg-K" << std::endl;
            output_stream << "Molecular Weight\t\t:\t" << getMolecularWeight() << "\tkg/mol" << std::endl;
            output_stream << "Heat Conductivity\t\t:\t" << getHeatConductivity() << "\tW/m-K" << std::endl;
            output_stream << "Standard Enthalpy of Formation\t:\t" << getStandardEnthalpyOfFormation() << "\tJ/kg" << std::endl;
        }

    private :

        /**
         * @brief Density of the substance
         */
        real_t _density;

        /**
         * @brief Heat capacity of the substance
         */
        real_t _heat_capacity;

        /**
         * @brief Molecular Weight of the substance
         */
        real_t _molecular_weight;

        /**
         * @brief Heat Conductivity of the substance
         */
        real_t _heat_conductivity;

        /**
         * @brief Standard Enthalpy of Formation of the substance
         */
        real_t _std_enthalpy_formation;
};

#endif