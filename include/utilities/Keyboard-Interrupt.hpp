/**
 * @file Keyboard_Interrupt.hpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief Contains functions to handle keyboard interrupts
 * @version 0.1
 * @date 2021-07-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#ifndef __KEYBOARD_INTERRUPT__
#define __KEYBOARD_INTERRUPT__

#include <exception>

class InterruptException : public std::exception
{
    public:
        
        InterruptException(int s) : S(s) {}
        
        int S;
};

void setUpKeyboardInterrupt();

#endif