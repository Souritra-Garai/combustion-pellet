/**
 * @file Keyboard_Interrupt.cpp
 * @author Souritra Garai (souritra.garai@iitgn.ac.in)
 * @brief This file defines the functions used for keyboard interrupt
 * @version 0.1
 * @date 2021-07-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "utilities/Keyboard_Interrupt.hpp"

#include <signal.h>

void throwInterruptException(int signal_code)
{
    throw InterruptException(signal_code);
}

struct sigaction sig_int_handler;

void setUpKeyboardInterrupt()
{
    sig_int_handler.sa_handler = throwInterruptException;
    sigemptyset(&sig_int_handler.sa_mask);
    sig_int_handler.sa_flags = 0;
    sigaction(SIGINT, &sig_int_handler, NULL);
}