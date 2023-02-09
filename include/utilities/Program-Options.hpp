#ifndef __PROGRAM_OPTIONS__
#define __PROGRAM_OPTIONS__

#include <iostream>

#include "utilities/ezOptionParser.hpp"

void setHelpOption(ez::ezOptionParser &opt);

void displayHelpOption(ez::ezOptionParser &opt);

void setTemperatureOption(ez::ezOptionParser &opt);

double getTemperatureOption(ez::ezOptionParser &opt, double default_value);

void setPhiOption(ez::ezOptionParser &opt);

double getPhiOption(ez::ezOptionParser &opt, double default_value);

void setIgnitionTemperatureOption(ez::ezOptionParser &opt);

double getIgnitionTemperatureOption(ez::ezOptionParser &opt, double default_value);

void setIgnitionLengthOption(ez::ezOptionParser &opt);

double getIgnitionLengthOption(ez::ezOptionParser &opt, double default_value);

#endif