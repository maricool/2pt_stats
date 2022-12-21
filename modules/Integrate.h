#ifndef INTEGRATE_H
#define INTEGRATE_H

// Guass-Legendre integration 
using namespace std;

#include "globaldef.h"
#include "function_cosebis.h"
#include <gsl/gsl_integration.h>

number gaussianIntegrate_gsl(function_cosebis&,number,number,int=20);

#endif