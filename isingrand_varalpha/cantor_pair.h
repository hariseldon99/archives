/*
// C++ Interface: using gsl polynomial functions to generate Cantor pairs
// Description:
// Use a Cantor Pairing function to create unique seeds: http://szudzik.com/ElegantPairing.pdf
// This is so that different mpi procs have different seeds and don't 
// generate the same pseudo random numbers
// Author: Analabha Roy <daneel@utexas.edu>, (C) 2013
//
//
*/
#include <gsl/gsl_math.h>
