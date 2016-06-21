/*
// C++ Interface: magnetization and correlations
//
// Description:
//
//
// Author: Analabha Roy <daneel@utexas.edu>, (C) 2013
//
//
*/


#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#define TOL 1E-8

gsl_complex gsl_determinant_complex (gsl_matrix_complex * A, int inPlace);
