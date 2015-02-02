#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_spline.h>

#include <glib.h>

#define INT_TOLERANCE 1E-2
#define DOUBLE_TOLERANCE 0.5

void prn_array (double *a, long size);
void begin_out ();
void begin_tout ();

double norm (const double phi[]);
void normalize (double phi[]);
double epsilon (double ampmax, double omega, double t, double offset);
double mu (double mu0, double muamp, double omega, double t, double offset);
double dmudt (double mu0, double muamp, double omega, double t,
	      double offset);
void data_out (double t, double *phi, double phi_norm);
