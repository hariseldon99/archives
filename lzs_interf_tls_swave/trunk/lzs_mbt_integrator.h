// C++ Interface: integrator
//
// Description:
//
//
// Author: Analabha Roy <daneel@bose.res.in>, (C) 2010
//
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>

/***********************************************ERROR TOLERANCES***********************************/
/*
The step-size adjustment procedure for this method begins by computing the desired error level D_i for each component,

D_i = ABSERROR + RELERROR * (YERROR |y_i| + YPRIMEERROR*DT |y'_i|)

and comparing it with the observed error E_i = |yerr_i|. If the observed error E exceeds the desired error level D by more than 10% for any component then the
method reduces the step-size by an appropriate factor
*/
#define DT 1e-9			/*better value 1E-6 but slower run */
#define ABSERROR 1e-11		/*Low Accuracy in order to run in a reasonable time, better value 1e-8 or 0 */
#define RELERROR 1e-11		/*Low Accuracy in order to run in a reasonable time, better value 1e-8or 1e-5 */
#define YERROR 0.5
#define YPRIMEERROR 0.2

#define NORMTOL 1e-5
#define DT_NONADAPTIVE 1e-6

void prn_array (double *a, long size);
void begin_out ();
void begin_tout ();

double norm (const double phi[]);
void normalize (double phi[]);
double epsilon (double ampmax, double omega, double t, double offset);
double mu (double mu0, double muamp, double omega, double t, double offset);
double dmudt (double mu0, double muamp, double omega, double t,
	      double offset);
int adiabatic_trans (double phi_adb[], double phi_dia[], double t,
		     void *params);
int adiabatic_trans_mbt (double phi_adb[], double phi_dia[], double t,
			 void *params);
void data_out (double t, double *phi, double phi_norm);
