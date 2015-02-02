/*
// C++ Interface: integrator
//
// Description:
//
//
// Author: Analabha Roy <daneel@utexas.edu>, (C) 2013
//
//
*/
/***********************************************ERROR TOLERANCES***********************************/
/*
The step-size adjustment procedure for this method begins by computing the desired error level D_i for each component,

D_i = ABSERROR + RELERROR * |y_i|  

and comparing it with the observed error E_i = |yerr_i|. If the observed error E exceeds the desired error level D by more than 10% for any component then the
method reduces the step-size by an appropriate factor
*/

#define STEP_TYPE gsl_odeiv2_step_rk8pd	/*Choice of integrator algorithm. For details, see GSL ref manual */
#define DT 1e-7			/*better value 1E-6 but slower run */
#define ABSERROR 1e-9		/*Low Accuracy in order to run in a reasonable time, better value 1e-8 or 0 */
#define RELERROR 1e-8		/*Low Accuracy in order to run in a reasonable time, better value 1e-8or 1e-5 */

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
