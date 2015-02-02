/*
// C++ Interface: integrator
//
// Description:
//
//
// Author: Analabha Roy <daneel@physics.utexas.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
*/
/***********************************************ERROR TOLERANCES***********************************/
/*
The step-size adjustment procedure for this method begins by computing the desired error level D_i for each component,

D_i = ABSERROR + RELERROR * (YERROR |y_i| + YPRIMEERROR*DT |y'_i|)

and comparing it with the observed error E_i = |yerr_i|. If the observed error E exceeds the desired error level D by more than 10% for any component then the
method reduces the step-size by an appropriate factor
*/
#define DT 1e-6			/*better value 1E-6 but slower run */
#define ABSERROR 1e-9		/*Low Accuracy in order to run in a reasonable time, better value 1e-8 or 0 */
#define RELERROR 0.0		/*Low Accuracy in order to run in a reasonable time, better value 1e-8or 1e-5 */
#define YERROR 0.5
#define YPRIMEERROR 0.2
double hamilt (void *param, int nt, int mt, double t);
double hamilt_deriv (void *param, int nt, int mt, double t);
