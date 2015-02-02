#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_odeiv.h>
#include "params.h"

double force (double x1, double x2, double v0, double u0);
double classical_pot (double x, double v0), in (double x1, double x2,
						double u0);

/*Force on particle x1*/
double
force (double x1, double x2, double v0, double u0)
{
  double wellforce, interac;
  interac = (x1 - x2) * (x1 - x2);
  interac = interac / (2.0 * SIGMA * SIGMA);
  interac = exp (-interac);
  interac = interac * (x1 - x2);
  interac = interac * u0 / (SIGMA * SIGMA * SIGMA * sqrt (2.0 * M_PI));
  wellforce = 4.0 * v0 * (x1 - x1 * x1 * x1);
  return (wellforce + interac);
}

/*This defines the function (RHS) of the dynamical system*/
int
func (double t, const double y[], double f[], void *params)
{
  double u0 = U0, v0 = V0;
  f[0] = 2.0 * y[2];
  f[1] = 2.0 * y[3];
  f[2] = force (y[0], y[1], v0, u0);
  f[3] = force (y[1], y[0], v0, u0);
  return GSL_SUCCESS;
}

/* This defines the Jacobian Matrix of the dynamical system     */
/*(Only used in Bulirsch Stoer Methods, not Runge-Kutta)       */
/* Since I will be using Runge-Kutta Prince-Dormand Method always, this is a dummy function */
int
jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  return GSL_SUCCESS;
}

/*This function actually carries out the NONADAPTIVE integration of the ODE for a set of ICs */
int
integrate (double *y, double initial, double final, double epsilon,
	   double *out)
{
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
  int count = 0;
  gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 4);
  gsl_odeiv_system sys = { func, jac, 4, &epsilon };
  double nadt = NADT, t = initial;
  double y_err[4], dydt_in[4], dydt_out[4];
  int status;
  GSL_ODEIV_FN_EVAL (&sys, t, y, dydt_in);
  while (t < final)
    {
      status =
	gsl_odeiv_step_apply (s, t, nadt, y, y_err, dydt_in, dydt_out, &sys);
      dydt_in[0] = dydt_out[0];
      dydt_in[1] = dydt_out[1];
      dydt_in[2] = dydt_out[2];
      dydt_in[3] = dydt_out[3];
      t += nadt;
      /*Poincare Condition(s) HERE */
      if ((fabs (y[1] - 1.0) <= STROBERR) && (y[3] >= 0.0))
	{
	  out[count] = y[0];
	  out[count + 1] = y[2];
	  count = count + 2;
	}
      if (status != GSL_SUCCESS)
	{
	  printf ("\n Errors occurred in the integration process!Bailing...");
	  fflush (stdout);
	  break;
	}
    }
  return (count);
}
