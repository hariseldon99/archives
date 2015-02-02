#include "params_ics.h"
#include "lzs_tls_integrator.h"


//Dynamics:
//H|psi> = i delt |psi>
//where:
//|psi> = {uk, vk}
//
//and 
//            |epsilon     Delta| 
//H = -(1/2)  |                 |
//            |Delta*   -epsilon|


int
func_tls (double t, const double phi[], double f[], void *params)
{
  paramspace_pt *p = (paramspace_pt *) params;
  double delta_real = p->delta_real;
  double delta_imag = p->delta_imag;
  double omega = p->omega;
  double offset = p->offset;
  double ampmax = p->ampmax;

  double eps = epsilon (ampmax, omega, t, offset);

  gsl_complex delta, uk, vk, rhs_up, rhs_dn;
  gsl_complex i2;

  GSL_SET_COMPLEX (&i2, 0, 0.5);

  GSL_SET_COMPLEX (&delta, delta_real, delta_imag);

  GSL_SET_COMPLEX (&uk, phi[0], phi[2]);
  GSL_SET_COMPLEX (&vk, phi[1], phi[3]);

  gsl_complex rhs_up1 = gsl_complex_mul_real (uk, eps);
  gsl_complex rhs_up2 = gsl_complex_mul (delta, vk);

  rhs_up = gsl_complex_add (rhs_up1, rhs_up2);
  rhs_up = gsl_complex_mul (i2, rhs_up);

  gsl_complex rhs_dn1 = gsl_complex_mul (gsl_complex_conjugate (delta), uk);
  gsl_complex rhs_dn2 = gsl_complex_mul_real (vk, eps);

  rhs_dn = gsl_complex_sub (rhs_dn1, rhs_dn2);
  rhs_dn = gsl_complex_mul (i2, rhs_dn);

  f[0] = GSL_REAL (rhs_up);
  f[1] = GSL_REAL (rhs_dn);
  f[2] = GSL_IMAG (rhs_up);
  f[3] = GSL_IMAG (rhs_dn);


  return GSL_SUCCESS;
}


int
jac_tls (double t, const double phi[], double *dfdphi, double dfdt[],
	 void *params)
{
  //Leaving blank for now. Not using Burlisch Stoer so no Jacobian needed
  return GSL_SUCCESS;
}


/*****This function actually runs the full integration of a particulat set of IC's from 0-t**********/
int
integrate (double *phi, void *initconds, void *params, GArray * data,
	   GArray * data_probs, GArray * data_probs_prime)
{
  int status = 0, i;
  double t_init, t_extra, t_periods, t_final;
  double dt, t, omega, phi_adiabatic[4], phi_norm;

  /*This sets up the gsl ODE system structure */
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rkck;
  /*Working algorithms: rk4imp(most accurate, slowest), rkf45(less accurate,faster),rkck(less accurate still, faster still) */
  gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 4);
  gsl_odeiv_control *c =
    gsl_odeiv_control_standard_new (ABSERROR, RELERROR, YERROR, YPRIMEERROR);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc (4);
  gsl_odeiv_system tls_sys = { func_tls, jac_tls, 4, params };
  dt = DT;

  //extract initial conditions
  paramspace_pt *p = (paramspace_pt *) params;
  omega = p->omega;

  initconds_ic *ic = (initconds_ic *) initconds;
  char verbose = ic->verbosity;

  for (i = 0; i < 4; i++)
    phi[i] = ic->phi_in[i];

  t_init = ic->t_init;
  t_periods = ic->t_periods;
  t_extra = ic->t_extra;

  t_final = t_init + t_periods + t_extra;

  t = t_init;
  if (verbose == 'y')
    begin_tout ();
  while (t < t_final)
    {
      status =
	gsl_odeiv_evolve_apply (e, c, s, &tls_sys, &t, t_final, &dt, phi);
      if (status != GSL_SUCCESS)
	{
	  status = GSL_FAILURE;
	  break;
	}

      phi_norm = norm (phi);
      //Output t, phi and norm
      g_array_append_val (data, t);
      int ierr = adiabatic_trans (phi_adiabatic, phi, t, params);
      status = ierr;
      g_array_append_vals (data, phi_adiabatic, 4);
      g_array_append_val (data, phi_norm);

      //Compute P_+(t) and append to data_probs
      g_array_append_val (data_probs, t);
      double phi_p_prob =
	phi_adiabatic[0] * phi_adiabatic[0] +
	phi_adiabatic[2] * phi_adiabatic[2];
      g_array_append_val (data_probs, phi_p_prob);

      //Compute P_+(t) in basis rotated by 45 degrees 
      //i.e. P'_+(t) and append to data_probs_prime
      g_array_append_val (data_probs_prime, t);
      double phi_p_prob_prime =
	phi_adiabatic[0] * phi_adiabatic[3] -
	phi_adiabatic[1] * phi_adiabatic[2];
      phi_p_prob_prime = 1.0 - 2.0 * phi_p_prob_prime;
      phi_p_prob_prime = phi_p_prob_prime / 2.0;
      g_array_append_val (data_probs_prime, phi_p_prob_prime);


      //If prog is running in verbose mode, then output to stdout
      if (verbose == 'y')
	{
	  data_out (t, phi_adiabatic, norm (phi_adiabatic));
	}


    }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  return status;
}
