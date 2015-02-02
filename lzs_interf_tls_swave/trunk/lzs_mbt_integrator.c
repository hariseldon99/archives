#include "params_ics.h"
#include "lzs_mbt_integrator.h"
//Dynamics:
//H|phi> = i delt |phi>
//where:
//|phi> = {uk, vk}
//
//and 
//          |epsilon     Delta| 
//H =   |                              |
//        |Delta*   -epsilon|


int
func_mbt (double t, const double phi[], double f[], void *params)
{
  paramspace_pt *p = (paramspace_pt *) params;
  double delta_real = p->delta_real;
  double delta_imag = p->delta_imag;
  double omega = p->omega;
  double offset = p->offset;
  double ampmax = p->ampmax;

  double mu0 = p->mu0;
  double muamp = p->muamp;

  double eps = ampmax - mu (mu0, muamp, omega, t, offset);
  gsl_complex delta, uk, vk, rhs_up, rhs_dn;
  gsl_complex i2;

  GSL_SET_COMPLEX (&i2, 0, -1.0);

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
jac_mbt (double t, const double phi[], double *dfdphi, double dfdt[],
	 void *params)
{
  paramspace_pt *p = (paramspace_pt *) params;
  double delta_real = p->delta_real;
  double delta_imag = p->delta_imag;
  double omega = p->omega;
  double offset = p->offset;
  double ampmax = p->ampmax;

  double mu0 = p->mu0;
  double muamp = p->muamp;

  double eps = ampmax - mu (mu0, muamp, omega, t, offset);
  double depsdt = -dmudt (mu0, muamp, omega, t, offset);
  gsl_matrix_view dfdphi_mat = gsl_matrix_view_array (dfdphi, 4, 4);
  gsl_matrix *jac = &dfdphi_mat.matrix;

  //"jac" is the real 4-dimensional jacobian of the system, taken from the array dfdphi 
  //Define the complex Jacobian as J (complexjac). Essentially the 2x2 complex matrix -iH
  //with H as the Hamiltonian defined above
  //If J = J_real + I J_imag
  //where J_real and J_imag are both real 2x2 matrices,
  //then the final Jacobian "jac" is the 4x4 matrix 
  //|J_real  -J_imag|
  //|J_imag   J_real|
  //where |phi> = |phi_real> + I |phi_imag>
  //|phi_real(imag)> = {uk_real(imag), vk_real(imag)}
  //J_real is the real part of iH and J_imag is the imaginary part of -iH

  //Set J_real and J_imag
  gsl_matrix *complexjac_real = gsl_matrix_alloc (2, 2);
  gsl_matrix *complexjac_imag = gsl_matrix_alloc (2, 2);

  gsl_matrix_set (complexjac_real, 0, 0, 0.0);
  gsl_matrix_set (complexjac_real, 0, 1, delta_imag);
  gsl_matrix_set (complexjac_real, 1, 0, -delta_imag);
  gsl_matrix_set (complexjac_real, 1, 1, 0.0);

  gsl_matrix_set (complexjac_imag, 0, 0, eps);
  gsl_matrix_set (complexjac_imag, 0, 1, -delta_real);
  gsl_matrix_set (complexjac_imag, 1, 0, -delta_real);
  gsl_matrix_set (complexjac_imag, 1, 1, -eps);

  //create submatrix views of "jac"
  gsl_matrix_view jac_topleft_view = gsl_matrix_submatrix (jac, 0, 0, 2, 2);
  gsl_matrix *jac_topleft = &jac_topleft_view.matrix;

  gsl_matrix_view jac_topright_view = gsl_matrix_submatrix (jac, 0, 2, 2, 2);
  gsl_matrix *jac_topright = &jac_topright_view.matrix;

  gsl_matrix_view jac_bottomleft_view =
    gsl_matrix_submatrix (jac, 2, 0, 2, 2);
  gsl_matrix *jac_bottomleft = &jac_bottomleft_view.matrix;

  gsl_matrix_view jac_bottomright_view =
    gsl_matrix_submatrix (jac, 2, 2, 2, 2);
  gsl_matrix *jac_bottomright = &jac_bottomright_view.matrix;

  //copy J_real and J_imag to these submatrices as shown in comments above
  gsl_matrix_memcpy (jac_topleft, complexjac_real);
  gsl_matrix_memcpy (jac_bottomleft, complexjac_imag);

  gsl_matrix_scale (complexjac_imag, -1.0);
  gsl_matrix_memcpy (jac_topright, complexjac_imag);
  gsl_matrix_memcpy (jac_bottomright, complexjac_real);


  //Now, set the explicit time derivative of the rhs to dfdt array
  //In 2 dimensional complex representation, this is just
  //-I*{\dot epsilon uk + Delta vk, Delta* uk-\dot epsilon vk}
  //In the 4 dimensional real repr. Split the above up to real part (top) and imaginary part (down)
  //this will be the dfdt

  gsl_complex delta, uk, vk, rhs_up, rhs_dn;
  gsl_complex i2;

  GSL_SET_COMPLEX (&i2, 0, -1.0);

  GSL_SET_COMPLEX (&delta, delta_real, delta_imag);

  GSL_SET_COMPLEX (&uk, phi[0], phi[2]);
  GSL_SET_COMPLEX (&vk, phi[1], phi[3]);

  gsl_complex rhs_up1 = gsl_complex_mul_real (uk, depsdt);
  gsl_complex rhs_up2 = gsl_complex_mul (delta, vk);
  rhs_up = gsl_complex_add (rhs_up1, rhs_up2);

  gsl_complex rhs_dn1 = gsl_complex_mul (gsl_complex_conjugate (delta), uk);
  gsl_complex rhs_dn2 = gsl_complex_mul_real (vk, depsdt);
  rhs_dn = gsl_complex_sub (rhs_dn1, rhs_dn2);

  //Real
  dfdt[0] = GSL_REAL (rhs_up);
  dfdt[1] = GSL_REAL (rhs_dn);

  //Imag
  dfdt[2] = GSL_IMAG (rhs_up);
  dfdt[3] = GSL_IMAG (rhs_dn);

  gsl_matrix_free (complexjac_real);
  gsl_matrix_free (complexjac_imag);

  return GSL_SUCCESS;
}

/*****This function actually runs the full integration of a particulat set of IC's from 0-t**********/
int
integrate_kvec_adiabatic (double *phi, double *phi_adiabatic, void *initconds,
			  void *params, double t_init, double t_final)
{
  int status = 0;

  double dt, t, omega, phi_norm;

  /*This sets up the gsl ODE system structure */
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rkck;
  /*Working algorithms: rk4imp(most accurate, slowest), rkf45(less accurate,faster),rkck(less accurate still, faster still) */
  gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 4);
  gsl_odeiv_control *c =
    gsl_odeiv_control_standard_new (ABSERROR, RELERROR, YERROR, YPRIMEERROR);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc (4);
  gsl_odeiv_system mbt_sys = { func_mbt, jac_mbt, 4, params };
  dt = DT;

  //extract initial conditions
  paramspace_pt *p = (paramspace_pt *) params;
  omega = p->omega;

  initconds_ic *ic = (initconds_ic *) initconds;
  char verbose = ic->verbosity;


  t = t_init;
  if (verbose == 'y')
    begin_tout ();
  while (t < t_final)
    {
      status =
	gsl_odeiv_evolve_apply (e, c, s, &mbt_sys, &t, t_final, &dt, phi);
      if (status != GSL_SUCCESS)
	{
	  status = GSL_FAILURE;
	  break;
	}

      int ierr = adiabatic_trans_mbt (phi_adiabatic, phi, t, params);
      status = ierr;

      //If prog is running in verbose mode, then output to stdout
      if (verbose == 'y')
	{
	  phi_norm = norm (phi);
	  data_out (t, phi, phi_norm);
	}


    }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  return status;
}

/*****This function actually runs the full integration of a particulat set of IC's from 0-t**********/
int
integrate_kvec_adiabatic_nonadaptive (double *phi, double *phi_adiabatic,
				      void *initconds, void *params,
				      double t_init, double t_final)
{
  int status = 0;

  double dt, t, omega, phi_norm;
  double dphidt_in[4], dphidt_out[4], phi_err[4];

  /*This sets up the gsl ODE system structure */
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rkck;
  /*Working algorithms: rk4imp(most accurate, slowest), rkf45(less accurate,faster),rkck(less accurate still, faster still) */
  gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 4);
  gsl_odeiv_system mbt_sys = { func_mbt, jac_mbt, 4, params };


  //extract initial conditions
  paramspace_pt *p = (paramspace_pt *) params;
  omega = p->omega;

  dt = DT_NONADAPTIVE;
  dt = dt * 2.0 * M_PI / omega;



  initconds_ic *ic = (initconds_ic *) initconds;
  char verbose = ic->verbosity;

  t = t_init;
  if (verbose == 'y')
    begin_tout ();

  /* initialise dydt_in from system parameters */
  GSL_ODEIV_FN_EVAL (&mbt_sys, t, phi, dphidt_in);

  while (t < t_final)
    {
      status =
	gsl_odeiv_step_apply (s, t, dt, phi, phi_err, dphidt_in, dphidt_out,
			      &mbt_sys);
      if (status != GSL_SUCCESS)
	{
	  status = GSL_FAILURE;
	  break;
	}

      normalize (phi);
      dphidt_in[0] = dphidt_out[0];
      dphidt_in[1] = dphidt_out[1];
      dphidt_in[2] = dphidt_out[2];
      dphidt_in[3] = dphidt_out[3];

      t += dt;

      int ierr = adiabatic_trans_mbt (phi_adiabatic, phi, t, params);
      status = ierr;


      //If prog is running in verbose mode, then output to stdout
      if (verbose == 'y')
	{
	  phi_norm = norm (phi);
	  data_out (t, phi, phi_norm);
	}


    }

  gsl_odeiv_step_free (s);
  return status;
}


/*****This function actually runs the full integration of a particulat set of IC's from 0-t**********/
int
integrate_kvec_diabatic (double *phi, void *initconds, void *params,
			 double t_init, double t_final)
{
  int status = 0;
  double dt, t, omega, phi_norm;

  /*This sets up the gsl ODE system structure */
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rkck;
  /*Working algorithms: rk4imp(most accurate, slowest), rkf45(less accurate,faster),rkck(less accurate still, faster still) */
  gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 4);
  gsl_odeiv_control *c =
    gsl_odeiv_control_standard_new (ABSERROR, RELERROR, YERROR, YPRIMEERROR);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc (4);
  gsl_odeiv_system mbt_sys = { func_mbt, jac_mbt, 4, params };
  dt = DT;

  //extract initial conditions
  paramspace_pt *p = (paramspace_pt *) params;
  omega = p->omega;

  initconds_ic *ic = (initconds_ic *) initconds;
  char verbose = ic->verbosity;

  t = t_init;
  if (verbose == 'y')
    begin_tout ();
  while (t < t_final)
    {
      status =
	gsl_odeiv_evolve_apply (e, c, s, &mbt_sys, &t, t_final, &dt, phi);
      if (status != GSL_SUCCESS)
	{
	  status = GSL_FAILURE;
	  break;
	}

      //If prog is running in verbose mode, then output to stdout
      if (verbose == 'y')
	{
	  phi_norm = norm (phi);
	  data_out (t, phi, phi_norm);
	}


    }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
  return status;
}

/*****This function actually runs the full integration of a particulat set of IC's from 0-t**********/
int
integrate_kvec_diabatic_nonadaptive (double *phi, void *initconds,
				     void *params, double t_init,
				     double t_final)
{
  int status = 0;

  double dt, t, omega, phi_norm;
  double dphidt_in[4], dphidt_out[4], phi_err[4];

  /*This sets up the gsl ODE system structure */
  const gsl_odeiv_step_type *T = gsl_odeiv_step_rkck;
  /*Working algorithms: rk4imp(most accurate, slowest), rkf45(less accurate,faster),rkck(less accurate still, faster still) */
  gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 4);
  gsl_odeiv_system mbt_sys = { func_mbt, jac_mbt, 4, params };


  //extract initial conditions
  paramspace_pt *p = (paramspace_pt *) params;
  omega = p->omega;

  dt = DT_NONADAPTIVE;
  dt = dt * 2.0 * M_PI / omega;


  initconds_ic *ic = (initconds_ic *) initconds;
  char verbose = ic->verbosity;


  t = t_init;
  if (verbose == 'y')
    begin_tout ();


  /* initialise dydt_in from system parameters */
  GSL_ODEIV_FN_EVAL (&mbt_sys, t, phi, dphidt_in);

  while (t < t_final)
    {
      status =
	gsl_odeiv_step_apply (s, t, dt, phi, phi_err, dphidt_in, dphidt_out,
			      &mbt_sys);
      if (status != GSL_SUCCESS)
	{
	  status = GSL_FAILURE;
	  break;
	}
      normalize (phi);

      dphidt_in[0] = dphidt_out[0];
      dphidt_in[1] = dphidt_out[1];
      dphidt_in[2] = dphidt_out[2];
      dphidt_in[3] = dphidt_out[3];

      t += dt;

      //If prog is running in verbose mode, then output to stdout
      if (verbose == 'y')
	{
	  phi_norm = norm (phi);
	  data_out (t, phi, phi_norm);
	}


    }

  gsl_odeiv_step_free (s);
  return status;
}

/*****This function actually runs the integration of a particulat set of IC's by one dt**********/
/*****This routine does so by internally updating the system by Euler's method to first order*****/
int
integrate_kvec_adiabatic_euler1 (double phi[], double phi_adiabatic[],
				 void *initconds, void *params, double t_init,
				 double t_final)
{

  int status = GSL_SUCCESS;

  double t, phi_norm;
  double dt = (t_final - t_init);

  gsl_complex idt = gsl_complex_rect (0, dt);

  gsl_complex uk, vk;
  gsl_complex uk_rhs, vk_rhs;
  gsl_complex uk_rhs1, vk_rhs1;
  gsl_complex uk_rhs2, vk_rhs2;
  gsl_complex uk_updated, vk_updated;

  //extract initial conditions
  paramspace_pt *p = (paramspace_pt *) params;
  double omega = p->omega;

  //Get delta
  double delta_real = p->delta_real;
  double delta_imag = p->delta_imag;
  gsl_complex delta = gsl_complex_rect (delta_real, delta_imag);

  double offset = p->offset;
  double ampmax = p->ampmax;

  double mu0 = p->mu0;
  double muamp = p->muamp;

  double eps = ampmax - mu (mu0, muamp, omega, t_init, offset);

  initconds_ic *ic = (initconds_ic *) initconds;
  char verbose = ic->verbosity;
  t = t_init;
  if (verbose == 'y')
    begin_tout ();
  while (t < t_final)
    {
      //Update phi here
      GSL_SET_COMPLEX (&uk, phi[0], phi[2]);
      GSL_SET_COMPLEX (&vk, phi[1], phi[3]);

      uk_rhs1 = gsl_complex_mul_real (uk, eps);
      uk_rhs2 = gsl_complex_mul (vk, delta);
      uk_rhs = gsl_complex_add (uk_rhs1, uk_rhs2);
      uk_rhs = gsl_complex_mul (uk_rhs, idt);
      uk_updated = gsl_complex_sub (uk, uk_rhs);

      vk_rhs1 = gsl_complex_mul (uk, gsl_complex_conjugate (delta));
      vk_rhs2 = gsl_complex_mul_real (vk, eps);
      vk_rhs = gsl_complex_sub (vk_rhs1, vk_rhs2);
      vk_rhs = gsl_complex_mul (vk_rhs, idt);
      vk_updated = gsl_complex_sub (vk, vk_rhs);

      phi[0] = GSL_REAL (uk_updated);
      phi[1] = GSL_REAL (vk_updated);
      phi[2] = GSL_IMAG (uk_updated);
      phi[3] = GSL_IMAG (vk_updated);
      //Finished updating phi

      if (status != GSL_SUCCESS)
	{
	  status = GSL_FAILURE;
	  break;
	}

      int ierr = adiabatic_trans_mbt (phi_adiabatic, phi, t, params);
      status = ierr;

      //If prog is running in verbose mode, then output to stdout
      if (verbose == 'y')
	{
	  phi_norm = norm (phi);
	  data_out (t, phi, phi_norm);
	}

      t += dt;
    }
  return status;
}
