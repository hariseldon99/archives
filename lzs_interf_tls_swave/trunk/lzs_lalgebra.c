#include "lzs_lalgebra.h"
#include "params_ics.h"
int
compare (double a, double b, double tol)
{

  double diff = fabs (a - b);
  if (diff <= tol)
    {
      return GSL_SUCCESS;
    }
  else
    return GSL_FAILURE;

}

//returns GSL_SUCCESS  if time t is an integer multiple of period
int
strobe (double t, double period)
{

  int result = GSL_FAILURE;
  double ratio = t / period;
  long int_ratio = (long) ratio;
  double difference = fabs (ratio - int_ratio);
  if (difference <= INT_TOLERANCE)
    {
      result = GSL_SUCCESS;
    }
  return result;
}


double
norm (const double phi[])
{
  int i;
  double ip = 0.0;
  for (i = 0; i < 4; i++)
    ip = ip + phi[i] * phi[i];
  return ip;
}


//Returns the inner product mod square between langle and rangle |<phi|psi>|^2
double
inner_product_modsq (const double langle[], const double rangle[])
{

  gsl_complex uk, vk;
  gsl_complex ukstar, vkstar;

  GSL_SET_COMPLEX (&uk, rangle[0], rangle[2]);
  GSL_SET_COMPLEX (&vk, rangle[1], rangle[3]);

  GSL_SET_COMPLEX (&ukstar, langle[0], -langle[2]);	//conjugate
  GSL_SET_COMPLEX (&vkstar, langle[1], -langle[3]);	//conjugate

  gsl_complex inner_product, inner_product1, inner_product2;
  double ipmodsq;

  inner_product1 = gsl_complex_mul (ukstar, uk);
  inner_product2 = gsl_complex_mul (vkstar, vk);

  inner_product = gsl_complex_add (inner_product1, inner_product2);

  ipmodsq = gsl_complex_abs2 (inner_product);

  return ipmodsq;
}



void
normalize (double phi[])
{
  int i;
  double normal = norm (phi);
  for (i = 0; i < 4; i++)
    phi[i] = phi[i] / normal;
}

//Translates to adiabatic basis. See arXiv:0911.1917 p5 for details
int
adiabatic_trans (double phi_adb[], double phi_dia[], double t, void *params)
{
  paramspace_pt *p = (paramspace_pt *) params;
  double delta_real = p->delta_real;
  double delta_imag = p->delta_imag;
  double omega = p->omega;
  double offset = p->offset;
  double ampmax = p->ampmax;
  double e = epsilon (ampmax, omega, t, offset);

  double gapsq = delta_real * delta_real + delta_imag * delta_imag;
  double omega_cap = e * e;
  omega_cap = omega_cap + gapsq;
  omega_cap = sqrt (omega_cap);

  double beta_p = omega_cap + e;
  beta_p = beta_p / (2.0 * omega_cap);
  beta_p = sqrt (beta_p);

  double beta_m = omega_cap - e;
  beta_m = beta_m / (2.0 * omega_cap);
  beta_m = sqrt (beta_m);

  gsl_complex uk, vk, phi_p, phi_m;

  GSL_SET_COMPLEX (&uk, phi_dia[0], phi_dia[2]);
  GSL_SET_COMPLEX (&vk, phi_dia[1], phi_dia[3]);

  gsl_complex rhs_p1, rhs_p2;
  gsl_complex rhs_m1, rhs_m2;

  rhs_p1 = gsl_complex_mul_real (uk, beta_m);
  rhs_p2 = gsl_complex_mul_real (vk, beta_p);
  phi_p = gsl_complex_sub (rhs_p1, rhs_p2);

  rhs_m1 = gsl_complex_mul_real (uk, beta_p);
  rhs_m2 = gsl_complex_mul_real (vk, beta_m);
  phi_m = gsl_complex_add (rhs_m1, rhs_m2);

  phi_adb[0] = GSL_REAL (phi_p);
  phi_adb[1] = GSL_REAL (phi_m);
  phi_adb[2] = GSL_IMAG (phi_p);
  phi_adb[3] = GSL_IMAG (phi_m);

  return GSL_SUCCESS;

}

//Translates to adiabatic basis for the BCS (BdG) Hamiltonian where chempot varies as time
//See adiabatic_trans_mbt_formula.jpg for details
//switch to gsl_complex at some time
int
adiabatic_trans_mbt (double phi_adb[], double phi_dia[], double t,
		     void *params)
{
  paramspace_pt *p = (paramspace_pt *) params;
  double delta_real = p->delta_real;
  double delta_imag = p->delta_imag;
  double gapsq = delta_real * delta_real + delta_imag * delta_imag;
  double omega = p->omega;
  double offset = p->offset;
  double mu0 = p->mu0;
  double muamp = p->muamp;
  double ampmax = p->ampmax;
  double epsilon = ampmax - mu (mu0, muamp, omega, t, offset);

  double ek = sqrt (epsilon * epsilon + gapsq);

  gsl_complex uk_adb, vk_adb, uk_dia, vk_dia;

  GSL_SET_COMPLEX (&uk_dia, phi_dia[0], phi_dia[2]);
  GSL_SET_COMPLEX (&vk_dia, phi_dia[1], phi_dia[3]);

  double beta_p, beta_m;

  beta_p = (ek + epsilon) / (2.0 * ek);
  beta_p = sqrt (beta_p);
  beta_m = (ek - epsilon) / (2.0 * ek);
  beta_m = sqrt (beta_m);

  gsl_complex rhs_up1, rhs_up2;
  gsl_complex rhs_dn1, rhs_dn2;

  rhs_up1 = gsl_complex_mul_real (uk_dia, beta_m);
  rhs_up2 = gsl_complex_mul_real (vk_dia, beta_p);
  uk_adb = gsl_complex_add (rhs_up1, rhs_up2);////

  rhs_dn1 = gsl_complex_mul_real (uk_dia, beta_p);
  rhs_dn2 = gsl_complex_mul_real (vk_dia, beta_m);
  vk_adb = gsl_complex_add (rhs_dn1, rhs_dn2);

  phi_adb[0] = GSL_REAL (uk_adb);
  phi_adb[1] = GSL_REAL (vk_adb);
  phi_adb[2] = GSL_IMAG (uk_adb);
  phi_adb[3] = GSL_IMAG (vk_adb);

  return GSL_SUCCESS;
}

void
set_initconds_diabatic_gnd (double phi_in[])
{


  double phi_gnd_dia[4] = { 0.0, 1.0, 0.0, 0.0 };
  phi_in[0] = phi_gnd_dia[0];
  phi_in[1] = phi_gnd_dia[1];
  phi_in[2] = phi_gnd_dia[2];
  phi_in[3] = phi_gnd_dia[3];

}

void
set_initconds_bcs_gnd (void *params, void *initconds, double phi_in[])
{

  //Get params
  paramspace_pt *p = (paramspace_pt *) params;
  double omega = p->omega;
  double offset = p->offset;
  double ampmax = p->ampmax;
  double mu0 = p->mu0;
  double muamp = p->muamp;
  double delta_real = p->delta_real;
  double delta_imag = p->delta_imag;
  double gapsq = delta_real * delta_real + delta_imag * delta_imag;

  //Get initconds
  initconds_ic *ic = (initconds_ic *) initconds;
  double t_init = ic->t_init;
  double epsilon = ampmax;	//Keeping fixed
  double chempot = mu (mu0, muamp, omega, t_init, offset);
  epsilon = epsilon - chempot;

  double ek = sqrt (epsilon * epsilon + gapsq);
  double epsilonoverek = epsilon / ek;
  double uk = 1.0 + epsilonoverek;
  uk = sqrt (uk);
  uk = M_SQRT1_2 * uk;

  double vk = 1.0 - epsilonoverek;
  vk = sqrt (vk);
  vk = M_SQRT1_2 * vk;

  phi_in[0] = uk;		//Re(uk aka uk)
  phi_in[1] = vk;		//Re(vk aka vk)
  phi_in[2] = 0.0;		//Im(uk aka uk)
  phi_in[3] = 0.0;		//Im (vk aka vk)

}

//Set semi-random uk(0), vk(0)
//See "random_number_ics.tex" for details
void
set_initconds_rand_gnd (void *params, void *initconds, double phi_in[], double rk)
{
  double uk0,vk0;
  //Get params
  paramspace_pt *p = (paramspace_pt *) params;
  double omega = p->omega;
  double offset = p->offset;
  double ampmax = p->ampmax;
  double mu0 = p->mu0;
  double muamp = p->muamp;
  double delta_real = p->delta_real;
  double delta_imag = p->delta_imag;
  double gapsq = delta_real * delta_real + delta_imag * delta_imag;
  //Get initconds
  initconds_ic *ic = (initconds_ic *) initconds;
  double t_init = ic->t_init;
  double fk = ampmax;	//Keeping fixed
  double chempot = mu (mu0, muamp, omega, t_init, offset);
  fk = fk - chempot;
  double ek = sqrt (fk * fk + gapsq);
  double rkfkratio=rk/fk;
 
   uk0 = 1.0 + (fk/ek);
   uk0  = uk0 +  rkfkratio;
   
   if((uk0>1.0)||(uk0<0.0)){
    uk0  = uk0 -rkfkratio;
    }
   
   uk0 = sqrt(uk0);
   uk0 = M_SQRT1_2 * uk0;
  
   vk0 = 1.0 - (fk/ek);
   vk0 = vk0 - rkfkratio;
   
   if((vk0>1.0)||(vk0<0.0)){
    vk0  = vk0 + rkfkratio;
    }
   
   vk0 = sqrt(vk0);
   vk0 = M_SQRT1_2 * vk0;

  phi_in[0] = uk0;		//Re(uk aka uk)
  phi_in[1] = vk0;		//Re(vk aka vk)
  phi_in[2] = 0.0;		//Im(uk aka uk)
  phi_in[3] = 0.0;		//Im (vk aka vk)
  
}

//gets the time average of probability
double
get_qfactor (GArray * a)
{
  double qfactor = 0.0;
  long probstride = 1;
  long datasize = a->len;
  while (probstride < datasize)
    {
      qfactor = qfactor + ((double *) (void *) a->data)[probstride];
      probstride = probstride + 2;
    }
  qfactor = 2.0 * qfactor / datasize;
  return qfactor;
}

//
//Calculates and returns <phi(t)|H|phi(t)> 
//          |epsilon     Delta| 
//H =       |                 |
//          |Delta*   -epsilon|
//|phi(t)> = {uk,vk}
//
double
energy_expectation (double t, const double phi[], void *params)
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

  gsl_complex delta;

  gsl_complex uk, vk;
  gsl_complex ukstar, vkstar;
  gsl_complex rhs_uk, rhs_vk;
  gsl_complex inner_product1, inner_product2;
  gsl_complex expectation;

  GSL_SET_COMPLEX (&delta, delta_real, delta_imag);

  GSL_SET_COMPLEX (&uk, phi[0], phi[2]);
  GSL_SET_COMPLEX (&vk, phi[1], phi[3]);

  gsl_complex rhs_uk1 = gsl_complex_mul_real (uk, eps);
  gsl_complex rhs_uk2 = gsl_complex_mul (delta, vk);

  rhs_uk = gsl_complex_add (rhs_uk1, rhs_uk2);

  gsl_complex rhs_vk1 = gsl_complex_mul (gsl_complex_conjugate (delta), uk);
  gsl_complex rhs_vk2 = gsl_complex_mul_real (vk, eps);

  rhs_vk = gsl_complex_sub (rhs_vk1, rhs_vk2);

  ukstar = gsl_complex_conjugate (uk);
  vkstar = gsl_complex_conjugate (vk);

  inner_product1 = gsl_complex_mul (ukstar, rhs_uk);
  inner_product2 = gsl_complex_mul (vkstar, rhs_vk);

  expectation = gsl_complex_add (inner_product1, inner_product2);

  return gsl_complex_abs(expectation);

}


//Compute the  derivative of data using cubic spline interpolation
void compute_deriv(GArray *time_quantity, GArray *time_deriv_quantity){
  
  long datasize = (time_quantity->len);
  long count;
  double t, deriv;
  GArray *time = g_array_new (FALSE, FALSE, sizeof (double));
  GArray *quantity = g_array_new (FALSE, FALSE, sizeof (double));
  
  //Split the data into time and quantity
  count=0;
  while(count<datasize){
    g_array_append_val(time,g_array_index(time_quantity, double, count));
    count++;
    g_array_append_val(quantity,g_array_index(time_quantity, double, count));
    count++;
  }
  //Convert GArray to reagular arrays
  double *time_arr = ((double *) (void *) time->data);  
  double *quantity_arr = ((double *) (void *) quantity->data);
  datasize = datasize/2;
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, datasize);
  
  gsl_spline_init (spline, time_arr, quantity_arr, datasize);
  
  for (count=0;count<datasize;count++){
    t = time_arr[count];
    deriv = gsl_spline_eval_deriv (spline, t, acc);
    g_array_append_val(time_deriv_quantity,t);
    g_array_append_val(time_deriv_quantity, deriv);
  }
  
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  g_array_free (time,TRUE);
  g_array_free (quantity,TRUE);
}

//Takes data from GArray as {{time, quantity}...}
//that is pre-sorted in time
//and returns the earliest time (scaled in units of period) when quantity is 'nonzero'
//ie the magnitude exceeds DOUBLE_TOLERANCE
double compute_nonzero (GArray *time_quantity, void *params){
  paramspace_pt *p = (paramspace_pt *) params;
  double omega = p->omega;
  
  long datasize = (time_quantity->len);
  long count=0;
  double t=0.0, scaled_time, quantity;
  
  while(count<datasize){
    t = g_array_index(time_quantity, double, count);
    count++;
    quantity = g_array_index(time_quantity, double, count);
    quantity = fabs(quantity);
    if(quantity>=DOUBLE_TOLERANCE) break;
    else count++;
  }
  scaled_time = omega * t/(2.0 * M_PI);
  return scaled_time;
}
