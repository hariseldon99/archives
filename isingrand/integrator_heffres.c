#include "integrator_heffres.h"
#include "params.h"

/*********************This evaluates the Differential Equation System*****************************/
int
func (double t, const double y[], double dydt[], void *param)
{
  double elem1, elem2;

  //Get the parameters
  paramspace_pt *p = (paramspace_pt *) param;


 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Now, build the u and v matrices by matrix views
   * The first N^2 elements of y are the real parts of u
   * The next N^2 elements of y are the imaginary parts of u
   * The next N^2 elements of y are the real parts of v
   * The final N^2 elements of y are the imaginary parts of v
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  gsl_matrix_const_view ureal_view =
    gsl_matrix_const_view_array (&y[0], N, N);
  const gsl_matrix *ureal = &ureal_view.matrix;
  gsl_matrix_const_view uimag_view =
    gsl_matrix_const_view_array (&y[N * N], N, N);
  const gsl_matrix *uimag = &uimag_view.matrix;
  gsl_matrix_const_view vreal_view =
    gsl_matrix_const_view_array (&y[2 * N * N], N, N);
  const gsl_matrix *vreal = &vreal_view.matrix;
  gsl_matrix_const_view vimag_view =
    gsl_matrix_const_view_array (&y[3 * N * N], N, N);
  const gsl_matrix *vimag = &vimag_view.matrix;

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Now, build the dudt and dvdt matrices by matrix views
   * The first N^2 elements of dydt are the real parts of dudt
   * The next N^2 elements of dydt are the imaginary parts of dudt
   * The next N^2 elements of dydt are the real parts of dvdt
   * The final N^2 elements of dydt are the imaginary parts of dvdt
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  gsl_matrix_view dudtreal_view = gsl_matrix_view_array (&dydt[0], N, N);
  gsl_matrix *dudtreal = &dudtreal_view.matrix;
  gsl_matrix_view dudtimag_view = gsl_matrix_view_array (&dydt[N * N], N, N);
  gsl_matrix *dudtimag = &dudtimag_view.matrix;
  gsl_matrix_view dvdtreal_view
    = gsl_matrix_view_array (&dydt[2 * N * N], N, N);
  gsl_matrix *dvdtreal = &dvdtreal_view.matrix;
  gsl_matrix_view dvdtimag_view
    = gsl_matrix_view_array (&dydt[3 * N * N], N, N);
  gsl_matrix *dvdtimag = &dvdtimag_view.matrix;

  gsl_matrix *a = gsl_matrix_calloc (N, N);	//set all elements to 0
  gsl_matrix *b = gsl_matrix_calloc (N, N);

  //Build the A and B matrices 
  int i;
  for (i = 0; i < N; i++)
    {
      elem1 = (p->hrand[i]);
      elem1 = -elem1;
      gsl_matrix_set (a, i, i, elem1);
      if (i < N - 1)
	{
	  elem2 = (p->jrand[i]);
	  elem2 = elem2 / 2.0;

	  gsl_matrix_set (b, i, i + 1, elem2);
	  gsl_matrix_set (b, i + 1, i, -elem2);
	}
    }

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * First, dudt = (2/i) (Au + Bv) using BLAS, i.e.
   * dudt_real =  2 ( A u_imag + B v_imag ) , 
   * dudt_imag = -2 ( A u_real + B v_real )
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, a, uimag, 0.0, dudtreal);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 2.0, b, vimag, 2.0, dudtreal);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, a, ureal, 0.0, dudtimag);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -2.0, b, vreal, -2.0, dudtimag);

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Then, dvdt = -(2/i) (Av + Bu) using BLAS, i.e.
   * dvdt_real = -2 ( B u_imag + A v_imag ),
   * dvdt_imag =  2 ( B u_real + A v_real )
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, b, uimag, 0.0, dvdtreal);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -2.0, a, vimag, -2.0, dvdtreal);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, b, ureal, 0.0, dvdtimag);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 2.0, a, vreal, 2.0, dvdtimag);

  gsl_matrix_free (a);
  gsl_matrix_free (b);

  return GSL_SUCCESS;
}

/***** IT'S ONLY NEEDED FOR BSIMP NOT NEEDED FOR RUNGE KUTTA METHODS******/
/****                         f=dydt                                *****/
/****   Currently just a dummy function. It does nothing            ****/
int
jac (double t, const double y[], double *dfdy, double dfdt[], void *param)
{

  return GSL_SUCCESS;
}


/*****This function actually runs the full integration of a particulat set of IC's from 0-T**********/
int
integrate (double *input, double initial, double final, void *param)
{
  int status;
  double dt = DT;
  gsl_odeiv2_system sys = { func, jac, 4 * N * N, param };
  gsl_odeiv2_driver *driver =
    gsl_odeiv2_driver_alloc_y_new (&sys, STEP_TYPE, dt, ABSERROR, RELERROR);

  status = gsl_odeiv2_driver_apply (driver, &initial, final, input);
  if (status != GSL_SUCCESS)
    {
      printf ("\n GSL execution of integration failed, bailing..");
    }
  gsl_odeiv2_driver_free (driver);
  return status;
}
