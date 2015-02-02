#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <stdlib.h>
#include "params.h"
#include "integrator.h"
/*********************This evaluates the Differential Equation System*****************************/
int
func (double t, const double y[], double dydt[], void *param)
{
  drive_and_tower *p = (drive_and_tower *) param;
  double y_re[DIMS], y_im[DIMS], dydt_re[DIMS], dydt_im[DIMS];
  double hamilt_mat[DIMS*DIMS];
  int nt, mt;
  int dim=p->dim;

  /*Build the hamiltonian matrix elements*/
  
  for(nt=0;nt<dim;nt++)
	  for(mt=0;mt<=nt;mt++){
	  hamilt_mat[nt+dim*mt]=hamilt (param,nt,mt,t);
	  hamilt_mat[mt+dim*nt]=hamilt_mat[nt+dim*mt];//Always hermitian
	  }
	  
 /*Write out the (complicated) expression for RHS & assign to dydt[i] */
  for (nt = 0; nt < dim; nt++)
    {
      y_re[nt] = y[nt];
      y_im[nt] = y[nt + dim];
    }

  {
    for (nt = 0; nt < dim; nt++)
      {
	dydt_re[nt] = 0.0;
	dydt_im[nt] = 0.0;
	for (mt = 0; mt < dim; mt++)
	  {
	    dydt_re[nt] = dydt_re[nt] + hamilt_mat[nt+dim*mt] * y_im[mt];
	    dydt_im[nt] = dydt_im[nt] - hamilt_mat[nt+dim*mt] * y_re[mt];
	  }
      }
  }

  for (nt = 0; nt < dim; nt++)
    {
      dydt[nt] = dydt_re[nt];
      dydt[dim + nt] = dydt_im[nt];
    }
  
  return GSL_SUCCESS;
}

/***** IT'S ONLY NEEDED FOR BSIMP NOT NEEDED FOR RUNGE KUTTA METHODS******/
int
jac (double t, const double y[], double *dfdy, double dfdt[], void *param)
{
  drive_and_tower *p = (drive_and_tower *) param;
  int dim=p->dim;
  int nt,mt;
  double hamilt_mat, hamilt_mat_deriv,jacobian_mat_deriv;
  
  /*view the array dfdy as a gsl_matrix and assign it to jacobian*/
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2 * dim, 2 * dim);
  gsl_matrix *jacobian = &dfdy_mat.matrix;
  gsl_matrix *jacobian_deriv =gsl_matrix_alloc (2*dim, 2*dim);
  
  
  /*Because I'm using poor man's complex numbers, if the Hamiltonian matrix is H, then the Jacobian J is: */
  /*					(0	H)
  					(-H	0)				*/
  
  /*Get the relevant submatrices of the jacobian...*/
  gsl_matrix_view submatrix1 = gsl_matrix_submatrix(jacobian,0,dim,dim,dim);
  gsl_matrix *upperright_submatrix=&submatrix1.matrix;
  
  gsl_matrix_view submatrix2 = gsl_matrix_submatrix(jacobian,dim,0,dim,dim);
  gsl_matrix *lowerleft_submatrix=&submatrix2.matrix;
  
  /* .. and it's derivative*/
  gsl_matrix_view submatrix1_deriv = gsl_matrix_submatrix(jacobian_deriv,0,dim,dim,dim);
  gsl_matrix *upperright_submatrix_deriv=&submatrix1_deriv.matrix;
  
  gsl_matrix_view submatrix2_deriv = gsl_matrix_submatrix(jacobian_deriv,dim,0,dim,dim);
  gsl_matrix *lowerleft_submatrix_deriv=&submatrix2_deriv.matrix;
  
  
  /*First, set jacobian and derivative to 0*/
  gsl_matrix_set_zero(jacobian);
  gsl_matrix_set_zero(jacobian_deriv);
  
  /*Then assign the upper right and lower left submatrices, which are hermitian.    *
   *The jacobian is of column size 2*dim. Therefore the submatrices (hamiltonian) *
   * are of column size dim                 									   */
  
  for(nt=0;nt<dim;nt++)
	  for(mt=0;mt<=nt;mt++){
	  	hamilt_mat=hamilt (param,nt,mt,t);
		hamilt_mat_deriv=hamilt_deriv(param,nt,mt,t);
		
		/*submatrices of the jacobian*/
  		gsl_matrix_set(upperright_submatrix,nt,mt,hamilt_mat);
  		gsl_matrix_set(upperright_submatrix,mt,nt,hamilt_mat);
				
		gsl_matrix_set(lowerleft_submatrix,nt,mt,-hamilt_mat);
		gsl_matrix_set(lowerleft_submatrix,mt,nt,-hamilt_mat);
		
		/*submatrices of the jacobian's derivative*/
		gsl_matrix_set(upperright_submatrix_deriv,nt,mt,hamilt_mat_deriv);
		gsl_matrix_set(upperright_submatrix_deriv,mt,nt,hamilt_mat_deriv);
		
		gsl_matrix_set(lowerleft_submatrix_deriv,nt,mt,-hamilt_mat_deriv);
		gsl_matrix_set(lowerleft_submatrix_deriv,mt,nt,-hamilt_mat_deriv);
		
		}

 /*Now, calculate dfdt using the jacobian's derivative*/
	for(nt=0;nt<2*dim;nt++){
		dfdt[nt]=0.0;
		for(mt=0;mt<2*dim;mt++){
		jacobian_mat_deriv=gsl_matrix_get(jacobian_deriv,nt,mt);
		dfdt[nt] = dfdt[nt] + jacobian_mat_deriv * y[mt];
		}
	}
  gsl_matrix_free (jacobian_deriv);
  return GSL_SUCCESS;
}



/*****This function actually runs the full integration of a particulat set of IC's from 0-T**********/
void
integrate (double *input, double initial, double final, void *param)
{
  drive_and_tower *p = (drive_and_tower *) param;
  int dim=p->dim;
  int status;
  double dt, t;
  /*This sets up the gsl ODE system structure */
  const gsl_odeiv_step_type *T = gsl_odeiv_step_bsimp;
  /*Working algorithms: rk4imp(most accurate, slowest), rkf45(less accurate,faster),rkck(less accurate still, faster still)*/
  gsl_odeiv_step *s = gsl_odeiv_step_alloc (T, 2 * dim);
  gsl_odeiv_control *c =
    gsl_odeiv_control_standard_new (ABSERROR, RELERROR, YERROR, YPRIMEERROR);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc (2 * dim);
  gsl_odeiv_system sys = { func, jac, 2 * dim, param };
  dt = DT;
  t = initial;
  while (t < final)
    {
      status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, final, &dt, input);
      if (status != GSL_SUCCESS)
	{
	  printf ("\n GSL execution of integration failed, bailing..");
	  break;
	}
    }
  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);
}
