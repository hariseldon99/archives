#include "initconds.h"
#include "params.h"

/*********************This sets the initial conditions to GHZ state*****************************/
int
set_initconds_ghz (double y[])
{
  int i, j;
  double elem = -0.5;

  gsl_matrix *H = gsl_matrix_calloc (2 * N, 2 * N);

  gsl_vector *eval = gsl_vector_alloc (2 * N);
  gsl_matrix *evec = gsl_matrix_alloc (2 * N, 2 * N);

    /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     *Now, build the Hamiltonian H as
     *(+A +B)
     *(-B -A)
     *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  //+A submatrix
  gsl_matrix_view p_amat = gsl_matrix_submatrix (H, 0, 0, N, N);
  gsl_matrix *plusa = &p_amat.matrix;

  //+B submatrix
  gsl_matrix_view p_bmat = gsl_matrix_submatrix (H, 0, N, N, N);
  gsl_matrix *plusb = &p_bmat.matrix;

  //-B submatrix
  gsl_matrix_view m_bmat = gsl_matrix_submatrix (H, N, 0, N, N);
  gsl_matrix *minusb = &m_bmat.matrix;

  //-A submatrix
  gsl_matrix_view m_amat = gsl_matrix_submatrix (H, N, N, N, N);
  gsl_matrix *minusa = &m_amat.matrix;

   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    *With \alpha = 1, h_0 = 0, h_i = 0 +-A is tridiagonal with 
    * no diagonal components
    *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  gsl_matrix_set_zero (plusa);
  gsl_matrix_set_zero (minusa);

  for (i = 0; i < N; i++)
    {
      if (i < N - 1)
	{
	  gsl_matrix_set (plusa, i, i + 1, elem);
	  gsl_matrix_set (plusa, i + 1, i, elem);

	  gsl_matrix_set (minusa, i, i + 1, -elem);
	  gsl_matrix_set (minusa, i + 1, i, -elem);

	}

    }
   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    *Now, build the +B matrix
    *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  for (i = 0; i < N; i++)
    {
      if (i < N - 1)
	{
	  gsl_matrix_set (plusb, i, i + 1, elem);
	  gsl_matrix_set (plusb, i + 1, i, -elem);
	}
    }

   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    *Now, build the -B matrix
    *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  for (i = 0; i < N; i++)
    {
      if (i < N - 1)
	{
	  gsl_matrix_set (minusb, i, i + 1, -elem);
	  gsl_matrix_set (minusb, i + 1, i, elem);
	}
    }

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   *Now, diagonalize H using gsl
   *and get the ground state
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (2 * N);
  gsl_eigen_symmv (H, eval, evec, w);

  //Sort from -\epsilon_\mu to +\epsilon_\mu
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  * The first N evecs constitute the GHZ state
  * After sorting, these evecs are the ones with nonpositive energies
  * so...
  * The first RIGHTMOST N X N submatrix of evecs are  u
  * The N X N block of elements below this block in evec are v 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  gsl_matrix_view vview = gsl_matrix_submatrix (evec, 0, N, N, N);
  gsl_matrix *v = &vview.matrix;

  gsl_matrix_view uview = gsl_matrix_submatrix (evec, N, N, N, N);
  gsl_matrix *u = &uview.matrix;

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  *The first N^2 elements of y[] are the real parts of u
  *The next N^2 elements of y[] are the imaginary parts of u 
  *The next N^2 elements of y[] are the real parts of v
  *The final N^2 elements of y[] are the imaginary parts of v
  * so...
  *The first N^2 elements of y[] are u
  *The next N^2 elements of y[] are 0 since u is real
  *The next N^2 elements of y[] are v
  *The final N^2 elements of y[] are 0 since v is real
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
	{
	  y[i * N + j] = gsl_matrix_get (u, i, j);
	  y[i * N + j + N * N] = 0.0;
	  y[i * N + j + 2 * N * N] = gsl_matrix_get (v, i, j);
	  y[i * N + j + 3 * N * N] = 0.0;
	}
    }

  gsl_eigen_symmv_free (w);
  gsl_matrix_free (H);
  gsl_vector_free (eval);
  gsl_matrix_free (evec);

  return GSL_SUCCESS;
}

/*********************This sets the init conds to the noninteracting ground state*****************************/
/****************************************i.e. all ising spins up**********************************************/

int
set_initconds_allup (double y[])
{
  int i, j;

  gsl_matrix *u = gsl_matrix_alloc (N, N);
  gsl_matrix_set_zero (u);

  gsl_matrix *v = gsl_matrix_alloc (N, N);
  gsl_matrix_set_identity (v);


  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  *The first N^2 elements of y[] are the real parts of u
  *The next N^2 elements of y[] are the imaginary parts of u 
  *The next N^2 elements of y[] are the real parts of v
  *The final N^2 elements of y[] are the imaginary parts of v
  * so...
  *The first N^2 elements of y[] are u
  *The next N^2 elements of y[] are 0 since u is real
  *The next N^2 elements of y[] are v
  *The final N^2 elements of y[] are 0 since v is real
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  for (i = 0; i < N; i++)
    {
      for (j = 0; j < N; j++)
	{
	  y[i * N + j] = gsl_matrix_get (u, i, j);
	  y[i * N + j + N * N] = 0.0;
	  y[i * N + j + 2 * N * N] = gsl_matrix_get (v, i, j);
	  y[i * N + j + 3 * N * N] = 0.0;
	}
    }

  gsl_matrix_free (u);
  gsl_matrix_free (v);

  return GSL_SUCCESS;
}
