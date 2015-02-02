#include "responses.h"
#include "params.h"

/************************This evaluates the Magnetization ************************************/
double
magnetization (const double y[])
{
  int i, j;
  double mag;
  gsl_matrix_complex *v = gsl_matrix_complex_alloc (N, N);
  gsl_complex velem;

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Now, build the v matrix
   * The first N^2 elements of y are the real parts of u
   * The next N^2 elements of y are the imaginary parts of u
   * The next N^2 elements of y are the real parts of v
   * The final N^2 elements of y are the imaginary parts of v
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      {
	GSL_SET_COMPLEX (&velem, y[i * N + j + 2 * N * N],
			 y[i * N + j + 3 * N * N]);
	gsl_matrix_complex_set (v, i, j, velem);
      }
  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Now, compute Tr (v^\dagger v) = 
   * 			\sum_i v^\dagger_{row i} \cdot v_{col i} 
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  gsl_complex trace_loc;
  gsl_vector_complex_view col;

  gsl_complex trace = GSL_COMPLEX_ZERO;
  for (i = 0; i < N; i++)
    {

      col = gsl_matrix_complex_column (v, i);
      gsl_blas_zdotc (&col.vector, &col.vector, &trace_loc);
      trace = gsl_complex_add (trace, trace_loc);
    }

  mag = gsl_complex_abs (trace);
  mag = (2.0 / N) * mag;	//Calculate mag per site
  mag = mag - 1.0;

  gsl_matrix_complex_free (v);
  return mag;
}

/***********************************G_{ij} from writeup***************************************/
gsl_complex
g_ij (int i, int j, gsl_matrix_complex * u, gsl_matrix_complex * v)
{
  gsl_complex gij;
  gsl_matrix_complex *vpu = gsl_matrix_complex_alloc (N, N);
  gsl_matrix_complex *vmu = gsl_matrix_complex_alloc (N, N);
  gsl_matrix_complex_memcpy (vpu, v);
  gsl_matrix_complex_memcpy (vmu, v);
  gsl_matrix_complex_add (vpu, u);
  gsl_matrix_complex_sub (vmu, u);

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * G = vmu \times vpu^\dagger 
   * Extract the ith row of vmu and jth ROW of vpu
   * Then do blas dot prod vpu^conj and vmu
   *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  gsl_vector_complex_view vmu_row = gsl_matrix_complex_row (vmu, i);
  gsl_vector_complex_view vpu_row = gsl_matrix_complex_row (vpu, j);
  gsl_blas_zdotc (&vpu_row.vector, &vmu_row.vector, &gij);

  gsl_matrix_complex_free (vpu);
  gsl_matrix_complex_free (vmu);

  return gij;
}


/************************Returns outprod = a \ocross b***************************/
gsl_matrix_complex *
kronecker_prod (gsl_matrix_complex * a, gsl_matrix_complex * b)
{
  int i, j, k, l;
  int m, p, n, q;
  m = a->size1;
  p = a->size2;
  n = b->size1;
  q = b->size2;

  gsl_matrix_complex *c = gsl_matrix_complex_alloc (m * n, p * q);
  gsl_complex da, db;

  for (i = 0; i < m; i++)
    {
      for (j = 0; j < p; j++)
	{
	  da = gsl_matrix_complex_get (a, i, j);
	  for (k = 0; k < n; k++)
	    {
	      for (l = 0; l < q; l++)
		{
		  db = gsl_matrix_complex_get (b, k, l);
		  gsl_matrix_complex_set (c, n * i + k, q * j + l,
					  gsl_complex_mul (da, db));
		}
	    }
	}
    }

  return c;
}



/*******This gets the density matrix of the (N/2-1)\otimes N/2 \otimes (N/2+1) bipartite*******/
gsl_matrix_complex *
midpoint_density_matrix_get (const double y[])
{
  gsl_complex temp;
  int i, j;
  int site = floor (N / 2);
  gsl_matrix_complex *rho = gsl_matrix_complex_alloc (4, 4);

  //See Guiseppe Santoro's notes page 12 for the meaning of this
  gsl_matrix_complex *u = gsl_matrix_complex_alloc (N, N);
  gsl_matrix_complex *v = gsl_matrix_complex_alloc (N, N);

  gsl_complex uelem, velem;

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Now, build the u and v matrices
   * The first N^2 elements of y are the real parts of u
   * The next N^2 elements of y are the imaginary parts of u
   * The next N^2 elements of y are the real parts of v
   * The final N^2 elements of y are the imaginary parts of v
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      {
	GSL_SET_COMPLEX (&uelem, y[i * N + j], y[i * N + j + N * N]);
	gsl_matrix_complex_set (u, i, j, uelem);

	GSL_SET_COMPLEX (&velem, y[i * N + j + 2 * N * N],
			 y[i * N + j + 3 * N * N]);
	gsl_matrix_complex_set (v, i, j, velem);
      }

  //Pauli matrices
  gsl_matrix_complex *id = gsl_matrix_complex_alloc (2, 2);
  gsl_matrix_complex *sx = gsl_matrix_complex_alloc (2, 2);
  gsl_matrix_complex *sy = gsl_matrix_complex_alloc (2, 2);
  gsl_matrix_complex *sz = gsl_matrix_complex_alloc (2, 2);
  gsl_matrix_complex_set_identity (id);
  gsl_matrix_complex_set_zero (sx);
  gsl_matrix_complex_set (sx, 0, 1, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set (sx, 1, 0, GSL_COMPLEX_ONE);
  gsl_matrix_complex_set_zero (sy);
  gsl_matrix_complex_set (sy, 0, 1, gsl_complex_rect (0, -1));
  gsl_matrix_complex_set (sy, 1, 0, gsl_complex_rect (0, 1));
  gsl_matrix_complex_set_identity (sz);
  gsl_matrix_complex_set (sz, 1, 1, GSL_COMPLEX_NEGONE);

  //Expectation Values in eq 13-15 of Hatano arXiv:0706.4162
  gsl_complex sz_expct_site = g_ij (site, site, u, v);
  gsl_complex sz_expct_sitep1 = g_ij (site + 1, site + 1, u, v);
  gsl_complex sxsx_expct = g_ij (site, site + 1, u, v);
  gsl_complex sysy_expct = g_ij (site + 1, site, u, v);
  gsl_complex szsz_expct = gsl_complex_mul (g_ij (site, site, u, v),
					    g_ij (site + 1, site + 1, u, v));
  temp =
    gsl_complex_mul (g_ij (site, site + 1, u, v),
		     g_ij (site + 1, site, u, v));
  szsz_expct = gsl_complex_sub (szsz_expct, temp);

  //Now calculate eq 12 in Hatano paper
  //These are the six terms in eq 12
  gsl_matrix_complex *id_full = gsl_matrix_complex_alloc (4, 4);
  gsl_matrix_complex_set_identity (id_full);

  gsl_matrix_complex *term2, *term3, *term4, *term5, *term6;
  term2 = kronecker_prod (sz, id);
  gsl_matrix_complex_scale (term2, sz_expct_site);

  term3 = kronecker_prod (id, sz);
  gsl_matrix_complex_scale (term3, sz_expct_sitep1);

  term4 = kronecker_prod (sx, sx);
  gsl_matrix_complex_scale (term4, sxsx_expct);

  term5 = kronecker_prod (sy, sy);
  gsl_matrix_complex_scale (term5, sysy_expct);

  term6 = kronecker_prod (sz, sz);
  gsl_matrix_complex_scale (term6, szsz_expct);

  //Now, add them and scale by 1/4
  gsl_matrix_complex_memcpy (rho, id_full);
  gsl_matrix_complex_add (rho, term2);
  gsl_matrix_complex_add (rho, term3);
  gsl_matrix_complex_add (rho, term4);
  gsl_matrix_complex_add (rho, term5);
  gsl_matrix_complex_add (rho, term6);
  gsl_matrix_complex_scale (rho, gsl_complex_rect (0.25, 0.0));

  gsl_matrix_complex_free (id_full);
  gsl_matrix_complex_free (term2);
  gsl_matrix_complex_free (term3);
  gsl_matrix_complex_free (term4);
  gsl_matrix_complex_free (term5);
  gsl_matrix_complex_free (term6);

  gsl_matrix_complex_free (id);
  gsl_matrix_complex_free (sx);
  gsl_matrix_complex_free (sy);
  gsl_matrix_complex_free (sz);

  gsl_matrix_complex_free (u);
  gsl_matrix_complex_free (v);

  return rho;
}


/***********This evaluates the von Neumann Entropy of a density matrix ***********************/
/******************The density matrix is destroyed in the process*****************************/
double
von_neumann_entropy (gsl_matrix_complex * rho)
{
  int i;
  double eval_i = 0.0;
  gsl_complex term = GSL_COMPLEX_ZERO;

  gsl_complex entropy = GSL_COMPLEX_ZERO;

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Diagonalize the hermitian density matrix
   * Get eigenvalues only. No need to sort
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  gsl_vector *eval = gsl_vector_alloc (rho->size1);
  gsl_eigen_herm_workspace *w = gsl_eigen_herm_alloc (rho->size1);
  gsl_eigen_herm (rho, eval, w);
  gsl_eigen_herm_free (w);

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Calculate entropy = \sum_i eval_i * ln (eval_i)
   * Remember to exclude zero evals
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  for (i = 0; i < rho->size1; i++)
    {
      eval_i = gsl_vector_get (eval, i);
      if (fabs (eval_i) >= TOL)
	{
	  term =
	    gsl_complex_mul_real (gsl_complex_log
				  (gsl_complex_rect (eval_i, 0.0)), eval_i);
	  entropy = gsl_complex_add (entropy, term);
	}
    }

  gsl_vector_free (eval);
  return -gsl_complex_abs (entropy);
}

/**************************This checks for unitary evolution *********************************/
double
check_unitarity (const double y[])
{
  int i, j;

  gsl_matrix_complex *u = gsl_matrix_complex_alloc (N, N);
  gsl_matrix_complex *v = gsl_matrix_complex_alloc (N, N);
  gsl_matrix_complex *udu_vdv_utv_vtu_mone = gsl_matrix_complex_alloc (N, N);

  gsl_matrix_complex *unity = gsl_matrix_complex_alloc (N, N);

  gsl_matrix_complex_set_identity (unity);

  gsl_complex uelem, velem;

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   * Now, build the u and v matrices
   * The first N^2 elements of y are the real parts of u
   * The next N^2 elements of y are the imaginary parts of u
   * The next N^2 elements of y are the real parts of v
   * The final N^2 elements of y are the imaginary parts of v
   * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      {
	GSL_SET_COMPLEX (&uelem, y[i * N + j], y[i * N + j + N * N]);
	GSL_SET_COMPLEX (&velem, y[i * N + j + 2 * N * N],
			 y[i * N + j + 3 * N * N]);
	gsl_matrix_complex_set (u, i, j, uelem);
	gsl_matrix_complex_set (v, i, j, velem);
      }

  //Compute udu + vdv
  gsl_blas_zgemm (CblasConjTrans, CblasNoTrans, GSL_COMPLEX_ONE, u, u,
		  GSL_COMPLEX_ZERO, udu_vdv_utv_vtu_mone);
  gsl_blas_zgemm (CblasConjTrans, CblasNoTrans, GSL_COMPLEX_ONE, v, v,
		  GSL_COMPLEX_ONE, udu_vdv_utv_vtu_mone);

  //Compute udu + vdv + utv + vtu
  gsl_blas_zgemm (CblasTrans, CblasNoTrans, GSL_COMPLEX_ONE, u, v,
		  GSL_COMPLEX_ONE, udu_vdv_utv_vtu_mone);
  gsl_blas_zgemm (CblasTrans, CblasNoTrans, GSL_COMPLEX_ONE, v, u,
		  GSL_COMPLEX_ONE, udu_vdv_utv_vtu_mone);

  //Compute udu + vdv + utv + vtu - 1
  gsl_matrix_complex_sub (udu_vdv_utv_vtu_mone, unity);

  //This should be 0. To check whether matrix A is 0, compute Tr(AdA)
  gsl_complex trace_loc;
  gsl_vector_complex_view col;

  gsl_complex trace = GSL_COMPLEX_ZERO;
  for (i = 0; i < N; i++)
    {
      col = gsl_matrix_complex_column (udu_vdv_utv_vtu_mone, i);
      gsl_blas_zdotc (&col.vector, &col.vector, &trace_loc);
      trace = gsl_complex_add (trace, trace_loc);
    }

  gsl_matrix_complex_free (u);
  gsl_matrix_complex_free (v);
  gsl_matrix_complex_free (udu_vdv_utv_vtu_mone);

  gsl_matrix_complex_free (unity);

  double tr = gsl_complex_abs (trace);
  tr = tr / N;

  return tr;

}
