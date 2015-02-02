#include "gsl_determinant_complex.h"
double
gsl_determinant_complex (gsl_matrix_complex * A, int inPlace)
{

/*
  inPlace = 1 => A is replaced with the LU decomposed copy.
  inPlace = 0 => A is retained, and a copy is used for LU.
*/
  gsl_complex det;
  int signum;
  gsl_permutation *p = gsl_permutation_alloc (A->size1);
  if (inPlace)
    {
      gsl_linalg_complex_LU_decomp (A, p, &signum);
      det = gsl_linalg_complex_LU_det (A, signum);
    }
  else
    {
      gsl_matrix_complex *tmpA =
	gsl_matrix_complex_alloc (A->size1, A->size2);
      gsl_matrix_complex_memcpy (tmpA, A);
      gsl_linalg_complex_LU_decomp (tmpA, p, &signum);
      det = gsl_linalg_complex_LU_det (tmpA, signum);
      gsl_matrix_complex_free (tmpA);
    }

  gsl_permutation_free (p);
  return gsl_complex_abs (det);
}
