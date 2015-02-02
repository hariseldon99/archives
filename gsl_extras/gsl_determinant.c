#include "gsl_determinant.h"
double
gsl_determinant (gsl_matrix * A, int inPlace)
{
/*
  inPlace = 1 => A is replaced with the LU decomposed copy.
  inPlace = 0 => A is retained, and a copy is used for LU.
*/
  double det;
  int signum;
  gsl_permutation *p = gsl_permutation_alloc (A->size1);
  if (inPlace)
    {
      gsl_linalg_LU_decomp (A, p, &signum);
      det = gsl_linalg_LU_det (A, signum);
    }

  else
    {
      gsl_matrix *tmpA = gsl_matrix_alloc (A->size1, A->size2);
      gsl_matrix_memcpy (tmpA, A);
      gsl_linalg_LU_decomp (tmpA, p, &signum);
      det = gsl_linalg_LU_det (tmpA, signum);
      gsl_matrix_free (tmpA);
    }

  gsl_permutation_free (p);
  return det;
}
