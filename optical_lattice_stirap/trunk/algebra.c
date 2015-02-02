/*
*  C Implementation: algebra
*
* Description:
*
*
* Author: Analabha Roy <daneel@physics.utexas.edu>, (C) 2008
*
* Copyright: See COPYING file that comes with this distribution
*
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "algebra.h"
#include "params.h"
/*****************************Defines Kronecker Delta**********************************************/
double
kdel (int i, int j)
{
  int result;
  if (i == j)
    {
      result = 1.0;
    }
  else
    result = 0.0;
  return (result);
}

int
arrcopy (double source[], double target[], int size)
{
  int i, err = 1;
  if (size < 0)
    err = 0;
  for (i = 0; i < size; i++)
    target[i] = source[i];
  return (err);
}

double
direct_prod (double a_real[], double a_imag[], double b_real[],
	     double b_imag[], int size)
{
  int k;
  double dot_absq, dot_real = 0.0, dot_imag = 0.0;
  double ar, ai, br, bi;
  for (k = 0; k < size; k++)
    {
      ar = a_real[k];
      ai = a_imag[k];
      br = b_real[k];
      bi = b_imag[k];
      dot_real = dot_real + ar * br + ai * bi;
      dot_imag = dot_imag + ar * bi - ai * br;
    }
  dot_absq = dot_real * dot_real + dot_imag * dot_imag;
  dot_absq = fabs (dot_absq);
  return (dot_absq);
}

/*Select the jth column of matrix and put into vector[]. Returns 1 if success, 0 if error*/
int
select_vec (double matrix[], int j, double vector[], int row_col_size)
{
  int i, err = 1;
  if (j < 0 || row_col_size < 0)
    err = 0;
  for (i = 0; i < row_col_size; i++)
    vector[i] = matrix[i + row_col_size * j];
  return (err);
}

/*Write vector[] into the jth column of a matrix. Returns 1 if success, 0 if error*/
int
write_vec (double vector[], int j, double matrix[], int row_col_size)
{
  int i, err = 1;
  if (j < 0 || row_col_size < 0)
    err = 0;
  for (i = 0; i < row_col_size; i++)
    matrix[i + row_col_size * j] = vector[i];
  return (err);
}


/*sorts floquet eigenvalues and eigenvetors in order of first iteration. Returns 1 if success, 0 if error*/
int
floquetsort (double evals_real[], double evals_imag[], double evecs_real[],
	     double evecs_imag[], double prev_real[], double prev_imag[],
	     double vecs_prev_real[], double vecs_prev_imag[], int size)
{
  int i, j, err;
  double dot;
  double vec_prev_real[STATENO], vec_prev_imag[STATENO];
  double vec_test_real[STATENO], vec_test_imag[STATENO];
  double evals_temp_real[STATENO], evals_temp_imag[STATENO];
  double evecs_temp_real[STATENO * STATENO],
    evecs_temp_imag[STATENO * STATENO];

  for (i = 0; i < size; i++)
    {
      err = select_vec (vecs_prev_real, i, vec_prev_real, size);
      err = select_vec (vecs_prev_imag, i, vec_prev_imag, size);
      for (j = 0; j < size; j++)
	{
	  err = select_vec (evecs_real, j, vec_test_real, size);
	  err = select_vec (evecs_imag, j, vec_test_imag, size);
	  dot =
	    direct_prod (vec_prev_real, vec_prev_imag, vec_test_real,
			 vec_test_imag, size);
	  if (dot >= TOLERANCE)
	    {
	      err = write_vec (vec_test_real, i, evecs_temp_real, size);
	      err = write_vec (vec_test_imag, i, evecs_temp_imag, size);
	      evals_temp_real[i] = evals_real[j];
	      evals_temp_imag[i] = evals_imag[j];
	      break;
	    }
	}
    }
  err = arrcopy (evals_temp_real, evals_real, size);
  err = arrcopy (evals_temp_imag, evals_imag, size);
  err = arrcopy (evecs_temp_real, evecs_real, size * size);
  err = arrcopy (evecs_temp_imag, evecs_imag, size * size);
  return (err);
}

void
print_arr (double matrix_real[], double matrix_imag[], int size)
{
  int i, j;
  for (i = 0; i < size; i++)
    {
      printf ("\n");
      for (j = 0; j < size; j++)
	printf ("%lf+I*%lf  ", matrix_real[i + size * j],
		matrix_imag[i + size * j]);
    }

}

/*Separates 4 integers into 2 pairs of equals. Returns 0 if error*/
int
getpair (int a1, int a2, int a3, int a4, int pair[])
{
  int error = 0;		//Failure by default
  if ((a1 == a2) && (a3 == a4))
    {
      error = 1;
      pair[0] = a1;
      pair[1] = a3;
    }
  if ((a1 == a3) && (a2 == a4))
    {
      error = 1;
      pair[0] = a1;
      pair[1] = a2;
    }
  if ((a1 == a4) && (a2 == a3))
    {
      error = 1;
      pair[0] = a1;
      pair[1] = a3;
    }
  return (error);
}

/*Checks to see in input matrix is unitary. returns # of non-unitary matrix elements*/
int
check_unitarity (double matrix_real[], double matrix_imag[], int size)
{
  int i, j, k, err = 0;
  double product_real, product_imag;
  for (i = 0; i < size; i++)
    {
      for (j = 0; j < size; j++)
	{
	  product_real = 0.0;
	  product_imag = 0.0;
	  for (k = 0; k < size; k++)
	    {
	      product_real =
		product_real + matrix_real[i + k * size] * matrix_real[j +
								       k *
								       size] +
		matrix_imag[i + k * size] * matrix_imag[j + k * size];
	      product_imag =
		product_imag + matrix_imag[i + k * size] * matrix_real[j +
								       k *
								       size] -
		matrix_real[i + k * size] * matrix_imag[j + k * size];
	    }
	  if (fabs (product_real - kdel (i, j)) > MNUM)
	    err = err + 1;
	  if (fabs (product_imag) > MNUM)
	    err = err + 1;
	}
    }
  return (err);
}

/*Folds unperturbed eigenstates into Floquet Brillouin Zone*/
int
fold_quasi (double array[], double omega, int size)
{
  int i, folded = 0, err = 1;
  double diff_upper, diff_lower;
  int check_foldedness (double array[], double omega, int size);
  while (folded == 0)
    {
      for (i = 0; i < size; i++)
	{
	  diff_upper = array[i] + (omega / 2.0);
	  diff_lower = array[i] - (omega / 2.0);
	  if (diff_upper <= MACHINENUM)
	    {
	      array[i] = array[i] + omega;
	    }
	  if (diff_lower >= MACHINENUM)
	    {
	      array[i] = array[i] - omega;
	    }
	}
      folded = check_foldedness (array, omega, size);
    }
  return (err);
}

int
check_foldedness (double array[], double omega, int size)
{
  int i, count = 0, folded = 0;
  double diff_upper, diff_lower;
  for (i = 0; i < size; i++)
    {
      diff_upper = array[i] + (omega / 2.0);
      diff_lower = array[i] - (omega / 2.0);
      if ((diff_lower <= MACHINENUM) && (diff_upper >= MACHINENUM))
	count++;
    }
  if (count == size)
    folded = 1;
  return (folded);
}

/*Tags floquet states by unperturbed energy*/
int
tag_quasi (int tagarray[], double quasi[], double unperturbed_quasi[],
	   int size)
{
  int i, j, err = 1;
  double quasienergy, unperturbed_quasienergy, diff;
  for (i = 0; i < size; i++)
    {
      quasienergy = quasi[i];
      for (j = 0; j < size; j++)
	{
	  unperturbed_quasienergy = unperturbed_quasi[j];
	  diff = quasienergy - unperturbed_quasienergy;
	  if (fabs (diff) <= MACHINENUM)
	    tagarray[i] = j;
	}
    }
  return (err);
}
