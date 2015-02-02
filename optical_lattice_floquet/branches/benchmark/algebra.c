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
