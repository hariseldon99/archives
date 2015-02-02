/*
*  C Implementation: quantum_hamilt_melements
*
* Description:
*
*
* Author: Analabha Roy <daneel@physics.utexas.edu>, (C) 2007
*
* Copyright: See COPYING file that comes with this distribution
*
*/

#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include<gsl/gsl_math.h>
#include "params.h"
double kdel (int, int);
/*Subroutines for calculating quantum hamiltonian matrix elements*/

/**************************** This is a matrix element of p^2***************************************/
double
psq_mat (int n)
{
  double value;
  n = n - N;			/*levels counted from -N to +N */
  value = n * n;
  return (value);
}

/**************************** This is a matrix element of cos(x)***************************************/
double
c_mat (int np, int n)
{
  double value = 0.0, k = 1;
  if ((n == k && np == 0) || (n == 0 && np == k))
    {
      value = M_SQRT1_2;
    }
  else if ((n > 0 && np > 0) || (n < 0 && np < 0))
    {
      value = 0.5 * (kdel (np, n + k) + kdel (np, n - k));
    }
  else
    value = 0.0;
  return (value);
}
