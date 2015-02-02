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
int getpair (int a1, int a2, int a3, int a4, int pair[]);

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
c_mat (int m, int n)
{
  double value = 0.0, k = 1;
  n = n - N;
  m = m - N;
  if ((n == k && m == 0) || (n == 0 && m == k))
    {
      value = M_SQRT1_2;
    }
  else if ((n > 0 && m > 0) || (n < 0 && m < 0))
    {
      value = 0.5 * (kdel (m, n + k) + kdel (m, n - k));
    }
  else
    value = 0.0;
  return (value);
}

/**************************** This is a matrix element of sin(kx)***************************************/
//THIS IS PROBABLY WRONG. DON'T USE IT FOR NOW (IE NO MOVING LATTICE
double
s_mat (int m, int n)
{
  double value = 0.0;
  int nonzero = 0;
  n = n - N;
  m = m - N;
  if (n == 0 && m < 0)
    {
      nonzero = m;
    }

  if (m == 0 && n < 0)
    {
      nonzero = n;
    }

  if (nonzero == -1)
    {
      value = M_SQRT1_2;
    }

  if ((m < 0 && n > 0) || (n < 0 && m > 0))
    {
      value = 0.5 * (kdel (n + m + 1, 0) - kdel (n + m - 1, 0));
    }
  return (value);
}

/**************************** This is the 2ptcl interaction matrix element ***************************************/
double
in_mat (int m1, int m2, int n1, int n2)
{

  double value = 0.0;
  int error = 0;
  int pair1 = 0, pair2 = 0, pair[2];
  m1 = m1 - N;
  m2 = m2 - N;
  n1 = n1 - N;
  n2 = n2 - N;
  //Separate the integers into 2 equal pairs pair1 and pair2
  error = getpair (m1, m2, n1, n2, pair);
  if (error != 0)
    {
      pair1 = pair[0];
      pair2 = pair[1];
      if ((pair1 == pair2) && (pair1 != 0))
	value = 0.75;
      if ((pair1 == pair2) && (pair2 == 0))
	value = 0.5;
      if (((pair1 > 0) && (pair2 == 0)) || ((pair2 > 0) && (pair1 == 0)))
	value = 0.5;
      if (((pair1 < 0) && (pair2 == 0)) || ((pair2 < 0) && (pair1 == 0)))
	value = 0.5;
      if (((pair1 > 0) && (pair2 < 0)) || ((pair1 > 0) && (pair2 < 0)))
	value = 0.5;
      if ((pair1 > 0) && (pair2 > 0) && (pair1 != pair2))
	value = 0.5;
      if ((pair1 < 0) && (pair2 < 0) && (pair1 != pair2))
	value = 0.5;
    }

  value = value / (LATTICESIZE * M_PI);
  value = 2.0 * value;		//Bosonization
  return (value);
}
