
#include <stdio.h>
#include <stdlib.h>
#include<math.h>

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

/*Subroutines for calculating quantum hamiltonian matrix elements*/

/*****************This is the lth single particle-in-a-box energy level*****************************/
double
En (int l)
{
  double energy;
  l = l + 1;
  energy = M_PI / 2.0;
  energy = energy * energy;
  energy = energy * l * l;
  energy = energy / (L * L);
  return (energy);
}

/**************************** This is a matrix element of x*****************************************/
double
f1 (int m, int n)
{
  double value;
  m = m + 1;			/*Arrays are counted from n=0, but levels from n=1 */
  n = n + 1;
  if ((m + n) % 2 == 1)
    {
      value = 16.0 * m * n;
      value = value / (M_PI * M_PI);
      value = value / (m * m - n * n);
      value = value / (m * m - n * n);
      value = value * L;
    }
  else
    {
      value = 0.0;
    }
  return (value);
}

/**************************** This is a matrix element of x^2**************************************/
double
f2 (int m, int n)
{
  double value, temp;
  m = m + 1;			/*Arrays are counted from n=0, but levels from n=1 */
  n = n + 1;
  if ((m + n) % 2 == 0)
    {
      if (m == n)
	{
	  value = ((1 / 3.0) - (2.0 / (n * n * M_PI * M_PI))) * L * L;
	}
      else
	{
	  value =
	    32.0 * m * n * L * L / (M_PI * M_PI * (m * m - n * n) *
				    (m * m - n * n));
	}
    }
  else
    {
      value = 0.0;
    }
  return (value);
}

/**************************** This is a matrix element of x^4***************************************/
double
f4 (int m, int n)
{
  double value;
  m = m + 1;			/*Arrays are counted from n=0, but levels from n=1 */
  n = n + 1;
  if ((m + n) % 2 == 0)
    {
      if (m == n)
	{
	  value =
	    L * L * L * L * ((1 / 5.0) - (4.0 / (n * n * M_PI * M_PI)) +
			     (24.0 /
			      (n * n * n * n * M_PI * M_PI * M_PI * M_PI)));
	}
      else
	{
	  value =
	    (64.0 * m * n / (M_PI * M_PI * M_PI * M_PI)) *
	    (pow (M_PI / (m * m - n * n), 2) -
	     (48.0 * (m * m + n * n) / pow (m * m - n * n, 4))) * L * L * L *
	    L;
	}
    }
  else
    {
      value = 0.0;
    }
  return (value);
}

/*Matrix elements of the Double Well Potential*/
double
pot (int n, int m)
{
  double f2 (int, int), f4 (int, int);
  double potelm;
  potelm = (-2.0 * f2 (n, m) + f4 (n, m));
  return (potelm);
}

/*Matrix Element of the Interaction*/
double
in (int n1, int n2, int n3, int n4)
{
  double term1 = 1.0 / (4.0 * L), term2 = 1.0 / (2.0 * L), term3 =
    3.0 / (4.0 * L);
  double interelement = 0.0;	/*Default Value of Interaction */
  n1 = n1 + 1;
  n2 = n2 + 1;
  n3 = n3 + 1;
  n4 = n4 + 1;

  if ((n1 - n2) == (n3 + n4))
    {
      if ((n1 - n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = -term1;
	}
    }

  if ((n1 + n2) == (n3 + n4))
    {
      if ((n1 + n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = term1;
	}
    }

  if ((n1 - n2) == (n4 - n3))
    {
      if ((n1 - n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = term1;
	}
    }

  if ((n1 + n2) == (n4 - n3))
    {
      if ((n1 + n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = -term1;
	}
    }

  if ((n1 - n2) == (n3 - n4))
    {
      if ((n1 - n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = term1;
	}
    }

  if ((n1 + n2) == (n3 - n4))
    {
      if ((n1 + n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = -term1;
	}
    }

  if ((n1 - n2) == (-n3 - n4))
    {
      if ((n1 - n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = -term1;
	}
    }

  if ((n1 + n2) == (-n3 - n4))
    {
      if ((n1 + n2) == 0)
	{
	  interelement = term2;
	}
      else
	{
	  interelement = term1;
	}
    }

  if ((n1 == n2) && (n2 == n3) && (n3 == n4))
    {
      interelement = term3;
    }
  /*Bosonization */
  interelement = 2.0 * interelement;
  return (interelement);
}
