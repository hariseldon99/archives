#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include <gsl/gsl_math.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "params.h"
double classical_int (double x1, double x2, double u0);
double classical_pot (double x, double v0);
double kdel (int i, int j);
/*Husimi Function Integrals*/
gsl_complex
gaussian_intg (double x, double p, double sigma, double a)
{
  gsl_complex result, exponent;
  exponent = gsl_complex_rect (-a * p * sigma * sigma, 2.0 * x);
  exponent = gsl_complex_mul_real (exponent, 0.5 * a * p);
  result = gsl_complex_exp (exponent);
  result = gsl_complex_mul_real (result, M_SQRTPI * M_SQRT2 * sigma);
  return (result);
}

gsl_complex
gaussian_prod (double x1, double x2, double p1, double p2, double s1,
	       double s2, double bm1m2, double am1, double am2)
{
  gsl_complex result;
  result = gsl_complex_polar (1.0, bm1m2);
  result = gsl_complex_mul (result, gaussian_intg (x1, p1, s1, am1));
  result = gsl_complex_mul (result, gaussian_intg (x2, p2, s2, am2));
  return (result);
}

gsl_complex
pbox_husimi01 (int m, double x1, double x2, double p1, double p2, double s1,
	       double s2, void *param)
{
  drive_and_tower *p = (drive_and_tower *) param;
  gsl_complex result;
  int m1 = p->tower[m][0];
  int m2 = p->tower[m][1];
  double am1, am2, bm1m2, denr;
  /*don;t forget to increment m1 & m2 coz they are counted from 0 in the tower of states */
  m1++;
  m2++;

  am1 = m1 * M_PI / (2.0 * L);
  am2 = m2 * M_PI / (2.0 * L);
  bm1m2 = -(m1 + m2) * M_PI / 2.0;
  result = gaussian_prod (x1, x2, p1, p2, s1, s2, bm1m2, am1, am2);

  am1 = am1;
  am2 = -am2;
  bm1m2 = (m2 - m1) * M_PI / 2.0;
  result =
    gsl_complex_sub (result,
		     gaussian_prod (x1, x2, p1, p2, s1, s2, bm1m2, am1, am2));

  am1 = -am1;
  am2 = -am2;
  bm1m2 = -bm1m2;
  result =
    gsl_complex_sub (result,
		     gaussian_prod (x1, x2, p1, p2, s1, s2, bm1m2, am1, am2));

  am1 = am1;
  am2 = -am2;
  bm1m2 = (m1 + m2) * M_PI / 2.0;
  result =
    gsl_complex_add (result,
		     gaussian_prod (x1, x2, p1, p2, s1, s2, bm1m2, am1, am2));

  denr = sqrt (2.0 * (1 + kdel (m1, m2)));
  result = gsl_complex_div_real (result, denr);
  return (result);
}

gsl_complex
pbox_husimi02 (int m, double x1, double x2, double p1, double p2, double s1,
	       double s2, void *param)
{
  drive_and_tower *p = (drive_and_tower *) param;
  gsl_complex result;
  /*m1 & m2 re interchanged from 01 function above bcoz of bosonization */
  int m1 = p->tower[m][1];
  int m2 = p->tower[m][0];
  double am1, am2, bm1m2, denr;
  /*don;t forget to increment m1 & m2 coz they are counted from 0 in the tower of states */
  m1++;
  m2++;

  am1 = m1 * M_PI / (2.0 * L);
  am2 = m2 * M_PI / (2.0 * L);
  bm1m2 = -(m1 + m2) * M_PI / 2.0;
  result = gaussian_prod (x1, x2, p1, p2, s1, s2, bm1m2, am1, am2);

  am1 = am1;
  am2 = -am2;
  bm1m2 = (m2 - m1) * M_PI / 2.0;
  result =
    gsl_complex_sub (result,
		     gaussian_prod (x1, x2, p1, p2, s1, s2, bm1m2, am1, am2));

  am1 = -am1;
  am2 = -am2;
  bm1m2 = -bm1m2;
  result =
    gsl_complex_sub (result,
		     gaussian_prod (x1, x2, p1, p2, s1, s2, bm1m2, am1, am2));

  am1 = am1;
  am2 = -am2;
  bm1m2 = (m1 + m2) * M_PI / 2.0;
  result =
    gsl_complex_add (result,
		     gaussian_prod (x1, x2, p1, p2, s1, s2, bm1m2, am1, am2));

  denr = sqrt (2.0 * (1 + kdel (m1, m2)));
  result = gsl_complex_div_real (result, denr);
  return (result);
}

gsl_complex
husimi (int n, double x1, double p1, void *param)
{
  drive_and_tower *p = (drive_and_tower *) param;
  double s1, s2, u0, v0;
  gsl_complex result = gsl_complex_rect (0.0, 0.0), temp;
  double x2 = STROBE, p2;
  /*For the first time, we extract ALL the parameters from param */
  double sortedcoeff[STATENO][STATENO], statetower[STATENO][2],
    energylevels[STATENO];
  int i, j, j1, j2, m;
  u0 = p->u0;
  v0 = p->v0;
  s1 = p->s1;
  s2 = p->s2;
  for (i = 0; i < STATENO; i++)
    {
      statetower[i][0] = p->tower[i][0];
      statetower[i][1] = p->tower[i][1];
      energylevels[i] = p->eigenvalues[i];
      for (j = 0; j < STATENO; j++)
	sortedcoeff[i][j] = p->eigenvectors[i][j];
    }
  /*test classical energy conservation */
  p2 =
    energylevels[n] - (p1 * p1) - classical_pot (x1, v0) - classical_pot (x2,
									  v0)
    - classical_int (x1, x2, u0);
  if (p2 >= 0.0)
    {
      /*sum over m */
      for (m = 0; m < STATENO; m++)
	{
	  temp =
	    gsl_complex_add (pbox_husimi01 (m, x1, x2, p1, p2, s1, s2, p),
			     pbox_husimi02 (m, x1, x2, p1, p2, s1, s2, p));
	  temp = gsl_complex_mul_real (temp, sortedcoeff[m][n]);
	  result = gsl_complex_add (result, temp);
	}
    }
  return (result);
}
