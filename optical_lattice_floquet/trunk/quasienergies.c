#include<stdio.h>
#include<math.h>
#include "quasienergies.h"
#include "params.h"
long
quasienergies (double evals_real[], double evals_imag[], int size,
	       double period, double quasi[])
{
  int i;
  long errorcount = 0;
  double shuldbeone;
  printf
    ("\n=========================================================================");
  printf
    ("\n=========================================================================");
  printf ("\nChecking for Unitarity of Floquet Matrix...\n");
  for (i = 0; i < size; i++)
    {
      quasi[i] = -atan2 (evals_imag[i], evals_real[i]) / period;
      shuldbeone =
	evals_real[i] * evals_real[i] + evals_imag[i] * evals_imag[i];
      if (fabs (shuldbeone - 1.0) >= MACHINENUM)
	{
	  errorcount = errorcount + 1;
	  printf
	    ("Eigenvalue %d =  %lf +I %lf , quasienergy = %lf has non-unit modulus = %lf\n",
	     i, evals_real[i], evals_imag[i], quasi[i], shuldbeone);
	}
    }
  printf ("\nQuasienergies Calculated. Number of quasienergy errors = %ld \n",
	  errorcount);
  printf
    ("\n=========================================================================");
  return (errorcount);
}
