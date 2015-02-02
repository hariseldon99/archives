
#include "params.h"
#include <gsl/gsl_math.h>
/*subroutines for calculating classical hamiltonian functions*/
double
classical_int (double x1, double x2, double u0)
{
  double interaction;
  double argument;
  argument = (x1 - x2) * (x1 - x2) / (2 * SIGMA * SIGMA);
  interaction = u0 / SIGMA;
  interaction = interaction / sqrt (2 * M_PI);
  interaction = interaction * exp (-argument);
  return (interaction);
}

double
classical_pot (double x, double v0)
{
  double term;
  term = v0 * (-2 * x * x + x * x * x * x);
  return (term);
}
