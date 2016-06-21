#include "cantor_pair.h"

double
cantor_pair (long x, long y)
{
  long p1 = 2 * y;
  long p2 = (y + 1) * y;
  double res = (x * x) + (3 + p1) * x + p2;
  return res / 2.0;
}
