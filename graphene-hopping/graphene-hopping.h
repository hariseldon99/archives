#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#include <mpi.h>
#include <glib.h>

typedef struct
{
  double vxbar;
  double vx;
  double vy;
} pot_amps;

typedef struct
{
  double min;
  double max;
  double inc;
  long size;
} grid;