/*
// C++ Interface: floquet_optical_lattice
//
// Description:
//
//
// Author: Analabha Roy <daneel@utexas.edu>, (C) 2013
//
// Copyright: See COPYING file that comes with this distribution
//
*/

//Include header files
#include <mpi.h>		/* for multinode jobs, MPI/OPENMP Hybridized */
#include <omp.h>		/* for parallelization using OpenMP */
#include <stdio.h>		/* for input/output */
#include <stdarg.h>
#include <stdlib.h>

#include <string.h>		/* for strings */

#include <math.h>		/* for basic math routines */
#include <gsl/gsl_math.h>

#include <time.h>		/* for profiling runtime */

#include <unistd.h>		/* for getopt */
#include <getopt.h>		/* for getopt */

#include <gsl/gsl_rng.h>	/* for random numbers */
#include <gsl/gsl_randist.h>

#include <gsl/gsl_sf_bessel.h>

#include <glib.h>


//Declare the integrator function
int integrate (double *input, double initial, double final, void *param);

//Declare the initial conditions setters
int set_initconds_ghz (double y[]);
int set_initconds_allup (double y[]);

//Declare the magnetization and unitarity check
double magnetization (const double y[]);
double check_unitarity (const double y[]);

#define BETATRUNC 50

struct globalArgs_t
{
  double gamma_a;		/* --drv-amp -a option */
  int besj_zero;		/* --zero_order -e option */
  double eta;
  double omega;
  double alpha;			/* --pert -p option */
  double finaltime;		/* --final-time -t option */
  int tsteps;			/* --time-steps -n option */
  char *outFileName_Mag;	/* --mag_output -m option */
  FILE *outFile_Mag;
  int check_unitarity;		/* --chk_unit -u option */
  int verbosity;		/* --verbose -v option */
  char **inputFiles;		/* input files */
  int numInputFiles;		/* # of input files */
} globalArgs;

static const char *optString = "a:e:p:t:n:m:uh:vh";


static const struct option longOpts[] = {
  {"drv-amp", required_argument, 0, 'a'},
  {"zero_order", required_argument, 0, 'e'},
  {"pert", required_argument, 0, 'p'},
  {"final-time", required_argument, 0, 't'},
  {"time-steps", required_argument, 0, 'n'},
  {"mag_output", required_argument, 0, 'm'},
  {"chk_unit", no_argument, 0, 'u'},
  {"verbose", no_argument, 0, 'v'},
  {0, 0, 0, 0}
};

const char options[] = {
  "OPTIONS:\tFLAGS:\t\tVALUE:\n --drv-amp\t-a\t\tdrive amplitude\n --zero_order\t-e\t\tthe order of zero of bessel function\n --pert\t\t-p\t\tperturbation amplitude\n --final-time\t-t\t\tfinal time\n --time-steps\t-n\t\tnumber of time steps\n --mag_output\t-m\t\t/path/to/output/file\n --chk_unit\t-u\t\tchecks for unitary evolution\n --verbose\t-v\t\tverbose output\n   none\t \t\t\tprint this help"
};

const char envv[] = {
  "ENVIRONMENT VARIABLES:\n GSL_RNG_TYPE\t\ttype of random number generator to be used\n GSL_RNG_SEED\t\tseed value of random number generator\n OMP_NUM_THREADS\tnumber of OpenMP threads"
};
