/*
// C++ Interface: lzs_tls
//
// Description:
//
//
// Author: Analabha Roy <daneel@bose.res.in>, (C) 2010
//
// Copyright: See COPYING file that comes with this distribution
//
*/
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <glib.h>

#define ARGNO 16		//16 arguments to lzs_tls bin

int arg_err (void);
void init_out (void *initconds);
void params_out (void *params);
int integrate (double *phi, void *initconds, void *params, GArray * data,
	       GArray * data_probs, GArray * data_probs_prime);
double norm (const double *phi);
void lzs_intgerr (int p);
void dump_out (GArray * a, FILE * filepointer);
void dump_out_probs (GArray * a, FILE * filepointer, void *params);

void norm_stats (void *normstats, GArray * a, double normexpect);
void final_out (time_t begin, time_t end, double normexpect, void *stats);
