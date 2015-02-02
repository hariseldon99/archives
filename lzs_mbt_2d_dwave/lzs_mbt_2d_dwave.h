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

#include <gsl/gsl_complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>

#include <mpi.h>
#include <glib.h>

#define ARGNO 22		//22 arguments to bin
#define ALPHA 0.0       //Anisotropy offset. Default is 0
int arg_err (int argc);
void init_out (void *initconds);
void begin_out ();
void params_out (void *params);
double delta_k (double kx, double ky, double lattice_size);
int copyvec (double source[], double target[], long size);
int fermisurface_2d (double kx, double ky, double ef, double spread);
void set_initconds_diabatic_gnd (double phi_in[]);
void set_initconds_bcs_gnd (void *params, void *initconds, double phi_in[]);
void set_initconds_rand_gnd (void *params, void *initconds, double phi_in[], double rk);
int adiabatic_trans_mbt (double phi_adb[], double phi_dia[], double t,
			 void *params);
void copytobuffer (GArray * a, double buffer[], long bufsize, long kstride);
long update_delta_loc (GArray * phik, void *params);

int integrate_kvec_diabatic (double *phi, void *initconds,
			     void *params, double t_init, double t_final);
int integrate_kvec_adiabatic (double *phi, double *phi_adiabatic,
			      void *initconds, void *params, double t_init,
			      double t_final);
double energy_expectation (double t, const double phi[], void *params);
double inner_product_modsq (const double langle[], const double rangle[]);
void copyfrombuffer (double buffer[], long bufsize, GArray * a, long kstride);
int strobe (double t, double period);
void dump_out_strobe (double scaled_time, double mag, double fidelity,
		      double delta, double defect_density, double resen,
		      double fidsus);
double norm (const double *phi);
int compare (double a, double b, double tol);
void lzs_intgerr (int p);
void compute_deriv(GArray *time_quantity, GArray *time_deriv_quantity);
void dump_out_probs (GArray * a, FILE * filepointer, void *params);
void dump_out_delta (GArray * a, FILE * filepointer, void *params);

void prn_array (double *a, long size);
double get_qfactor (GArray * a);
double get_magavg (GArray * a);

double compute_nonzero (GArray *time_quantity, void *params);
void final_out_mbt (time_t begin, time_t end, double gridsize,
		    void *qfactors);
void dump(const char quantity[], double value);
void output_stdev (GArray *time_fidelity_deriv_mag);
