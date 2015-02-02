#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "params.h"
#include "stirap_optical_lattice.h"
#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char **argv)
{
  time_t begin, end, rawtime;
  struct tm *timeinfo;
  paramspace_pt param;

  int  n, m, argv_count = 1;
  int m1, m2, n1, n2, count;
  //double nf, ns, error;
  double cosmatrix, sinmatrix, matelement;
  //double y[2 * STATENO], y_real[STATENO], y_imag[STATENO];
  double sortedenergies[STATENO];
  //double period;
  double u0;
  //double prob;
  double kappa;
  double parameter_begin, parameter_end, parameter_inc;

  gsl_matrix *h = gsl_matrix_alloc (STATENO, STATENO);	//2 ptcl hamiltonian
  gsl_vector *energies = gsl_vector_alloc (STATENO);
  gsl_matrix *evecs = gsl_matrix_alloc (STATENO, STATENO);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (STATENO);

  double eigenvectors[STATENO * STATENO], nullvec[STATENO];

  FILE *energy;

  energy = fopen (argv[argv_count++], "w");
  //state = fopen (argv[argv_count++], "w");
  parameter_begin = atof (argv[argv_count++]);
  parameter_end = atof (argv[argv_count++]);
  //parameter_tot = atof (argv[argv_count++]);
  parameter_inc = atof (argv[argv_count++]);
  //lambda0 = atof (argv[argv_count++]);

  //wf = atof (argv[argv_count++]);
  //ws = atof (argv[argv_count++]);

  //nf = atof (argv[argv_count++]);
  //ns = atof (argv[argv_count++]);

  u0 = atof (argv[argv_count++]);
  //kappa = atof (argv[argv_count++]);
  //statechosen = atoi (argv[argv_count++]);

  begin = time (NULL);

  //error = (wf / ws) - (nf / ns);
  //error = fabs (error);

  //if (error >= RATTOL)
    //{
      //printf ("Error. Frequencies are not commensurate, exiting...");
      //exit (1);
    //}
  //else
    //{
      //printf ("Commensurate frequencies found %lf / %lf = %lf / %lf ", wf, ws,
	      //nf, ns);
    //}
  //Calculate commensurate period
  //period = M_PI * ((nf / wf) + (ns / ws));
  //omega = 2.0 * M_PI / period;
  //param.kappa = kappa;

  for(kappa=parameter_begin;kappa<=parameter_end;kappa=kappa+parameter_inc){
  
  param.kappa = kappa;
  param.lambda0 = 0.0;
  param.wf = 0.0;
  param.ws = 0.0;
  param.u0 = u0;

  param.tf = 1.0;
  param.ts = 1.0;
  param.td = 1.0;

  /*create a null vec */
  for (n = 0; n < STATENO; n++)
    nullvec[n] = 0.0;

  /*Calculate the cosine(x) Matrix elements in the Hamiltonian Representation
   *And put them in struct param */
  for (n = 0; n < DIM; n++)
    {
      param.psq[n] = psq_mat (n);
      for (m = 0; m <= n; m++)
	{
	  cosmatrix = c_mat (n, m);
	  param.cosmatrix[n + m * DIM] = cosmatrix;
	  param.cosmatrix[m + n * DIM] = cosmatrix;	//Hermitian

	  sinmatrix = s_mat (n, m);
	  param.sinmatrix[n + m * DIM] = sinmatrix;
	  param.sinmatrix[m + n * DIM] = sinmatrix;	//Hermitian
	}
    }

  /*Build the 2-particle tower of states */
  count = 0;
  for (n = 0; n < DIM; n++)
    {
      for (m = 0; m <= n; m++)
	{
	  param.towerofstates[count][0] = n;
	  param.towerofstates[count][1] = m;
	  count++;
	}
    }

  /*Calculate the interaction matrix elements and put them in struct param */
  for (m = 0; m < STATENO; m++)
    {
      for (n = 0; n <= m; n++)
	{
	  m1 = param.towerofstates[m][0];
	  m2 = param.towerofstates[m][1];
	  n1 = param.towerofstates[n][0];
	  n2 = param.towerofstates[n][1];
	  matelement = u0 * in_mat (m1, m2, n1, n2);
	  param.interaction[m + STATENO * n] = matelement;
	  param.interaction[n + STATENO * m] = matelement;
	}
    }
  //t_win = 0.0;
  /*Generate the 2 ptcl hamiltonian */
  for (n = 0; n < STATENO; n++)
    {
      for (m = 0; m <= n; m++)
	{
	  matelement = hamilt (&param, n, m, 0.0);
	  gsl_matrix_set (h, m, n, matelement);
	  gsl_matrix_set (h, n, m, matelement);	//Hermitian
	  //printf("\n %lf : (%d, %d) = %lf", kappa,n,m,matelement);
	}
    }

  /*diagonalize the  2 ptcl hamiltonian */
  gsl_eigen_symmv (h, energies, evecs, w);
  gsl_eigen_symmv_sort (energies, evecs, GSL_EIGEN_SORT_VAL_ASC);

  /*write out the evecs to eigenvectors array */
  for (n = 0; n < STATENO; n++)
    {
      for (m = 0; m < STATENO; m++)
	{
	  eigenvectors[n + STATENO * m] = gsl_matrix_get (evecs, n, m);
	}
    }

  /*tell root process to write out parameter and the eigenvalues */
  fprintf (energy, "\n %lf ", kappa);
  for (n = 0; n < STATENO; n++)
    {
      double energies_n = gsl_vector_get (energies, n);
      sortedenergies[n] = energies_n;
      fprintf (energy, "%lf ", energies_n);
    }
  }
  end = time (NULL);

  final_output (0.0, begin, end, 1);

  /*This is not needed if slepc is called already coz slepc will do the finalizing */
  gsl_eigen_symmv_free (w);
  gsl_matrix_free (evecs);
  gsl_matrix_free (h);
  gsl_vector_free (energies);

  fclose (energy);

  time (&rawtime);
  timeinfo = localtime (&rawtime);
  return 0;
}
