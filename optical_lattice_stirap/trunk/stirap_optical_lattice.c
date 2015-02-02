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
  time_t begin, end, now, rawtime;
  struct tm *timeinfo;
  paramspace_pt param;
  int i, j, n, m, argv_count = 1;
  int m1, m2, n1, n2, count, statechosen;
  double nf, ns, error;
  double cosmatrix, sinmatrix, matelement;
  double y[2 * STATENO], y_real[STATENO], y_imag[STATENO];
  double sortedenergies[STATENO];
  double period;
  double u0;
  double prob;
  double kappa, lambda0, wf, ws, omega, t_win;
  double parameter_begin, parameter_end, parameter_tot, parameter_inc;

  gsl_matrix *h = gsl_matrix_alloc (STATENO, STATENO);	//2 ptcl hamiltonian
  gsl_vector *energies = gsl_vector_alloc (STATENO);
  gsl_matrix *evecs = gsl_matrix_alloc (STATENO, STATENO);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (STATENO);

  double eigenvectors[STATENO * STATENO], nullvec[STATENO],
    estatevec[STATENO];

  FILE *state, *energy;

  energy = fopen (argv[argv_count++], "w");
  state = fopen (argv[argv_count++], "w");
  parameter_begin = atof (argv[argv_count++]);
  parameter_end = atof (argv[argv_count++]);
  parameter_tot = atof (argv[argv_count++]);
  parameter_inc = atof (argv[argv_count++]);
  lambda0 = atof (argv[argv_count++]);

  wf = atof (argv[argv_count++]);
  ws = atof (argv[argv_count++]);

  nf = atof (argv[argv_count++]);
  ns = atof (argv[argv_count++]);

  u0 = atof (argv[argv_count++]);
  kappa = atof (argv[argv_count++]);
  statechosen = atoi (argv[argv_count++]);

  begin = time (NULL);

  error = (wf / ws) - (nf / ns);
  error = fabs (error);

  if (error >= RATTOL)
    {
      printf ("Error. Frequencies are not commensurate, exiting...");
      exit (1);
    }
  else
    {
      printf ("Commensurate frequencies found %lf / %lf = %lf / %lf ", wf, ws,
	      nf, ns);
    }
  //Calculate commensurate period
  period = M_PI * ((nf / wf) + (ns / ws));
  omega = 2.0 * M_PI / period;
  //param.kappa = kappa;
  param.kappa = kappa;
  param.lambda0 = lambda0;
  param.wf = wf;
  param.ws = ws;
  param.u0 = u0;

  param.tf = (1 / 3.0) * parameter_tot;
  param.ts = (2.0 / 3.0) * parameter_tot;
  param.td = (1 / 14.0) * parameter_tot;

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
  t_win = 0.0;
  /*Generate the 2 ptcl hamiltonian */
  for (n = 0; n < STATENO; n++)
    {
      for (m = 0; m <= n; m++)
	{
	  matelement = hamilt (&param, n, m, 0.0);
	  gsl_matrix_set (h, m, n, matelement);
	  gsl_matrix_set (h, n, m, matelement);	//Hermitian
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
  fprintf (energy, "\n %lf ", t_win);
  for (n = 0; n < STATENO; n++)
    {
      double energies_n = gsl_vector_get (energies, n);
      sortedenergies[n] = energies_n;
      fprintf (energy, "%lf ", energies_n);
    }

  j = statechosen;

  /* Initial condition=chosen state */
  for (i = 0; i < 2 * STATENO; i++)
    {
      if (i < STATENO)
	{
	  y[i] = gsl_matrix_get (evecs, i, statechosen);
	}
      else
	{
	  y[i] = 0.0;
	}
    }

  /*parameter  loop begins here */
  for (t_win = parameter_begin; t_win < parameter_end;
       t_win = t_win + parameter_inc)
    {
      param.t_win = t_win;

      now = time (NULL);
      printf ("\nELAPSED TIME = %ld sec", now - begin);
      printf ("\n\nEVOLVING SYSTEM WITH T_WIN =  %lf\n", t_win);


      /*DONE. Now evolve the initial state */
      integrate (y, t_win, t_win + parameter_inc, &param);

      /*Split into real and imaginary parts */
      for (i = 0; i < STATENO; i++)
	{
	  y_real[i] = y[i];
	  y_imag[i] = y[i + STATENO];
	}

      /*output the state elements in the hamiltonian representation */
      fprintf (state, "%lf ", t_win + parameter_inc);
      for (i = 0; i < STATENO; i++)
	{
	  /*select the ith eigenstate */
	  select_vec (eigenvectors, i, estatevec, STATENO);
	  /*Compute the probability of y in the ith eigenstate */
	  prob = direct_prod (estatevec, nullvec, y_real, y_imag, STATENO);
	  fprintf (state, "%lf ", prob);
	}
      fprintf (state, "\n");
    }
  /*parameter loop ends here */

  end = time (NULL);

  final_output (period, begin, end, 1);

  /*This is not needed if slepc is called already coz slepc will do the finalizing */
  gsl_eigen_symmv_free (w);
  gsl_matrix_free (evecs);
  gsl_matrix_free (h);
  gsl_vector_free (energies);

  fclose (state);
  fclose (energy);

  time (&rawtime);
  timeinfo = localtime (&rawtime);
  return 0;
}
