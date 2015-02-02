/*C Program for calculating 2-particle husimi functions after p2 is integrated and plotting it's cross section*/
/*Here, the full 4-dimentional husimi function is calculated*/

#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#include <mpi.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>

#include "params.h"
double classical_int (double x1, double x2, double u0);
double classical_pot (double x, double v0);
double En (int l);
double f1 (int m, int n);
double f2 (int m, int n);
double f3 (int m, int n);
double pot (int m, int n);
double in (int n1, int n2, int n3, int n4);
double undriven_hamilt (void *param, int nt, int mt);
gsl_complex gaussian_intg (double x, double p, double sigma, double a);
gsl_complex gaussian_prod (double x1, double x2, double p1, double s1,
			   double s2, double bm1m2, double am1, double am2);
gsl_complex pbox_husimi01 (int m, double x1, double x2, double p1, double s1,
			   double s2, void *param);
gsl_complex pbox_husimi02 (int m, double x1, double x2, double p1, double s1,
			   double s2, void *param);
gsl_complex husimi (int n, double x1, double p1, void *param);

int
main (int argc, char *argv[])
{
  time_t begin, end, rawtime;
  struct tm *timeinfo;
  drive_and_tower param;
  int i, j, k, l, n, m, n1, n2, m1, m2, towerofstates[STATENO][2];
  int np, pid, statechosen;
  long errorcount = 0, plotpoints;
  double v0, u0, s1, s2, element, x, p, energy, radius;
  gsl_complex h;
  int id, nprocs, dist_size, quotient, remainder, rank;
  double xinitglobal = XINIT, xfinalglobal = XFINAL;
  double xinit, xfinal, pinit, pfinal;
  int *datasizes, *offsets;
  double points[3 * BUFSIZE];
  int count = 0, total = 0;
  FILE *output;

  output = fopen (argv[1], "w");
  statechosen = atof (argv[2]);
  s1 = atof (argv[3]);
  s2 = atof (argv[4]);
  plotpoints = atoi (argv[5]);
  plotpoints = sqrt (plotpoints);

  gsl_matrix *hamiltonian = gsl_matrix_alloc (STATENO, STATENO);
  gsl_matrix *eigenbasis = gsl_matrix_alloc (STATENO, STATENO);
  gsl_vector *eigenvalues = gsl_vector_alloc (STATENO);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (STATENO);
  MPI_Init (&argc, &argv);

  /*Barrier synchrinization so that all processes continue at the same time  and time offsets are limited */
  MPI_Barrier (MPI_COMM_WORLD);

  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  begin = time (NULL);
  datasizes = (int *) malloc (nprocs * sizeof (int));
  offsets = (int *) malloc (nprocs * sizeof (int));
  /*Build the 2-particle tower of states */
  count = 0;
  for (i = 0; i < ALPHA; i++)
    {
      for (j = 0; j <= i; j++)
	{
	  towerofstates[count][0] = i;
	  towerofstates[count][1] = j;
	  count++;
	}
    }
  count = 0;
  u0 = U0;
  v0 = V0;

  param.v0 = v0;
  param.u0 = u0;
  param.s1 = s1;
  param.s2 = s2;
  for (i = 0; i < STATENO; i++)
    {
      param.tower[i][0] = towerofstates[i][0];
      param.tower[i][1] = towerofstates[i][1];
    }
  /*Formulate and Diagonalize the Unperturbed Hamiltonian here */
  for (i = 0; i < STATENO; i++)
    {
      for (j = 0; j < STATENO; j++)
	{
	  element = undriven_hamilt (&param, i, j);
	  gsl_matrix_set (hamiltonian, i, j, element);
	}
    }
  /*Diagonalize the hamiltonian */
  gsl_eigen_symmv (hamiltonian, eigenvalues, eigenbasis, w);
  gsl_eigen_symmv_free (w);

  /*Sort eigenvectors in order of increasing eigenvalues */
  gsl_eigen_symmv_sort (eigenvalues, eigenbasis, GSL_EIGEN_SORT_VAL_ASC);

  /*Put eigenbasis in an array & in struct param */
  for (i = 0; i < STATENO; i++)
    {
      param.eigenvalues[i] = gsl_vector_get (eigenvalues, i);
      for (j = 0; j < STATENO; j++)
	{
	  param.eigenvectors[i][j] = gsl_matrix_get (eigenbasis, i, j);
	}
    }
  /*Get chosen state from argv[] */
  k = statechosen;
  energy = gsl_vector_get (eigenvalues, k);
  quotient = plotpoints / nprocs;
  remainder = plotpoints % nprocs;
  rank = id + 1;
  if (rank <= remainder)
    {
      dist_size = quotient + 1;
    }
  else
    {
      dist_size = quotient;
    }
  /*parallelized in x only, not in p for simplicity */
  xinit = (xfinalglobal - xinitglobal) * (id) / nprocs;
  xinit = xinit + xinitglobal;

  xfinal = (xfinalglobal - xinitglobal) * (id + 1) / nprocs;
  xfinal = xfinal + xinitglobal;


  for (i = 0; i < dist_size; i++)
    {
      x = xinit + i * ((xfinal - xinit) / dist_size);
      radius =
	classical_pot (x, v0) + classical_pot (STROBE, v0) + classical_int (x,
									    STROBE,
									    u0);
      radius = energy - radius;
      if (radius >= 0.0)
	{
	  radius = sqrt (radius);
	  pinit = -radius;
	  pfinal = radius;
	  for (j = 0; j < plotpoints; j++)
	    {
	      p = pinit + j * ((pfinal - pinit) / plotpoints);

	      h = husimi (k, x, p, &param);
	      points[count] = x;
	      points[count + 1] = p;
	      points[count + 2] = gsl_complex_abs2 (h);
	      count = count + 3;
	    }
	}
    }

  /*Use collective communication to send the array size to root */
  MPI_Gather (&count, 1, MPI_INT, datasizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
  /*Then, build an array of offsets */
  for (i = 0; i < nprocs; i++)
    {
      offsets[i] = 0;
      if (datasizes[i] != 0)
	{
	  for (j = 1; j < i; j++)
	    offsets[i] = offsets[i] + datasizes[j];
	}
    }
  /*Now, use collective communication to gather  all the data  to root */
  /*All data is contiguous with size=datasize, root pid is 0 */
  MPI_Gatherv (points, count, MPI_DOUBLE, points, datasizes, offsets,
	       MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Finalize ();
  /*Total Data Size */
  for (i = 0; i < nprocs; i++)
    total = total + datasizes[i];
  end = time (NULL);
  if (id == 0)
    {
      for (i = 0; i < total; i = i + 3)
	fprintf (output, "%lf %lf %lf\n", points[i], points[i + 1],
		 points[i + 2]);
      fprintf (output, "\n");	/*One last linebreak for gri contours */
      printf ("\n#######################################################");
      printf ("\n#xrange =   %10.3lf - %10.3lf", xinitglobal, xfinalglobal);
      printf ("\n#Total time taken =     %ld ", end - begin);
      printf
	("\n#Program successfully completed. Output is in dat file. Plot it.");
      printf ("\n#######################################################\n");
      fflush (stdout);
    }
  fclose (output);
}
