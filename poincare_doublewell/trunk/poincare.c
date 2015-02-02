/*This program uses gsl to compute the poincare section of 2 particles in a Double Well*/
/*and parallelizes it using MPI*/

#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include "mpi.h"
#include<gsl/gsl_errno.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_odeiv.h>
#include "params.h"

double classical_int (double x1, double x2, double u0);
double classical_pot (double x, double v0);
double force (double x1, double x2, double v0, double u0);
int func (double t, const double y[], double f[], void *params);
int jac (double t, const double y[], double *dfdy, double dfdt[],
	 void *params);
int integrate (double *y, double initial, double final, double epsilon,
	       double *out);

int
main (int argc, char *argv[])
{
  double x1init, p1init, x2init, p2init;
  double y[4], temp, energy;
  double v0 = V0, u0 = U0;
  double x1initglobal, p1initglobal, x1finalglobal, p1finalglobal;	/*Range of initial conditions input from file */
  double points[BUFSIZE];	/*This is where all the data goes */
  int id, p, numprocs;
  int datasize, *datasizes, *offsets;
  int count = 0, energyconservation = 1, trajno;
  int i, j;
  double walltime;
  long total = 0;
  double period = TIME, epsilon = 0.0;
  FILE *output;
  output = fopen (argv[1], "w");	/*Output data file */

  MPI_Init (&argc, &argv);

  /*Barrier synchrinization so that all processes continue at the same time  and time offsets are limited */
  MPI_Barrier (MPI_COMM_WORLD);

  /*Start walltime counter */
  walltime = -MPI_Wtime ();

  MPI_Comm_rank (MPI_COMM_WORLD, &id);
  MPI_Comm_size (MPI_COMM_WORLD, &p);
  /*assign global range of initial conditions */
  x1initglobal = atof (argv[2]);
  p1initglobal = atof (argv[3]);
  x1finalglobal = atof (argv[4]);
  p1finalglobal = atof (argv[5]);
  energy = atof (argv[6]);
  numprocs = p;
  datasizes = (int *) malloc (numprocs * sizeof (int));
  offsets = (int *) malloc (numprocs * sizeof (int));
  /*Set Initial Conditions according to the conservation of energy */
  x1init = (x1finalglobal - x1initglobal) * (id + 1) / p;
  x1init = x1init + x1initglobal;
  p1init = (p1finalglobal - p1initglobal) * (id + 1) / p;
  p1init = p1init + p1initglobal;
  x2init = 1.0;
  temp =
    p1init * p1init + classical_pot (x1init, v0) + classical_pot (x2init,
								  v0) +
    classical_int (x1init, x2init, u0);
  temp = energy - temp;
  if (temp < 0.0)
    {
      /*Whoops!Energy is not conserved. Nothing to see here, move aint... */
      energyconservation = 0;
    }
  else
    {
      p2init = sqrt (temp);
      y[0] = x1init;
      y[1] = x2init;
      y[2] = p1init;
      y[3] = p2init;

      /*Do the integration, store it in points and return the number of poitns in count */
      count = integrate (y, 0.0, period, epsilon, points);
    }

  /*OK, Now you're done calculating the trajectories. */
  /*Data size=count*energyconservation bcoz no data means that energy is not conserved */
  datasize = energyconservation * count;

  /*Use collective communication to send the array size to root */
  MPI_Gather (&datasize, 1, MPI_INT, datasizes, 1, MPI_INT, 0,
	      MPI_COMM_WORLD);

  /*Then, build an array of offsets */
  for (i = 0; i < numprocs; i++)
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
  MPI_Gatherv (points, datasize, MPI_DOUBLE, points, datasizes, offsets,
	       MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /*Calculates the number of trajectories by sum reduction */
  MPI_Reduce (&energyconservation, &trajno, 1, MPI_INT, MPI_SUM, 0,
	      MPI_COMM_WORLD);

  walltime += MPI_Wtime ();
  MPI_Finalize ();

  /*Total Data Size */
  for (i = 0; i < numprocs; i++)
    total = total + datasizes[i];

  /*Dump points to output file */
  if (id == 0)
    {
      fprintf (output, "%lf %lf \n", x1init, p1init);
      for (i = 0; i < total; i = i + 2)
	fprintf (output, "%lf %lf\n", points[i], points[i + 1]);

      printf ("\n#######################################################");
      printf ("\n#xrange =   %10.3lf - %10.3lf", x1initglobal, x1finalglobal);
      printf ("\n#prange =     %10.3lf - %10.3lf", p1initglobal,
	      p1finalglobal);
      printf ("\n#energy =     %10.3lf", energy);
      printf ("\n#Number of trajectories =  %d", trajno);
      printf ("\n#Total time taken =     %lf sec ", walltime);
      printf
	("\n#Program successfully completed. Output is in dat file. Plot it.");
      printf ("\n#######################################################\n");
      fflush (stdout);
      fclose (output);
      return 0;
    }
}
