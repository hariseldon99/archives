#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "params.h"
#include "benchmark.h"
#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char **argv)
{
  double local_time,global_time;
  double runtime_avg;
  drive_and_tower param;
  int i, j, n, m, argv_count = 1;
  int pid, nprocs;
  double cosmatrix;
  double y[2*DIMS], y_real[DIMS], y_imag[DIMS];
  double period;
  double kappa, lambda, omega;
  int parameter_begin, parameter_end, probsize,dim;
  FILE *benchmark;
  benchmark = fopen (argv[argv_count++], "w");

  MPI_Init (&argc, &argv);
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Comm_rank (MPI_COMM_WORLD, &pid);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  
  local_time = -gtod_timer();

 parameter_begin = atof (argv[argv_count++]);
 parameter_end = atof (argv[argv_count++]);
 kappa=atof (argv[argv_count++]);
 lambda = atof (argv[argv_count++]);
 omega = atof (argv[argv_count++]);

  period = 2.0 * M_PI / omega;
  param.kappa = kappa; 
  param.lambda = lambda;
  param.omega = omega;

  /*Calculate the cosine(x) Matrix elements in the Hamiltonian Representation 
   *And put them in struct param */
  for (n = 0; n < DIMS; n++)
    {
      param.psq[n] = psq_mat (n);
      for (m = 0; m <=n; m++)
	{
	  cosmatrix = c_mat (n, m);
	  param.cosmatrix[n + DIMS * m] = cosmatrix;
	  param.cosmatrix[m + DIMS * n] = cosmatrix;
	}
    }
    local_time+=gtod_timer();

  /*parameter  loop begins here */
  for (probsize = parameter_begin; probsize < parameter_end;probsize=probsize+5)
    {
	    local_time = -gtod_timer();
	    dim=2*probsize+1;
	    param.dim=dim;
	    if(pid==0) printf ("\nintegrating with problem size =  %d, dimensionality= %d\n", probsize,dim);
    
	    j = pid;
	/*This is the jth column of the matrix */
	/* Initial conditions, */
	for (i = 0; i < 2 * dim; i++)
	  {
	    if (i < dim)
	      {
		y[i] = kdel (i, j);
	      }
	    else
	      {
		y[i] = 0.0;
	      }
	  }

	/*DONE. Now evolve the jth column of the floquet matrix (stored in array y[j]) across time period */
	integrate (y, 0.0, period, &param);
	 /*DONE*/

	  /*Split into real and imaginary parts */
	  for (i = 0; i < dim; i++)
	  {
	    y_real[i] = y[i];
	    y_imag[i] = y[i + dim];
	  }
	  local_time+=gtod_timer();
	  MPI_Reduce (&local_time,&global_time,1, MPI_DOUBLE, MPI_SUM,0,MPI_COMM_WORLD);
	  if(pid==0){
		  runtime_avg=global_time/nprocs;
	 	  fprintf(benchmark,"\n %d %lf",probsize,runtime_avg);
	  	}
	if(pid==0) printf ("\nELAPSED TIME = %lf\n", local_time);
    }
  /*parameter loop ends here */
  MPI_Finalize(); 
  fclose (benchmark);

  return 0;
}
