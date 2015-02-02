#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <slepceps.h>
#include "params.h"
#include "floquet_optical_lattice.h"
#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char **argv)
{
  time_t begin, end, now, rawtime;
  struct tm *timeinfo;
  paramspace_pt param;
  int i, j, n, m, err = 1, argv_count = 1;
  int m1, m2, n1, n2, count;
  double nf, ns, error;
  int nprocs, pid, ret, start = 1;
  int tagarray[STATENO];
  PetscErrorCode ierr;
  long errorcount = 0;
  double cosmatrix, sinmatrix, matelement;
  double y[2 * STATENO], y_real[STATENO], y_imag[STATENO];
  double quasi[STATENO], sortedenergies[STATENO];
  double period;
  double u0;
  double kappa, lambda0, wf, ws, omega, tfix;
  double parameter_begin, parameter_end, parameter_tot, parameter_inc;

  double floquetmatrix_real[STATENO * STATENO],
    floquetmatrix_imag[STATENO * STATENO];
  double floquetevals_real[STATENO], floquetevals_imag[STATENO];
  double floquetprev_real[STATENO], floquetprev_imag[STATENO];
  double floquetevecs_real[STATENO * STATENO],
    floquetevecs_imag[STATENO * STATENO];
  double floquetvecs_prev_real[STATENO * STATENO],
    floquetvecs_prev_imag[STATENO * STATENO];

  gsl_matrix *h = gsl_matrix_alloc (STATENO, STATENO);	//2 ptcl hamiltonian
  gsl_vector *energies = gsl_vector_alloc (STATENO);
  gsl_matrix *evecs = gsl_matrix_alloc (STATENO, STATENO);
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (STATENO);

  double eigenvectors[STATENO * STATENO], eigenvector_i[STATENO],
    nullvec[STATENO];

  int taga, tagb, tagc, tagd;
  double proba, probb, probc, probd;
  double taga_state_real[STATENO], taga_state_imag[STATENO];
  double tagb_state_real[STATENO], tagb_state_imag[STATENO];
  double tagc_state_real[STATENO], tagc_state_imag[STATENO];
  double tagd_state_real[STATENO], tagd_state_imag[STATENO];

  FILE *floqueteigenval, *energy;
  FILE *floquetstate_taga, *floquetstate_tagb, *floquetstate_tagc,
    *floquetstate_tagd;
  //FILE *floquetfile_real, *floquetfile_imag;

  energy = fopen (argv[argv_count++], "w");
  floqueteigenval = fopen (argv[argv_count++], "w");
  parameter_begin = atof (argv[argv_count++]);
  parameter_end = atof (argv[argv_count++]);
  parameter_tot=atof (argv[argv_count++]);
  parameter_inc = atof (argv[argv_count++]);
  lambda0 = atof (argv[argv_count++]);

  wf = atof (argv[argv_count++]);
  ws = atof (argv[argv_count++]);

  nf = atof (argv[argv_count++]);
  ns = atof (argv[argv_count++]);

  u0 = atof (argv[argv_count++]);
  kappa = atof (argv[argv_count++]);

  taga = atoi (argv[argv_count++]);
  tagb = atoi (argv[argv_count++]);
  tagc = atoi (argv[argv_count++]);
  tagd = atoi (argv[argv_count++]);

  floquetstate_taga = fopen (argv[argv_count++], "w");
  floquetstate_tagb = fopen (argv[argv_count++], "w");
  floquetstate_tagc = fopen (argv[argv_count++], "w");
  floquetstate_tagd = fopen (argv[argv_count++], "w");


  //floquetfile_real=fopen("result/floquetmatrix_real.dat","w");
  //floquetfile_imag=fopen("result/floquetmatrix_imag.dat","w");

  MPI_Init (&argc, &argv);
  MPI_Barrier (MPI_COMM_WORLD);
  MPI_Comm_rank (MPI_COMM_WORLD, &pid);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
  SlepcInitialize (&argc, &argv, (char *) 0, help);

  begin = time (NULL);

  error = (wf / ws) - (nf / ns);
  error = fabs (error);

  if (error >= RATTOL && pid == 0)
    {
      printf ("Error. Frequencies are not commensurate, exiting...");
      exit (1);
    }
  else if (pid == 0)
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
  tfix = 0.0;
  if (pid == 0)
    {

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
      fprintf (energy, "\n %lf ", tfix);
      for (n = 0; n < STATENO; n++)
	{
	  double energies_n = gsl_vector_get (energies, n);
	  sortedenergies[n] = energies_n;
	  fprintf (energy, "%lf ", energies_n);
	}
      /*Now, fold sortedenergies into quasienergies */
      err = fold_quasi (sortedenergies, omega, STATENO);
    }

  /*parameter  loop begins here */
  for (tfix = parameter_begin; tfix < parameter_end;
       tfix = tfix + parameter_inc)
    {
      param.tfix = tfix;

		/****************HERE IS WHERE THE PARALLEL FUN BEGINS********************************/
      {

	if (pid == 0)
	  {
	    now = time (NULL);
	    printf ("\nELAPSED TIME = %ld sec", now - begin);
	    printf ("\n\nEVOLVING FLOQUET MATRIX WITH TFIX =  %lf\n", tfix);
	  }

	/*The jth processor with pid of "j" initializes a row "y" with jth element=1 & all else 0 */
	j = pid;

	/*This is the jth column of the matrix */
	/* Initial conditions, */
	for (i = 0; i < 2 * STATENO; i++)
	  {
	    if (i < STATENO)
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
	  for (i = 0; i < STATENO; i++)
	  {
	    y_real[i] = y[i];
	    y_imag[i] = y[i + STATENO];
	  }

	/*DONE, Now aggregate the final values of floquet matrices (real and imaginary) for that array into
	 *the floquetmatrix arrays using MPI */

	MPI_Allgather (y_real, STATENO, MPI_DOUBLE, floquetmatrix_real,
		       STATENO, MPI_DOUBLE, MPI_COMM_WORLD);
	MPI_Allgather (y_imag, STATENO, MPI_DOUBLE, floquetmatrix_imag,
		       STATENO, MPI_DOUBLE, MPI_COMM_WORLD);
	/*End */
	/*All Processes must have floquetmatrix in order to use slepc */
      }
		/****************HERE IS WHERE THE PARALLEL FUN ENDS********************************/

      /*Check for the unitarity of the Floquet Matrix */
      if (pid == 0)
	{
	  err =
	    check_unitarity (floquetmatrix_real, floquetmatrix_imag, STATENO);
	  if (err != 0)
	    printf
	      ("***************ERROR FLOQUET MATRIX NOT UNITARY FOR KAPPA=%lf!!!! Has %d wrong elements***************",
	       kappa, err);
	}
      /*We want to diagonalize the floquet matrix now */

      ret =
	slepc_eval_nhep (argc, argv, floquetmatrix_real, floquetmatrix_imag,
			 STATENO, floquetevals_real, floquetevals_imag,
			 floquetevecs_real, floquetevecs_imag,
			 MPI_COMM_WORLD);
      /*ret=# of  converged eigenpairs if more than an aceptable minimum, else ret=0 */

      /*End Diagonalization */
      if (pid == 0)
	{
	  if (start == 0)
	    {
	      err =
		floquetsort (floquetevals_real, floquetevals_imag,
			     floquetevecs_real, floquetevecs_imag,
			     floquetprev_real, floquetprev_imag,
			     floquetvecs_prev_real, floquetvecs_prev_imag,
			     STATENO);
	      if (err == 0)
		printf
		  ("***************ERROR IN SORTING FLOQUET STATES!!!!***************");
	    }
	  err = arrcopy (floquetevals_real, floquetprev_real, STATENO);
	  err = arrcopy (floquetevals_imag, floquetprev_imag, STATENO);
	  err =
	    arrcopy (floquetevecs_real, floquetvecs_prev_real,
		     STATENO * STATENO);
	  err =
	    arrcopy (floquetevecs_imag, floquetvecs_prev_imag,
		     STATENO * STATENO);
	  if (err == 0)
	    printf
	      ("***************ERROR IN COPYING FLOQUET STATES!!!!***************");

	  if (ret != 0)
	    {
	      /*Calculate quasienergies */
	      errorcount =
		quasienergies (floquetevals_real, floquetevals_imag, STATENO,
			       period, quasi);
	      /*tag and print the quasienergies if this is the first iteration */
	      if (start == 1)
		{
		  err = tag_quasi (tagarray, quasi, sortedenergies, STATENO);
		  if (errorcount <= QUASIMAX && err == 1)
		    {
		      fprintf (floqueteigenval, "Param ");
		      for (i = 0; i < ret; i++)
			fprintf (floqueteigenval, "%10d ", tagarray[i]);
		    }
		}
	      if (errorcount <= QUASIMAX && err == 1)
		{
		  fprintf (floqueteigenval, "\n%lf", tfix);
		  fprintf (floquetstate_taga, "\n %lf", tfix);
		  fprintf (floquetstate_tagb, "\n %lf", tfix);
		  fprintf (floquetstate_tagc, "\n %lf", tfix);
		  fprintf (floquetstate_tagd, "\n %lf", tfix);

		  /*select the tagged floquet states */

		  err =
		    select_vec (floquetevecs_real, taga, taga_state_real,
				STATENO);
		  err =
		    select_vec (floquetevecs_imag, taga, taga_state_imag,
				STATENO);

		  err =
		    select_vec (floquetevecs_real, tagb, tagb_state_real,
				STATENO);
		  err =
		    select_vec (floquetevecs_imag, tagb, tagb_state_imag,
				STATENO);

		  err =
		    select_vec (floquetevecs_real, tagc, tagc_state_real,
				STATENO);
		  err =
		    select_vec (floquetevecs_imag, tagc, tagc_state_imag,
				STATENO);

		  err =
		    select_vec (floquetevecs_real, tagd, tagd_state_real,
				STATENO);
		  err =
		    select_vec (floquetevecs_imag, tagd, tagd_state_imag,
				STATENO);

		  /*Now, calculate dotproduct with ith eigenstate and dump it as p_i */
		  for (i = 0; i < STATENO; i++)
		    {
		      /*select the ith eigenstate */
		      err =
			select_vec (eigenvectors, i, eigenvector_i, STATENO);
		      /*compute the dotproducts */
		      proba =
			direct_prod (eigenvector_i, nullvec, taga_state_real,
				     taga_state_imag, STATENO);
		      probb =
			direct_prod (eigenvector_i, nullvec, tagb_state_real,
				     tagb_state_imag, STATENO);
		      probc =
			direct_prod (eigenvector_i, nullvec, tagc_state_real,
				     tagc_state_imag, STATENO);
		      probd =
			direct_prod (eigenvector_i, nullvec, tagd_state_real,
				     tagd_state_imag, STATENO);

		      fprintf (floquetstate_taga, " %10.6lf", proba);
		      fprintf (floquetstate_tagb, " %10.6lf", probb);
		      fprintf (floquetstate_tagc, " %10.6lf", probc);
		      fprintf (floquetstate_tagd, " %10.6lf", probd);
		    }


		  for (i = 0; i < ret; i++)
		    fprintf (floqueteigenval, " %10.6lf", quasi[i]);
		}
	    }

	}
      start = 0;
    }
  /*parameter loop ends here */

  end = time (NULL);

  if (pid == 0)
    final_output (period, begin, end, nprocs);

  /*This is not needed if slepc is called already coz slepc will do the finalizing */
  /*MPI_Finalize(); */
  gsl_eigen_symmv_free (w);
  gsl_matrix_free (evecs);
  gsl_matrix_free (h);
  gsl_vector_free (energies);
  ierr = SlepcFinalize ();
  CHKERRQ (ierr);
  fclose (floqueteigenval);
  fclose (energy);
  fclose (floquetstate_taga);
  fclose (floquetstate_tagb);
  fclose (floquetstate_tagc);
  fclose (floquetstate_tagd);

  //fclose(floquetfile_real);
  //fclose(floquetfile_imag);

  time (&rawtime);
  timeinfo = localtime (&rawtime);
  return 0;
}
