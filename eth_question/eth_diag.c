/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Rigol Lattice: 
      Diagonalization of the lattice in
      Rigol et al, Nature 452, 854 (2008), doi:10.1038/nature06838
      for any given lattice size.
   Copyright (c) 2013 Analabha Roy (daneel@utexas.edu)

   This is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include "eth_diag.h"
#include "filenames.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char **argv)
{
  double begin, end;		//MPI timers 
  double point1, point2;	//More MPI timers
  Mat H, H_init;		// Hamiltonian matrix 
  Mat U_parallel;		//Unitary matrix of eigenvalues of H in row order
  EPS eps, eps_init;		// eigenproblem solver context 
  EPSType type;
  PetscReal error, tol, re, im;
  PetscScalar kr, ki;
  PetscScalar e_gnd_re, e_gnd_im;
  Vec evals, evals_err;
  Vec xr, xi;
  Vec fullstate_re, fullstate_im;
  PetscScalar *xrarray;		//Will point this to the vec xr 
  PetscInt i, dim, nev, maxit, its, nconv;
  PetscErrorCode ierr;
  PetscInt count;
  PetscBool willdraw = PETSC_FALSE;
  //Local indices and array of indices for the parallel eigenvector
  PetscInt loc_first, loc_last, loc_size;

  int mpirank;

  SlepcInitialize (&argc, &argv, (char *) 0, help);

  begin = MPI_Wtime ();

  /*********************************Input Block**********************************/


  ierr = PetscOptionsGetBool (NULL, "-draw_esys", &willdraw, NULL);
  CHKERRQ (ierr);

  /*******************************End Input Block********************************/


  MPI_Comm_rank (PETSC_COMM_WORLD, &mpirank);

  point1 = MPI_Wtime ();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Hx=kx
     Also get the root process to setup the sequential matrix U of the x'es
     This matrix needs to be row ordered.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*******************************Mandatory Input*********************************/
  PetscViewer likhbo;

  ierr = MatCreate (PETSC_COMM_WORLD, &H);
  CHKERRQ (ierr);
  ierr = MatSetType (H, MATMPIAIJ);
  CHKERRQ (ierr);

  ierr = MatCreate (PETSC_COMM_WORLD, &H_init);
  CHKERRQ (ierr);
  ierr = MatSetType (H_init, MATMPIAIJ);
  CHKERRQ (ierr);

  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 " Reading hamiltonian matrix in binary from %s ...",
		 hamilt_fnamesbin[0]);
  CHKERRQ (ierr);
  ierr =
    PetscViewerBinaryOpen (PETSC_COMM_WORLD, hamilt_fnamesbin[0],
			   FILE_MODE_READ, &likhbo);

  ierr = MatLoad (H, likhbo);
  CHKERRQ (ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
  CHKERRQ (ierr);

  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 " Reading lower hamiltonian matrix in binary from %s ...",
		 hamilt_fnamesbin[1]);
  CHKERRQ (ierr);
  ierr =
    PetscViewerBinaryOpen (PETSC_COMM_WORLD, hamilt_fnamesbin[1],
			   FILE_MODE_READ, &likhbo);

  ierr = MatLoad (H_init, likhbo);
  CHKERRQ (ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
  CHKERRQ (ierr);
  /*****************************End Mandatory Input*******************************/

  ierr = MatGetSize (H, &dim, NULL);
  CHKERRQ (ierr);

  //Set up evals vector
  ierr = MatGetVecs (H, NULL, &evals);
  CHKERRQ (ierr);
  ierr = MatGetVecs (H, NULL, &evals_err);
  CHKERRQ (ierr);

  //Set up evecs
  ierr = MatGetVecs (H, NULL, &xr);
  CHKERRQ (ierr);
  ierr = MatGetVecs (H, NULL, &xi);
  CHKERRQ (ierr);
  //Set up evec matrix
  ierr = MatDuplicate (H, MAT_DO_NOT_COPY_VALUES, &U_parallel);
  CHKERRQ (ierr);

  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "   ----------------- ------------------\n");
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
   */

  point1 = MPI_Wtime ();
  ierr = EPSCreate (PETSC_COMM_WORLD, &eps);
  CHKERRQ (ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
   */
  ierr = EPSSetOperators (eps, H, NULL);
  CHKERRQ (ierr);
  ierr = EPSSetProblemType (eps, EPS_HEP);
  CHKERRQ (ierr);
  ierr = EPSSetType (eps, PROBLEMTYPE);
  CHKERRQ (ierr);
  /*
     Set solver parameters at runtime
   */
  ierr = EPSSetFromOptions (eps);
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  nev = dim;			//Number of requested eigenvalues = all of them!

  ierr = EPSSolve (eps);
  CHKERRQ (ierr);
  /*
     Optional: Get some information from the solver and display it
   */
  point2 = MPI_Wtime ();
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "\n Diagonalized Hamiltonian in %f secs\n\n",
		 point2 - point1);
  CHKERRQ (ierr);
  ierr = EPSGetIterationNumber (eps, &its);
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 " Number of iterations of the method: %D\n", its);
  CHKERRQ (ierr);
  ierr = EPSGetType (eps, &type);
  CHKERRQ (ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, " Solution method: %s\n\n", type);
  CHKERRQ (ierr);
  ierr = EPSSetDimensions (eps, nev, PETSC_DECIDE, PETSC_DECIDE);
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n",
		 nev);
  CHKERRQ (ierr);
  ierr = EPSGetTolerances (eps, &tol, &maxit);
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 " Stopping condition: tol=%.4G, maxit=%D\n", tol, maxit);
  CHKERRQ (ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Get number of converged approximate eigenpairs
   */
  ierr = EPSGetConverged (eps, &nconv);
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD, " Number of converged eigenpairs: %D\n\n",
		 nconv);
  CHKERRQ (ierr);

  PetscInt *loc_idx = malloc (dim * sizeof (PetscInt));

  if (nconv == dim)
    {
      point1 = MPI_Wtime ();

      /*
         Evaluate eigenvalues, eigenvectors and relative errors
       */
      for (i = 0; i < nconv; i++)
	{
	  /*
	     Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
	     ki (imaginary part)
	   */
	  ierr = EPSGetEigenpair (eps, i, &kr, &ki, xr, xi);
	  CHKERRQ (ierr);
	  /*
	     Compute the relative error associated to each eigenpair
	   */
	  ierr = EPSComputeRelativeError (eps, i, &error);
	  CHKERRQ (ierr);

#if defined(PETSC_USE_COMPLEX)
	  re = PetscRealPart (kr);
	  im = PetscImaginaryPart (kr);
#else
	  re = kr;
	  im = ki;
#endif
	  if (im != 0.0)
	    {
	      ierr = VecSetValue (evals, i, kr, INSERT_VALUES);
	      CHKERRQ (ierr);
	    }
	  else
	    {
	      ierr = VecSetValue (evals, i, re, INSERT_VALUES);
	      CHKERRQ (ierr);
	    }
	  ierr = VecSetValue (evals_err, i, error, INSERT_VALUES);
	  CHKERRQ (ierr);

	  /*******************Build the Eigenvector Matrix********************/

	  /*
	     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	     xr is a distributed vector. Get the local components and ranges
	     and insert them rowwise into U_parallel
	     IGNORE xi FOR NOW. ASSUME IT IS ALL REAL 
	     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	   */

	  ierr = VecGetOwnershipRange (xr, &loc_first, &loc_last);
	  CHKERRQ (ierr);
	  loc_size = loc_last - loc_first;
	  ierr = VecGetArray (xr, &xrarray);
	  CHKERRQ (ierr);

	  for (count = 0; count < loc_size; count++)
	    loc_idx[count] = loc_first + count;

	  //Now, store these data in the ith row of the U matrix

	  ierr =
	    MatSetValues (U_parallel, 1, &i, loc_size, loc_idx, xrarray,
			  INSERT_VALUES);
	  CHKERRQ (ierr);

	  ierr = VecRestoreArray (xr, &xrarray);	//Undo VecGetArray calls
	  CHKERRQ (ierr);

	  /*******************Built the Eigenvector Matrix********************/

	}

      //Assemble the U matrix, the matrix of eigenvectors in root process 
      ierr = MatAssemblyBegin (U_parallel, MAT_FINAL_ASSEMBLY);
      CHKERRQ (ierr);
      ierr = MatAssemblyEnd (U_parallel, MAT_FINAL_ASSEMBLY);
      CHKERRQ (ierr);
      point2 = MPI_Wtime ();
      ierr =
	PetscPrintf (PETSC_COMM_WORLD,
		     "\n Gathered eigensystem in %lf secs\n",
		     point2 - point1);
      CHKERRQ (ierr);
    }
  else
    {
      ierr =
	PetscPrintf (PETSC_COMM_WORLD,
		     "Number of converged eigenpairs too small: %d\n", nconv);
      SETERRQ (PETSC_COMM_WORLD, ierr,
	       "ERROR: Number of converged eigenpairs too small:\n");
      CHKERRQ (ierr);
    }

  //Assemble the evals vector. Parallel eigenvalues and errors
  ierr = VecAssemblyBegin (evals);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (evals);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (evals_err);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (evals_err);
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   *     Now, diagonalize H_init and find the lowest eigenvalue and eigenvector
   *     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  point1 = MPI_Wtime ();
  //First, prepare the state to receive the eigenvector
  ierr = MatGetVecs (H_init, NULL, &fullstate_re);
  CHKERRQ (ierr);
  ierr = MatGetVecs (H_init, NULL, &fullstate_im);
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   *                Create the eigensolver and set various options
   *     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
   *     Create eigensolver context
   */
  ierr = EPSCreate (PETSC_COMM_WORLD, &eps_init);
  CHKERRQ (ierr);

  /*
   *     Set operators. In this case also, it is a standard eigenvalue problem
   */
  ierr = EPSSetOperators (eps_init, H_init, NULL);
  CHKERRQ (ierr);
  ierr = EPSSetProblemType (eps_init, EPS_HEP);
  CHKERRQ (ierr);
  ierr = EPSSetType (eps_init, PROBLEMTYPE);
  CHKERRQ (ierr);
  /*
   *     Set solver parameters at runtime
   */
  ierr = EPSSetFromOptions (eps_init);
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   *                      Solve the eigensystem
   *     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve (eps_init);
  CHKERRQ (ierr);
  point2 = MPI_Wtime ();

  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "\n Diagonalized lower part of Hamiltonian and \n set initial condition in %f secs\n\n",
		 point2 - point1);
  CHKERRQ (ierr);
  /*
   *     Optional: Get some information from the solver and display it
   */
  ierr = EPSGetIterationNumber (eps_init, &its);
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 " Number of iterations of the method: %D\n", its);
  CHKERRQ (ierr);
  ierr = EPSGetType (eps_init, &type);
  CHKERRQ (ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, " Solution method: %s\n\n", type);
  CHKERRQ (ierr);
  ierr = EPSGetDimensions (eps_init, &nev, NULL, NULL);
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n",
		 nev);
  CHKERRQ (ierr);
  ierr = EPSGetTolerances (eps_init, &tol, &maxit);
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 " Stopping condition: tol=%.4G, maxit=%D\n", tol, maxit);
  CHKERRQ (ierr);

  /*
     Get number of converged approximate eigenpairs                                                  *
   */

  ierr = EPSGetConverged (eps_init, &nconv);
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD, " Number of converged eigenpairs: %D\n\n",
		 nconv);
  CHKERRQ (ierr);
  if (nconv > 0)
    {

      /*
         Get the smallest eigenpair and equate the eigenvector to the full state at time t=0.
         According to slepc docs, the eigenpairs are sorted by decreasing modulus with positive first
         so second value is the highest negative value ie the lowest one   
       */
      ierr =
	EPSGetEigenpair (eps_init, 1, &e_gnd_re, &e_gnd_im,
			 fullstate_re, fullstate_im);
      CHKERRQ (ierr);

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart (e_gnd_re);
      im = PetscImaginaryPart (e_gnd_im);
#else
      re = e_gnd_re;
      im = e_gnd_im;
#endif
      ierr = PetscPrintf (PETSC_COMM_WORLD,
			  " Ground state energy: %F + I %F\n", re, im);
    }
  else
    {
      ierr =
	PetscPrintf (PETSC_COMM_WORLD,
		     "None of the eigenpairs converged: %d\n", nconv);
      SETERRQ (PETSC_COMM_WORLD, ierr,
	       "ERROR: Number of converged eigenpairs too small:\n");
      CHKERRQ (ierr);
    }
  //This is the initial condition in the diagonal basis.

  /***********************End Set initial condition*******************************/


  /**********************************End time*************************************/

  end = MPI_Wtime ();
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "   ----------------- ------------------\n");
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD, "\n Total runtime: %f\n", end - begin);
  CHKERRQ (ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, "\n");
  CHKERRQ (ierr);

  /********************************End time done***********************************/


  /*******************************Optional Output**********************************/

  /*Graphical representation of matrices and vectors */
  if (willdraw)
    {
      PetscViewer dekhbo;	//for viewing vectors and matrices
      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0, "Hamiltonian Matrix",
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     PETSC_DECIDE, &dekhbo);
      CHKERRQ (ierr);
      ierr = MatView (H, dekhbo);	//This prints the matrix to stdout
      CHKERRQ (ierr);
      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0, "Eigenvalues", PETSC_DECIDE,
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     &dekhbo);
      CHKERRQ (ierr);
      ierr = PetscViewerSetFormat (dekhbo, PETSC_VIEWER_DRAW_BASIC);
      CHKERRQ (ierr);
      ierr = VecView (evals, dekhbo);
      CHKERRQ (ierr);

      ierr =
	PetscPrintf (PETSC_COMM_WORLD,
		     "\nEigenspectrum of full Hamiltonian:");
      CHKERRQ (ierr);
      ierr = EPSPrintSolution (eps, NULL);
      CHKERRQ (ierr);

      ierr =
	PetscPrintf (PETSC_COMM_WORLD,
		     "\nEigenspectrum of lower part of Hamiltonian:");
      CHKERRQ (ierr);
      ierr = EPSPrintSolution (eps_init, NULL);
      CHKERRQ (ierr);
      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0, "Eigenvalue Errors",
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     PETSC_DECIDE, &dekhbo);
      CHKERRQ (ierr);
      ierr = VecView (evals_err, dekhbo);
      CHKERRQ (ierr);
      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0, "Eigenvector Matrix",
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     PETSC_DECIDE, &dekhbo);
      ierr = MatView (U_parallel, dekhbo);
      CHKERRQ (ierr);
      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0,
			     "Lower part of the Hamiltonian Matrix",
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     PETSC_DECIDE, &dekhbo);
      CHKERRQ (ierr);
      ierr = MatView (H_init, dekhbo);	//This prints the matrix to stdout
      CHKERRQ (ierr);

      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0,
			     "Ground state of lower part", PETSC_DECIDE,
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     &dekhbo);
      CHKERRQ (ierr);
      ierr = VecView (fullstate_re, dekhbo);
      CHKERRQ (ierr);

      ierr = PetscViewerDestroy (&dekhbo);
      CHKERRQ (ierr);
    }

  /*****************************End Optional Output********************************/


  /*******************************Mandatory Output*********************************/

  PetscViewer outview;
  {
    ierr =
      PetscPrintf (PETSC_COMM_WORLD,
		   "\n Writing Eigenvalues in bin to %s\n ...",
		   filenames_bin[0]);
    CHKERRQ (ierr);
    ierr =
      PetscViewerBinaryOpen (PETSC_COMM_WORLD, filenames_bin[0],
			     FILE_MODE_WRITE, &outview);
    CHKERRQ (ierr);
    ierr = VecView (evals, outview);
    CHKERRQ (ierr);
    ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
    CHKERRQ (ierr);
    ierr =
      PetscPrintf (PETSC_COMM_WORLD,
		   "\n Writing Eigenvectors in bin to %s\n ...",
		   filenames_bin[1]);
    CHKERRQ (ierr);


    ierr =
      PetscViewerBinaryOpen (PETSC_COMM_WORLD, filenames_bin[1],
			     FILE_MODE_WRITE, &outview);
    CHKERRQ (ierr);
    ierr = MatView (U_parallel, outview);
    CHKERRQ (ierr);

    ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
    CHKERRQ (ierr);
    ierr =
      PetscPrintf (PETSC_COMM_WORLD,
		   "\n Writing initial state in bin to %s\n ...",
		   filenames_bin[2]);
    CHKERRQ (ierr);
    ierr =
      PetscViewerBinaryOpen (PETSC_COMM_WORLD, filenames_bin[2],
			     FILE_MODE_WRITE, &outview);
    CHKERRQ (ierr);
    ierr = VecView (fullstate_re, outview);
    CHKERRQ (ierr);
    ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
    CHKERRQ (ierr);
    ierr =
      PetscPrintf (PETSC_COMM_WORLD,
		   "   ----------------- ------------------\n");
    CHKERRQ (ierr);
    ierr = PetscViewerDestroy (&outview);
    CHKERRQ (ierr);
  }


  /*****************************End Mandatory Output*******************************/


  /*******************************Free Work Space*********************************/

  ierr = EPSDestroy (&eps_init);
  CHKERRQ (ierr);
  ierr = EPSDestroy (&eps);
  CHKERRQ (ierr);
  ierr = MatDestroy (&H);
  CHKERRQ (ierr);
  ierr = MatDestroy (&H_init);
  CHKERRQ (ierr);
  ierr = MatDestroy (&U_parallel);
  CHKERRQ (ierr);
  ierr = VecDestroy (&fullstate_re);
  CHKERRQ (ierr);
  ierr = VecDestroy (&fullstate_im);
  CHKERRQ (ierr);
  ierr = VecDestroy (&xr);
  CHKERRQ (ierr);
  ierr = VecDestroy (&xi);
  CHKERRQ (ierr);
  ierr = VecDestroy (&evals);
  CHKERRQ (ierr);
  ierr = VecDestroy (&evals_err);
  CHKERRQ (ierr);
  ierr = SlepcFinalize ();

  free (loc_idx);


  /******************************Freed Work Space*********************************/


  return 0;
}
