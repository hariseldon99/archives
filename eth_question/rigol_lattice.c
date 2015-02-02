/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Rigol Lattice: 
      Formulation of the lattice Hamiltonian in
      Rigol et al, Nature 452, 854 (2008), doi:10.1038/nature06838
      for any given lattice size.
   Copyright (c) 2013 Analabha Roy (daneel@utexas.edu)

   This is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include "rigol_lattice.h"
#include "filenames.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char **argv)
{
  double begin, end;		//MPI timers 
  double point1, point2;	//More MPI timers
  Mat H, H_init;		// Hamiltonian matrix 
  PetscScalar value;
  PetscReal u = DEFAULTREP;
  PetscInt latsize = DEFAULTLATSIZE, vecsize = DEFAULTSIZE, Istart, Iend;
  PetscErrorCode ierr;
  PetscInt count, istride;
  PetscInt alpha, beta;		//These will stride over alpha and beta ie matrix elements
  PetscInt sitestride1, sitestride2;	//These will stride over all lattice sites
  PetscBool willdraw = PETSC_FALSE;
  //Local indices and array of indices for the parallel eigenvector

  int mpirank;

  PetscInitialize (&argc, &argv, (char *) 0, help);

  begin = MPI_Wtime ();

  /*********************************Input Block**********************************/

  ierr = PetscOptionsGetInt (NULL, "-lattice_size", &latsize, NULL);
  ierr = PetscOptionsGetInt (NULL, "-vector_size", &vecsize, NULL);
  if (latsize % 2 == 0)
    SETERRQ (PETSC_COMM_WORLD, ierr, "Lattice size must be odd\n");
  CHKERRQ (ierr);
  if (vecsize >= (latsize * latsize))
    SETERRQ (PETSC_COMM_WORLD, ierr,
	     "vector size (# of hard-core bosons) can't be more than the # of sites\n");
  CHKERRQ (ierr);

  ierr = PetscOptionsGetReal (NULL, "-repulsion", &u, NULL);
  CHKERRQ (ierr);

  ierr = PetscOptionsGetBool (NULL, "-draw_mats", &willdraw, NULL);
  CHKERRQ (ierr);

  ierr = PetscPrintf (PETSC_COMM_WORLD, "\n");
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD, " Lattice size = %D X %D\n", latsize,
		 latsize);
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD, " Number of hard sphere bosons = %D\n",
		 vecsize);
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD, " Nearest-Neighbor repulsion = %F\n", u);
  CHKERRQ (ierr);

  /*******************************End Input Block********************************/


  MPI_Comm_rank (PETSC_COMM_WORLD, &mpirank);

  point1 = MPI_Wtime ();
  /***************************Generate the basis states***************************/

  gsl_combination *c = gsl_combination_calloc (latsize * latsize, vecsize);
  long **basis;
  long dim = 0;

  //Evaluate the dimensionality of the Hilbert space (# of allowed combinations)
  do
    {
      /*
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         c is the ith combination of N integers from M^2 integers
         Get all the non forbidden elements
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       */
      if (is_forbidden (c->data, vecsize, latsize) == GSL_FAILURE)
	{			//If the site is NOT forbidden, append dim
	  dim++;
	}
    }
  while (gsl_combination_next (c) == GSL_SUCCESS);

  ierr =
    PetscPrintf (PETSC_COMM_WORLD, " Hilbert Space Dimensions = %D\n\n", dim);
  CHKERRQ (ierr);

  /*
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Now, generate the actual basis
     basis is an array of arrays that will contain the 
     selected basis vectors. Each basis vector is an array of N integers
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   */
  basis = (long **) malloc (dim * sizeof (long *));
  for (count = 0; count < dim; count++)
    basis[count] = (long *) malloc (vecsize * sizeof (long));

  //Reset the combination to the one that is lexicographically first
  gsl_combination_init_first (c);
  count = 0;
  //Now, append all allowed combinations to the basis
  do
    {
      /*
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         c is the ith combination of N integers from M^2 integers
         Get all the non forbidden elements
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       */
      if (is_forbidden (c->data, vecsize, latsize) == GSL_FAILURE)
	{			//If the site is NOT forbidden, append data in c to basis
	  for (istride = 0; istride < vecsize; istride++)
	    basis[count][istride] = c->data[istride];
	  count++;
	}
    }
  while (gsl_combination_next (c) == GSL_SUCCESS);

  point2 = MPI_Wtime ();
  ierr =
    PetscPrintf (PETSC_COMM_WORLD, "\n Assembled Basis in %f secs\n",
		 point2 - point1);
  CHKERRQ (ierr);
  /**************************Generated the basis states***************************/


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Hx=kx
     Also get the root process to setup the sequential matrix U of the x'es
     This matrix needs to be row ordered.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  point1 = MPI_Wtime ();

  //This is just an array {0,1,2,3 ... dim-1}
  PetscInt *idx;
  idx = malloc (dim * sizeof (PetscInt));
  for (count = 0; count < dim; count++)
    idx[count] = count;

  ierr = MatCreate (PETSC_COMM_WORLD, &H);
  CHKERRQ (ierr);
  ierr = MatSetSizes (H, PETSC_DECIDE, PETSC_DECIDE, dim, dim);
  CHKERRQ (ierr);
  ierr = MatSetFromOptions (H);
  CHKERRQ (ierr);
  ierr = MatSetUp (H);
  CHKERRQ (ierr);

  //Set all elements to zero by default
  ierr = MatZeroEntries (H);

  ierr = MatCreate (PETSC_COMM_WORLD, &H_init);
  CHKERRQ (ierr);
  ierr = MatSetSizes (H_init, PETSC_DECIDE, PETSC_DECIDE, dim, dim);
  CHKERRQ (ierr);
  ierr = MatSetFromOptions (H_init);
  CHKERRQ (ierr);
  ierr = MatSetUp (H_init);
  CHKERRQ (ierr);

  ierr = MatGetOwnershipRange (H, &Istart, &Iend);
  CHKERRQ (ierr);
  //Create loop local basis vectors
  long *alphaloc = malloc (vecsize * sizeof (long));
  long *betaloc = malloc (vecsize * sizeof (long));
  //This array stores the values of a particular row of the Hamiltonian
  PetscScalar *colvalues = malloc (dim * sizeof (PetscScalar));

  //Build the matrix elements of the Hamiltonian H
  for (alpha = Istart; alpha < Iend; alpha++)
    {
      for (beta = 0; beta < (alpha + 1); beta++)
	{
	  value = 0.0;
	  for (sitestride1 = 0; sitestride1 < latsize * latsize;
	       sitestride1++)
	    for (sitestride2 = 0; sitestride2 < latsize * latsize;
		 sitestride2++)
	      {
		if (are_nn (sitestride1, sitestride2, latsize) == GSL_SUCCESS)
		  {
		    //Create loop local copies of basis vectors
		    copyvec (basis[alpha], alphaloc, vecsize);
		    copyvec (basis[beta], betaloc, vecsize);
		    //Kinetic Energy
		    bdaggerb (sitestride1, sitestride2, betaloc, vecsize,
			      latsize);
		    value = value - can_sdot (alphaloc, betaloc, vecsize);

		    //Create loop local copies of basis vectors
		    copyvec (basis[alpha], alphaloc, vecsize);
		    copyvec (basis[beta], betaloc, vecsize);
		    //Potential energy
		    nsquared (sitestride1, sitestride2, betaloc, vecsize);
		    value =
		      value + (u / 2.0) * can_sdot (alphaloc, betaloc,
						    vecsize);
		  }
	      }
	  colvalues[beta] = value;
	  CHKERRQ (ierr);
	}
      //Insert the array of 'beta' values in 'colvalues' to the alpha^th row of H   
      ierr =
	MatSetValues (H, 1, &alpha, (alpha + 1), idx, colvalues,
		      INSERT_VALUES);
      CHKERRQ (ierr);
    }

  ierr = MatAssemblyBegin (H, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (H, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - 
   * Now, compute the hermitian part of H = (H + H^\dagger)/2
   * THIS IS FASTER THAN SETTING THE UPPER TRIANGULAR PART MANUALLY
   * - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - 
   */

  Mat Htrans;
  ierr = MatTranspose (H, MAT_INITIAL_MATRIX, &Htrans);
  CHKERRQ (ierr);
  ierr = MatAXPY (H, 1.0, Htrans, DIFFERENT_NONZERO_PATTERN);
  CHKERRQ (ierr);
  ierr = MatScale (H, 0.5);
  CHKERRQ (ierr);
  ierr = MatDestroy (&Htrans);
  CHKERRQ (ierr);

  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "   ----------------- ------------------\n");
  CHKERRQ (ierr);
  point2 = MPI_Wtime ();
  ierr =
    PetscPrintf (PETSC_COMM_WORLD, "\n Assembled Hamiltonian in %f secs",
		 point2 - point1);
  CHKERRQ (ierr);



  /*************************Set initial condition*********************************/
  /*
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Currently sets initial condition to the one from the Rigol paper
     i.e. the ground state of N bosons in the lower right quarter of the lattice  
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   */
  point1 = MPI_Wtime ();
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "\n   ----------------- ------------------\n");
  //First Copy H to H_init;
  ierr = MatCopy (H, H_init, DIFFERENT_NONZERO_PATTERN);
  CHKERRQ (ierr);
  ierr = MatAssemblyBegin (H_init, MAT_FLUSH_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (H_init, MAT_FLUSH_ASSEMBLY);
  CHKERRQ (ierr);
  point2 = MPI_Wtime ();

  ierr =
    PetscPrintf (PETSC_COMM_WORLD, "\n Copied H to H_init in %lf secs",
		 point2 - point1);
  CHKERRQ (ierr);
  /*
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Now, if the alpha, beta matrix element of H_init corresponds to a state outside
     the lower right quarter of the lattice, then set that to a unit matrix elements
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   */
  point1 = MPI_Wtime ();
  ierr = MatGetOwnershipRange (H_init, &Istart, &Iend);
  CHKERRQ (ierr);
  for (alpha = Istart; alpha < Iend; alpha++)
    for (beta = 0; beta < dim; beta++)
      {
	//if either the \alpha or \beta state is disallowed, the submatrix is identity
	if ((is_outside_lower_right_quarter (basis[alpha], vecsize, latsize)
	     == GSL_SUCCESS)
	    || (is_outside_lower_right_quarter (basis[beta], vecsize, latsize)
		== GSL_SUCCESS))
	  {
	    if (alpha == beta)
	      {
		value = 1.0;
		//Change the value of (H_init)_{\alpha\beta} to 'value'
		ierr =
		  MatSetValue (H_init, alpha, beta, value, INSERT_VALUES);
		CHKERRQ (ierr);
	      }
	    else
	      {
		value = 0.0;
		//Change the value of (H_init)_{\alpha\beta} to 'value'
		ierr =
		  MatSetValue (H_init, alpha, beta, value, INSERT_VALUES);
		CHKERRQ (ierr);
	      }
	  }
	//If both \alpha and \beta states are allowed, leave the Hamiltonian invariant  
      }

  //Re-assemble H_init
  ierr = MatAssemblyBegin (H_init, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (H_init, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  point2 = MPI_Wtime ();
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "\n Reset H_init to lower part in %lf secs",
		 point2 - point1);
  CHKERRQ (ierr);

  //This is the initial condition in the diagonal basis.

  /***********************End Set initial condition*******************************/


  /**********************************End time*************************************/

  end = MPI_Wtime ();
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "\n   ----------------- ------------------\n");
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
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0,
			     "Lower part of the Hamiltonian Matrix",
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     PETSC_DECIDE, &dekhbo);
      CHKERRQ (ierr);
      ierr = MatView (H_init, dekhbo);	//This prints the matrix to stdout
      CHKERRQ (ierr);

      ierr = PetscViewerDestroy (&dekhbo);
      CHKERRQ (ierr);
    }

  /*****************************End Optional Output********************************/


  /*******************************Mandatory Output*********************************/

  PetscViewer likhbo;
  {
    ierr =
      PetscPrintf (PETSC_COMM_WORLD,
		   "\n Writing full Hamiltonian in bin to %s\n ...",
		   hamilt_fnamesbin[0]);
    CHKERRQ (ierr);

    ierr =
      PetscViewerBinaryOpen (PETSC_COMM_WORLD, hamilt_fnamesbin[0],
			     FILE_MODE_WRITE, &likhbo);
    CHKERRQ (ierr);
    ierr = MatView (H, likhbo);
    CHKERRQ (ierr);

    ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
    CHKERRQ (ierr);

    ierr =
      PetscPrintf (PETSC_COMM_WORLD,
		   "\n Writing lower part Hamiltonian in bin to %s\n ...",
		   hamilt_fnamesbin[1]);
    CHKERRQ (ierr);

    ierr =
      PetscViewerBinaryOpen (PETSC_COMM_WORLD, hamilt_fnamesbin[1],
			     FILE_MODE_WRITE, &likhbo);
    CHKERRQ (ierr);
    ierr = MatView (H_init, likhbo);
    CHKERRQ (ierr);

    ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
    CHKERRQ (ierr);


    ierr =
      PetscPrintf (PETSC_COMM_WORLD,
		   "   ----------------- ------------------\n");
    CHKERRQ (ierr);
    ierr = PetscViewerDestroy (&likhbo);
    CHKERRQ (ierr);
  }


  /*****************************End Mandatory Output*******************************/


  /*******************************Free Work Space*********************************/

  ierr = MatDestroy (&H);
  CHKERRQ (ierr);
  ierr = MatDestroy (&H_init);
  CHKERRQ (ierr);
  ierr = PetscFinalize ();

  free (idx);
  free (colvalues);
  free (alphaloc);
  free (betaloc);
  free (basis);

  gsl_combination_free (c);

  /******************************Freed Work Space*********************************/


  return 0;
}
