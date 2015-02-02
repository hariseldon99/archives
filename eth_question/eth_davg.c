/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Eigenstate Thermalization Hypothesis (eth): 
      Computation of diagonal averages of the lattice in
      Rigol et al, Nature 452, 854 (2008), doi:10.1038/nature06838
      for any given lattice size. Lattice eigensystem obtained from rigol_lattice.c
   Copyright (c) 2013 Analabha Roy (daneel@utexas.edu)

   This is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include "eth_davg.h"
#include "filenames.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char **argv)
{
  double begin, end;		//MPI timers 
  double point1, point2;	//More MPI timers
  Mat U_parallel, U_seq;	//Unitary matrix of eigenvalues of H in column order
  Vec evals;
  PetscReal u = DEFAULTREP;
  Vec fullstate;
  PetscInt latsize = DEFAULTLATSIZE, vecsize = DEFAULTSIZE, Istart, Iend;
  PetscErrorCode ierr;
  PetscInt count, istride;
  PetscInt sitestride1, sitestride2;	//These will stride over all lattice sites
  PetscBool willdraw = PETSC_FALSE;
  int mpisize, mpirank;

  PetscInitialize (&argc, &argv, (char *) 0, help);

  begin = MPI_Wtime ();

  MPI_Comm_rank (PETSC_COMM_WORLD, &mpirank);

  MPI_Comm_size (PETSC_COMM_WORLD, &mpisize);

  point1 = MPI_Wtime ();
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

  ierr = PetscOptionsGetBool (NULL, "-draw_out", &willdraw, NULL);
  CHKERRQ (ierr);

  PetscInt matsize = latsize * latsize;
  //Get size and range of momentum space. If no input then use the default
  //values below
  PetscInt k_matsize = matsize;
  ierr = PetscOptionsGetInt (NULL, "-ksize", &k_matsize, NULL);
  CHKERRQ (ierr);
  PetscScalar kxinit = -M_PI;
  PetscScalar kxfinal = M_PI;
  PetscScalar kyinit = -M_PI;
  PetscScalar kyfinal = M_PI;

  ierr = PetscOptionsGetScalar (NULL, "-krange", &kxfinal, NULL);
  kxinit = -kxfinal;
  kyinit = kxinit;
  kyfinal = kxfinal;

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


  /**************************Generate the basis states***************************/

  gsl_combination *c = gsl_combination_calloc (latsize * latsize, vecsize);
  long **basis;
  long dim = 0;

  //Evaluate the dimensionality of the Hilbert space (# of allowed combinations)
  do
    {
      /*
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         c is the ith combination of N integers from M^2 integers
         Get all the non forbidden elements
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Now, generate the actual basis
     basis is an array of arrays that will contain the 
     selected basis vectors. Each basis vector is an array of N integers
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         c is the ith combination of N integers from M^2 integers
         Get all the non forbidden elements
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

  /*************************Generated the basis states***************************/


  /*******************************Mandatory Input*********************************/


  PetscViewer likhbo;

  ierr = VecCreate (PETSC_COMM_WORLD, &evals);
  CHKERRQ (ierr);

  ierr = MatCreate (PETSC_COMM_WORLD, &U_parallel);
  CHKERRQ (ierr);
  ierr = MatSetType (U_parallel, MATMPIAIJ);
  CHKERRQ (ierr);

  ierr = VecCreate (PETSC_COMM_WORLD, &fullstate);
  CHKERRQ (ierr);

  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "\n Reading evals in binary from %s ...", filenames_bin[0]);
  CHKERRQ (ierr);
  ierr =
    PetscViewerBinaryOpen (PETSC_COMM_WORLD, filenames_bin[0],
			   FILE_MODE_READ, &likhbo);
  ierr = VecLoad (evals, likhbo);
  CHKERRQ (ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 " Reading evec matrix in binary from %s ...",
		 filenames_bin[1]);
  CHKERRQ (ierr);
  ierr =
    PetscViewerBinaryOpen (PETSC_COMM_WORLD, filenames_bin[1],
			   FILE_MODE_READ, &likhbo);

  ierr = MatLoad (U_parallel, likhbo);
  CHKERRQ (ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 " Reading init state in binary from %s ...",
		 filenames_bin[2]);
  CHKERRQ (ierr);
  ierr =
    PetscViewerBinaryOpen (PETSC_COMM_WORLD, filenames_bin[2],
			   FILE_MODE_READ, &likhbo);
  ierr = VecLoad (fullstate, likhbo);
  CHKERRQ (ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
  CHKERRQ (ierr);
  ierr = PetscViewerDestroy (&likhbo);
  CHKERRQ (ierr);


  /*****************************End Mandatory Input*******************************/


  /******************************Diagonal Average*********************************/

  point1 = MPI_Wtime ();

  //First, rotate the initial state to diagonal basis
  Vec D;
  ierr = VecDuplicate (fullstate, &D);
  CHKERRQ (ierr);
  ierr = MatMult (U_parallel, fullstate, D);
  CHKERRQ (ierr);
  ierr = VecAssemblyBegin (D);
  CHKERRQ (ierr);
  ierr = VecAssemblyEnd (D);
  CHKERRQ (ierr);
  //These are the components D_\alpha = < \alpha | \psi(o) >
  //\alpha counts the eigenstates

  /*
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Now, square D, so it stores |D_\alpha|^2
     Note that D_\alpha is the |C_\alpha|^2 from the rigol paper
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   */

  ierr = VecPointwiseMult (D, D, D);
  CHKERRQ (ierr);

  //Now, scatter D to each processor since each i,j will need ALL values of D
  Vec D_seq;
  VecScatter ctx;		//Scatter context
  ierr = VecScatterCreateToAll (D, &ctx, &D_seq);
  CHKERRQ (ierr);
  ierr = VecScatterBegin (ctx, D, D_seq, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRQ (ierr);
  ierr = VecScatterEnd (ctx, D, D_seq, INSERT_VALUES, SCATTER_FORWARD);
  CHKERRQ (ierr);

  /*
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create a sequential matrix and copy U_parallel to it
     so that each processor has complete local copy of U
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   */
  ierr =
    MatGetRedundantMatrix (U_parallel, mpisize, PETSC_COMM_SELF, dim,
			   MAT_INITIAL_MATRIX, &U_seq);
  CHKERRQ (ierr);

  ierr = MatAssemblyBegin (U_seq, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (U_seq, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);


  ierr = PetscPrintf (PETSC_COMM_WORLD, "\n Computing diagonal average ...");
  CHKERRQ (ierr);

  /*
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Now, have to compute (b^\dagger_i b_j )_{\alpha\alpha} for each i,j on lattice site
     Note that | \alpha > is the \alpha-th evec or the \alpha-th row of U_parallel
     In canonical representation,  | \alpha > = \sum_l (U)_{\alpha l} | l >
     where | l > is the lth lexicographical combination, i.e. basis[l]
     So, for each i,j: 
     (b^\dagger_i b_j )_{\alpha\alpha} = \sum_{lm}
     U_{\alpha m} U_{\alpha l} can_sdot(<m|, bdaggerb(i,j,|l>,vecsize.latsize) )
     Code below does this for each \alpha and assemble a vector (b^\dagger_i b_j )_{\alpha\alpha} 
     Now, for each i,j, avg(b^\dagger_i b_j ) = 
     \sum_\alpha |Dseq_\alpha|^2 (b^\dagger_i b_j )_{\alpha\alpha}
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   */

  //Now, assemble a matrix for storing avg <b^\dagger_i b_j>

  Mat AVG_BDIBJ;

  ierr = MatCreate (PETSC_COMM_WORLD, &AVG_BDIBJ);
  CHKERRQ (ierr);
  ierr =
    MatSetSizes (AVG_BDIBJ, PETSC_DECIDE, PETSC_DECIDE, matsize, matsize);
  CHKERRQ (ierr);
  ierr = MatSetFromOptions (AVG_BDIBJ);
  CHKERRQ (ierr);
  ierr = MatSetUp (AVG_BDIBJ);
  CHKERRQ (ierr);

  ierr = MatGetOwnershipRange (AVG_BDIBJ, &Istart, &Iend);
  CHKERRQ (ierr);

  //This is just an array {0,1,2,3 ... matsize-1}
  PetscInt *idx;
  idx = malloc (matsize * sizeof (PetscInt));
  for (count = 0; count < matsize; count++)
    idx[count] = count;

  //This array stores the values of a particular row of the matrix
  PetscScalar *colvalues = malloc (matsize * sizeof (PetscScalar));

  PetscInt alpha, lambda, mu;
  //Create loop local basis vectors
  long *lambda_loc = malloc (vecsize * sizeof (long));
  long *mu_loc = malloc (vecsize * sizeof (long));
  PetscScalar U_al, U_am;
  size_t sitestride1_vec[1], sitestride2_vec[2];

  //This will store a particular row of U_seq
  const PetscScalar *u_row_arr;

  PetscScalar dasq, element, value;
  for (sitestride1 = Istart; sitestride1 < Iend; sitestride1++)
    {
      sitestride1_vec[0] = sitestride1;

      for (sitestride2 = 0; sitestride2 < matsize; sitestride2++)
	{
	  sitestride2_vec[0] = sitestride2;
	  //If either site is forbidden
	  if ((is_forbidden (sitestride1_vec, 1, latsize) == GSL_SUCCESS)
	      || (is_forbidden (sitestride2_vec, 1, latsize) == GSL_SUCCESS))
	    colvalues[sitestride2] = 0.0;
	  else
	    {
	      element = 0.0;
	      for (alpha = 0; alpha < dim; alpha++)
		{

		  //Get the alpha^th row of U_seq
		  ierr = MatGetRow (U_seq, alpha, NULL, NULL, &u_row_arr);
		  CHKERRQ (ierr);
		  //Get the alpha^th element of |D_\alpha|^2
		  ierr = VecGetValues (D_seq, 1, &alpha, &dasq);
		  CHKERRQ (ierr);
		  for (lambda = 0; lambda < dim; lambda++)
		    {
		      U_al = u_row_arr[lambda];

		      for (mu = 0; mu < dim; mu++)
			{
			  U_am = u_row_arr[mu];

			  //Create loop local copy of basis vectors
			  copyvec (basis[lambda], lambda_loc, vecsize);
			  copyvec (basis[mu], mu_loc, vecsize);
			  //bdaggerb on |l>
			  bdaggerb (sitestride1, sitestride2, lambda_loc,
				    vecsize, latsize);
			  //<m| on the result above
			  value = can_sdot (mu_loc, lambda_loc, vecsize);
			  element = element + value * U_al * U_am * dasq;
			}
		    }
		  ierr = MatRestoreRow (U_seq, alpha, NULL, NULL, &u_row_arr);
		  CHKERRQ (ierr);
		}
	      colvalues[sitestride2] = element;
	    }
	}

      //Insert the array of colvalues to the sitestride1^th row of H   
      ierr =
	MatSetValues (AVG_BDIBJ, 1, &sitestride1, matsize, idx, colvalues,
		      INSERT_VALUES);
      CHKERRQ (ierr);
    }

  ierr = MatAssemblyBegin (AVG_BDIBJ, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (AVG_BDIBJ, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);

  //Scatter this matrix to all processors

  Mat AVG_BDIBJ_seq;
  ierr =
    MatGetRedundantMatrix (AVG_BDIBJ, mpisize, PETSC_COMM_SELF, matsize,
			   MAT_INITIAL_MATRIX, &AVG_BDIBJ_seq);
  CHKERRQ (ierr);
  ierr = MatAssemblyBegin (AVG_BDIBJ_seq, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (AVG_BDIBJ_seq, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);

  /*
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Now, have to compute n(kx,ky)
     n_{l,m} = (1/latsize^2) \sum_{i,j} P(l,m)_{i,j} AVG_BDIBJ_{i,j}
     where P(l,m)_{i,j} = \exp{-I 2\pi (k_lm,k_m) \cdot ( (x_i,y_i) - (x_j,y_j) ) }
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   */

  //Brillouin Zone
  PetscScalar kx, ky;
  PetscInt xi, xj, yi, yj;
  PetscInt k_sitestride1, k_sitestride2;

  PetscScalar kxinc = (kxfinal - kxinit) / (k_matsize - 1);
  PetscScalar kyinc = (kyfinal - kyinit) / (k_matsize - 1);

  //This array stores the values of a particular row of the matrix
  PetscScalar *nk_colvalues = malloc (k_matsize * sizeof (PetscScalar));
  //This is just an array {0,1,2,3 ... k_matsize-1}
  PetscInt *k_idx;
  k_idx = malloc (k_matsize * sizeof (PetscInt));
  for (count = 0; count < k_matsize; count++)
    k_idx[count] = count;

  Mat Nk;
  PetscScalar Nk_elem, Phase_elem, bb_elem;
  ierr = MatCreate (PETSC_COMM_WORLD, &Nk);
  CHKERRQ (ierr);
  ierr = MatSetSizes (Nk, PETSC_DECIDE, PETSC_DECIDE, k_matsize, k_matsize);
  CHKERRQ (ierr);
  ierr = MatSetFromOptions (Nk);
  CHKERRQ (ierr);
  ierr = MatSetUp (Nk);
  CHKERRQ (ierr);
  PetscInt Nkrowstart, Nkrowend;
  ierr = MatGetOwnershipRange (Nk, &Nkrowstart, &Nkrowend);
  CHKERRQ (ierr);

  const PetscScalar *matrow;
  for (k_sitestride1 = Nkrowstart; k_sitestride1 < Nkrowend; k_sitestride1++)
    {
      kx = kxinit + k_sitestride1 * kxinc;
      for (k_sitestride2 = 0; k_sitestride2 < k_matsize; k_sitestride2++)
	{
	  ky = kyinit + k_sitestride2 * kyinc;

	  Nk_elem = 0.0;
	  for (sitestride1 = 0; sitestride1 < matsize; sitestride1++)
	    {
	      ierr =
		MatGetRow (AVG_BDIBJ_seq, sitestride1, NULL, NULL, &matrow);
	      CHKERRQ (ierr);
	      xi = sitestride1 / latsize;
	      yi = sitestride1 % latsize;
	      for (sitestride2 = 0; sitestride2 < matsize; sitestride2++)
		{
		  xj = sitestride2 / latsize;
		  yj = sitestride2 % latsize;
		  Phase_elem = kx * (xi - xj) + ky * (yi - yj);
		  //Complex math is from <complex.h>
		  Phase_elem = -I * Phase_elem * 2 * M_PI;
		  Phase_elem = cexp (Phase_elem);
		  bb_elem = matrow[sitestride2];
		  Nk_elem = Nk_elem + (Phase_elem * bb_elem);
		}
	      ierr =
		MatRestoreRow (AVG_BDIBJ_seq, sitestride1, NULL, NULL,
			       &matrow);
	      CHKERRQ (ierr);
	    }
	  //Evaluate the Nk for a particular k
	  nk_colvalues[k_sitestride2] = cabs (Nk_elem) / (latsize * latsize);
	}
      //Insert the array of colvalues to the k_sitestride1^th row of Nk   
      ierr =
	MatSetValues (Nk, 1, &k_sitestride1, k_matsize, k_idx, nk_colvalues,
		      INSERT_VALUES);
      CHKERRQ (ierr);
    }

  ierr = MatAssemblyBegin (Nk, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (Nk, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);

  point2 = MPI_Wtime ();
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "\n Computed diagonal average in %f secs\n",
		 point2 - point1);
  CHKERRQ (ierr);
  /**************************** End Diagonal Average*******************************/


  /*******************************Optional Output**********************************/

  /*Graphical representation of matrices and vectors */
  if (willdraw)
    {
      PetscViewer dekhbo;	//for viewing vectors and matrices

      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0,
			     "Initial State", PETSC_DECIDE,
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     &dekhbo);
      CHKERRQ (ierr);
      ierr = VecView (fullstate, dekhbo);
      CHKERRQ (ierr);

      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0,
			     "Eigenspace distribution",
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     PETSC_DECIDE, &dekhbo);
      CHKERRQ (ierr);
      ierr = VecView (D, dekhbo);
      CHKERRQ (ierr);

      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0, "Hopping Matrix Diag Avg",
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     PETSC_DECIDE, &dekhbo);
      ierr = MatView (AVG_BDIBJ, dekhbo);
      CHKERRQ (ierr);

      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0, "Momentum Distribution",
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
			     PETSC_DECIDE, &dekhbo);
      ierr = MatView (Nk, dekhbo);
      CHKERRQ (ierr);


      ierr = PetscViewerDestroy (&dekhbo);
      CHKERRQ (ierr);
    }

  /*****************************End Optional Output********************************/


  /******************************Mandatory Output**********************************/

  PetscViewer outview;		//for viewing vectors and matrices
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "\n Writing diagonal average hopping matrix in binary to %s ...",
		 filenames_eth_bin_davg[0]);
  CHKERRQ (ierr);
  ierr =
    PetscViewerBinaryOpen (PETSC_COMM_WORLD, filenames_eth_bin_davg[0],
			   FILE_MODE_WRITE, &outview);
  CHKERRQ (ierr);
  ierr = MatView (AVG_BDIBJ, outview);
  CHKERRQ (ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
  CHKERRQ (ierr);
  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 " Writing diagonal average momentum distribution in binary to %s ...",
		 filenames_eth_bin_davg[1]);
  CHKERRQ (ierr);
  ierr =
    PetscViewerBinaryOpen (PETSC_COMM_WORLD, filenames_eth_bin_davg[1],
			   FILE_MODE_WRITE, &outview);
  CHKERRQ (ierr);
  ierr = MatView (Nk, outview);
  CHKERRQ (ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
  CHKERRQ (ierr);

  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 " Writing initial state in rotated basis to %s ...",
		 filenames_eth_bin_davg[2]);
  CHKERRQ (ierr);

  ierr =
    PetscViewerBinaryOpen (PETSC_COMM_WORLD, filenames_eth_bin_davg[2],
			   FILE_MODE_WRITE, &outview);
  CHKERRQ (ierr);
  ierr = VecView (D, outview);
  CHKERRQ (ierr);

  ierr = PetscPrintf (PETSC_COMM_WORLD, " Done!\n");
  CHKERRQ (ierr);

  ierr = PetscViewerDestroy (&outview);
  CHKERRQ (ierr);

  /*****************************End Mandatory Output*******************************/


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


  /*******************************Free Work Space*********************************/
  ierr = VecScatterDestroy (&ctx);
  CHKERRQ (ierr);
  ierr = VecDestroy (&D);
  CHKERRQ (ierr);
  ierr = VecDestroy (&D_seq);
  CHKERRQ (ierr);
  ierr = VecDestroy (&evals);
  CHKERRQ (ierr);
  ierr = MatDestroy (&AVG_BDIBJ);
  CHKERRQ (ierr);
  ierr = MatDestroy (&AVG_BDIBJ_seq);
  CHKERRQ (ierr);
  ierr = MatDestroy (&Nk);
  CHKERRQ (ierr);
  ierr = MatDestroy (&U_parallel);
  CHKERRQ (ierr);
  ierr = MatDestroy (&U_seq);
  CHKERRQ (ierr);
  ierr = VecDestroy (&fullstate);
  CHKERRQ (ierr);

  free (idx);
  free (k_idx);
  free (basis);
  free (lambda_loc);
  free (mu_loc);
  free (colvalues);
  free (nk_colvalues);
  gsl_combination_free (c);

  /******************************Freed Work Space*********************************/

  PetscFinalize ();
  return 0;
}
