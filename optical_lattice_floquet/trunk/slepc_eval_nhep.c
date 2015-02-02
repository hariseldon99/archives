/*Don't forget to transclude this in the file which contains the function that executes the one below:
************************************************************************************************************
slepc_eval_nhep (int argc, char **argv, double matrixreal[], double matriximag[], int size, double eigenreal[], double eigenimag[], double ev_real[], double ev_im[]);
************************************************************************************************************
Notes:
1. Matrix A_[i,j] is passed as 2 1-d arrays A_real[i+n*j] and A_imag[i+n*j] (Real and imaginary parts of each element)
*/
//There's a memory leak in here somewhere...
#include <stdlib.h>
#include <slepceps.h>
#include <mpi.h>
#include "slepc_eval_nhep.h"

#undef __FUNCT__
#define __FUNCT__ "slepc_eval_nhep"
int
slepc_eval_nhep (int argc, char **argv, double matrix_real[],
		 double matrix_imag[], int size, double evals_real[],
		 double evals_imag[], double evecs_real[],
		 double evecs_imag[], MPI_Comm Comm)
{
  Mat A;			/* operator matrix */
  EPS eps;			/* eigenproblem solver context */
  //EPSType type;
  PetscReal error, re, im;
  PetscScalar kr, ki;
  PetscScalar *temp_global_real, *temp_global_imag;
  Vec xr, xi, xrseq, xiseq;
  VecScatter scatter_real, scatter_imag;
  IS toandfrom;
  int *idxfromto;
  PetscErrorCode ierr;
  PetscInt n = size, i, j, Istart, Iend;
  int its, nconv, ret;
  //PetscTruth FirstBlock = PETSC_FALSE, LastBlock = PETSC_FALSE;
  PetscScalar value, value_real, value_imag;

  ierr = PetscOptionsGetInt (PETSC_NULL, "-n", &n, PETSC_NULL);
  CHKERRQ (ierr);
  /*Get local pid */

  /*array containing the component index list for the eigenvector(s) */
  idxfromto = (int *) malloc (size * sizeof (int));
  for (i = 0; i < size; i++)
    idxfromto[i] = i;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate (Comm, &A);
  CHKERRQ (ierr);
  ierr = MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE, n, n);
  CHKERRQ (ierr);
  ierr = MatSetFromOptions (A);
  CHKERRQ (ierr);

  ierr = MatGetOwnershipRange (A, &Istart, &Iend);
  CHKERRQ (ierr);
  for (i = Istart; i < Iend; i++)
    {
      for (j = 0; j < n; j++)
	{
	  value_real = matrix_real[i + n * j];
	  value_imag = matrix_imag[i + n * j];
	  value = value_real + PETSC_i * value_imag;
	  ierr = MatSetValue (A, i, j, value, INSERT_VALUES);
	}
    }

  ierr = MatAssemblyBegin (A, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (A, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);

  ierr = MatGetVecs (A, PETSC_NULL, &xr);
  CHKERRQ (ierr);
  ierr = MatGetVecs (A, PETSC_NULL, &xi);
  CHKERRQ (ierr);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
   */
  ierr = EPSCreate (Comm, &eps);
  CHKERRQ (ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
   */
  ierr = EPSSetOperators (eps, A, PETSC_NULL);
  CHKERRQ (ierr);
  ierr = EPSSetProblemType (eps, EPS_NHEP);
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

  ierr = EPSSolve (eps);
  CHKERRQ (ierr);
  ierr = EPSGetIterationNumber (eps, &its);
  CHKERRQ (ierr);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Get number of converged approximate eigenpairs
   */
  ierr = EPSGetConverged (eps, &nconv);
  CHKERRQ (ierr);

  if (nconv == n)
    {
      /*
         Display eigenvalues and eigenmatrix
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

#ifdef PETSC_USE_COMPLEX
	  re = PetscRealPart (kr);
	  im = PetscImaginaryPart (kr);
#else
	  re = kr;
	  im = ki;
#endif
	  evals_real[i] = re;
	  evals_imag[i] = im;

	  /*xr and xi are distributed vectors. Need to scatter them into sequential vectors */
	  /*create the sequential eigenvectors */
	  VecCreateSeq (PETSC_COMM_SELF, size, &xrseq);
	  VecCreateSeq (PETSC_COMM_SELF, size, &xiseq);
	  /*create the index sets */
	  ISCreateGeneral (PETSC_COMM_SELF, size, idxfromto, &toandfrom);
	  /*create the scatter context */
	  VecScatterCreate (xr, toandfrom, xrseq, toandfrom, &scatter_real);
	  VecScatterCreate (xi, toandfrom, xiseq, toandfrom, &scatter_imag);
	  /*scatter the vectors */
	  VecScatterBegin (scatter_real, xr, xrseq, INSERT_VALUES,
			   SCATTER_FORWARD);
	  VecScatterEnd (scatter_real, xr, xrseq, INSERT_VALUES,
			 SCATTER_FORWARD);
	  VecScatterBegin (scatter_imag, xi, xiseq, INSERT_VALUES,
			   SCATTER_FORWARD);
	  VecScatterEnd (scatter_imag, xi, xiseq, INSERT_VALUES,
			 SCATTER_FORWARD);
	  /*Cleanup */
	  ISDestroy (toandfrom);
	  VecScatterDestroy (scatter_real);
	  VecScatterDestroy (scatter_imag);
	  /*Done */

	  VecGetArray (xrseq, &temp_global_real);
	  VecGetArray (xiseq, &temp_global_imag);

	  for (j = 0; j < nconv; j++)
	    {
	      /*jth row, ith column.Evecs are column ordered */
	      evecs_real[j + n * i] = temp_global_real[j];	/*xr{i} */
	      evecs_imag[j + n * i] = temp_global_imag[j];	/*xi{i} */
	    }

	  VecRestoreArray (xrseq, &temp_global_real);
	  VecRestoreArray (xiseq, &temp_global_imag);
	}

      ret = nconv;

    }
  else
    {
      ierr =
	PetscPrintf (Comm,
		     "ERROR: Number of converged eigenpairs too small: %d\n",
		     nconv);
      CHKERRQ (ierr);
      ret = 0;
    }

  /*
     Free work space
   */
  ierr = EPSDestroy (eps);
  CHKERRQ (ierr);
  ierr = MatDestroy (A);
  CHKERRQ (ierr);
  ierr = VecDestroy (xr);
  CHKERRQ (ierr);
  ierr = VecDestroy (xi);
  CHKERRQ (ierr);
  ierr = VecDestroy (xrseq);
  CHKERRQ (ierr);
  ierr = VecDestroy (xiseq);
  CHKERRQ (ierr);


  free (idxfromto);

  return ret;
}
