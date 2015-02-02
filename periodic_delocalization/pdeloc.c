/*TODO: Add periodic measurement*/
static char help[] =
  "Solves a simple time-dependent Schroedinger equation for a TB model with periodic drive.\n\
Apart from PetSc's routine options, input parameters include:\n\
  -m <pnts>,   where <pnts>   = lattice size\n\
  -a <ampl>,   where <ampl>   = drive amplitude\n\
  -t <freq>,   where <freq>   = measurement frequency\n\
  -w <freq>,   where <freq>   = drive frequency\n\
  -debug              : Activate debugging printouts\n\n";

/*
   Default problem parameters. These are overwritten by argument vectors
*/
#define DEFAULT_SIZE 100
#define DEFAULT_T 3.0
#define DEFAULT_AMPL 3.0
#define DEFAULT_W 3.0
#define TIME_TOTAL_MAX 100.0	/* default max total time */
#define TIME_STEPS_MAX 1E5	/* default max timesteps */
/*
   Default problem parameters. These are always hard coded and never 
   overwritten
*/
#define INIT_STEP 1E-4
#define INIT_TIME 0.0
#define ABSERROR 1E-9
#define RELERROR 1E-9

/* ------------------------------------------------------------------------

   This program solves the one-dimensional Schroedinger equation in a 
   tightly bound lattice with a periodic time dependence
       \partial_t u_m(t) = \sum_m J_{lm}(t) u_m(t),
   on the domain 0 <= m < L, with the boundary conditions
       u_m(0) = 0 for all m except m=[L/2], where u_m(0) = 1
   This is a set of linear, first-order ODEs
   The time evolution using the various TS methods can be run by
   running the program via
       mpiexec -n <procs> ex3 -ts_type <timestepping solver>
  ------------------------------------------------------------------------- */

/*
   Include "petscdmda.h" so that we can use distributed arrays (DMDAs) to manage
   the parallel grid.  Include "petscts.h" so that we can use TS solvers.
   Note that this file automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h  - vectors
     petscmat.h  - matrices
     petscis.h     - index sets            petscksp.h  - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h   - preconditioners
     petscksp.h   - linear solvers        petscsnes.h - nonlinear solvers
*/


#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscdraw.h>
#include <petscksp.h>

/*
   User-defined application context - contains data needed by the
   application-provided call-back routines.
*/
typedef struct
{
  MPI_Comm comm;		/* communicator */
  DM da;			/* distributed array data structure */
  Vec localwork;		/* local ghosted work vector */
  Vec u_local;			/* local ghosted approximate solution vector */
  PetscInt m;			/* total number of grid points */
  PetscBool debug;		/* flag (1 indicates activation of debugging printouts) */
  PetscViewer viewer1;		/* viewer for the solution  */
  PetscReal norm_2;		/* wavefunction norm */
  PetscReal measurement_freq;	/* measurement frequency */
  PetscReal h0;			/* drive amplitude */
  PetscReal w;			/* drive frequency */
  PetscMPIInt rank;
  PetscMPIInt size;
} AppCtx;

/*
   User-defined routines
*/
extern PetscInt kdel (PetscInt, PetscInt);
extern PetscErrorCode InitialConditionsLocalized (Vec, AppCtx *);
extern PetscErrorCode RHSMatrixSchrodinger (TS, PetscReal, Vec, Mat, Mat,
					    void *);
extern PetscErrorCode Monitor (TS, PetscInt, PetscReal, Vec, void *);

#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char **argv)
{
  AppCtx appctx;		/* user-defined application context */
  TS ts;			/* timestepping context */
  Mat A;			/* matrix data structure */
  Vec u;			/* approximate solution vector */
  PetscReal time_total_max = TIME_TOTAL_MAX;	/* default max total time */
  PetscInt time_steps_max = TIME_STEPS_MAX;	/* default max timesteps */
  PetscDraw draw;		/* drawing context */
  PetscErrorCode ierr;
  PetscInt steps, m;
  PetscMPIInt size, rank;
  PetscReal dt, ftime;
  PetscReal measurement_freq;	/* Measurement frequency */
  PetscReal h0;			/* Drive amplitude */
  PetscReal w;			/* drive frequency */
  PetscBool flg;

  double init_runtime, final_runtime;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program and set problem parameters
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = PetscInitialize (&argc, &argv, (char *) 0, help);
  CHKERRQ (ierr);

  init_runtime = MPI_Wtime ();

  appctx.comm = PETSC_COMM_WORLD;

  m = DEFAULT_SIZE;
  ierr = PetscOptionsGetInt (NULL, "-m", &m, NULL);
  CHKERRQ (ierr);

  h0 = DEFAULT_AMPL;
  ierr = PetscOptionsGetReal (NULL, "-a", &h0, NULL);
  CHKERRQ (ierr);

  measurement_freq = DEFAULT_T;
  ierr = PetscOptionsGetReal (NULL, "-t", &measurement_freq, NULL);
  CHKERRQ (ierr);

  w = DEFAULT_W;
  ierr = PetscOptionsGetReal (NULL, "-w", &w, NULL);
  CHKERRQ (ierr);

  ierr = PetscOptionsHasName (NULL, "-debug", &appctx.debug);
  CHKERRQ (ierr);

  /*If this flag is true, then do not evaluate anything and just dump out the help */
  flg = PETSC_FALSE;
  ierr = PetscOptionsGetBool (NULL, "-help", &flg, NULL);
  CHKERRQ (ierr);

  appctx.m = m;
  appctx.measurement_freq = measurement_freq;
  appctx.h0 = h0;
  appctx.w = w;

  appctx.norm_2 = 1.0;

  ierr = MPI_Comm_size (PETSC_COMM_WORLD, &size);
  CHKERRQ (ierr);
  ierr = MPI_Comm_rank (PETSC_COMM_WORLD, &rank);
  CHKERRQ (ierr);

  appctx.size = size;
  appctx.rank = rank;

  ierr =
    PetscPrintf (PETSC_COMM_WORLD,
		 "Solving a linear TS problem, number of processors = %d\n",
		 size);
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create vector data structures
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create distributed array (DMDA) to manage parallel grid and vectors
     and to set up the ghost point communication pattern.  There are M
     total grid values spread equally among all the processors.
   */

  ierr =
    DMDACreate1d (PETSC_COMM_WORLD, DM_BOUNDARY_NONE, m, 1, 1, NULL,
		  &appctx.da);
  CHKERRQ (ierr);

  /*
     Extract global and local vectors from DMDA; we use these to store the
     approximate solution.  Then duplicate these for remaining vectors that
     have the same types.
   */
  ierr = DMCreateGlobalVector (appctx.da, &u);
  CHKERRQ (ierr);
  ierr = DMCreateLocalVector (appctx.da, &appctx.u_local);
  CHKERRQ (ierr);

  /*
     Create local work vector for use in evaluating right-hand-side function;
   */
  ierr = VecDuplicate (appctx.u_local, &appctx.localwork);
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set up displays to show graphs of the solution 
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  if (appctx.debug)
    {
      ierr =
	PetscViewerDrawOpen (PETSC_COMM_WORLD, 0, "Probability Distribution",
			     PETSC_DECIDE, PETSC_DECIDE, PETSC_DRAW_HALF_SIZE,
			     PETSC_DRAW_HALF_SIZE, &appctx.viewer1);
      CHKERRQ (ierr);
      ierr = PetscViewerDrawGetDraw (appctx.viewer1, 0, &draw);
      CHKERRQ (ierr);
      ierr = PetscDrawSetDoubleBuffer (draw);
      CHKERRQ (ierr);
    }
  else if (!flg)
    {
      ierr =
	PetscPrintf (PETSC_COMM_WORLD,
		     "\n##########################################");
      ierr = PetscPrintf (PETSC_COMM_WORLD, "\n# time\t\t# variance\t");
      ierr =
	PetscPrintf (PETSC_COMM_WORLD,
		     "\n##########################################");
      CHKERRQ (ierr);

    }


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = TSCreate (PETSC_COMM_WORLD, &ts);
  CHKERRQ (ierr);
  /* In this case, the dynamics is always linear */
  TSSetProblemType (ts, TS_LINEAR);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set user-defined monitoring routine
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSMonitorSet (ts, Monitor, &appctx, NULL);
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix data structure; set matrix evaluation routine.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate (PETSC_COMM_WORLD, &A);
  CHKERRQ (ierr);
  ierr = MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE, m, m);
  CHKERRQ (ierr);
  ierr = MatSetFromOptions (A);
  CHKERRQ (ierr);
  ierr = MatSetUp (A);
  CHKERRQ (ierr);



  /*
     For linear problems with a time-dependent f(u,t) in the equation
     u_t = f(u,t), the user provides the discretized right-hand-side
     as a time-dependent matrix.
   */
  ierr = TSSetRHSFunction (ts, NULL, TSComputeRHSFunctionLinear, &appctx);
  CHKERRQ (ierr);
  ierr = TSSetRHSJacobian (ts, A, A, RHSMatrixSchrodinger, &appctx);
  CHKERRQ (ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set solution vector and initial timestep
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  dt = INIT_STEP;
  ierr = TSSetInitialTimeStep (ts, INIT_TIME, dt);
  CHKERRQ (ierr);
  ierr = TSSetSolution (ts, u);
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize timestepping solver:
     - Set the solution method to be the Backward Euler method.
     - Set timestepping duration info
     - Set default tolerances
     Then set runtime options, which can override these defaults.
     For example,
     -ts_max_steps <maxsteps> -ts_final_time <maxtime>
     to override the defaults set by TSSetDuration().
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = TSSetType (ts, TSRK);
  CHKERRQ (ierr);
  /*Default solver is Runge Kutta (5) Prince-Dormand (4) */
  ierr = TSRKSetType (ts, TSRK5DP);
  CHKERRQ (ierr);
  ierr = TSSetDuration (ts, time_steps_max, time_total_max);
  CHKERRQ (ierr);
  ierr = TSSetExactFinalTime (ts, TS_EXACTFINALTIME_INTERPOLATE);
  CHKERRQ (ierr);

  ierr = TSSetTolerances (ts, ABSERROR, NULL, RELERROR, NULL);
  CHKERRQ (ierr);

  ierr = TSSetFromOptions (ts);
  CHKERRQ (ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve the problem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Evaluate initial conditions
   */
  ierr = InitialConditionsLocalized (u, &appctx);
  CHKERRQ (ierr);

  /*
     Run the timestepping solver if help flag is not set
   */
  if (!flg)
    {
      ierr = TSSolve (ts, u);
      CHKERRQ (ierr);
      ierr = TSGetSolveTime (ts, &ftime);
      CHKERRQ (ierr);
      ierr = TSGetTimeStepNumber (ts, &steps);
      CHKERRQ (ierr);

      final_runtime = MPI_Wtime ();
      ierr =
	PetscPrintf (PETSC_COMM_WORLD,
		     "\n##########################################");
      /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         View timestepping solver info
         - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
      ierr =
	PetscPrintf (PETSC_COMM_WORLD,
		     "\nTotal timesteps %D, Final time %g\n", steps,
		     (double) ftime);
      CHKERRQ (ierr);
      ierr =
	PetscPrintf (PETSC_COMM_WORLD, "Final. (2 norm) = %g \n",
		     (double) (appctx.norm_2));
      CHKERRQ (ierr);
      ierr =
	PetscPrintf (PETSC_COMM_WORLD, "Total runtime (s) = %g \n",
		     (double) (final_runtime - init_runtime));
      CHKERRQ (ierr);
    }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = TSDestroy (&ts);
  CHKERRQ (ierr);
  ierr = MatDestroy (&A);
  CHKERRQ (ierr);
  ierr = VecDestroy (&u);
  CHKERRQ (ierr);

  if (appctx.debug)
    {
      ierr = PetscViewerDestroy (&appctx.viewer1);
      CHKERRQ (ierr);
    }
  ierr = VecDestroy (&appctx.localwork);
  CHKERRQ (ierr);
  ierr = VecDestroy (&appctx.u_local);
  CHKERRQ (ierr);
  ierr = DMDestroy (&appctx.da);
  CHKERRQ (ierr);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
     - finalizes the PETSc libraries as well as MPI
     - provides summary and diagnostic information if certain runtime
     options are chosen (e.g., -log_summary).
   */
  ierr = PetscFinalize ();
  return 0;
}

/* --------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "InitialConditionsLocalized"
/*
   InitialConditionsLocalized - Computes the solution at the initial time.

   Input Parameter:
   u - uninitialized solution vector (global)
   appctx - user-defined application context

   Output Parameter:
   u - vector with solution at initial time (global)
*/
PetscErrorCode
InitialConditionsLocalized (Vec u, AppCtx * appctx)
{
  PetscScalar *u_localptr;
  PetscInt i, mybase, myend;
  PetscErrorCode ierr;
  PetscScalar zero = 0.0, one = 1.0;

  /*
     Determine starting point of each processor's range of
     grid values.
   */
  ierr = VecGetOwnershipRange (u, &mybase, &myend);
  CHKERRQ (ierr);

  /*
     Get a pointer to vector data.
     - For default PETSc vectors, VecGetArray() returns a pointer to
     the data array.  Otherwise, the routine is implementation dependent.
     - You MUST call VecRestoreArray() when you no longer need access to
     the array.
   */
  ierr = VecGetArray (u, &u_localptr);
  CHKERRQ (ierr);

  /*
     We initialize the solution array by simply writing the solution
     directly into the array locations.
     This will set the midpoint to unity
   */
  for (i = mybase; i < myend; i++)
    if (i == floor ((appctx->m) / 2))
      u_localptr[i - mybase] = one;
    else
      u_localptr[i - mybase] = zero;

  /*
     Restore vector
   */
  ierr = VecRestoreArray (u, &u_localptr);
  CHKERRQ (ierr);

  /*
     Print debugging information if desired
   */
  if (appctx->debug)
    {
      ierr =
	PetscPrintf (appctx->comm,
		     "\n########## Initial Vector ##########\n");
      CHKERRQ (ierr);
      ierr = VecView (u, PETSC_VIEWER_STDOUT_WORLD);
      CHKERRQ (ierr);
      ierr =
	PetscPrintf (appctx->comm,
		     "\n####################################\n");
    }

  return 0;
}

/* --------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Monitor"
/*
   Monitor - User-provided routine to monitor the solution computed at
   each timestep.  This example plots the solution and returns the norm

   Input Parameters:
   ts     - the timestep context
   step   - the count of the current step (with 0 meaning the
             initial condition)
   time   - the current time
   u      - the solution at this timestep
   ctx    - the user-provided context for this monitoring routine.
            In this case we use the application context which contains
            information about the problem size, workspace and the exact
            solution.
*/
PetscErrorCode
Monitor (TS ts, PetscInt step, PetscReal time, Vec u, void *ctx)
{
  AppCtx *appctx = (AppCtx *) ctx;	/* user-defined application context */
  PetscErrorCode ierr;

  /*
     If not debugging, print out regular output 
     - Regular output consists of variance, mean and edge mobility  
     - First get matrices for operators x, x^2 and diag(0,0,...1)
     - Then compute the three quantitites above using linear algebra
   */
  if (!appctx->debug)
    {
      Mat X, Xsq;
      Vec sites, temp;
      PetscInt u_loc_size;
      PetscScalar xbar, xsqbar, variance;

      /*Duplicate the parallel structure of u into temp. */
      ierr = VecDuplicate (u, &temp);
      CHKERRQ (ierr);

      /*
         Calculation of mean and variance  
       */

      /*Duplicate the parallel structure of u into sites, then set the sites */
      ierr = VecDuplicate (u, &sites);
      CHKERRQ (ierr);
      PetscInt j, start, end;
      PetscScalar jval;
      ierr = VecGetOwnershipRange (sites, &start, &end);
      for (j = start; j < end; j++)
	{
	  jval = j;
	  ierr = VecSetValues (sites, 1, &j, &jval, INSERT_VALUES);
	  CHKERRQ (ierr);
	}
      ierr = VecAssemblyBegin (sites);
      CHKERRQ (ierr);
      ierr = VecAssemblyEnd (sites);
      CHKERRQ (ierr);

      /*Get local size of u */
      ierr = VecGetLocalSize (u, &u_loc_size);
      CHKERRQ (ierr);

      /* Build position (site) operator */
      ierr = MatCreate (MPI_COMM_WORLD, &X);
      CHKERRQ (ierr);
      /*The local block of the matrix must multiply with local part of u */
      ierr = MatSetSizes (X, u_loc_size, u_loc_size, appctx->m, appctx->m);
      CHKERRQ (ierr);
      ierr = MatSetUp (X);
      CHKERRQ (ierr);
      ierr = MatDiagonalSet (X, sites, INSERT_VALUES);
      CHKERRQ (ierr);

      ierr = MatMult (X, u, temp);	/* Computes |temp> = X |u> */
      CHKERRQ (ierr);

      /*Computes xbar = <u|temp> = <u|X|u> */
      ierr = VecDot (temp, u, &xbar);
      CHKERRQ (ierr);

      /*Build position^2 operator trivially. Set Xsq = X^2 */
      ierr = MatMatMult (X, X, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Xsq);
      CHKERRQ (ierr);
      ierr = MatMult (Xsq, u, temp);	/* Computes |temp> = X^2 |u> */
      CHKERRQ (ierr);
      /*Computes xsqbar = <u|temp> = <u|X^2|u> */
      ierr = VecDot (temp, u, &xsqbar);
      CHKERRQ (ierr);

      variance = xsqbar - (xbar * xbar);
      variance = variance / ((appctx->m) * (appctx->m));	/*Normalization */

      xbar = xbar / (appctx->m);
      /*
         Dump output to stdout  
       */
      ierr = MPI_Barrier (PETSC_COMM_WORLD);
      CHKERRQ (ierr);
      ierr =
	PetscPrintf (PETSC_COMM_WORLD, "\n# %2.4lf\t# %1.8lf\t", time,
		     variance);
      CHKERRQ (ierr);

      ierr = VecDestroy (&temp);
      CHKERRQ (ierr);
      ierr = VecDestroy (&sites);
      CHKERRQ (ierr);
      ierr = MatDestroy (&X);
      CHKERRQ (ierr);
      ierr = MatDestroy (&Xsq);
      CHKERRQ (ierr);
    }

  /*
     Print debugging information if desired
   */
  if (appctx->debug)
    {
      PetscReal norm_2;

      /*
         Compute and store norm of wavefunction in norm_2 
       */
      ierr = VecNorm (u, NORM_2, &norm_2);
      CHKERRQ (ierr);
      appctx->norm_2 = norm_2;

      /* Build the lattice sites */
      Vec uabs;
      ierr = VecDuplicate (u, &uabs);
      CHKERRQ (ierr);

      ierr = VecCopy (u, uabs);
      CHKERRQ (ierr);
      ierr = VecAbs (uabs);
      CHKERRQ (ierr);

      ierr = VecPointwiseMult (uabs, uabs, uabs);
      CHKERRQ (ierr);
      /*
         View a graph of the current iterate
       */
      ierr =
	PetscPrintf (appctx->comm, "time = %lf\t norm = %lf\n", time, norm_2);

      ierr = VecView (uabs, appctx->viewer1);
      CHKERRQ (ierr);

      ierr = VecDestroy (&uabs);
      CHKERRQ (ierr);

    }
  return 0;
}

/* --------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "RHSMatrixSchrodinger"
/*
   RHSMatrixSchrodinger - User-provided routine to compute the right-hand-side
   matrix for the heat equation.

   Input Parameters:
   ts - the TS context
   t - current time
   global_in - global input vector
   dummy - optional user-defined context, as set by TSetRHSJacobian()

   Output Parameters:
   AA - Jacobian matrix
   BB - optionally different preconditioning matrix
   str - flag indicating matrix structure

  RHSMatrixSchrodinger computes entries for the locally owned part of the system.
   - Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.
   - Each processor needs to insert only elements that it owns
     locally (but any non-local elements will be sent to the
     appropriate processor during matrix assembly).
   - Always specify global row and columns of matrix entries when
     using MatSetValues(); we could alternatively use MatSetValuesLocal().
   - Here, we set all entries for a particular row at once.
   - Note that MatSetValues() uses 0-based row and column numbers
     in Fortran as well as in C.
*/
PetscErrorCode
RHSMatrixSchrodinger (TS ts, PetscReal t, Vec X, Mat AA, Mat BB, void *ctx)
{
  Mat A = AA;			/* Jacobian matrix */
  AppCtx *appctx = (AppCtx *) ctx;	/* user-defined application context */
  PetscErrorCode ierr;
  PetscInt m = appctx->m;
  PetscScalar drive;
  PetscReal h0 = appctx->h0;
  PetscReal w = appctx->w;


  drive = w * t;
  drive = PetscSinReal (drive);
  drive = 2.0 * PETSC_i * h0 * drive;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute entries for the locally owned part of the matrix
     The Jacobian Matrix in our case is blocked
     - The values are A_{lm} 
     =  i\left(\delta_{lm+1}+\delta_{lm-1}\right) + 2i h(t) \delta_{lm}
     - Thus the superdiagonal and subdiagonal elements are just i, and 
     diagonal elements are the drive
     - The elements are set rowwise.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInt rowstart, rowend, rowcount;
  PetscScalar edge[2], middle[3];
  PetscInt idxedge[2], idxmiddle[3];
  ierr = MatGetOwnershipRange (A, &rowstart, &rowend);
  CHKERRQ (ierr);
  for (rowcount = rowstart; rowcount < rowend; rowcount++)
    {
      if (rowcount == 0)	//If it is the first row of the whole matrix
	{
	  edge[0] = drive;
	  edge[1] = PETSC_i;
	  idxedge[0] = 0;
	  idxedge[1] = 1;
	  ierr =
	    MatSetValues (A, 1, &rowcount, 2, idxedge, edge, INSERT_VALUES);
	  CHKERRQ (ierr);
	}
      else if (rowcount == m - 1)	//If it is the last row of the whole matrix
	{
	  edge[0] = PETSC_i;
	  edge[1] = drive;
	  idxedge[0] = m - 2;
	  idxedge[1] = m - 1;
	  ierr =
	    MatSetValues (A, 1, &rowcount, 2, idxedge, edge, INSERT_VALUES);
	  CHKERRQ (ierr);
	}
      else
	{
	  middle[0] = PETSC_i;
	  middle[1] = drive;
	  middle[2] = PETSC_i;
	  idxmiddle[0] = rowcount - 1;
	  idxmiddle[1] = rowcount;
	  idxmiddle[2] = rowcount + 1;
	  ierr =
	    MatSetValues (A, 1, &rowcount, 3, idxmiddle, middle,
			  INSERT_VALUES);
	  CHKERRQ (ierr);
	}
    }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Complete the matrix assembly process and set some options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Assemble matrix, using the 2-step process:
     MatAssemblyBegin(), MatAssemblyEnd()
     Computations can be done while messages are in transition
     by placing code between these two statements.
   */
  ierr = MatAssemblyBegin (A, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);
  ierr = MatAssemblyEnd (A, MAT_FINAL_ASSEMBLY);
  CHKERRQ (ierr);

  /*
     Set and option to indicate that we will never add a new nonzero location
     to the matrix. If we do, it will generate an error.
   */
  ierr = MatSetOption (A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  CHKERRQ (ierr);
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "kdel"
PetscInt
kdel (PetscInt i, PetscInt j)
{
  PetscInt value = 0;
  if (i == j)
    value = 1;
  return value;
}
