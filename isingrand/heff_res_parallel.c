#include "heff_res_parallel.h"
#include "params.h"

/* Display program usage, and exit.
 */
void
display_usage (int argc, char *argv[])
{
  printf
    ("\nUsage: %s [--option1 value1] [--option2 value2] [-flag1 value1] [-flag2 value2] ...\n",
     argv[0]);
  puts
    ("Integrates the time independent renormalized Hamiltonian of a driven disordered Ising model.");
  printf
    ("Current hardcoded lattice size is %d. To change, edit file 'params.h' and recompile\n",
     N);
  printf
    ("Current hardcoded # of random instances per MPI process is %d. To change, edit file 'params.h' and recompile\n",
     NRAND);
#if defined(FIELD_DISORDER_OFF)
  puts ("\nField disorder is OFF!\n");
#else
  puts ("\nField disorder is ON!\n");
#endif
#if defined(INITCONDS_GHZ)
  puts
    ("Current hardcoded initial state is the GHZ state. To change, edit file 'params.h' and recompile\n");
#else
  puts
    ("Current hardcoded initial state is the ordered ground state (all spins up). To change, edit file 'params.h' and recompile\n");
#endif
  printf
    ("Current hardcoded integrator method is given at runtime. To change, edit file 'integrator.h' and recompile\n\n");
  puts (options);
  printf ("\n");
  puts (envv);
  printf ("\nNOTE:\n");
  printf
    (" For environment variable options, please consult the GNU Scientific Library Reference Manual.\n");
  exit (EXIT_FAILURE);
}

//Function for loading static array of GArrays in memory
void
load_arrays (GArray * garrays[NRAND])
{
  gint j;
  // put the load arrays stuff here.
  for (j = 0; j < NRAND; j++)
    {
      garrays[j] = g_array_new (FALSE, FALSE, sizeof (gdouble));
    }
}

/* Progress bar function for verbose mode.
   Process has done x out of n rounds,
   and we want a bar of width w and resolution r */
static inline void
loadbar (int x, int n, int r, int w)
{
  // Only update r times.
  if (x % (n / r) != 0)
    return;

  // Calculuate the ratio of complete-to-incomplete.
  float ratio = x / (float) n;
  int c = ratio * w;

  // Show the percentage complete.
  printf ("%3d%% [", (int) (ratio * 100));

  // Show the load bar.
  for (x = 0; x < c; x++)
    printf ("=");

  for (x = c; x < w; x++)
    printf (" ");

  // ANSI Control codes to go back to the
  // previous line and clear it.
  printf ("]\n\033[F\033[J");
}

//Calculates beta function
double
betaval (double eta)
{
  int i, n;
  double beta = 0.0;

  for (i = 0; i < BETATRUNC; i++)
    {
      n = 2 * i + 1;
      beta = beta + (gsl_sf_bessel_Jn (n, eta) / n);
    }
  beta = 2.0 * beta;
  return beta;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char *argv[])
{
  int mpi_numprocs, mpi_rank;
  time_t begin, end;
  double wall_begin, wall_end, wall_diff, wall_res;
  paramspace_pt param;
  double random_no;
  int tid = omp_get_thread_num ();
  int nthreads;
  //Create output block
  double rand_arrh[NRAND][N], rand_arrj[NRAND][N];
  double data[2];
  GArray *magrand_out[NRAND];

  //Initiate MPI environment. This is for one process per node
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  //Now, each of these processes runs the doce below.
  //Each MPI process spawns $OMP_NUM_THREADS number of OMP threads

  load_arrays (magrand_out);

  //Create array of NRAND files. One for each instance of randomness
  FILE *outfiles[NRAND];

  char appender[100], tempstring[100];

  /************************************  Input Block  ************************************/

  //Initiate defaults
  int opt = 0;
  int noentry = 1;
  int longIndex = 0;
  globalArgs.gamma_a = 0.0;
  globalArgs.omega = 0.0;
  globalArgs.alpha = 0.0;
  globalArgs.finaltime = 1.0;
  globalArgs.tsteps = 100;
  globalArgs.outFileName_Mag = "mags.dat";
  globalArgs.check_unitarity = 0;
  globalArgs.verbosity = 0;
  globalArgs.inputFiles = NULL;
  globalArgs.numInputFiles = 0;

  //Process the arguments with getopt_long(), then populate globalArgs. 
  opt = getopt_long (argc, argv, optString, longOpts, &longIndex);
  while (opt != -1)
    {
      switch (opt)
	{

	case 'a':
	  globalArgs.gamma_a = atof (optarg);
	  param.gamma_a = globalArgs.gamma_a;
	  noentry = 0;
	  break;

	case 'e':
	  globalArgs.besj_zero = atoi (optarg);
	  globalArgs.eta = gsl_sf_bessel_zero_J0 (globalArgs.besj_zero);
	  globalArgs.omega = 0.0;	//No drive frequency here since this is a simple quench problem.
	  param.eta = globalArgs.eta;
	  param.omega = globalArgs.omega;
	  noentry = 0;
	  break;

	case 'p':
	  globalArgs.alpha = atof (optarg);
	  param.alpha = globalArgs.alpha;
	  noentry = 0;
	  break;

	case 't':
	  globalArgs.finaltime = atof (optarg);
	  noentry = 0;
	  break;

	case 'n':
	  globalArgs.tsteps = atoi (optarg);
	  noentry = 0;
	  break;

	case 'm':
	  globalArgs.outFileName_Mag = optarg;
	  noentry = 0;
	  break;

	case 'u':
	  globalArgs.check_unitarity++;
	  noentry = 0;

	case 'v':
	  globalArgs.verbosity++;
	  noentry = 0;
	  break;

	default:
	  display_usage (argc, argv);
	  exit (EXIT_FAILURE);
	}

      opt = getopt_long (argc, argv, optString, longOpts, &longIndex);

    }
  if ((noentry != 0) && (mpi_rank == 0) && (tid == 0))
    {
      display_usage (argc, argv);
      exit (EXIT_FAILURE);
    }

  globalArgs.inputFiles = argv + optind;
  globalArgs.numInputFiles = argc - optind;


  /**********************************  End Input Block  **********************************/

  (void) time (&begin);
  wall_begin = omp_get_wtime ();
  int i, j;

  double y[DIM], y_init[DIM];	//This is the full state
  int status;
  double t_dyn, t_inc = globalArgs.finaltime / (globalArgs.tsteps - 1);

  int x = 0, res = BAR_RES, w = BAR_RES;	//loadbar point, resolution and width

/*******************************  Random Numbers Block  *******************************/
  gsl_rng *r;
  const gsl_rng_type *rtype;
  gsl_rng_env_setup ();

  rtype = gsl_rng_default;
  r = gsl_rng_alloc (rtype);

  //This is so that different mpi procs have different seeds and don't 
  //generate the same pseudo random numbers
  gsl_rng_default_seed = gsl_rng_default_seed * (mpi_rank + 1);
  gsl_rng_set (r, gsl_rng_default_seed);

  if ((globalArgs.verbosity != 0) && (mpi_rank == 0) && (tid == 0))
    {
      printf ("\n");
      puts ("Random Number Generator Info:");
      printf ("=============================================\n");
      printf ("Generator type \t=\t%s\n", gsl_rng_name (r));
      printf ("Seed value \t=\t%lu\n", gsl_rng_default_seed);
      printf ("First value \t=\t%lf\n", gsl_rng_uniform (r));
      printf ("=============================================\n");
    }

  int randcount;
  double res_elem;
  //Evaluate betaval(eta)
  double beta = betaval (param.eta);

  //Initialize NRAND * N random numbers in the range (-1,1) to use as disorder in lattice
  for (i = 0; i < NRAND; i++)
    for (j = 0; j < N; j++)
      {
	random_no = gsl_rng_uniform (r);	//This generates a rand number in the interval [0,1) with uniform probability
	random_no = 2.0 * random_no - 1.0;	//Scale it to the interval (-1,1)
	rand_arrh[i][j] = random_no;
#if defined(FIELD_DISORDER_OFF)
	rand_arrh[i][j] = 0.0;
#endif

	random_no = gsl_rng_uniform (r);	//This generates a random number in the interval [0,1) with uniform probability
	random_no = 2.0 * random_no - 1.0;	//Scale it to the interval (-1,1)
	res_elem = param.alpha * random_no;
	res_elem = res_elem * param.eta * beta;
	res_elem = res_elem * rand_arrh[i][j] / param.gamma_a;
	rand_arrj[i][j] = res_elem;
      }

  /******************************  End Random Numbers Block  *****************************/

  //This loop initializes N instances of random arrays and uses them in N integrations
  //This loop is parallelized using OpenMP
#pragma omp parallel private(status, randcount, i, random_no, y, x, t_dyn, data) firstprivate(param, y_init) shared(rand_arrh, rand_arrj, magrand_out, globalArgs)
  {
#pragma omp for schedule(guided, CHUNKSIZE)
    for (randcount = 0; randcount < NRAND; randcount++)
      {
	for (i = 0; i < N; i++)
	  {
	    //Copy from random array of h and j for each instance of randcount
#pragma omp flush (rand_arrh, rand_arrj)
	    {
	      param.hrand[i] = param.alpha * rand_arrh[randcount][i];	//Scale random number by the perturbation
	      param.jrand[i] = param.alpha * rand_arrj[randcount][i];
	    }
	  }
   /******************************  Initial Conditions Block  *****************************/
#if defined(INITCONDS_GHZ)
	//This sets the initial conditions to a random but normalized state
	set_initconds_ghz (y_init);
#else
	//This sets the initial conditions to all ising spins up ie the ordered ground state
	set_initconds_allup (y_init);
#endif
  /******************************  End Initial Conditions Block  *****************************/

  /******************************  Set Initial Conditions Block  *****************************/
	for (i = 0; i < DIM; i++)
	  {
	    y[i] = y_init[i];
	  }

  /******************************  End Set Initial Conditions Block  *****************************/

#pragma omp flush (globalArgs)
  /********************************    Integration  Block    ********************************/
	if ((globalArgs.verbosity != 0) && (mpi_rank == 0) && (tid == 0))
	  {
	    printf ("\n");
	    printf ("Now beginning integration with magnetization = %lf\n",
		    magnetization (y));
	    printf ("Number of threads = %d\n", omp_get_num_threads ());
	    printf ("=============================================\n");
	  }
	x = 0;			//x counts the number of executions of integrator  
	status = GSL_SUCCESS;
	for (t_dyn = 0.0; t_dyn < globalArgs.finaltime; t_dyn = t_dyn + t_inc)
	  {

	    if (tid == 0)
	      {
		if (status != GSL_SUCCESS)
		  printf ("\n Error. Integrator failed at t=%lf", t_dyn);
		else
		  {
		    data[0] = t_dyn;
		    data[1] = magnetization (y);

		    if (globalArgs.check_unitarity != 0)
		      printf
			("\n Nonunitary parameter at t = %f is = %f",
			 t_dyn, check_unitarity (y));

#pragma omp critical
		    g_array_append_vals (magrand_out[randcount], &data, 2);

		  }
#pragma omp flush (globalArgs)
		{
		  if ((globalArgs.verbosity != 0) && (mpi_rank == 0)
		      && (tid == 0))
		    loadbar (x, globalArgs.tsteps, res, w);
		}
	      }
	    status = integrate (y, t_dyn, t_dyn + t_inc, &param);
	    x++;		//x counts the number of executions of integrator  
	  }

#pragma omp flush (globalArgs)
	if ((globalArgs.verbosity != 0) && (mpi_rank == 0) && (tid == 0))
	  {
	    printf ("\n");
	    printf ("=============================================\n");
	    printf ("End integration with magnetization = %lf",
		    magnetization (y));
	  }
  /********************************  End Integration  Block  ********************************/
	nthreads = omp_get_num_threads ();
      }
  }
  (void) time (&end);
  wall_end = omp_get_wtime ();
  wall_diff = wall_end - wall_begin;
  wall_res = 1.0 / omp_get_wtick ();


  /************************************  Output Block  ************************************/
  if (tid == 0)
    {

      //Output magrand_out. magrand_out[i] is the time, mag data of the ith random number iteration
      //First, open NRAND files with NRAND filenames
      for (i = 0; i < NRAND; i++)
	{
	  sprintf (appender, ".%d_seed_%ld", i, gsl_rng_default_seed);
	  strcpy (tempstring, globalArgs.outFileName_Mag);
	  outfiles[i] =
	    fopen (strcat (globalArgs.outFileName_Mag, appender), "w");
	  strcpy (globalArgs.outFileName_Mag, tempstring);
	}
      //Then, write output to each filenames
      for (i = 0; i < NRAND; i++)
	{
	  for (j = 0; j < magrand_out[i]->len; j = j + 2)
	    {
	      fprintf (outfiles[i], "%lf %lf\n",
		       g_array_index (magrand_out[i], gdouble, j),
		       g_array_index (magrand_out[i], gdouble, j + 1));
	    }
	}
      printf ("\n");
      //Print entered values
      printf ("\n");
      puts ("Runtime Parameters:");
      printf ("=============================================\n");
      printf ("Drive amplitude \t=\t%lf\n", param.gamma_a);
      printf ("eta \t\t\t=\t%lf\n", param.eta);
      printf ("beta \t\t\t=\t%lf\n", beta);
      printf ("Perturbation amplitude\t=\t%lf\n", param.alpha);
      printf ("Final time\t\t=\t%lf\n", globalArgs.finaltime);
      printf ("Time steps\t\t=\t%d\n", globalArgs.tsteps);
      printf ("# of rands\t\t=\t%d\n", NRAND * mpi_numprocs);
      printf ("Lattice size\t\t=\t%d\n", N);
      printf ("=============================================\n");
      puts ("Final Output:");
      printf ("=============================================\n");
      printf ("Total Runtime (root) \t=\t%lf\n", difftime (end, begin));
      printf ("Total Walltime \t\t=\t%lf\n", wall_diff);
      printf ("Clock resolution \t=\t%lf\n", wall_res);
      printf ("Number of threads \t=\t%d\n", nthreads);
      printf ("=============================================\n");
    }
  /**********************************  End Output Block  **********************************/

  //fclose (globalArgs.outFile);
  gsl_rng_free (r);

  for (i = 0; i < NRAND; i++)
    {
      fclose (outfiles[i]);
      g_array_free (magrand_out[i], TRUE);
    }
  MPI_Finalize ();
  return EXIT_SUCCESS;
}
