#include "params_ics.h"
#include "lzs_mbt_2d_dwave.h"
#undef __FUNCT__
#define __FUNCT__ "main"
int
main (int argc, char **argv)
{
  time_t begin, end;
  paramspace_pt params;
  initconds_ic initconds;
  qfactors qf;
  int argv_count = 1, result;
  char stateout, gndstate;
  char reltime[]="Relaxation time for fidelity";
  char resenstr[]="Initial Energy";

  int nprocs, pid;
  int pinit_loc, pfinal_loc;
  int strobecount = 0;
  long local_kcount, global_kcount;

  double tau;
  double t, t_init, t_periods, scaled_time;
  double t_extra, t_final, t_inc;
  double periods, t_extra_frac, wt;
  double phi[4], phi_init[4], phi_adb[4];
  double omega, time_period;
  double kx, ky, ktot, strobeoff;
  double kxinit, kxfinal, kxinc;
  double kyinit, kyfinal, kyinc;

  long p, gridsize, t_gridsize, gridsize_local;

  GArray *phik = g_array_new (FALSE, FALSE, sizeof (double));	//phi and momentum updated as (kx+Nky,phi[0],phi[1],phi[2],phi[3]...)
  GArray *time_phik = g_array_new (FALSE, FALSE, sizeof (double));	//time, phi and momentum updated as (t,kx,ky,phi[0-4],defect_density,fidelity,...)
  GArray *time_magnetization = g_array_new (FALSE, FALSE, sizeof (double));
  GArray *time_fidelity = g_array_new (FALSE, FALSE, sizeof (double));
  GArray *time_delta = g_array_new (FALSE, FALSE, sizeof (double));

  GArray *time_defectdensity = g_array_new (FALSE, FALSE, sizeof (double));
  GArray *time_residualenergy = g_array_new (FALSE, FALSE, sizeof (double));
  GArray *time_fidsuscpt = g_array_new (FALSE, FALSE, sizeof (double));
  GArray *time_adbmag = g_array_new (FALSE, FALSE, sizeof (double));
  
  GArray *time_fidelity_deriv_mag = g_array_new (FALSE, FALSE, sizeof (double));
  GArray *time_gfunc = g_array_new (FALSE, FALSE, sizeof (double));
  
  MPI_File tphikfile;

  //argument error checking
  if (argc != ARGNO + 1)
    {
      int aerr = arg_err (argc);
      exit (aerr);
    }
  //start mpi
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &pid);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs);

  if (pid == 0)
    begin_out ();

  /***************************************Input block****************************************/

  strobeoff = atof (argv[argv_count++]);	//offset error for strobing at time period
  FILE *tmagfile = fopen (argv[argv_count++], "w");
  FILE *tfidfile = fopen (argv[argv_count++], "w");
  FILE *deltafile = fopen (argv[argv_count++], "w");
  FILE *ndfile = fopen (argv[argv_count++], "w");
  FILE *tresenfile = fopen (argv[argv_count++], "w");
  FILE *tfidsfile = fopen (argv[argv_count++], "w");
  FILE *tmagadbfile = fopen (argv[argv_count++], "w");
  FILE *tfidderivfile = fopen (argv[argv_count++], "w");
  FILE *tgfuncfile = fopen (argv[argv_count++], "w");


  //Parallel IO file for t-phik. Output buffer
  int ierr = MPI_File_open (MPI_COMM_WORLD, argv[argv_count++],
			    MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_EXCL,
			    MPI_INFO_NULL, &tphikfile);

  //Inout number of period sweeps
  periods = atof (argv[argv_count++]);
  //Input Extra time after n sweeps in fractions of the period
  t_extra_frac = atof (argv[argv_count++]);

  //Input system parameters
  params.omega = atof (argv[argv_count++]);	//Drive frequency
  params.delta_real = atof (argv[argv_count++]);
  params.delta_imag = atof (argv[argv_count++]);

  params.mu0 = atof (argv[argv_count++]);	//Chempot at thermal eqm
  params.muamp = atof (argv[argv_count++]);	//amplitude of chempot

  gridsize = atoi (argv[argv_count++]);	//\vec{k} grid size
  t_gridsize = atoi (argv[argv_count++]);	//t grid size

  //Determine if full state output is desired
  if (strcmp (argv[argv_count++], "y") == 0)
    {
      stateout = 'y';
      params.stateout = stateout;
    }
  else
    {
      stateout = 'n';
      params.stateout = stateout;
    }

  if (strcmp (argv[argv_count], "d") == 0)
    {				//whether t=0 gnd state is diabatic , bcs, or random. Answer d, s, or r
      gndstate = 'd';
      params.gndstate = gndstate;
    }
  else if (strcmp (argv[argv_count], "s") == 0)
    {
      gndstate = 's';
      params.gndstate = gndstate;
    }
   else
   {
     gndstate = 'r';
     params.gndstate = gndstate;
  } 

  /************************************End Input block***************************************/

  initconds.verbosity = NOT_VERBOSE;

  params.offset = 0.0;		//Drive amplitude offset is hard coded to 0
  omega = params.omega;
  time_period = 2.0 * M_PI / omega;

  //If random initial conditions is set, then initialize the random number generator
  //and generate gridsize * gridsize random numbers between -1 and 1
   double *rkarray;
   long rkcount;
   double rk;
   long fullgridsize = gridsize * gridsize ;
  rkarray = malloc (fullgridsize * sizeof (*rkarray));
  if(params.gndstate=='r'){
  const gsl_rng_type * rangen_type;
  gsl_rng * random_number_object;
  gsl_rng_env_setup();
  rangen_type = gsl_rng_default;
  random_number_object = gsl_rng_alloc (rangen_type);
   for(rkcount=0;rkcount<fullgridsize;rkcount++){
     rk = gsl_rng_uniform_pos (random_number_object);
     rk = 2.0 * rk - 1.0;
     rkarray[rkcount]=rk;
   }
   gsl_rng_free(random_number_object);
  }
 
  //Calculate the \vec{k} grid
  kxinit = -M_PI;		//FBZ within delta
  kyinit = -M_PI;
  kxfinal = M_PI;		//FBZ
  kyfinal = M_PI;
  kxinc = (kxfinal - kxinit) / (gridsize - 1);
  kyinc = (kyfinal - kyinit) / (gridsize - 1);
  t_init = 0.0;			//Initial time is now hard coded to 0
  initconds.t_init = t_init;

  //Calculate t_extra
  t_extra = 2 * t_extra_frac * M_PI / omega;
  initconds.t_extra = t_extra;

  //Calculate t_periods
  t_periods = 2.0 * periods * M_PI / omega;
  initconds.t_periods = t_periods;
  t_final = t_init + t_periods + t_extra;
  t_inc = (t_final - t_init) / (t_gridsize - 1);

  MPI_Barrier (MPI_COMM_WORLD);
  (void) time (&begin);

  if (nprocs > gridsize)
    {

      puts ("\n Error, number of processes is bigger than row gridsize...\n");
      exit (1);

    }
    
  //Divide the momemtum grid rowwise
  gridsize_local = gridsize / nprocs;
  long extra_rows = gridsize % nprocs;
  long gridsize_root = gridsize_local + extra_rows;
  //create local kx index ranges for each process
  if (pid == 0)
    {
      gridsize_local = gridsize_root;
      pinit_loc = 0;
      pfinal_loc = gridsize_local - 1;
    }
  else
    {
      pinit_loc = gridsize_root;
      pinit_loc = pinit_loc + (pid - 1) * gridsize_local;
      pfinal_loc = pinit_loc + gridsize_local - 1;

    }
    
  double probloc, prob, probavg;
  gsl_complex ukinit, vkinit;
  double phi_unit[4],phi_adbmag[4];
  gsl_complex ukadb,vkadb;
  double vkadb_abs;
  gsl_complex uk, vk;
  gsl_complex uv,vu,uvdiff;

  gsl_complex veff_cmplx;
  gsl_complex delta = gsl_complex_rect (params.delta_real, params.delta_imag);
  gsl_complex delta_rhs_local = gsl_complex_rect (0.0, 0.0);
  gsl_complex delta_rhs;
  gsl_complex temp;
  double delta_rhs_loc_real, delta_rhs_loc_imag;
  double delta_rhs_real, delta_rhs_imag;

  // For fermisurface integrations only. Disabling it for now
  //double spread = SPREAD_FRAC * gsl_complex_abs (delta);

  double uksq, vksq;
  double defect_density, defect_density_local, defect_density_global;
  double norm_var, norm_local, norm_global;
  double ekinit_local, ekinit_global, energy_k, ek_local, ek_global;
  double fidelity_susceptibility_local, fidelity_susceptibility_global;
  double probloc_adb, fidk_dia;
  double fidelity_local, fidelity_global;
  double adbmag_local, adbmag_global;
  double fk;
  gsl_complex gk,pk;
  gsl_complex gfunc_local;
  double gfunc_loc_real,gfunc_loc_imag;
  double gfunc_glob_real, gfunc_glob_imag;

  long probcount_loc, probcount;
  long kcount_loc, kcount;

  double deltaloc_real, deltaloc_imag;
  double deltaglob_real, deltaglob_imag;
  double deltaglob_abs;

  long kstride = 0;
  long  kycount;
  //Set initial conditions for phik
  //Also, calculate veff
  //Input initial conditions
  //THIS IS IN THE INITIAL CONDITIONS IN THE DIABATIC BASIS .
  kcount_loc = 0;
  ekinit_local = 0.0;
    
  for (p = pinit_loc; p <= pfinal_loc; p++)
    {
      kx = kxinit + p * kxinc;
      kycount = 0;
      for (ky = kyinit; ky <= kyfinal; ky = ky + kyinc)
	{
	  params.ampmax = cos (kx) + cos (ky);	//set the lattice kinetic energy
	  params.pk = cos (kx) - cos (ky)-ALPHA; //set the d-wave anisotropy factor

	  if (gndstate == 'd')
	    {
	      set_initconds_diabatic_gnd (phi_init);	//set initconds for a particular k; diabatic ground state
	    }
	  else if(gndstate == 's')
	    {
	      set_initconds_bcs_gnd (&params, &initconds, phi_init);	//set initconds for a particular k; BCS ground state
	    }
	  else
	  {
	    rk = rkarray[p + gridsize * kycount];
	    set_initconds_rand_gnd (&params, &initconds, phi_init, rk);	//set initconds for a particular k; random ground state
	  }  

	  //calculate local rhs of veff expression
	  GSL_SET_COMPLEX (&ukinit, phi_init[0], phi_init[2]);
	  GSL_SET_COMPLEX (&vkinit, phi_init[1], phi_init[3]);
	  temp = gsl_complex_mul (gsl_complex_conjugate (ukinit), vkinit);
	  delta_rhs_local = gsl_complex_add (delta_rhs_local, temp);

	  //Calculate initial energy
	  ekinit_local =
	    ekinit_local + energy_expectation (t_init, phi_init, &params);
	  //Rotate phi to adiabatic basis
	  adiabatic_trans_mbt (phi_adb, phi_init, t_init, &params);

	  //Copy phi_init to phi
	  int status = copyvec (phi_init, phi, 4);
	  if (status == GSL_FAILURE)
	    {
	      puts ("\n Copy error, size mismatch");
	      exit (9);
	    }
	  //copy phi to phiks
	  ktot = kx + gridsize * ky;
	  g_array_append_val (phik, ktot);
	  g_array_append_vals (phik, phi, 4);

	  kcount_loc++;
	  kycount++;
	}
    }
  MPI_Barrier (MPI_COMM_WORLD);
  //MPI_reduce to root process: The residual energy 
  MPI_Reduce (&ekinit_local, &ekinit_global, 1, MPI_DOUBLE, MPI_SUM, 0,
	      MPI_COMM_WORLD);
  //Split delta_rhs_local into real and imaginary parts
  delta_rhs_loc_real = GSL_REAL (delta_rhs_local);
  delta_rhs_loc_imag = GSL_IMAG (delta_rhs_local);

  //MPI_Allreduce to ALL processes: The local delta_rhs_loc datatypes
  MPI_Allreduce (&delta_rhs_loc_real, &delta_rhs_real, 1, MPI_DOUBLE, MPI_SUM,
		 MPI_COMM_WORLD);
  MPI_Allreduce (&delta_rhs_loc_imag, &delta_rhs_imag, 1, MPI_DOUBLE, MPI_SUM,
		 MPI_COMM_WORLD);

  //MPI_Allreduce to all processes: The momenta counted
  MPI_Allreduce (&kcount_loc, &kcount, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

  //average out the delta_rhs datatypes
  delta_rhs_real = delta_rhs_real / kcount;
  delta_rhs_imag = delta_rhs_imag / kcount;

  //Now, calculate veff
  GSL_SET_COMPLEX (&delta_rhs, delta_rhs_real, delta_rhs_imag);
  veff_cmplx = gsl_complex_div (delta, delta_rhs);
  params.veff = GSL_REAL (veff_cmplx);

  if (pid == 0)
    {
      //Report parameters inputted/computed
      params_out (&params);
      ekinit_global = ekinit_global / kcount;
      dump(resenstr,ekinit_global);
    }

  for (t = t_init; t <= t_final; t = t + t_inc)
    {
      probloc = 0.0;
      fidelity_local = 0.0;
      adbmag_local=0.0;
      probcount_loc = 0.0;
      defect_density_local = 0.0;
      norm_local = 0.0;
      ek_local = 0;
      fidelity_susceptibility_local = 0.0;
      gfunc_local = gsl_complex_rect(0.0,0.0);
      kstride = 0;

      global_kcount = 0;
      /************Many body correction that updates the gap 'delta' at time t.***************/
      local_kcount = update_delta_loc (phik, &params);
      //MPI allreduce the local deltas from local instances of 'params'
      deltaloc_real = params.delta_real;
      deltaloc_imag = params.delta_imag;
      MPI_Barrier (MPI_COMM_WORLD);
      MPI_Allreduce (&deltaloc_real, &deltaglob_real, 1, MPI_DOUBLE, MPI_SUM,
		     MPI_COMM_WORLD);
      MPI_Allreduce (&deltaloc_imag, &deltaglob_imag, 1, MPI_DOUBLE, MPI_SUM,
		     MPI_COMM_WORLD);
      MPI_Allreduce (&local_kcount, &global_kcount, 1, MPI_LONG, MPI_SUM,
		     MPI_COMM_WORLD);
      //Overwrite the local deltas in local 'params' with the global sum
      params.delta_real = deltaglob_real / global_kcount;	//correction for average sum
      params.delta_imag = deltaglob_imag / global_kcount;	//correction for average sum
      /***************Comment the above lines out to remove many body effects*****************/
      //Update the logarithm of delta
      delta = gsl_complex_rect(params.delta_real, params.delta_imag);
      for (p = pinit_loc; p <= pfinal_loc; p++)
	{			//kx is being parallelized
	  kx = kxinit + p * kxinc;
	  kycount = 0;
	  for (ky = kyinit; ky <= kyfinal; ky = ky + kyinc)
	    {
	      //copy phik to buffer phi
	      copytobuffer (phik, phi, 4, kstride);

	      params.ampmax = cos (kx) + cos (ky);	//set the lattice kinetic energy
	      params.pk = cos (kx) - cos (ky)-ALPHA; //anisotropy term

	      //integrate the probs for a particular value of \vec{k}= (kx, ky)
	      result =
		integrate_kvec_adiabatic (phi, phi_adb, &initconds, &params,
					  t, t + t_inc);

	      //Calculate local defect density
	      GSL_SET_COMPLEX (&uk, phi[0], phi[2]);
	      GSL_SET_COMPLEX (&vk, phi[1], phi[3]);

	      if (gndstate == 'd')
		{
		  set_initconds_diabatic_gnd (phi_init);	//set initconds for a particular k; diabatic ground state
		  //phi_unit[]=phi_init[];
		  phi_unit[0] = phi_init[0];
		  phi_unit[1] = phi_init[1];
		  phi_unit[2] = phi_init[2];
		  phi_unit[3] = phi_init[3];
		}
	      else if(gndstate == 's')
		{
	      set_initconds_bcs_gnd (&params, &initconds, phi_init);	//set initconds for a particular k; BCS ground state
		  //phi_unit[]={1.0,0.0,0.0,0.0};
		  phi_unit[0] = 0.0;
		  phi_unit[1] = 1.0;
		  phi_unit[2] = 0.0;
		  phi_unit[3] = 0.0;	
		}
	      else
	      {
	      rk = rkarray[p + gridsize * kycount];
	      set_initconds_rand_gnd (&params, &initconds, phi_init, rk);	//set initconds for a particular k; random ground state
	      
		  //phi_unit[]=phi_init[];
		  phi_unit[0] = phi_init[0];
		  phi_unit[1] = phi_init[1];
		  phi_unit[2] = phi_init[2];
		  phi_unit[3] = phi_init[3];	      
	      }
	      
	      GSL_SET_COMPLEX (&vkinit, phi_init[1], phi_init[3]);
	      uksq = gsl_complex_abs2 (uk);
	      vksq = gsl_complex_abs2 (vk);
	     
	      adiabatic_trans_mbt(phi_adbmag,phi_unit,t,&params);
	      GSL_SET_COMPLEX(&ukadb,phi_adbmag[0],phi_adbmag[2]);
	      GSL_SET_COMPLEX(&vkadb,phi_adbmag[1],phi_adbmag[3]);
	      vkadb_abs = gsl_complex_abs2(vkadb);
	     
	      uv = gsl_complex_mul(ukadb,vk);
	      vu = gsl_complex_mul(vkadb,uk);
	      uvdiff = gsl_complex_sub(uv,vu);
	      
	      defect_density = gsl_complex_abs2(uvdiff);
	      defect_density_local = defect_density_local + defect_density;
	      
	      //Set adiabatic vk for this particular delta
	      adbmag_local = adbmag_local + vkadb_abs;

	      norm_var = uksq + vksq;
	      norm_local = norm_local + norm_var;

	      energy_k = energy_expectation (t, phi, &params);
	      energy_k = energy_k - energy_expectation (t, phi_adbmag, &params);
	      //Calculate local energy expectation
	      ek_local = ek_local + energy_k;
	      //Calculate local fidelity susceptibility
	      fidk_dia = inner_product_modsq (phi_init, phi);
	      fidk_dia = log (fidk_dia);
	      fidelity_susceptibility_local =
		fidelity_susceptibility_local + fidk_dia;

         //THIS MAY NEED TO BE UPDATED WITH ANISOTROPY pk
	     //Evaluate locally the function \sum_k uk* x vk x exp(i x fk x t)
	     //and update to gfunc
		fk = params.ampmax - params.mu0;
		pk = gsl_complex_polar(1.0,2.0 * fk * t);
		gk = gsl_complex_conjugate(uk);
		gk = gsl_complex_mul(gk,vk);
		gk = gsl_complex_mul(gk,pk);
		gfunc_local = gsl_complex_add(gfunc_local,gk);
		
	      //copy new phi to old phik
	      copyfrombuffer (phi, 4, phik, kstride);

	      if (result == GSL_FAILURE)
		{
		  lzs_intgerr (result);
		  break;
		}
	      probloc = probloc + phi[1] * phi[1] + phi[3] * phi[3];
	      probloc_adb = phi_adb[1] * phi_adb[1] + phi_adb[3] * phi_adb[3];
	      fidelity_local = fidelity_local + log (probloc_adb);

	      //append t, k and phik to time t in output buffer
	      if (stateout == 'y')
		{
		  if((kx>=0.0) && (ky>=0.0))
		    if(fabs(fk)<=fabs(params.omega))
		      {
			wt = omega * t / (2.0 * M_PI);
			g_array_append_val (time_phik, wt);
			g_array_append_val (time_phik, kx);
			g_array_append_val (time_phik, ky);
			g_array_append_vals (time_phik, phi, 4);
			g_array_append_val (time_phik, defect_density);
			g_array_append_val (time_phik, probloc_adb);
		      }
		}
	      probcount_loc++;
	      kstride = kstride + 5;
	      kycount++;
	    }
	}
      //Now, probloc is the sum of all probabilities at time t for this process
      //probcount_loc is the number of k's taken at time t for this process
      //So evaluate the total probability and total number by MPI_Reduce and dum out the average
      prob = 0.0;
      probcount = 0;

      gfunc_loc_real = GSL_REAL (gfunc_local);
      gfunc_loc_imag = GSL_IMAG (gfunc_local);
      
      MPI_Barrier (MPI_COMM_WORLD);
      MPI_Reduce (&probloc, &prob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce (&ek_local, &ek_global, 1, MPI_DOUBLE, MPI_SUM, 0,
		  MPI_COMM_WORLD);
      MPI_Reduce (&defect_density_local, &defect_density_global, 1,
		  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce (&norm_local, &norm_global, 1,
		  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce (&fidelity_local, &fidelity_global, 1, MPI_DOUBLE, MPI_SUM,
		  0, MPI_COMM_WORLD);
      MPI_Reduce (&fidelity_susceptibility_local,
		  &fidelity_susceptibility_global, 1, MPI_DOUBLE, MPI_SUM, 0,
		  MPI_COMM_WORLD);
      MPI_Reduce(&adbmag_local,&adbmag_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&gfunc_loc_real,&gfunc_glob_real,1,MPI_DOUBLE,MPI_SUM,0,
		MPI_COMM_WORLD);
      MPI_Reduce(&gfunc_loc_imag,&gfunc_glob_imag,1,MPI_DOUBLE,MPI_SUM,0,
		MPI_COMM_WORLD);      
      MPI_Reduce (&probcount_loc, &probcount, 1, MPI_LONG, MPI_SUM, 0,
		  MPI_COMM_WORLD);
      if (pid == 0)
	{
	  g_array_append_val (time_magnetization, t);
	  g_array_append_val (time_fidelity, t);
	  g_array_append_val (time_delta, t);
	  g_array_append_val (time_defectdensity, t);
	  g_array_append_val (time_residualenergy, t);
	  g_array_append_val (time_fidsuscpt, t);
	  g_array_append_val (time_adbmag,t);
	  g_array_append_val (time_gfunc,t);
	  
	  adbmag_global = adbmag_global/probcount;
	  probavg = prob / probcount;

	  //Calculate defect density
	  defect_density_global = defect_density_global / probcount;
	  norm_global = norm_global / probcount;
	  g_array_append_val (time_defectdensity, defect_density_global);
	  

	  //Calculate magnetizations
	  adbmag_global = 2.0 * adbmag_global -1.0;
	  adbmag_global = -adbmag_global;
	  probavg = 2.0 * probavg - 1.0;

	  //Calculate fidelity
	  fidelity_global = fidelity_global / probcount;
	  fidelity_global = exp (fidelity_global);

	  //calculate gap
	  deltaglob_abs =
	    params.delta_real * params.delta_real +
	    params.delta_imag * params.delta_imag;
	  deltaglob_abs = sqrt (deltaglob_abs);

	  //Calculate residual energy
	  ek_global = ek_global / probcount;
	  //ek_global = ek_global - ekinit_global;
	  //This is for the old defn of residual energy
	  
	  ek_global = fabs (ek_global);	//Just to keep it positive 

	  //calculate fidelity susceptibility
	  fidelity_susceptibility_global =
	    fidelity_susceptibility_global / probcount;
	  fidelity_susceptibility_global =
	    exp (fidelity_susceptibility_global);

	  gfunc_glob_real = params.veff * (gfunc_glob_real/probcount);
	  gfunc_glob_imag = params.veff * (gfunc_glob_imag/probcount);

	  //Store it all in memory
	  g_array_append_val (time_magnetization, probavg);
	  g_array_append_val (time_fidelity, fidelity_global);
	  
	  g_array_append_val (time_delta, deltaglob_abs);

	  g_array_append_val (time_residualenergy, ek_global);
	  g_array_append_val (time_fidsuscpt, fidelity_susceptibility_global);
	  g_array_append_val (time_adbmag,adbmag_global);
	  g_array_append_val (time_gfunc, gfunc_glob_real);
	  g_array_append_val (time_gfunc, gfunc_glob_imag);
	}

      scaled_time = omega * t / (2.0 * M_PI);
      if ((compare (scaled_time, 1.0, strobeoff) == GSL_SUCCESS)
	  && strobecount == 0)
	{
	  if (pid == 0)
	    {
	      dump_out_strobe (scaled_time, probavg, fidelity_global,
			       deltaglob_abs, defect_density_global,
			       ek_global, fidelity_susceptibility_global);
	      strobecount = 1;
	    }
	}
    }
  (void) time (&end);
  MPI_Barrier (MPI_COMM_WORLD);

  //Dump time_phik in parallel io to putput buffer in binary form
  //Get GArray size
  long time_phik_size = time_phik->len;
  //Let each process broadcast it's data size to all other processes. Store in array datasizes
  long *datasizes;
  datasizes = malloc (nprocs * sizeof (*datasizes));
  if (stateout == 'y')
    {
      ierr =
	MPI_Allgather (&time_phik_size, 1, MPI_LONG, datasizes, 1, MPI_LONG,
		       MPI_COMM_WORLD);
    }
  //Calculate displacement of each process in file based on displacements of all other processes 
  //of lesser rank
  long displacement = 0, cnt;
  if (pid != 0)
    for (cnt = 0; cnt <= pid; cnt++)
      {
	displacement = displacement + datasizes[cnt];
      }
  //Create file view for each process now
  if (stateout == 'y')
    {
      MPI_File_set_view (tphikfile, displacement * sizeof (double),
			 MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
      //Each prcess now writes to file
      MPI_File_write (tphikfile, ((double *) (void *) time_phik->data),
		      time_phik_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
    }
  MPI_Barrier (MPI_COMM_WORLD);
  if (pid == 0)
    {
      //Compute the derivative of fidelity from GArrays time_fidelity
      compute_deriv(time_fidelity, time_fidelity_deriv_mag);
      //Final output
      dump_out_probs (time_magnetization, tmagfile, &params);
      dump_out_probs (time_fidelity, tfidfile, &params);
      dump_out_probs (time_delta, deltafile, &params);
      dump_out_probs (time_defectdensity, ndfile, &params);
      dump_out_probs (time_residualenergy, tresenfile, &params);
      dump_out_probs (time_fidsuscpt, tfidsfile, &params);
      dump_out_probs (time_adbmag,tmagadbfile,&params);
      dump_out_probs (time_fidelity_deriv_mag, tfidderivfile, &params);
      dump_out_delta (time_gfunc,tgfuncfile,&params);
      tau = compute_nonzero (time_fidelity_deriv_mag,&params);
      dump (reltime,tau);

      qf.probavg = get_qfactor (time_magnetization);
      qf.fidavg = get_qfactor (time_fidelity);
      qf.deltavg = get_qfactor (time_delta);
      qf.ndavg = get_qfactor (time_defectdensity);
      qf.ekavg = get_qfactor (time_residualenergy);
      output_stdev (time_fidelity_deriv_mag);

      final_out_mbt (begin, end, gridsize, &qf);
    }

  g_array_free (phik, TRUE);
  g_array_free (time_phik, TRUE);
  g_array_free (time_magnetization, TRUE);
  g_array_free (time_fidelity, TRUE);
  g_array_free (time_delta, TRUE);
  g_array_free (time_defectdensity, TRUE);
  g_array_free (time_residualenergy, TRUE);
  g_array_free (time_fidsuscpt, TRUE);
  g_array_free (time_adbmag,TRUE);
  g_array_free (time_fidelity_deriv_mag, TRUE);
  g_array_free (time_gfunc, TRUE);

  free (rkarray);
  free (datasizes);

  fclose (tmagfile);
  fclose (tfidfile);
  fclose (deltafile);
  fclose (ndfile);
  fclose (tresenfile);
  fclose (tfidsfile);
  fclose (tmagadbfile);
  fclose (tfidderivfile);
  fclose (tgfuncfile);
  
  MPI_File_close (&tphikfile);

  MPI_Finalize ();

  return 0;
}
