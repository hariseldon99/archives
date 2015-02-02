#include "params_ics.h"
#include "lzs_tls.h"
#undef __FUNCT__
#define __FUNCT__ "main"

int
main (int argc, char **argv)
{
  time_t begin, end;
  paramspace_pt params;
  initconds_ic initconds;
  stats normstats;
  int i, argv_count = 1, result;

  double t_init, t_periods, t_extra;
  double t_init_frac, periods, t_extra_frac;
  double phi[4], omega, phi_norm;
  double normexpect = 1.0;	//Expected norm
  GArray *data = g_array_new (FALSE, FALSE, sizeof (double));
  GArray *data_probs = g_array_new (FALSE, FALSE, sizeof (double));
  GArray *data_probs_prime = g_array_new (FALSE, FALSE, sizeof (double));


  FILE *lzs_tls_phi, *lzs_tls_phi_probs, *lzs_tls_phi_probs_prime;

  //argument error checking
  if (argc != ARGNO + 1)
    {
      int aerr = arg_err ();
      exit (aerr);
    }

  lzs_tls_phi = fopen (argv[argv_count++], "w");	//input filenames first
  lzs_tls_phi_probs = fopen (argv[argv_count++], "w");
  lzs_tls_phi_probs_prime = fopen (argv[argv_count++], "w");

  //Input initial conditions
  for (i = 0; i < 4; i++)
    {
      initconds.phi_in[i] = atof (argv[argv_count++]);	//input initial phi
      phi[i] = initconds.phi_in[i];
    }

  t_init_frac = atof (argv[argv_count++]);	//initial time_t in fractions of the period
  periods = atof (argv[argv_count++]);	//Number of period sweeps

  t_extra_frac = atof (argv[argv_count++]);	//Extra time after n sweeps in fractions of the period

  //Input system parameters
  params.delta_real = atof (argv[argv_count++]);
  params.delta_imag = atof (argv[argv_count++]);

  params.omega = atof (argv[argv_count++]);
  params.offset = atof (argv[argv_count++]);
  params.ampmax = atof (argv[argv_count++]);;
  if (strcmp (argv[argv_count++], "y") == 0)
    {				//whether output is verbose or not, answer y or n
      initconds.verbosity = 'y';
    }
  else
    {
      initconds.verbosity = 'n';
    }

  omega = params.omega;
  params.gndstate = 'b';



  //Calculate t_init

  t_init = 2 * t_init_frac * M_PI / omega;
  initconds.t_init = t_init;

  //Calculate t_extra
  t_extra = 2 * t_extra_frac * M_PI / omega;
  initconds.t_extra = t_extra;

  //Calculate t_periods
  t_periods = 2.0 * periods * M_PI / omega;
  initconds.t_periods = t_periods;

  //Initialize phi data
  t_init = initconds.t_init;
  g_array_append_val (data, t_init);
  g_array_append_vals (data, &phi, 4);
  phi_norm = norm (phi);
  g_array_append_val (data, phi_norm);

  //Report initial conditions
  init_out (&initconds);
  //Report parameters inputted
  params_out (&params);

  //Call subroutine for integrating the system
  (void) time (&begin);
  result =
    integrate (phi, &initconds, &params, data, data_probs, data_probs_prime);
  if (result == GSL_FAILURE)
    {
      lzs_intgerr (result);
    }
  (void) time (&end);

  //Dump output to file
  dump_out (data, lzs_tls_phi);
  dump_out_probs (data_probs, lzs_tls_phi_probs, &params);
  dump_out_probs (data_probs_prime, lzs_tls_phi_probs_prime, &params);
  //Compute norm statistics
  norm_stats (&normstats, data, normexpect);

  //Report final output and time
  final_out (begin, end, normexpect, &normstats);

  fclose (lzs_tls_phi);
  fclose (lzs_tls_phi_probs);
  fclose (lzs_tls_phi_probs_prime);

  g_array_free (data, FALSE);
  g_array_free (data_probs, FALSE);
  g_array_free (data_probs_prime, FALSE);
  return 0;
}
