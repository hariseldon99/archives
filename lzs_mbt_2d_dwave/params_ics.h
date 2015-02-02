/*
// C++ Interface: 
//
// Description:
//
//
// Author: Analabha Roy <daneel@bose.res.in>, (C) 2010
//
//
*/

//paramspace_pt is a datatype that stores the system parameters, namely the gap 'delta' and the driving frequency 'omega'
typedef struct
{
  double delta_real;
  double delta_imag;
  double omega;
  double offset;
  double ampmax;
  double pk;
  double mu0;
  double muamp;
  double veff;
  char stateout;
  char gndstate;
} paramspace_pt;


//initconds_ic is a datatype that stores the initial conditions
typedef struct
{
  double phi_in[4];		//The state is phi=(uk_real, vk_real,uk_imag,vk_imag)
  double t_init;
  double t_periods;
  double t_extra;		//total time is calculated as (2*periods*pi/omega)+t_extra
  char verbosity;
} initconds_ic;

typedef struct
{
  double mean;
  double stddev;
  double max;
  double min;
} stats;

typedef struct
{
  double probavg;		//average magnetization
  double fidavg;		//average fidelity
  double deltavg;		//average gap
  double ndavg;			//average defect density
  double ekavg;			//average residual energy
} qfactors;

//Non-verbose output. Verbose output needed only for lzs_tls
#define NOT_VERBOSE 'n'
#define VERBOSE 'y'

//Define sread fractn of delta about 2d Fermi surface
#define SPREAD_FRAC 10.0
