/*
// C++ Interface: params
//
// Description:
//
//
// Author: Analabha Roy <daneel@utexas.edu>, (C) 2013
//
// Copyright: See COPYING file that comes with this distribution
//
*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * Hardcode the lattice size here
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#define N 10

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * Hardcode the number of random instances PER MPI PROCESS here
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#define NRAND 1

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * Set this macro if you want to set initial conditions to the GHZ 
 * (Greenberger-Horne-Zeilinger) state. The default is the classical 
 * ground state i.e. all spins up
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//#define INITCONDS_GHZ

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * Set this macro if you want to switch off the field disorder
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
//#define FIELD_DISORDER_OFF

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * Set this macro if you want to set periodic boundary conditions. The default is
 * open boundary conditions
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#define PERIODIC_BC



#define DIM 4*N*N

/* This is the kronecker delta*/
#define kdel(i,j) ((i==j)?1.0:0.0)

/*Global numerical tolerance*/
#define MACHINENUM 1E-5

/*Progress bar width and resolution*/
#define BAR_RES 30

typedef struct
{
  double gamma_a;		//Amplitude of Gamma(t)
  int zero_order;		//Order of zero of bessel function
  double eta;			//Zero of Bessel Function at order zero_order
  double omega;			//Frequency of Gamma(t)
  double alpha;			//Perturbation
  double hrand[N];		//Array of N random numbers for the magnetic field
  double jrand[N];		//Array of N random numbers for the hopping
} paramspace_pt;
