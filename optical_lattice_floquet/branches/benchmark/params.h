/*
// C++ Interface: params
//
// Description:
//
//
// Author: Analabha Roy <daneel@physics.utexas.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
*/

#define N 200
#define DIMS (2*N+1)
/* Truncsation of the hamiltonian basis is N, then DIM=2N+1 (-N to +N). Problem is embarassingly parallelized     						//across DIM processes so DON'T OVERSUBSCRIBE YOUR NODES*/

/*Global numerical tolerance*/
#define MACHINENUM 1E-5
/*Maximum number of nonunitary quasienergies*/
#define QUASIMAX 2
/********The Parameter epsilon & the tower of 2-particle states are stored in this datatype********/
typedef struct
{
  double kappa;
  double lambda;
  double omega;
  int dim;
  double psq[DIMS];
  double cosmatrix[DIMS * DIMS];
} drive_and_tower;
