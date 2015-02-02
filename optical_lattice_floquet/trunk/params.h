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
/* Truncation of the s.p. hamiltonian basis is N, then DIM=2N+1 (-N to +N). Problem is embarassingly parallelized     						//across DIM processes so DON'T OVERSUBSCRIBE YOUR NODES*/
#define N 1
#define DIM (2*N+1)
#define STATENO (DIM*(DIM+1))/2

/*Period size of lattice*/
#define LATTICESIZE 2
/*Global numerical tolerance*/
#define MACHINENUM 1E-5
/*Maximum number of nonunitary quasienergies*/
#define QUASIMAX 8

/********The Parameter epsilon & the tower of 2-particle states are stored in this "parameter space point" datatype********/
typedef struct {
        double kappa;
        double lambda0;
        double wf;
        double ws;
        double u0;
        double tfix;
        double tf;
        double ts;
        double td;
        double psq[DIM];
        double cosmatrix[DIM * DIM];
        double sinmatrix[DIM * DIM];
        double interaction[STATENO * STATENO];
        int towerofstates[STATENO][2];

} paramspace_pt;
