/*
*  C Implementation: hamiltonian
*
* Description:
*
*
* Author: Analabha Roy <daneel@physics.utexas.edu>, (C) 2007
*
* Copyright: See COPYING file that comes with this distribution
*
*/
/*This uses the above subroutines to evaluate the hamiltonian matrix elements*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "params.h"
#include "hamiltonian.h"

double
hamilt (void *param, int nt, int mt, double t)
{
  double matrixelm;
  drive_and_tower *p = (drive_and_tower *) param;
  /*Get epsilon & towerofstates from structure param */
  double kappa = p->kappa;
  double lambda = p->lambda;
  double omega = p->omega;
  int dim=p->dim;
  double psqmatrix, cosmatrix;
  psqmatrix = p->psq[nt];
  cosmatrix = p->cosmatrix[nt + dim * mt];

  /*Write out the matrix element in the hamiltonian representation */
  matrixelm = 2.0 * lambda * cos (omega * t);
  matrixelm = matrixelm + kappa;
  matrixelm = matrixelm * cosmatrix;
  matrixelm = matrixelm + psqmatrix * kdel (mt, nt);
  return (matrixelm);
}

/*Time derivative of hamiltonian*/
double
hamilt_deriv (void *param, int nt, int mt, double t)
{
 double matrixelm;
 drive_and_tower *p = (drive_and_tower *) param;
 /*Get epsilon & towerofstates from structure param */
 double kappa = p->kappa;
 double lambda = p->lambda;
 double omega = p->omega;
 int dim=p->dim;
 double psqmatrix, cosmatrix;
 psqmatrix = p->psq[nt];
 cosmatrix = p->cosmatrix[nt + dim * mt];
	
 /*Write out the matrix element in the hamiltonian representation */
 matrixelm = -2.0 * lambda * omega * sin (omega * t);
 matrixelm = matrixelm + kappa;
 matrixelm = matrixelm * cosmatrix;
 matrixelm = matrixelm + psqmatrix * kdel (mt, nt);
 return (matrixelm);
}
