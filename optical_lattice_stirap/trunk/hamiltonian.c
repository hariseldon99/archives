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
hamilt_sp (void *param, int nt, int mt, double t)
{
  double matrixelm;
  paramspace_pt *p = (paramspace_pt *) param;
  /*Get epsilon & towerofstates from structure param */
  double kappa = p->kappa;
  double lambda0 = p->lambda0;
  double wf = p->wf;
  double ws = p->ws;
  double tf = p->tf;
  double ts = p->ts;
  double td = p->td;
  double psqmatrix, cosmatrix, sinmatrix;
  double t_diff, t_diff_sq, td_sq, ratio;
  double lambda_f, lambda_s;
  psqmatrix = p->psq[nt];
  cosmatrix = p->cosmatrix[nt + DIM * mt];
  sinmatrix = p->sinmatrix[nt + DIM * mt];

  /*Calculate STIRAP Amplitudes */

  t_diff = (t - tf);
  t_diff_sq = t_diff * t_diff;
  td_sq = 4.0 * td * td;
  ratio = t_diff_sq / td_sq;
  lambda_f = exp (-ratio);


  t_diff = (t - ts);
  t_diff_sq = t_diff * t_diff;
  td_sq = 4.0 * td * td;
  ratio = t_diff_sq / td_sq;
  lambda_s = exp (-ratio);
  /*Write out the matrix element in the hamiltonian representation */
  matrixelm = cosmatrix;
  matrixelm = matrixelm * kappa;

  matrixelm =
    matrixelm + 2.0 * lambda0 * cosmatrix * (lambda_f * cos (wf * t) +
					     lambda_s * cos (ws * t));

  matrixelm = matrixelm + psqmatrix * kdel (mt, nt);

  return (matrixelm);
}

/*Time derivative of hamiltonian*/
double
hamilt_deriv_sp (void *param, int nt, int mt, double t)
{
  double matrixelm;
  paramspace_pt *p = (paramspace_pt *) param;
  /*Get epsilon & towerofstates from structure param */
  double kappa = p->kappa;
  double lambda0 = p->lambda0;
  double wf = p->wf;
  double ws = p->ws;
  double tf = p->tf;
  double ts = p->ts;
  double td = p->td;
  double psqmatrix, cosmatrix, sinmatrix;
  double t_diff, t_diff_sq, td_sq, ratio;
  double lambda_f, lambda_s;
  psqmatrix = p->psq[nt];
  cosmatrix = p->cosmatrix[nt + DIM * mt];
  sinmatrix = p->sinmatrix[nt + DIM * mt];


  /*Calculate STIRAP Amplitudes */

  t_diff = (t - tf);
  t_diff_sq = t_diff * t_diff;
  td_sq = 4.0 * td * td;
  ratio = t_diff_sq / td_sq;
  lambda_f = exp (-ratio);


  t_diff = (t - ts);
  t_diff_sq = t_diff * t_diff;
  td_sq = 4.0 * td * td;
  ratio = t_diff_sq / td_sq;
  lambda_s = exp (-ratio);

  /*Write out the matrix element in the hamiltonian representation */
  matrixelm = cosmatrix;
  matrixelm = matrixelm * kappa;
  matrixelm =
    matrixelm - 2.0 * lambda0 * cosmatrix * (lambda_f * wf * sin (wf * t) +
					     lambda_s * ws * sin (ws * t));

  matrixelm = matrixelm + psqmatrix * kdel (mt, nt);

  return (matrixelm);
}

double
hamilt (void *param, int nt, int mt, double t)
{
  double matrixelm;
  paramspace_pt *p = (paramspace_pt *) param;
  int n1, n2, m1, m2;

  /*Get epsilon & towerofstates from structure param */
  n1 = p->towerofstates[nt][0];
  n2 = p->towerofstates[nt][1];
  m1 = p->towerofstates[mt][0];
  m2 = p->towerofstates[mt][1];
  matrixelm = p->interaction[nt + STATENO * mt];
  matrixelm = matrixelm + hamilt_sp (p, n1, m1, t) * kdel (n2, m2);
  matrixelm = matrixelm + hamilt_sp (p, n1, m2, t) * kdel (n2, m1);
  matrixelm = matrixelm + hamilt_sp (p, n2, m1, t) * kdel (n1, m2);
  matrixelm = matrixelm + hamilt_sp (p, n2, m2, t) * kdel (n1, m1);
  //Normalization
  matrixelm = matrixelm / sqrt (1 + kdel (n1, n2));
  matrixelm = matrixelm / sqrt (1 + kdel (m1, m2));

  return (matrixelm);
}

double
hamilt_deriv (void *param, int nt, int mt, double t)
{
  double matrixelm;
  paramspace_pt *p = (paramspace_pt *) param;
  int n1, n2, m1, m2;
  /*Get epsilon & towerofstates from structure param */
  n1 = p->towerofstates[nt][0];
  n2 = p->towerofstates[nt][1];
  m1 = p->towerofstates[mt][0];
  m2 = p->towerofstates[mt][1];
  matrixelm = p->interaction[mt + STATENO * nt];
  matrixelm = matrixelm + hamilt_deriv_sp (p, n1, m1, t) * kdel (n2, m2);
  matrixelm = matrixelm + hamilt_deriv_sp (p, n1, m2, t) * kdel (n2, m1);
  matrixelm = matrixelm + hamilt_deriv_sp (p, n2, m1, t) * kdel (n1, m2);
  matrixelm = matrixelm + hamilt_deriv_sp (p, n2, m2, t) * kdel (n1, m1);
  //Normalization
  matrixelm = matrixelm / sqrt (1 + kdel (n1, n2));
  matrixelm = matrixelm / sqrt (1 + kdel (m1, m2));
  return (matrixelm);
}
