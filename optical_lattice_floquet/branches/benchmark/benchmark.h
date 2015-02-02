/*
// C++ Interface: floquet_optical_lattice
//
// Description:
//
//
// Author: Analabha Roy <daneel@physics.utexas.edu>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//
*/
double gtod_timer();
double c_mat (int m, int n), psq_mat (int n), kdel (int i, int j);
int integrate (double *y, double initial, double final, void *out);
double hamilt (void *param, int nt, int mt, double t);

