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

//paramspace_pt_serial is a datatype that stores the system parameters, namely the gap 'delta' and the driving frequency 'omega'
typedef struct
{
  double gridsize;
  double omega;
  double mu0;
  double muamp;
  double delta_real;
  double delta_imag;
  double veff;
} paramspace_pt_serial;
