/*
 T his *program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor,
 Boston, MA  02110-1301, USA.

 ---
 Copyright (C) 2011, Analabha Roy <daneel@physics.utexas.edu>
 */

#include "lattice.h"
#include "params_ics.h"
//Time dependence
//Energy
double
epsilon (double ampmax, double omega, double t, double offset)
{
  double amplitude;
  amplitude = ampmax * sin (omega * t);
  amplitude = amplitude + offset;
  return amplitude;

}

//Chemical potential
double
mu (double mu0, double muamp, double omega, double t, double offset)
{
  double amplitude;
  amplitude = mu0 - muamp * sin (omega * t) + offset;
  return amplitude;
}

//time derivative of above
double
dmudt (double mu0, double muamp, double omega, double t, double offset)
{
  double amplitude;
  offset = 0.0;			//redundant
  amplitude = -muamp * omega * cos (omega * t);
  return amplitude;
}

//Delta = veff * sum_k uk*   vk
//Returns the number of ks in array

long
update_delta_loc (GArray * phik, void *params)
{

  paramspace_pt *p = (paramspace_pt *) params;
  double veff = p->veff;
  //Here, phik is stored as (k,phi[0], phi[1],phi[2],phi[3]...)
  long kstride = 0, kcount = 1;
  long phiksize = phik->len;
  double phi[4];
  gsl_complex delta_loc = gsl_complex_rect (0.0, 0.0);
  gsl_complex uk, vk, rhs;

  GSL_SET_COMPLEX (&delta_loc, 0.0, 0.0);

  while (kstride < phiksize)
    {
      kstride++;
      phi[0] = g_array_index (phik, double, kstride);
      kstride++;
      phi[1] = g_array_index (phik, double, kstride);
      kstride++;
      phi[2] = g_array_index (phik, double, kstride);
      kstride++;
      phi[3] = g_array_index (phik, double, kstride);
      kstride++;

      GSL_SET_COMPLEX (&uk, phi[0], phi[2]);
      GSL_SET_COMPLEX (&vk, phi[1], phi[3]);

      rhs = gsl_complex_conjugate (uk);
      rhs = gsl_complex_mul (rhs, vk);

      delta_loc = gsl_complex_add (delta_loc, rhs);
      if (kstride < phiksize)
	{
	  kcount++;
	}

    }
  p->delta_real = veff * GSL_REAL (delta_loc);
  p->delta_imag = veff * GSL_IMAG (delta_loc);

  return kcount;

}



int
fermisurface_2d (double kx, double ky, double ef, double spread)
{

  int result = GSL_FAILURE;
  double energy = cos (kx) + cos (ky);
  double off = fabs (energy - ef);
  if (off <= spread)
    {
      result = GSL_SUCCESS;
    }
  return result;
}
