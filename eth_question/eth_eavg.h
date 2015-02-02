/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Eigenstate Thermalization Hypothesis (eth): 
      Computation of ensemble averages of the lattice in
      Rigol et al, Nature 452, 854 (2008), doi:10.1038/nature06838
      for any given lattice size. Lattice eigensystem obtained from rigol_lattice.c
   Copyright (c) 2013 Analabha Roy (daneel@utexas.edu)

   This is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
/*
 Note: NULL pointer means null vector
 */

static char help[] =
  "Computation of diagonal averages of lattice in Nature doi:10.1038/nature06838\n"
  "First, run the program that diagonalizes the Hamiltonian with parameters entered as below.\n"
  "Then, save the output in binary format ONLY and run this program with the options below:\n"
  "\t-lattice_size <n>\t<n>  =  ODD number of sites (defaults to 5). Lattice dimension is n^2.\n"
  "\t-vector_size <n>\t<n>  =  Size of vector i.e. number of hard sphere bosons (defaults to 2).\n"
  "\t-repulsion <U>\t\t<U>  =  Repulsive energy (defaults to 0.1) between occupied neighboring sites.\n"
  "\t-init_e <E>\t\t<E>  =  Energy of init state (defaults to the value from system defaults).\n"
  "\t-delta_e <dE>\t\t<dE> =  Energy shell (defaults to 0.1) around the init state.\n"
  "\t-ksize <nk>\t\t<nk> =  Number of momenta to compute in any given axis (defaults to hopmatrix size).\n"
  "\t-draw_out\t\t\tShows matrix plots of the Hamiltonian, its eigenvectors and a listplot of eigenvalues.\n"
  "\t-left_avg\t\t\tPerforms left average instead of the default right average.\n"
  "\t-h\t\t\t\tPrints this help\n\n";

//Default values of lattice size and repulsive energy
#define DEFAULTLATSIZE 5
#define DEFAULTSIZE 2
#define DEFAULTREP 0.1


#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscviewerhdf5.h>	//Will need this for output

#include <gsl/gsl_combination.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort_long.h>
#include <string.h>		/* for strings */

#define DEFAULTSHELL 0.1
#define INITENERGYDEFAULT -5.360044

void nullify (long vec[], const long vecsize);
int isnull (const long vec[], const long vecsize);
void copyvec (const long src[], long dest[], long vecsize);
long find_index (long vec[], long vecsize, long value);
int
is_forbidden (const size_t data[], const long vecsize, const long latsize);
int
is_outside_lower_right_quarter (const long vec[], const long vecsize,
				const long latsize);
double can_sdot (const long vec1[], const long vec2[], long vecsize);
int are_nn (const long i, const long j, const long latsize);
void
bdaggerb (const long i, const long j, long vec[], const long vecsize,
	  const long latsize);
void nsquared (const long i, const long j, long vec[], const long vecsize);
