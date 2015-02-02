/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Rigol Lattice: 
      Diagonalization of the lattice in
      Rigol et al, Nature 452, 854 (2008), doi:10.1038/nature06838
      for any given lattice size.
   Copyright (c) 2013 Analabha Roy (daneel@utexas.edu)

   This is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

static char help[] =
  "Formulation of the lattice Hamiltonian in Nature doi:10.1038/nature06838\n"
  "\tFirst, it builds the full Hilbert space canonical basis excluding forbidden cells,\n"
  "\tThen it formulates the Hamiltonian matrix elements in parallel using petsc and dumps to binary"
  "The command line options are:\n"
  "\t-lattice_size <n>\t<n> = ODD number of sites (defaults to 5). Lattice dimension is n^2.\n"
  "\t-vector_size <n>\t<n> = Size of vector i.e. number of hard sphere bosons (defaults to 2).\n"
  "\t-repulsion <U>\t\t<U> = Repulsive energy (defaults to 0.1) between occupied neighboring sites.\n"
  "\t-draw_mats\t\t      Shows matrix plots of the Hamiltonian, its eigenvectors and a listplot of eigenvalues.\n"
  "\t-h\t\t\t      Prints this help\n\n";

//Default values of lattice size and repulsive energy
#define DEFAULTLATSIZE 5
#define DEFAULTSIZE 2
#define DEFAULTREP 0.1

//Chooses which default eigensolver to use for matrix diagonalization. 
//Best to keep this at lapack
//Can be overridden at runtime. Consult README for details
#define PROBLEMTYPE EPSLAPACK


#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscviewerhdf5.h>	//Will need this for output

#include <gsl/gsl_combination.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort_long.h>
#include <string.h>		/* for strings */

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
