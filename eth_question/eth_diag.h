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
  "Diagonalization of the lattice in Nature doi:10.1038/nature06838\n"
  "\tFirst, it inputs the Hamiltonian froma  PEtsc binary output,\n"
  "\tthen performs symmetric eigenproblem using slepc.\n"
  "\tAfter that, it diagonalizes in parallel using slepc and dumps output in binary using petsc.\n\n"
  "The command line options are:\n"
  "\t-draw_esys\t\t      Shows matrix plots of the Hamiltonian, its eigenvectors and a listplot of eigenvalues.\n"
  "\t-h\t\t\t      Prints this help\n\n";

//Chooses which default eigensolver to use for matrix diagonalization. 
//Best to keep this at lapack
//Can be overridden at runtime. Consult README for details
#define PROBLEMTYPE EPSLAPACK

#include <slepceps.h>
#include <string.h>		/* for strings */
