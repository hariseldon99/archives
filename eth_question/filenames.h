/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Rigol Lattice: filenames
      Diagonalization of the lattice in
      Rigol et al, Nature 452, 854 (2008), doi:10.1038/nature06838
      for any given lattice size.
   Copyright (c) 2013 Analabha Roy (daneel@utexas.edu)

   This is free software: you can redistribute it and/or modify it under  the
   terms of version 3 of the GNU Lesser General Public License as published by
   the Free Software Foundation.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/
static char hamilt_fnamesbin[2][PETSC_MAX_PATH_LEN] =
  { "hamilt.dat", "hamilt_lower.dat" };
static char filenames_bin[3][PETSC_MAX_PATH_LEN] =
  { "evals.dat", "evecs.dat", "initstate.dat" };
static char filenames_eth_bin_davg[3][PETSC_MAX_PATH_LEN] =
  { "hopmat_davg.dat", "nk_davg.dat", "C_a-sq.dat" };
static char filenames_eth_bin_eavg[2][PETSC_MAX_PATH_LEN] =
  { "hopmat_eavg.dat", "nk_eavg.dat" };
