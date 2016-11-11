#!/usr/bin/env python
from __future__ import division, print_function

"""
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * Copyright (c) 2015 Analabha Roy (daneel@utexas.edu)
 *
 * This is free software: you can redistribute it and/or modify it under
 * the terms of version 2 of the GNU Lesser General Public License
 * as published by the Free Software Foundation.
 * Notes:
 * 1. The initial state is currently hard coded to be the classical ground
 *     state
 * 2. Primary references are
 *    Anatoli: Ann. Phys 325 (2010) 1790-1852
 *    Schachemmayer: arXiv:1408.4441
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""
desc = """Comparison of quantum analytical with Discrete Truncated Wigner
          Approximation for 1d Ising model (sxsx) with long range
          interactions and no fields"""

import numpy as np
import argparse, sys
from pprint import pprint
from scipy.sparse import *

# Hardcoded input Parameters. Default values
beta = 1.0
N = 10 # Lattice Size
t_init = 0.0 # Initial time
n_cycles = 60.0 # Fractions of the drive period for running
n_steps = 1000 # Number of time steps

def input():
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-omx', '--output_magx',\
      help="sx (magnetization) output file", default="sx_outfile.txt")

    parser.add_argument('-ox', '--output_sxvar',\
      help="sx variance output file", default="sxvar_outfile.txt")

    parser.add_argument("-v", '--verbose', action='store_true', \
      help="increase output verbosity")

    parser.add_argument('-pbc', '--periodic', \
      help="switch to periodic boundary conditions, default is open",\
        action="store_true")
    parser.add_argument('-n', '--nonorm', \
      help="do not normalize the energy per particle, default is yes",\
        action="store_true")

    parser.add_argument('-b', '--beta', type=np.float64, \
      help="power law index for long range interactions", default=beta)

    return parser.parse_args()


#This is the Jmn hopping matrix with power law decay for open boundary
#conditions.
#This is the Jmn hopping matrix with power law decay for open boundary
#conditions.
def get_jmat_obc(args):
    J = dia_matrix((N, N))
    for i in xrange(1,N):
        elem = pow(i, -args.beta)
        J.setdiag(elem, k=i)
        J.setdiag(elem, k=-i)
    return J.toarray()

def get_jmat_pbc(args):
    J = dia_matrix((N, N))
    mid_diag = np.floor(N/2).astype(int)
    for i in xrange(1,mid_diag+1):
        elem = pow(i, -args.beta)
        J.setdiag(elem, k=i)
        J.setdiag(elem, k=-i)
    for i in xrange(mid_diag+1, N):
        elem = pow(N-i, -args.beta)
        J.setdiag(elem, k=i)
        J.setdiag(elem, k=-i)
    return J.toarray()

def precalc(params):
    # calc all the 'outer' sums/diffs, and zero the k=i,j terms
    ii = np.arange(params.jmat.shape[0])
    JM1 = params.jmat[:, :, None] + params.jmat[:, None, :]
    JM2 = params.jmat[:, :, None] - params.jmat[:, None, :]
    JM1[ii, ii, :] = 0
    JM2[ii, ii, :] = 0
    JM1[ii, :, ii] = 0
    JM2[ii, :, ii] = 0
    JM1 = JM1.transpose([1, 2, 0])
    JM2 = JM2.transpose([1, 2, 0])
    return JM1, JM2

class Precalc_objects:
    description = """Class to store conditions and precalculated objects"""

    def __init__(self, args):
        #Copy arguments from parser to this class
        self.__dict__.update(args.__dict__)
        self.lattice_size = N
        if args.periodic:
            self.periodic_boundary_conditions = 1
            self.open_boundary_conditions = 0
            self.jmat = get_jmat_pbc(args)
            mid = np.floor(N/2).astype(int)
            self.norm =\
              2.0 * np.sum(1/(pow(np.arange(1, mid+1), args.beta).astype(float)))
        else:
            self.periodic_boundary_conditions = 0
            self.open_boundary_conditions = 1
            self.jmat = get_jmat_obc(args)
            self.norm =\
              np.sum(1/(pow(np.arange(1, N+1), args.beta).astype(float)))
        if args.nonorm:
            self.norm = 1.0
        self.jmat = self.jmat/self.norm
        self.j_p, self.j_m = precalc(self)

def corr_time(t, JM1, JM2):
    return np.prod(np.cos(2*JM1*t), axis=-1)+np.prod(np.cos(2*JM2*t), axis=-1)

def dtwa_ising_analytical_longrange(params):
    if params.verbose:
        pprint("# Run parameters:")
        pprint(vars(params), depth=1)
    else:
        pprint("# Starting run ...")
    t_final = t_init + n_cycles
    dt = (t_final-t_init)/(n_steps-1.0)
    t_output = np.arange(t_init, t_final, dt)
    #Initial condition is all spins up, s^+_k has expectation of 1 for all k
    sx_init = np.ones(N)
    sx_corrs_in = np.ones(N*N).reshape(N, N)
    # These are the matrices P_{ki} = Cos(2* J_{ki} * t) for each t
    tmats = np.cos([2.0*params.jmat*t for t in t_output])
    pprint("# Building sx expectations now .. ")
    splus_totals = np.array([np.prod(mat, axis=1) \
      for mat in tmats]).dot(sx_init)/N
    pprint("# Building sx correlations now .. ")
    ssq_totals = np.array([np.sum(sx_corrs_in * \
      corr_time(t, params.j_p, params.j_m)) for t in t_output])/(2.0 * N**2)
    ssq_totals = ssq_totals - (splus_totals.real)**2
    pprint('# Finished building. Now dumping to file')
    np.savetxt(params.output_magx, np.vstack((np.abs(t_output), \
      splus_totals.real)).T, delimiter=' ')
    np.savetxt(params.output_sxvar, np.vstack((np.abs(t_output), \
      ssq_totals)).T, delimiter=' ')
    pprint('# Done!')

if __name__ == '__main__':
    args = input()
    params = Precalc_objects(args)
    dtwa_ising_analytical_longrange(params)
