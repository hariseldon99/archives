#!/usr/bin/python

"""
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * Rigol Lattice: postprocessing of Petsc Matrix data
 * Copyright (c) 2013 Analabha Roy (daneel@utexas.edu)
 * 
 * This is free software: you can redistribute it and/or modify it under  the
 * terms of version 3 of the GNU Lesser General Public License as published by
 * the Free Software Foundation.
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""

"""
Python program to 
read a petsc binary matrix M using numpy
and tests to see whether it is heritian
by computing the norm of the eigenvalues
of (M-M^\dagger) and comparing with
TOLERANCE (i.e the trace distance)

Usage: mat_hermtest.py <PetSc binary matrix>
"""

TOLERANCE=1e-9

import numpy as np
import scipy.linalg
import sys,os
sys.path.append(os.path.join(os.environ['PETSC_DIR'],'bin','pythonscripts'))
import PetscBinaryIO

if __name__ == '__main__':
    petsc = PetscBinaryIO.PetscBinaryIO()
    fullmat =  petsc.readBinaryFile(sys.argv[1],mattype='dense')
    fullmat =  np.matrix(fullmat[0])
    fullmat_T = np.transpose(np.conjugate(fullmat))
    diffmat = fullmat - fullmat_T
    prodmat = diffmat * np.transpose(np.conjugate(diffmat))
    trace_vals = scipy.linalg.eigvals(prodmat)
    dist = scipy.linalg.norm(trace_vals)
    if(dist > TOLERANCE):
        print "!!! Matrix is NOT HERMITIAN !!!"
    else:
        print "@@@ Matrix is Hermitian to within tolerance", TOLERANCE ," @@@"