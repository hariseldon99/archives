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
read a petsc binary matrix and convert it to ascii
in stdout
Usage: mat_bin2ascii.py <PetSc binary matrix>
"""

import numpy as np
import sys,os
sys.path.append(os.path.join(os.environ['PETSC_DIR'],'bin','pythonscripts'))
import PetscBinaryIO

#Float tolerance
tol = np.finfo(np.float).eps

if __name__ == '__main__':
    petsc = PetscBinaryIO.PetscBinaryIO()
    fullmat =  petsc.readBinaryFile(sys.argv[1],mattype='sparse')
    (M,N), (I,J,V)  = fullmat[0][0], fullmat[0][1]
    for i in xrange(len(I)-1):
        start, end = I[i], I[i+1]
        colidx = J[start:end]
        values = V[start:end]
        values.real[abs(values.real) < tol] = 0.0
        values.imag[abs(values.imag) < tol] = 0.0
        print 'row %d:' % i ,  zip(colidx, values)
