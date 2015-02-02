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
read a petsc binary matrix using numpy
and matplots it

Usage: mat_binplot.py <PetSc binary matrix>
"""
import numpy as np
import matplotlib.pylab as pylab
import sys,os
sys.path.append(os.path.join(os.environ['PETSC_DIR'],'bin','pythonscripts'))
import PetscBinaryIO

if __name__ == '__main__':
    petsc = PetscBinaryIO.PetscBinaryIO()
    fullmat =  petsc.readBinaryFile(sys.argv[1],mattype='dense')
    fullmat =  np.matrix(fullmat[0].real)
    pylab.matshow(fullmat,cmap=pylab.cm.coolwarm_r)
    pylab.suptitle('matrix plot', fontsize=12)
    pylab.show()
    # Uncomment for remote systems
    # pylab.savefig('mat.svg')
    # print "Output in svg file"
