#!/usr/bin/python

"""
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * Rigol Lattice: postprocessing of Petsc Vector data
 * Copyright (c) 2013 Analabha Roy (daneel@utexas.edu)
 * 
 * This is free software: you can redistribute it and/or modify it under  the
 * terms of version 3 of the GNU Lesser General Public License as published by
 * the Free Software Foundation.
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""

"""
Python program to 
read a petsc binary vector using numpy
and plot by index

Usage: vec_binplot.py <PetSc binary vector>
"""
import numpy as np
import matplotlib.pyplot as plt
import sys,os
sys.path.append(os.path.join(os.environ['PETSC_DIR'],'bin','pythonscripts'))
import PetscBinaryIO

if __name__ == '__main__':
    petsc = PetscBinaryIO.PetscBinaryIO()
    vec = petsc.readBinaryFile(sys.argv[1])
    vec = list(vec)[0]
    x = range(len(vec))
    #np.ndarray.sort(vec)
    print zip(x,vec)
