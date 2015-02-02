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
and average the columns.

Thus, if the matrix is a discretized
2D function f(x,y) with x as rows and y 
as columns, then this program returns 
f(x) = (1/R) \int dy f(x,y), where R is 
the range of the x axis.

Also, it returns the closest approximation 
to the number f(0) where x,y are assumed to lie
from -pi, pi in equal steps

Usage: mat_bin21d.py <PetSc binary matrix>
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import sys,os
sys.path.append(os.path.join(os.environ['PETSC_DIR'],'bin','pythonscripts'))
import PetscBinaryIO

if __name__ == '__main__':
    petsc = PetscBinaryIO.PetscBinaryIO()
    fullmat =  petsc.readBinaryFile(sys.argv[1],mattype='dense')
    fullmat =  np.matrix(fullmat[0].real)
    (M,N) = fullmat.shape
    col_avg = fullmat.mean(axis = 0)
    col_avg = np.array(col_avg).reshape(-1,).tolist()
    #
    # Brillouin Zone
    #
    idx = range(N)
    kxinit = -math.pi
    kxfinal = math.pi
    kxinc = (kxfinal - kxinit)/(len(idx)-1)
    bzone = [kxinit + i * kxinc for i in idx]
    #
    # Plot col average of matrix as a function of kx
    #
    plt.plot(bzone,col_avg)
    plt.xlim((kxinit,kxfinal))
    plt.show()
    # Uncomment for remote systems
    # plt.savefig('nx.svg')
    # print "Output in svg file"
    mid_idx = 1 + int(len(col_avg)/2)
    print 'Column average for kx = 0 : %f' % col_avg [mid_idx]
    #print zip(idx,col_avg)
