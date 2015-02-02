#!/usr/bin/python

"""
Created on Tues Mar  4 2014

@author: Analabha Roy (daneel@utexas.edu)

Usage : ./isingrand_esys.py -h
"""

import numpy as np
import numpy.random as nprand
import numpy.linalg as nplinalg
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys, getopt, time

"""
Default Parameters are entered here
"""
#Lattice size
L = 7

#Transverse field
alpha = 10.0

#Random number seed
seed = 7

out_fname_evals = "isingrand_evals.npy"
out_fname_evecs = "isingrand_evecs.npy"
out_fname_mags =  'mags_ebasis_tfield_%f_latsize_%d.txt'
out_fname_cc = 'initstate_ebasis_tfield_%f_latsize_%d.txt'

helpstring = """Usage:\n\n ./isingrand_esys.py 
                    -l <lattice size> 
                    -s <random no generator seed> 
                    -f <transverse field> 
                    -e </path/to/evals/outfile> 
                    -u </path/to/evecs/outfile> 
                    -m </path/to/mags/outfile> 
                    -v <verbose>"""

#Pauli matrices
sx,sy,sz = np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),\
np.array([[1,0],[0,-1]])                    

def kemat(i):
    #Left Hand Side    
    if(i == 1):
        ke = np.kron(sx,sx)
    else:
        dim = 2**(i-1)
        ke = np.kron(np.eye(dim,dim),sx)
        ke = np.kron(ke,sx)
    #Right Hand Side    
    if(i < L-1):
        dim = 2**(L-i-1)
        ke = np.kron(ke,np.eye(dim,dim))
    return(ke)  

def nummat(i):
    #Left Hand Side
    if(i==1):
        num = sz 
    else:
        dim = 2**(i-1)
        num = np.kron(np.eye(dim,dim),sz)
    #Right Hand Side    
    if(i < L):
        dim = 2**(L-i)
        num = np.kron(num,np.eye(dim,dim))
    return(num)

def diagonalize(H):
    try:
        evals, U = nplinalg.eigh(H)
        idx = evals.argsort()
        evals = evals[idx]
        U = U[:,idx]
    except nplinalg.linalg.LinAlgError:
        evals, U = 0.0,0.0 
    return(evals, U)    

if __name__ == '__main__':
    verbosity = 0
    try:
        opts, args = getopt.getopt(sys.argv[1:] , "hl:s:f:e:u:m:v")
    except getopt.GetoptError:
        print helpstring
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print helpstring
            sys.exit(1)
        elif opt == '-l':
            L = int(arg)
        elif opt == '-s':
            seed = int(arg)
        elif opt == '-f':
            alpha = float(arg)
        elif opt == '-e':
            out_fname_evals = arg
        elif opt == '-u':
            out_fname_evecs = arg
        elif opt == '-m':
            out_fname_mags = arg
        elif opt == '-v':
            verbosity = 1
    print "\n\nFor usage help, run ./isingrand_esys.py -h\n\n"
    print "Executing diagonalization with parameters:"
    print "Lattice Size = " , L 
    print "Trasnverse field = " , alpha 
    print "\n\n ... Please Wait ...\n\n"
    print "-----------------------------------"
    print "Tfield |", "Gndst |" , "Mag avg |" , "time (s)"
    print "-----------------------------------"
    numstates = np.delete(np.arange(L+1),0)
    limits = (-2**L/100,2**L)
    start_time = time.time()
    nprand.seed(seed)
    hdata = np.array(nprand.rand(L))
    jdata = np.array(nprand.rand(L))
    H = np.zeros((2**L,2**L))
    for i in numstates:
        H = H - 2.0 * alpha * hdata[i-1] *  nummat(i)
        if(i < L-1):
            H = H - jdata[i-1] * kemat(i)
    evals, evecs = diagonalize(H)
    if evals.all != 0.0 and evecs.all != 0.0 :
        #Compute magnetization diagonal components m_{aa}|<apsi_0>|^2
        magop = np.zeros((2**L,2**L))
        for i in numstates:
            magop = magop + nummat(i)
        #mag_aa is the array of m_aa = <a|M|a> for eigenvector |a>    
        mag_aa = np.diagonal(np.dot(np.transpose(evecs),np.dot(magop,evecs)))
        #Assume that psi_0 is all spins up ie (1 0 0 0 ... 0)
        #Then <a|psi_0> is just the complex conjugate of the first element of 
        #|a>
        #So, the array of | <a|psi_0>|^2 is just the first row of evecs mod - 
        #squared
        diag_cc = np.abs(evecs[0,:])**2
        mags = np.multiply(mag_aa,diag_cc)
        np.set_printoptions(precision=3)
        print "%2.3f  | %2.3f | %2.3f   | %2.3f" % (alpha, \
        np.min(np.abs(evals.real)), np.sum(mags)/L ,time.time() - start_time)
        print "-----------------------------------"
        if(verbosity == 1):
            #Plot the eigenvalues and eigenvectors and magnetization components
            plt.plot(evals.real/L)
            plt.title('Eigenvalues - Sorted')
            plt.matshow(np.abs(evecs)**2,interpolation='nearest',\
            cmap=cm.coolwarm)
            plt.title('Eigenvectors - Sorted' ) 
            plt.colorbar()
            fig, ax = plt.subplots()
            plt.bar(np.delete(np.arange(2**L+1),0),diag_cc,edgecolor='blue')                
            plt.xlim(limits)
            plt.text(L,-np.max(diag_cc)/10.0,'Transverse field = %lf'%alpha,\
            horizontalalignment='center')
            #plt.title('Transverse field = %lf'%alpha )
            plt.xscale('log')
            ax.xaxis.tick_top()
            ax.set_xlabel('Lexicographical order of eigenstate')
            ax.xaxis.set_label_position('top') 
            ax.set_ylabel('Initial state probability wrt eigenstate')            
            plt.show()
            print "\nDumping data to file" ,  out_fname_evals , " and " , \
            out_fname_evecs , "..."
            np.save(out_fname_evals,evals)
            np.save(out_fname_evecs, evecs)
        out_fname_mags = out_fname_mags % (alpha,L)    
        print "\nDumping outputs to files" , out_fname_mags , " and " ,\
        out_fname_cc ,"..."
        outdat = np.vstack((np.delete(np.arange(2**L+1),0),mags)).T
        np.savetxt(out_fname_mags,outdat,delimiter=' ')
        out_fname_cc = out_fname_cc % (alpha,L)
        outdat = np.vstack((np.delete(np.arange(2**L+1),0),diag_cc)).T      
        np.savetxt(out_fname_cc,outdat,delimiter=' ')
        print 'Done'
    else:
        print "Error! Eigenvalues did not converge for this parameter,\
        skipping ..."