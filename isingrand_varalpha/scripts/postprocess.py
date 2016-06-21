#!/usr/bin/python

"""
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * isingrand: postprocessing of isingrand data
 * Copyright (c) 2013 Analabha Roy (daneel@utexas.edu)
 * 
 * This is free software: you can redistribute it and/or modify it under  the
 * terms of version 3 of the GNU Lesser General Public License as published by
 * the Free Software Foundation.
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""

"""
Python program to 
read files of quantities as functions of time
and dump their averages as functions of time
and do fft of the smoothed difference

Usage: disorder_average.py <files>
"""
import numpy as np
import scipy.signal 
import scipy.optimize as opt
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import sys,os.path,glob
from matplotlib import rc
rc('text', usetex=True)

gauss_winsize = 5
gauss_stdev = gauss_winsize/4

out_fname = "avg.dat"
out_fname_smooth = "avg_smooth.dat"
out_fname_osc = "avg_osc.dat"
out_fname_fft = "fft.dat"

def fit_func(t,mbar,tau):
    c = mbar - 1
    return 1 + c * (1 - np.exp(-t/tau))

def processFile(filename):
    fileHandle = open(filename, "r")
    out = list()
    for line in fileHandle:
        # do some processing
        line=line.strip()
        for number in line.split():
            out.append(float(number))
    fileHandle.close()
    return out

def processFiles(args):
    input_filemask = "log"
    directory = args[1]
    shape = (-1,2)
    listofdata = list()
    xvals_avg = list()
    if os.path.isdir(directory):
        print "processing a directory"
        list_of_files = glob.glob('%s/*.%s' % (directory, input_filemask))
    else:
        print "processing a list of files"
        list_of_files = sys.argv[1:]
    for file_name in list_of_files:
        print file_name
        data = np.array(processFile(file_name))
        data.shape = shape
        listofdata.append(data)
    listofdata = np.array(listofdata)
    listofdata.shape = shape
    #Select common times,
    xvals = np.unique(listofdata[:,0])
    for x in xvals:
        datasubset = listofdata[listofdata[:,0] == x]
        #append avg to xvals 
        xvals_avg.append(np.mean(datasubset[:,1]))
    return xvals, xvals_avg

if __name__ == '__main__':
    if (len(sys.argv) > 1):
        x,avg = processFiles(sys.argv)
    else:
        print 'Usage: disorder_average.py <files> or <directory>'
    windows = scipy.signal.gaussian(gauss_winsize,gauss_stdev)    
    avg_smooth = ndimage.filters.convolve1d(avg,windows/windows.sum())
    params = opt.curve_fit(fit_func,x,avg_smooth)
    [mbar,tau] = params[0]  
    [mberr, tauerr] = np.diag(params[1])
    fitdata = fit_func(x,mbar,tau)
    res_osc = avg- fitdata
    freqs = np.fft.fftfreq(res_osc.size, d=x[1]-x[0])
    spectrum = np.fft.fft(res_osc)
    s = np.array_split(spectrum,2)[0]
    s = s/np.amax(np.abs(s))
    f = np.array_split(freqs,2)[0] #Nyquist rate
    
    plt.figure()
    ax1 = plt.subplot(221)
    plt.ylabel("Mag.")
    plt.plot(x,avg,'b',label="Disorder average")
    plt.plot(x,fitdata,'r',label="Disorder average - Smoothed")
    plt.legend()
    ax2 = plt.subplot(223)
    plt.plot(x,res_osc, label="Difference of the two")
    plt.xlabel("Time")
    plt.legend()

    ax3 = plt.subplot(222)
    plt.plot(2.0 * np.pi * np.abs(f),np.abs(s), label="Spectra (FFT)")
    plt.legend()

    ax4 = plt.subplot(224)
    plt.plot(2.0 * np.pi * np.abs(f),np.abs(s), label="Spectra (FFT - logscale)")
    plt.xlabel(r"$\omega$")
    ax4.set_yscale('log')
    plt.legend()

    plt.show()
    print "\nDumping averages to file" ,  out_fname , "..."
    x,avg = np.array(x),np.array(avg)
    outdat = np.vstack((x, avg)).T
    np.savetxt(out_fname,outdat,delimiter=' ')
    print "Done!"
    print "\nDumping smoothed averages to file" ,  out_fname_smooth , "..."
    x,avg_smooth = np.array(x),np.array(avg_smooth)
    outdat = np.vstack((x, avg_smooth)).T
    np.savetxt(out_fname_smooth,outdat,delimiter=' ')
    print "Done!"    
    print "\nDumping residual oscillations to file" ,  out_fname_osc , "..."
    outdat = np.vstack((x,res_osc)).T
    np.savetxt(out_fname_osc,outdat,delimiter=' ')
    print "\n Dumping FFT to file", out_fname_fft , "..."
    outdat = np.vstack((2.0 * np.pi * np.abs(f),np.abs(s))).T
    np.savetxt(out_fname_fft,outdat,delimiter=' ')
    print "Done!"
