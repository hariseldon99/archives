#!/usr/bin/python
from __future__ import division, print_function

"""
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 * Postprocessing of DTWA
 * Copyright (c) 2014 Analabha Roy (daneel@utexas.edu)
 * 
 * This is free software: you can redistribute it and/or modify it under  the
 * terms of version 3 of the GNU Lesser General Public License as published by
 * the Free Software Foundation.
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""

"""
Python program to 
read files of single particle observables 
as functions of time obtained from the DTWA
and estimate the relaxation time, if any.

Usage: relaxation_time.py <input file>
"""
import numpy as np
from scipy.fftpack import rfft, rfftfreq
import sys,os.path,glob
import pylab as plt
import peakdetect as pkd

pass_band = (np.pi/30, 10.0)

def processFile(filename):
    fileHandle = open(filename, "r")
    out = list()
    for line in fileHandle:
        # do some processing
        line=line.strip()
        out.append(np.array([float(d) for d in line.split()]))
    fileHandle.close()
    return np.array(out)

if __name__ == '__main__':
    if (len(sys.argv) > 1):
        filedata = processFile(sys.argv[1])
        if (len(sys.argv) > 2):
	  pass_band = float(sys.argv[2]), float(sys.argv[3])
        print("pass-band =", pass_band)
        t = filedata[:,0]
        data = filedata[:,1]
    else:
       print("Usage: relaxation_time.py <data-file> pass-band. Default= "\
	      +format(pass_band))
       sys.exit()
    Delta = t[1]-t[0];   
    FFT = rfft(data)
    omega = rfftfreq(data.size, d=Delta)
    #Removed aliased frequencies
    nonaliased = FFT.size/2
    FFT = FFT[0:nonaliased]
    omega = omega[0:nonaliased]
    
    #Apply band pass filter
    FFT_cut = FFT.copy()
    FFT_cut[(omega<=pass_band[0])] = 0.0
    FFT_cut[(omega>=pass_band[1])] = 0.0
    
    #Estimate the HWHM by from local peaks
    #of FFT 
    indx = pkd.detect_peaks(np.absolute(FFT_cut), valley=False)
    absfft = np.absolute(FFT_cut[indx])
    wpeaks = omega[indx]
    w1 = wpeaks[absfft.argmax()]
    w2 = wpeaks[::-1][np.abs(absfft[::-1]-np.max(absfft)/2.0).argmin()]
    hwhm = w2-w1
    tau_fft = np.pi/hwhm
    
    #Identify the peaks in the original data
    modded_data = np.abs(data-np.mean(data))
    peak_inds = pkd.detect_peaks(modded_data)
    peaktimes = t[peak_inds]
    peakvals = modded_data[peak_inds]
    #Add the initial value to this data
    peaktimes = np.insert(peaktimes,0,t[0])
    peakvals = np.insert(peakvals,0,modded_data[0])
    #Decay point where peakval is close to 1/e of the max
    dp = np.max(peakvals)/10.0
    #Decay time is where this value is close
    tau = peaktimes[np.abs(peakvals-dp).argmin()]
    data_near_tau = peakvals[np.abs(peakvals-dp).argmin()]
    
    f = plt.figure()
    
    ax=plt.subplot(221)
    ax.set_title('Data ')
    ax.plot(t,data)
    
    ax=plt.subplot(222)
    ax.set_title('FFT of data')
    ax.plot(omega, np.abs(FFT))
    ax.set_xlim(-0.3,omega[-1])
    
    ax=plt.subplot(223)
    ax.set_title('Modded data w tau estimate')
    ax.scatter(peaktimes, np.abs(peakvals), color='red', linewidth=2)
    ax.plot(t, np.abs(modded_data))
    pointloc=(tau_fft, np.abs(data_near_tau))
    textloc=(tau, 1.5 * np.abs(data_near_tau))
    ax.annotate("tau = "+"{:1.3f}".format(tau), xy=pointloc,\
		xytext=textloc, textcoords='data')
    
    ax=plt.subplot(224)
    ax.set_title('FFT with filter, tau from HWHM= '+\
		  '{:1.3f}'.format(tau_fft))
    ax.plot(omega, np.absolute(FFT_cut), color='blue', linewidth=0.5)
    shadelow = np.where(omega==w1)[0][0]
    shadehigh = np.where(omega==w2)[0][-1]
    bulkfreqs =  omega[shadelow:shadehigh]
    bulkamps = np.absolute(FFT_cut[shadelow:shadehigh])
    
    ax.fill_between(bulkfreqs,bulkamps, facecolor='cyan', edgecolor='cyan', alpha=0.8)
    ax.set_xlim((-0.3,omega[-1]))

    plt.show()