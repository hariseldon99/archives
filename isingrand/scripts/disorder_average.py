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

Usage: disorder_average.py <files>
"""
import numpy as np
import scipy as sp
import scipy.signal 
import scipy.ndimage
import scipy.optimize as opt
import matplotlib.pyplot as plt
import sys,os.path,glob

gauss_winsize = 10
gauss_stdev = gauss_winsize/4
#Fraction of the total interval to be used for time average
frac = 1000

out_fname = "avg.dat"
out_fname_smooth = "avg_smooth.dat"

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
    avg_smooth = sp.ndimage.filters.convolve1d(avg,windows/windows.sum())
    params = opt.curve_fit(fit_func,x,avg_smooth)
    [mbar,tau] = params[0]  
    [mberr, tauerr] = np.diag(params[1])      
    plt.gca().set_color_cycle(['blue', 'red', 'green', 'yellow'])    
    plt.plot(x,avg)
    plt.plot(x,fit_func(x,mbar,tau))
    plt.xlim((0,x[-1]))
    plt.xlabel('Time')
    plt.ylabel('Average')
    plt.legend(['Disorder average','Disorder average - smooth fit'],loc='upper right')	
    #Uncomment below for file dump
    #plt.savefig('avg.png')
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
    lowerlim = outdat[0,0]
    upperlim = outdat[-1,0]
    diff = upperlim-lowerlim
    lowerlim = lowerlim+(diff/frac)
    m0 = avg_smooth[0]
    minfty = avg_smooth[-1]
    mfall = (m0+minfty)/2.0
    print
    print "Time avg of disorder avg from t = ",lowerlim,"- ",upperlim,"is:"
    dat_tavg = np.mean(outdat[outdat[:,0]>=lowerlim][:,1])     
    print dat_tavg
    print "Time it takes for data to fall by 1/2 of maximum:"
    print "t_1/2 = ", x[np.argmin(np.abs(avg_smooth-mfall))]
    print
    print "Timescale from curve fit:"
    print "tau = ", tau , "+/-" , tauerr
