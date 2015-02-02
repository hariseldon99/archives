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
and plot their fluctuations as function of time

Usage: disorder_stdev.py <files>
"""
import numpy as np
import scipy as sp
import scipy.signal 
import scipy.ndimage
import matplotlib.pyplot as plt
import sys,os.path,glob

gauss_winsize = 40
gauss_stdev = gauss_winsize/4
#Fraction of the total interval to be used for time average
frac = 5

out_fname = "stdev.dat"
out_fname_smooth = "stdev_smooth.dat"

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
    xvals_stdev = list()
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
        #append stdev to xvals 
        xvals_stdev.append(np.std(datasubset[:,1]))
    return xvals, xvals_stdev

if __name__ == '__main__':
    if (len(sys.argv) > 1):
        x,stdev = processFiles(sys.argv)
    else:
        print 'Usage: disorder_stdev.py <files> or <directory>'
    windows = scipy.signal.gaussian(gauss_winsize,gauss_stdev)    
    stdev_smooth = sp.ndimage.filters.convolve1d(stdev,windows/windows.sum())
    plt.gca().set_color_cycle(['blue', 'red', 'green', 'yellow'])    
    plt.plot(x,stdev)
    plt.plot(x,stdev_smooth)
    plt.xlim((0,x[-1]))
    plt.xlabel('Time')
    plt.ylabel('Fluctuations')
    plt.legend(['Disorder fluctuations','Disorder fluctuations - smoothed'],loc='upper right')	    
    #Uncomment below for file dump
    #plt.savefig('stdev.png')    
    plt.show()
    print "\nDumping stdev to file" ,  out_fname , "..."
    x,stdev = np.array(x),np.array(stdev)
    outdat = np.vstack((x, stdev)).T
    np.savetxt(out_fname,outdat,delimiter=' ')
    print "Done!"
    print "\nDumping stdev smooth to file" ,  out_fname_smooth , "..."
    x,stdev_smooth = np.array(x),np.array(stdev_smooth)
    outdat = np.vstack((x, stdev_smooth)).T
    print "Time avg of disorder avg = ",np.mean(outdat[:,1])
    np.savetxt(out_fname_smooth,outdat,delimiter=' ')
    print "Done!"    
    lowerlim = outdat[0,0]
    upperlim = outdat[-1,0]
    diff = upperlim-lowerlim
    lowerlim = lowerlim+(diff/frac)
    print
    print "Time avg of disorder stdev from t = ",lowerlim,"- ",upperlim,"is:"
    print np.mean(outdat[outdat[:,0]>=lowerlim][:,1])
    print