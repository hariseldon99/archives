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
and dump their time averages

Usage: time_average.py <files>
"""
import numpy as np
import sys,os.path,glob

#Fraction of the total interval to be used for time average
frac = 5

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
    timeavgs=list()
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
        lowerlim = data[0,0]
        upperlim = data[-1,0]
        diff = upperlim-lowerlim
        lowerlim = lowerlim+(diff/frac)
        timeavgs.append(np.mean(data[data[:,0]>=lowerlim][:,1]))
    return timeavgs

if __name__ == '__main__':
    if (len(sys.argv) > 1):
        avgs = processFiles(sys.argv)
    else:
        print 'Usage: time_average.py <files> or <directory>'
    print "Time averages:"    
    avgs = np.array(avgs)
    print np.vstack(avgs)