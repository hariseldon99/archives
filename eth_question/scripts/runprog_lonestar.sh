#!/bin/bash
#$ -V                     # Inherit the submission environment 
#$ -cwd                   # Start job in  submission directory
#$ -N name                # Job Name
#$ -j y                   # combine stderr & stdout into stdout  
#$ -o $JOB_NAME.o$JOB_ID  # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 36           # Requests 12 cores/node, 24 cores total
#$ -q normal              # Queue name
#$ -l h_rt=01:30:00       # Run time (hh:mm:ss) - 1.5 hours

export LATSIZE=5
export NUMPTCLS=2
export REPULSION=0.1
export DELTA_E=0.1
export RUNPROG="./rigol_lattice"
	      
ibrun $RUNPROG -lattice_size $LATSIZE -vector_size $NUMPTCLS -repulsion $REPULSION

