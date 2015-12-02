#!/bin/bash
################################################################################################################################
########################################Lonestar Sun Grid Engine options.########################################################
################################################################################################################################
#$ -V                     	# Inherit the submission environment 
#$ -cwd                   	# Start job in  submission directory
#$ -N test         		# Job Name
#$ -j y                   	# Combine stderr & stdout into stdout  
#$ -o $JOB_NAME.o$JOB_ID  	# Name of the output file (eg. myMPI.oJobID)
#$ -pe 8way 12                   # Ex: -pe 1way 24.
#$ -q normal              	# Queue name ("serial" for single node run with openmp)
#$ -l h_rt=2:00:00       	# Run time (hh:mm:ss) 
################################################################################################################################
################################################################################################################################
SCRIPT="./run_dtwa.py"

# make sure I'm the only one that can read my output
umask 0077
# Load python openmpi, mpi4py
module load mkl python/2.7.3-epd-7.3.2

##Now run my prog
BEGINTIME=$(date +"%s")
ibrun python $SCRIPT
ENDTIME=$(date +"%s")
ELAPSED_TIME=$(($ENDTIME-$BEGINTIME))

echo "#Runtime: $(($ELAPSED_TIME / 60)) minutes and $(($ELAPSED_TIME % 60)) seconds."
