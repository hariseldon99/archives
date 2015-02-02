#!/bin/bash

################################################################################################################################
########################################Lonestar Sun Grid Engine options.########################################################
################################################################################################################################
#$ -V                     	# Inherit the submission environment 
#$ -cwd                   	# Start job in  submission directory
#$ -N pdeloc         		# Job Name
#$ -j y                   	# Combine stderr & stdout into stdout  
#$ -o $JOB_NAME.o$JOB_ID  	# Name of the output file (eg. myMPI.oJobID)
#$ -M daneel@utexas.edu		# Email for job notification
#$ -m be			# Email at Begin and End of job
#$ -pe 12way 12                  # Ex: -pe 1way 24. Requests number of cores (24). 
#$ -q normal              	# Queue name ("serial" for single node run with openmp)
#$ -l h_rt=6:00:00       	# Run time (hh:mm:ss) 
################################################################################################################################
##############

#Lattice Size
export LATSIZE=1000

#Varying parameters, default values
export AMPLITUDE=0.5
export OMEGA=3.0
export MFREQ=3.0 


#Fixed parameters
export FINALTIME=100.0
export TIMESTEPS=500

export RUNPROG="./pdeloc"
#Random number seed
export RANDSEED=3

#Scale for floating point arithmetic in this script
float_scale=8
#string comparators

if [ ! -e "$RUNPROG" ]; then
  echo 
  echo "Error! runtime binary not found ..."
  echo 
  exit
fi

# Run the MPI executable. Use TACC's affinity utility to position one task on each socket.
# See http://www.tacc.utexas.edu/user-services/user-guides/lonestar-user-guide#running:numa:socket
ibrun tacc_affinity $RUNPROG -m $LATSIZE -a $AMPLITUDE -w $OMEGA -t $MFREQ -ts_max_steps $TIMESTEPS -ts_final_time $FINALTIME -random_seed $RANDSEED
