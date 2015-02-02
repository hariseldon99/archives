#!/bin/bash

export MPI_NUM_PROCS=4

################################################################################################################################
########################################Default Sun Grid Engine options.########################################################
################################################################################################################################
#$ -S /bin/bash
#$ -V                           # Inherit the submission environment
#$ -cwd                         # Start job in  submission directory
#$ -N pdeloc                    # Job Name
#$ -j y                         # Combine stderr & stdout into stdout  
#$ -pe make 12                   # Requests number of SLOTS
				# The 'make' pe launches in round robin fashion one slot per node
#$ -o $JOB_NAME.o$JOB_ID        # Name of the output file (eg. myMPI.oJobID)
#$ -q all.q                     # Queue name
################################################################################################################################
################################################################################################################################

#Lattice Size
export LATSIZE=100

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

# Run the MPI executable.
mpirun -np $MPI_NUM_PROCS $RUNPROG -m $LATSIZE -a $AMPLITUDE -w $OMEGA -t $MFREQ -ts_max_steps $TIMESTEPS -ts_final_time $FINALTIME -random_seed $RANDSEED
