#!/bin/bash

#Number of parallel threads. Set to 1 if you want serial computing
export OMP_NUM_THREADS=12

################################################################################################################################
########################################Lonestar Sun Grid Engine options.########################################################
################################################################################################################################
#$ -V                     	# Inherit the submission environment 
#$ -cwd                   	# Start job in  submission directory
#$ -N test         		# Job Name
#$ -j y                   	# Combine stderr & stdout into stdout  
#$ -o $JOB_NAME.o$JOB_ID  	# Name of the output file (eg. myMPI.oJobID)
#$ -pe 1way 24                  # Ex: -pe 1way 24. Requests number of cores (24). Kept this '1way' or equivalent for 1 task/node
#$ -q normal              	# Queue name ("serial" for single node run with openmp)
#$ -l h_rt=02:30:00       	# Run time (hh:mm:ss) 
################################################################################################################################
##############

#Freq range
export W_INIT=10.0
export W_FINAL=30.0
export W_NOS=10.0

#Varying parameters, default values
export AMPLITUDE=12.02413 
export OMEGA=20.0 
export PERT_AMP=0.1 

#First zero of the Bessel function
#Source: http://wwwal.kuicr.kyoto-u.ac.jp/www/accelerator/a4/besselroot.htmlx
export BESJZERO=2.40482555769577 

#Fixed parameters
export FINALTIME=100.0
export TIMESTEPS=500
export MAG_FILE="mags.dat"
export ENT_FILE="ents.dat"
export ERROR_FILE=/dev/null
#Answer "y" or "n" for verbose output
export VERBOSITY="n"

#GSL random number generator parameters. See GNU Scientific Library Manual for options
export GSL_RNG_TYPE="ranlxs2" 
export GSL_RNG_SEED=3
#Actual run command.
export RUNPROG="./isingrand_parallel"


################################################################################################################################
########################################Do not edit this code block#############################################################
################################################################################################################################
if  [ $VERBOSITY = "y" ]; then
    VERBOSITY="-v"
else
    VERBOSITY=" "
fi

#Scale for floating point arithmetic in this script
float_scale=8
#string comparators

export W_INC=$(echo "scale=$float_scale ;($W_FINAL-$W_INIT)/$W_NOS" | bc -q 2>/dev/null)

if [ ! -e "$RUNPROG" ]; then
  echo 
  echo "Error! runtime binary not found ..."
  echo 
  exit
fi

for i in $(seq $W_NOS)
do
  #Amplitude
  export OMEGA=$(echo "scale=$float_scale ;$W_INIT + ($i-1)*$W_INC" | bc -q 2>/dev/null)
  #Output file path and name
  OUTDIR="result_omega_"
  OUTDIR=$OUTDIR$OMEGA
  export OUTDIR
  if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
  fi
  #Evaluate amplitude
  export AMPLITUDE=$(echo "scale=$float_scale ;($BESJZERO * $OMEGA)/4.0" | bc -q 2>/dev/null)
  # Run the MPI executable. Use TACC's affinity utility to position one task on each socket.
  # See http://www.tacc.utexas.edu/user-services/user-guides/lonestar-user-guide#running:numa:socket
  ibrun tacc_affinity $RUNPROG $VERBOSITY -a $AMPLITUDE -f $OMEGA -p $PERT_AMP -t $FINALTIME -n $TIMESTEPS -m "${OUTDIR}/${MAG_FILE}" -e "${OUTDIR}/${ENT_FILE}" 2>$ERROR_FILE
done
