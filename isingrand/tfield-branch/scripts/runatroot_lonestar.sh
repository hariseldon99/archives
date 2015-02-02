#!/bin/bash

#Number of parallel threads. Set to 1 if you want serial computing
export OMP_NUM_THREADS=12

################################################################################################################################
########################################Lonestar Sun Grid Engine options.########################################################
################################################################################################################################
#$ -V                     	# Inherit the submission environment 
#$ -cwd                   	# Start job in  submission directory
#$ -N betaroot         		# Job Name
#$ -j y                   	# Combine stderr & stdout into stdout  
#$ -o $JOB_NAME.o$JOB_ID  	# Name of the output file (eg. myMPI.oJobID)
#$ -pe 1way 24                  # Ex: -pe 1way 24. Requests number of cores (24). Kept this '1way' or equivalent for 1 task/node
#$ -q normal              	# Queue name ("serial" for single node run with openmp)
#$ -l h_rt=12:00:00       	# Run time (hh:mm:ss) 
################################################################################################################################
##############

#Freq range
export ETA_INIT=6.7
export ETA_FINAL=6.9
export ETA_NOS=10.0

#Varying parameters, default values
export HOPPING=5.0
export OMEGA=10.0

#Second zero of the Beta function
export BETAZERO=6.78102763986

#Fixed parameters
export FINALTIME=100.0
export TIMESTEPS=1000
export MAG_FILE="mags.dat"
export ENT_FILE="ents.dat"
export ERROR_FILE=/dev/null
#Answer "y" or "n" for verbose output
export VERBOSITY="n"

#GSL random number generator parameters. See GNU Scientific Library Manual for options
export GSL_RNG_TYPE="ranlxs2" 
export GSL_RNG_SEED=3
#Actual run command.
export RUNPROG="./isingrand_tfield"


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

#Calculate AMPLITUDES
export AMP_INIT=$(echo "scale=$float_scale ;$ETA_INIT * $OMEGA/4.0" | bc -q 2>/dev/null)
export AMP_FINAL=$(echo "scale=$float_scale ;$ETA_FINAL * $OMEGA/4.0" | bc -q 2>/dev/null)
export AMP_INC=$(echo "scale=$float_scale ;($AMP_FINAL-$AMP_INIT)/$ETA_NOS" | bc -q 2>/dev/null)

if [ ! -e "$RUNPROG" ]; then
  echo 
  echo "Error! runtime binary not found ..."
  echo 
  exit
fi

for i in $(seq $ETA_NOS)
do
  export AMPLITUDE=$(echo "scale=$float_scale ;$AMP_INIT + ($i-1)*$AMP_INC" | bc -q 2>/dev/null)
  #Output file path and name
  OUTDIR="result_amplitude_"
  OUTDIR=$OUTDIR$AMPLITUDE
  export OUTDIR
  if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
  fi
  # Run the MPI executable. Use TACC's affinity utility to position one task on each socket.
  # See http://www.tacc.utexas.edu/user-services/user-guides/lonestar-user-guide#running:numa:socket
  #ibrun tacc_affinity $RUNPROG $VERBOSITY -a $AMPLITUDE -f $OMEGA -p $HOPPING -t $FINALTIME -n $TIMESTEPS -m "${OUTDIR}/${MAG_FILE}" -e "${OUTDIR}/${ENT_FILE}" 2>$ERROR_FILE
  #Serial run
  mpirun -np $OMP_NUM_THREADS $RUNPROG $VERBOSITY -a $AMPLITUDE -f $OMEGA -p $HOPPING -t $FINALTIME -n $TIMESTEPS -m "${OUTDIR}/${MAG_FILE}" -e "${OUTDIR}/${ENT_FILE}" 2>$ERROR_FILE  
done
