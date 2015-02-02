#!/bin/bash

#Number of parallel threads. Set to 1 if you want serial computing
export OMP_NUM_THREADS=2
export MPI_NUM_PROCS=5
################################################################################################################################
########################################Default Sun Grid Engine options.########################################################
################################################################################################################################
#$ -S /bin/bash
#$ -V                           # Inherit the submission environment
#$ -cwd                         # Start job in  submission directory
#$ -N isingrand_root             # Job Name
#$ -j y                         # Combine stderr & stdout into stdout  
#$ -pe make 5                   # Requests number of SLOTS
				# The 'make' pe launches in round robin fashion one slot per node
#$ -o $JOB_NAME.o$JOB_ID        # Name of the output file (eg. myMPI.oJobID)
#$ -q all.q                     # Queue name
################################################################################################################################
################################################################################################################################

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
  mpirun -np $MPI_NUM_PROCS --bynode $RUNPROG $VERBOSITY -a $AMPLITUDE -f $OMEGA -p $HOPPING -t $FINALTIME -n $TIMESTEPS -m "${OUTDIR}/${MAG_FILE}" -e "${OUTDIR}/${ENT_FILE}" 2>$ERROR_FILE  
done
