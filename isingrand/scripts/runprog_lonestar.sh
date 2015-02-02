#!/bin/bash

#Number of parallel threads. Set to 1 if you want serial computing
export OMP_NUM_THREADS=15

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
################################################################################################################################

export VARYING_PARAM='alpha' #omega for frequency, ampl for amplitude, alpha for perturbation amplitude
export PARAM_INIT=0.1
export PARAM_FINAL=1.0
export PARAM_NOS=4.0

#Varying parameters, default values
export AMPLITUDE=12.02413 
export OMEGA=20.0 
export PERT_AMP=0.5 

#Fixed parameters
export FINALTIME=10.0
export TIMESTEPS=300
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
export S1='omega'
export S2='ampl'
export S3='alpha'
export PARAM_INC=$(echo "scale=$float_scale ;($PARAM_FINAL-$PARAM_INIT)/$PARAM_NOS" | bc -q 2>/dev/null)

if [ ! -e "$RUNPROG" ]; then
  echo 
  echo "Error! runtime binary not found ..."
  echo 
  exit
fi

for i in $(seq $PARAM_NOS)
do
  
  case "$VARYING_PARAM" in
    "$S1")
    #Frequency
    export OMEGA=$(echo "scale=$float_scale ;$PARAM_INIT + $i*$PARAM_INC" | bc -q 2>/dev/null)
    #Output file path and name
    OUTDIR="result_"
    OUTDIR=$OUTDIR:"$VARYING_PARAM"
    OUTDIR=$OUTDIR:"_"
    OUTDIR=$OUTDIR:"$OMEGA"
    OUTDIR=$OUTDIR:"_seed"
    OUTDIR=$OUTDIR:"$GSL_RNG_SEED"
    export OUTDIR
    if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
    fi
    ;;  
    "$S2")
    #Amplitude
    export AMPLITUDE=$(echo "scale=$float_scale ;$PARAM_INIT + $i*$PARAM_INC" | bc -q 2>/dev/null)
    #Output file path and name
    OUTDIR="result_"
    OUTDIR=$OUTDIR:"$VARYING_PARAM"
    OUTDIR=$OUTDIR:"_"
    OUTDIR=$OUTDIR:"$AMPLITUDE"
    OUTDIR=$OUTDIR:"_seed"
    OUTDIR=$OUTDIR:"$GSL_RNG_SEED"
    export OUTDIR
    if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
    fi
    ;;
    "$S3")
    export PERT_AMP=$(echo "scale=$float_scale ;$PARAM_INIT + $i*$PARAM_INC" | bc -q 2>/dev/null)
    #Output file path and name
    OUTDIR="result_"
    OUTDIR=$OUTDIR:"$VARYING_PARAM"
    OUTDIR=$OUTDIR:"_"
    OUTDIR=$OUTDIR:"$PERT_AMP"
    OUTDIR=$OUTDIR:"_seed"
    OUTDIR=$OUTDIR:"$GSL_RNG_SEED"
    export OUTDIR
    if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
    fi
    ;;
    *)
    echo 
    echo "Error! Please edit the script and select the varying parameter: 'omega', 'ampl' or 'alpha' ..."
    echo
    exit
    ;;
  esac

  # Run the MPI executable. Use TACC's affinity utility to position one task on each socket.
  # See http://www.tacc.utexas.edu/user-services/user-guides/lonestar-user-guide#running:numa:socket
  ibrun tacc_affinity $RUNPROG $VERBOSITY -a $AMPLITUDE -f $OMEGA -p $PERT_AMP -t $FINALTIME -n $TIMESTEPS -m "${OUTDIR}/${MAG_FILE}" -e "${OUTDIR}/${ENT_FILE}" 2>$ERROR_FILE
done
