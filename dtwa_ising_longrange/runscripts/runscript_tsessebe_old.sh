#!/bin/bash
#########################################################################
## Name of my job
#PBS -N test
#PBS -l walltime=1:00:00
#########################################################################
##Export all PBS environment variables
#PBS -V
#########################################################################
##Output file. Combine stdout and stderr into one
#PBS -o stdout.dat
#PBS -e stderr.dat
#PBS -j oe 
#########################################################################
##Number of nodes and procs per node.
##See docs at http://wiki.chpc.ac.za/howto:pbs-pro_job_submission_examples 
#PBS -l select=5:ncpus=8:mpiprocs=8

#########################################################################
##Send me email when my job aborts, begins, or ends
#PBS -m ea
#PBS -M daneel@sun.ac.za
#########################################################################

# Make sure I'm the only one that can read my output
umask 0077
# Load the module system
source /etc/profile.d/modules.sh
#Load relevant modules. Load them with THESE TWO LINES, NOT FROM ONE LINE
module load dot intel
module load gcc/4.9.1 Anaconda/2.1.0

cd $PBS_O_WORKDIR

#########################################################################
# How many cores total do we have?
NO_OF_CORES=$(cat $PBS_NODEFILE | wc -l)
#########################################################################
##Parameter values
export AMPL=1.0
export OMEGA=0.0
export BETA=1.0
export HX=0.2
export HY=0.0
export HZ=0.0
export JX=0.0
export JY=0.0
export JZ=1.0

##Output filenames
export OUTFILE_MAGX=sx_time_for_beta_${BETA}_hx_${HX}_hy_${HY}_hz_${HZ}_jx_${JX}_jy_${JY}_jz_${JZ}.txt
export OUTFILE_MAGY=sy_time_for_beta_${BETA}_hx_${HX}_hy_${HY}_hz_${HZ}_jx_${JX}_jy_${JY}_jz_${JZ}.txt
export OUTFILE_MAGZ=sz_time_for_beta_${BETA}_hx_${HX}_hy_${HY}_hz_${HZ}_jx_${JX}_jy_${JY}_jz_${JZ}.txt

export OUTFILE_SXVAR=sxvar_time_for_beta_${BETA}_hx_${HX}_hy_${HY}_hz_${HZ}_jx_${JX}_jy_${JY}_jz_${JZ}.txt
export OUTFILE_SYVAR=syvar_time_for_beta_${BETA}_hx_${HX}_hy_${HY}_hz_${HZ}_jx_${JX}_jy_${JY}_jz_${JZ}.txt
export OUTFILE_SZVAR=szvar_time_for_beta_${BETA}_hx_${HX}_hy_${HY}_hz_${HZ}_jx_${JX}_jy_${JY}_jz_${JZ}.txt

export OUTFILE_SXYVAR=sxyvar_time_for_beta_${BETA}_hx_${HX}_hy_${HY}_hz_${HZ}_jx_${JX}_jy_${JY}_jz_${JZ}.txt
export OUTFILE_SXZVAR=sxzvar_time_for_beta_${BETA}_hx_${HX}_hy_${HY}_hz_${HZ}_jx_${JX}_jy_${JY}_jz_${JZ}.txt
export OUTFILE_SYZVAR=syzvar_time_for_beta_${BETA}_hx_${HX}_hy_${HY}_hz_${HZ}_jx_${JX}_jy_${JY}_jz_${JZ}.txt

export NITER=2000

export SCRIPT="./dtwa_ising_longrange.py"
#########################################################################
##Now, run the code
BEGINTIME=$(date +"%s")
mpirun -np $NO_OF_CORES -machinefile $PBS_NODEFILE  python $SCRIPT \
    -v -pbc -omx $OUTFILE_MAGX -omy $OUTFILE_MAGY -omz $OUTFILE_MAGZ \
	 -ox $OUTFILE_SXVAR -oy $OUTFILE_SYVAR -oz $OUTFILE_SZVAR \
		-oxy $OUTFILE_SXYVAR -oxz $OUTFILE_SXZVAR -oyz $OUTFILE_SYZVAR \
		 	-a $AMPL -w $OMEGA -x $HX -y $HY -z $HZ -jx $JX -jy $JY -jz $JZ -t $NITER -b $BETA
ENDTIME=$(date +"%s")
ELAPSED_TIME=$(($ENDTIME-$BEGINTIME))

echo "#Runtime: $(($ELAPSED_TIME / 60)) minutes and $(($ELAPSED_TIME % 60)) seconds."
