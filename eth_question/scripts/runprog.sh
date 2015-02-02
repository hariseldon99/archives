################################################################################################################################
########################################Default Sun Grid Engine options.########################################################
################################################################################################################################
#$ -S /bin/bash
#$ -V                           # Inherit the submission environment
#$ -cwd                         # Start job in  submission directory
#$ -N rigol_lattice             # Job Name
#$ -j y                         # Combine stderr & stdout into stdout  
#$ -pe make 30                  # Requests number of SLOTS
				# The 'make' pe launches in round robin fashion one slot per node
#$ -o $JOB_NAME.o$JOB_ID        # Name of the output file (eg. myMPI.oJobID)
#$ -q all.q                     # Queue name
################################################################################################################################
################################################################################################################################

export LATSIZE=5
export NUMPTCLS=2
export REPULSION=0.1
export DELTA_E=0.1
export RUNPROG="./rigol_lattice"
export NUMPROCS=30

mpirun -np $NUMPROCS $RUNPROG -lattice_size $LATSIZE -vector_size $NUMPTCLS -repulsion $REPULSION
