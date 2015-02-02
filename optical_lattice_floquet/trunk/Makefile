include ${SLEPC_DIR}/conf/slepc_common

CC=/usr/bin/mpicc
CFLAGS= -Wall -O3 -c
LD_FLAGS=-lgsl -lgslcblas -lm
BIN=floquet_oplattice 
SCRIPT=runfloquet_landau
INSTALL=/usr/bin/install
WORK=${PWD}
RESULT_DIR=result
CC_SLEPC=$(CC) $(CFLAGS) $< ${SLEPC_INCLUDE} ${PETSC_INCLUDE}

#Dependencies
OBJECTS=quasienergies.o integrator.o algebra.o 																																					quantum_hamilt_melements.o hamiltonian.o 																																						slepc_eval_nhep.o final_output.o

%.o:%.c %.h params.h

floquet_oplattice: floquet_oplattice.o $(OBJECTS) params.h
	${CLINKER} -o $@ floquet_oplattice.o $(OBJECTS) ${SLEPC_LIB} $(LD_FLAGS)

floquet_oplattice.o: floquet_optical_lattice.c floquet_optical_lattice.h params.h
	$(CC_SLEPC)
	${MV} floquet_optical_lattice.o $@

slepc_eval_nhep.o: slepc_eval_nhep.c slepc_eval_nhep.h params.h
	$(CC_SLEPC)

install: 
	$(INSTALL) $(BIN) $(SCRIPT) $(WORK)
	mkdir $(WORK)/$(RESULT_DIR)	

uninstall:
	${RM} $(WORK)/$(BIN) $(WORK)/$(SCRIPT)

cleanup:
	${RM} $(BIN) *~ *.o

distclean:
	${RM} $(BIN) *~ *.o
	${RM} $(WORK)/$(RESULT_DIR)/*
