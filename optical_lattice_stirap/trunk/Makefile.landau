CC=/usr/bin/gcc
CFLAGS= -Wall -O3 -c
LD_FLAGS=-lgsl -lgslcblas -lm
BIN=stirap_oplattice 
SCRIPT=runstirap_landau
INSTALL=/usr/bin/install
WORK=${PWD}
RESULT_DIR=result
CC_COMPILE=$(CC) $(CFLAGS) $< 

#Dependencies
OBJECTS=quasienergies.o integrator.o algebra.o 																																					quantum_hamilt_melements.o hamiltonian.o 																																						final_output.o

%.o:%.c %.h params.h

stirap_oplattice: stirap_oplattice.o $(OBJECTS) params.h
	$(CC) stirap_oplattice.o $(OBJECTS) $(LD_FLAGS) -o $@

stirap_oplattice.o: stirap_optical_lattice.c stirap_optical_lattice.h params.h
	$(CC_COMPILE)
	mv stirap_optical_lattice.o $@

install: 
	$(INSTALL) $(BIN) $(SCRIPT) $(WORK)
	mkdir $(WORK)/$(RESULT_DIR)	

uninstall:
	${RM} $(WORK)/$(BIN) $(WORK)/$(SCRIPT)

clean:
	${RM} $(BIN) *~ *.o *out*.dat

distclean:
	${RM} $(BIN) *~ *.o *out*.dat
	${RM} $(WORK)/$(RESULT_DIR)/*
