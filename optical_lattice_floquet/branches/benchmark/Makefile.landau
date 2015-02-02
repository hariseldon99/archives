CC=/usr/bin/mpicc
CFLAGS= -Wall -O3 -c
LD_FLAGS=-lgsl -lgslcblas -lm
BIN=benchmark
MV=/bin/mv
SCRIPT=runprog_landau
INSTALL=/usr/bin/install
WORK=${PWD}
RESULT_DIR=result

#Dependencies
OBJECTS=algebra.o timer.o benchmark.o hamiltonian.o integrator.o quantum_hamilt_melements.o

%.o:%.c %.h params.h

benchmark: $(OBJECTS) params.h
	${CC} -o $@ $(OBJECTS) $(LD_FLAGS)

install: 
	$(INSTALL) $(BIN) $(SCRIPT) $(WORK)
	mkdir $(WORK)/$(RESULT_DIR)	

uninstall:
	${RM} $(WORK)/$(BIN) $(WORK)/$(SCRIPT)

clean:
	${RM} $(BIN) *~ *.o
