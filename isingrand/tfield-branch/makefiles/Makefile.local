CC=/usr/bin/mpicc
CFLAGS= -Wall -O3 -c
LD_FLAGS=-lgsl -lgslcblas -lm
GLIB_LOC = pkg-config --cflags --libs glib-2.0
GLIB_FLAGS = pkg-config --libs glib-2.0

#Note: To disable the use of ATLAS for BLAS operations, comment the line below
#LD_FLAGS=-lgsl -lcblas -latlas -lm

WORK=${PWD}
CC_COMPILE=$(CC) `$(GLIB_LOC)` $(CFLAGS) -fopenmp $< 
TAR=/bin/tar
NAME=isingrand
BIN_NAMES=isingrand_tfield

#Dependencies
OBJECTS=initconds.o responses.o gsl_determinant_complex.o

.DEFAULT_GOAL := all

all: $(BIN_NAMES)

indent:
	indent *.c *.h
	${RM} *~

%.o:%.c %.h params.h
	$(CC_COMPILE)

isingrand_tfield: isingrand_tfield.o integrator_tfield.o $(OBJECTS) params.h
	$(CC) isingrand_tfield.o integrator_tfield.o $(OBJECTS) `$(GLIB_FLAGS)` $(LD_FLAGS) -fopenmp -o $@

allclean: clean resclean

clean:
	${RM} $(BIN_NAMES) result/*.dat* *~ *.o *out *.dat* *.tgz

resclean:
	${RM} *result/*
	${RM} -r *result_*
