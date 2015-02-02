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
BIN_NAMES=isingrand_parallel heff_res_parallel
DIST=Makefile isingrand.kdev4 *.c *.h makefiles scripts

#Dependencies
OBJECTS=initconds.o responses.o gsl_determinant_complex.o

.DEFAULT_GOAL := all

all: $(BIN_NAMES)

indent:
	indent *.c *.h
	${RM} *~

%.o:%.c %.h params.h
	$(CC_COMPILE)

isingrand_parallel: isingrand_parallel.o integrator_isingrand.o $(OBJECTS) params.h
	$(CC) isingrand_parallel.o integrator_isingrand.o $(OBJECTS) `$(GLIB_FLAGS)` $(LD_FLAGS) -fopenmp -o $@

heff_res_parallel: heff_res_parallel.o integrator_heffres.o $(OBJECTS) params.h
	$(CC) heff_res_parallel.o integrator_heffres.o $(OBJECTS) `$(GLIB_FLAGS)` $(LD_FLAGS) -fopenmp -o $@

dist:
	$(TAR) czvf $(NAME).tar.gz $(DIST)
	
allclean: clean resclean

clean:
	${RM} $(BIN_NAMES) result/*.dat* *~ *.o *out *.dat* *.tgz

resclean:
	${RM} *result/*
	${RM} -r *result_*
