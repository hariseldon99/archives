CC=`which mpicc`
GCC=`which gcc`
GLIB_LOC = pkg-config --cflags --libs glib-2.0
GSL_LOC = gsl-config --cflags
GLIB_FLAGS = pkg-config --libs glib-2.0
CFLAGS= -Wall -O3 -c
GSL_FLAGS= -lgsl -lcblas -latlas -lm
#Uncomment for using GSL BLAS implementation
#GSL_FLAGS=-lgsl -lgslcblas -lm
OUTDIR="result"

%.o:%.c %.h
	$(CC) `$(GLIB_LOC)` `$(GSL_LOC)` $(CFLAGS) $<

graphenehop: graphene-hopping.o
	-$(CC) -o $@ graphene-hopping.o `$(GLIB_FLAGS)` $(GSL_FLAGS) 

libclean:
	${RM} *~ graphene-hopping.o

clean: libclean
	${RM} graphenehop *.dat *.txt $(OUTDIR)/*