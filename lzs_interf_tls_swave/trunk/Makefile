CC=`which mpicc`
GCC=`which gcc`
GLIB_LOC = pkg-config --cflags --libs glib-2.0
GSL_LOC = gsl-config --cflags
GLIB_FLAGS = pkg-config --libs glib-2.0
CFLAGS= -Wall -O3 -c
GSL_FLAGS= gsl-config --libs
#Binary for simple tls
BIN_TLS=lzs_tls
#binary for the lzs problem in a mbt system
BIN_MBT_2D=lzs_mbt_2d_swave
SCRIPT_TLS=runlzs_tls
SCRIPT_MBT_2D=runlzs_mbt_2d_local
INSTALL=/usr/bin/install
WORK=${PWD}
RESULT_DIR=result

#Dependencies
OBJECTS= lzs_lalgebra.o lzs_io.o lzs_normstats.o lattice.o

%.o:%.c %.h params_ics.h
	$(CC) `$(GLIB_LOC)` `$(GSL_LOC)` $(CFLAGS) $<

all: $(BIN_TLS) $(BIN_MBT_2D)

$(BIN_TLS): lzs_tls.o lzs_tls_integrator.o $(OBJECTS) params_ics.h
	$(CC) lzs_tls.o lzs_tls_integrator.o $(OBJECTS) `$(GLIB_FLAGS)` `$(GSL_FLAGS)` -o $@

$(BIN_MBT_2D): lzs_mbt_2d_swave.o lzs_mbt_integrator.o $(OBJECTS) params_ics.h
	$(CC) lzs_mbt_2d_swave.o lzs_mbt_integrator.o $(OBJECTS) `$(GLIB_FLAGS)` `$(GSL_FLAGS)` -o $@

libclean:
	${RM} *~ *.o

binclean:
	${RM} $(BIN_TLS) $(BIN_MBT_2D)

resclean:
	${RM} $(WORK)/$(RESULT_DIR)/*.dat $(WORK)/$(RESULT_DIR)/*.h5 $(WORK)/$(RESULT_DIR)/*.txt*

allclean: libclean binclean resclean
