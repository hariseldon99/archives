WORK=${PWD}
TAR=/bin/tar
USE_PDFLATEX=true
BIBTEX=true
VIEWPDF=okular
NAME=notes
TARNAME=periodic_delocalization
EXTRA_DIST+= *.c *mk Makefile
CLEAN_FILES+=*.backup *converted-to.pdf *.bbl

include latex.gmk

CFLAGS=

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

.DEFAULT_GOAL := pdeloc

indent: 
	indent *.c
	rm -f *~

pdeloc: pdeloc.o  chkopts
	-${CLINKER} -o pdeloc pdeloc.o  ${PETSC_TS_LIB}  
	${RM} pdeloc.o

tarball: dist
	mv $(NAME).tar.gz $(TARNAME).tar.gz

allclean: clean
	rm -f pdeloc ex4 *.tmp

include ${PETSC_DIR}/conf/test
