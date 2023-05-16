FC = gfortran
CC = gcc
FFLAGS = -c -fPIC
CFLAGS = -c -fPIC
f2pyFFLAGS = "-O4 -fdefault-real-8 -fdefault-double-8"

.PHONY: doc clean clean_all

GENERATED_FILES = scm_oce.F90

all : scm_gls.o scm_kpp.o scm_tke.o scm_oce.so
	@echo
	@echo ===================================================
	@echo       Fortran SCM tools OK
	@echo ===================================================
	@echo

scm_oce.so : scm_oce.F90
	f2py -c scm_oce.F90 scm_gls.o scm_kpp.o scm_tke.o --fcompiler=$(FC) --f90flags=$(f2pyFFLAGS) -m scm_oce

scm_gls.o : scm_gls.F90
	$(FC) -O4 -fdefault-real-8 -fdefault-double-8 -o scm_gls.o -c scm_gls.F90

scm_kpp.o : scm_kpp.F90
	$(FC) -O4 -fdefault-real-8 -fdefault-double-8 -o scm_kpp.o -c scm_kpp.F90

scm_tke.o : scm_tke.F90
	$(FC) -O4 -fdefault-real-8 -fdefault-double-8 -o scm_tke.o -c scm_tke.F90

clean:
	/bin/rm -f *.o *.mod *.so *.pyc *.msg* *~
clean_all:
	/bin/rm -f *.o *.mod *.so *.pyc *.msg* *~
