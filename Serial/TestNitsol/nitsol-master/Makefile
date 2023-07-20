ALL: bratu cavity pormed

# Macro definitions.

FC 		= ifort
FFLAGS 		= -qmkl -g -traceback -check all -fp-stack-check
FLINKER 	= ifort

NITSOL = ./Nitsol/libnitsol.a


# Default compilation rules.

.f.o:
	$(FC) -c $(FFLAGS) $*.f

# Rules to build libraries.

nitsol_lib:
	cd Nitsol; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

# lapack_lib:
# 	cd Lapack; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

# blas_lib:
# 	cd Blas; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

# Rules to build test programs.

bratu_obj:
	cd Bratu; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

bratu:  nitsol_lib bratu_obj
	$(FLINKER) $(FFLAGS) -o bratu Bratu/*.o $(NITSOL)

cavity_obj:
	cd Cavity; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

cavity:  nitsol_lib cavity_obj
	$(FLINKER) $(FFLAGS) -o cavity Cavity/*.o $(NITSOL)

pormed_obj:
	cd Pormed; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

pormed:  nitsol_lib pormed_obj
	$(FLINKER) $(FFLAGS) -o pormed Pormed/*.o $(NITSOL)

clean:
	- /bin/rm -f bratu cavity pormed core
	cd Bratu; $(MAKE) clean
	cd Cavity; $(MAKE) clean
	cd Pormed; $(MAKE) clean

veryclean:
	make clean
	cd Nitsol; $(MAKE) clean
	# cd Lapack; $(MAKE) clean
	# cd Blas; $(MAKE) clean

remake:
	make clean
	make

help:
	cat make.help
