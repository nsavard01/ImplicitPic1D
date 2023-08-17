FC = ifort
CC=mpiifort
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -ffloat-store -xHost#-xCORE-AVX512 -mt_mpi 
DEBUG = -qmkl -g -traceback -check all -fp-stack-check -warn all #-traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= generateModel.o Nitsol/libnitsol.a mod_BasicFunctions.o constants.o  # mod_collisions.o
RM = rm -f

all: cleanNitsol BoundPlasmaExample

debug: cleanNitsol setDebug BoundPlasmaExample

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -qopenmp -g -traceback -check all -fp-stack-check -warn all)
	
BoundPlasmaExample: $(OBJ)
	$(FC) $(BASE) -o $@ $(OBJ)


generateModel.o: generateModel.f90 Nitsol/libnitsol.a mod_BasicFunctions.o constants.o
	$(FC) $(BASE) -c  -o $@ $<

mod_BasicFunctions.o: ../OpenMP/SharedModules/src/mod_BasicFunctions.f90 constants.o
	$(FC) $(BASE) -c  -o $@ $<

constants.o: ../OpenMP/SharedModules/src/constants.f90
	$(FC) $(BASE) -c  -o $@ $<


Nitsol/libnitsol.a:
	cd Nitsol; $(MAKE) "FFLAGS=$(BASE)" "FC=$(FC)"

cleanNitsol:
	cd Nitsol; $(MAKE) clean





