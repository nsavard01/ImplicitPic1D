FC = ifort
CC=mpiifort
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -ffloat-store -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG = -qmkl -g -traceback -check all -fp-stack-check -warn all #-traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= constants.o mod_BasicFunctions.o mod_domain.o TestCombination.o 
RM = rm -f

all: clean TestCombination


setDebug: 
	$(eval BASE = -qmkl -g -traceback -check all -fp-stack-check -warn all)
	

TestCombination: $(OBJ)
	$(FC) $(DEBUG) -o $@ $(OBJ)

TestCombination.o: TestCombination.f90 constants.o mod_BasicFunctions.o mod_domain.o
	$(FC) $(DEBUG) -c  -o $@ $<

# mod_potentialSolver.o: mod_potentialSolver.f90 mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o
# 	$(FC) $(DEBUG) -c  -o $@ $<

mod_domain.o: mod_domain.f90 constants.o mod_BasicFunctions.o
	$(FC) $(DEBUG) -c  -o $@ $<

mod_BasicFunctions.o: mod_BasicFunctions.f90
	$(FC) $(DEBUG) -c  -o $@ $<

constants.o: constants.f90
	$(FC) $(DEBUG) -c  -o $@ $<

clean:
	$(RM) TestCombination *.o *.mod	



