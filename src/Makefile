FC = ifort
CC=mpiifort
OPT= $(BASE)
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -xHost -qmkl # -mtune=skylake -xCORE-AVX512 -mt_mpi -heap-arrays 8192
DEBUG= $(BASE) -traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OBJ= constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o mod_potentialSolver.o mod_collisions.o mod_simulation.o BoundPlasmaExample.o
RM = rm -f

all: cleanData clean BoundPlasmaExample

BoundPlasmaExample: $(OBJ)
	$(FC) $(DEBUG) -o $@ $(OBJ)

BoundPlasmaExample.o: BoundPlasmaExample.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o mod_potentialSolver.o mod_collisions.o mod_simulation.o
	$(FC) $(DEBUG) -c  -o $@ $<

mod_simulation.o: mod_simulation.f90 mod_collisions.o mod_potentialSolver.o mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(DEBUG) -c  -o $@ $<

mod_collisions.o: mod_collisions.f90 mod_potentialSolver.o mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(DEBUG) -c  -o $@ $<

mod_potentialSolver.o: mod_potentialSolver.f90 mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(DEBUG) -c  -o $@ $<

mod_domain.o: mod_domain.f90 constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(DEBUG) -c  -o $@ $<

mod_particle.o: mod_particle.f90 constants.o mod_BasicFunctions.o
	$(FC) $(DEBUG) -c  -o $@ $<

mod_BasicFunctions.o: mod_BasicFunctions.f90 constants.o
	$(FC) $(DEBUG) -c  -o $@ $<

constants.o: constants.f90
	$(FC) $(DEBUG) -c  -o $@ $<



clean:
	$(RM) BoundPlasmaExample *.o *.mod

cleanData:
	$(RM) ../Data/*


