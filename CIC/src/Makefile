FC = ifort
CC=mpiifort
OPT= $(BASE)
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -ffloat-store -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG= -traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OBJ= constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o mod_potentialSolver.o mod_collisions.o mod_simulation.o BoundPlasmaExample.o
RM = rm -f

all: cleanData clean BoundPlasmaExample

BoundPlasmaExample: $(OBJ)
	$(FC) $(OPT) -o $@ $(OBJ)

BoundPlasmaExample.o: BoundPlasmaExample.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o mod_potentialSolver.o mod_collisions.o mod_simulation.o
	$(FC) $(OPT) -c  -o $@ $<

mod_simulation.o: mod_simulation.f90 mod_collisions.o mod_potentialSolver.o mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(OPT) -c  -o $@ $<

mod_collisions.o: mod_collisions.f90 mod_potentialSolver.o mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(OPT) -c  -o $@ $<

mod_potentialSolver.o: mod_potentialSolver.f90 mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(OPT) -c  -o $@ $<

mod_particle.o: mod_particle.f90 constants.o mod_BasicFunctions.o mod_domain.o
	$(FC) $(OPT) -c  -o $@ $<

mod_domain.o: mod_domain.f90 constants.o mod_BasicFunctions.o
	$(FC) $(OPT) -c  -o $@ $<

mod_BasicFunctions.o: mod_BasicFunctions.f90 constants.o
	$(FC) $(OPT) -c  -o $@ $<

constants.o: constants.f90
	$(FC) $(OPT) -c  -o $@ $<



clean:
	$(RM) BoundPlasmaExample *.o *.mod

cleanData:
	$(RM) ../Data/*.dat
	$(RM) ../Data/Density/*.dat
	$(RM) ../Data/Phi/*.dat
	$(RM) ../Data/PhaseSpace/*.dat
	$(RM) ../Data/ElectronTemperature/*.dat


