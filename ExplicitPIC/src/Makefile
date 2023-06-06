FC = ifort
CC=mpiifort
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -ffloat-store -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG= -traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= TwoStreamInstability.o mod_simulation.o mod_potentialSolver.o mod_collisions.o mod_domain.o mod_particle.o mod_BasicFunctions.o constants.o     
RM = rm -f

all: cleanData clean TwoStreamInstability

debug: cleanData clean setDebug TwoStreamInstability

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -g -traceback -check all -fp-stack-check -warn all)

TwoStreamInstability: $(OBJ)
	$(FC) $(BASE) -o $@ $(OBJ)

TwoStreamInstability.o: TwoStreamInstability.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_collisions.o mod_domain.o mod_potentialSolver.o mod_simulation.o
	$(FC) $(BASE) -c  -o $@ $<

mod_simulation.o: mod_simulation.f90 mod_collisions.o mod_potentialSolver.o constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(BASE) -c  -o $@ $<

mod_potentialSolver.o: mod_potentialSolver.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o
	$(FC) $(BASE) -c  -o $@ $<

mod_collisions.o: ../../SharedModules/src/mod_collisions.f90 constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(BASE) -c  -o $@ $<

mod_particle.o: ../../SharedModules/src/mod_particle.f90 constants.o mod_BasicFunctions.o
	$(FC) $(BASE) -c  -o $@ $<

mod_domain.o: mod_domain.f90 constants.o mod_BasicFunctions.o
	$(FC) $(BASE) -c  -o $@ $<

mod_BasicFunctions.o: ../../SharedModules/src/mod_BasicFunctions.f90 constants.o
	$(FC) $(BASE) -c  -o $@ $<

constants.o: ../../SharedModules/src/constants.f90
	$(FC) $(BASE) -c  -o $@ $<



clean:
	$(RM) TwoStreamInstability *.o *.mod

cleanData:
	$(RM) ../Data/*.dat
	$(RM) ../Data/Density/*.dat
	$(RM) ../Data/Phi/*.dat
	$(RM) ../Data/PhaseSpace/*.dat
	$(RM) ../Data/ElectronTemperature/*.dat


