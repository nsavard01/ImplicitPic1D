FC = ifort
CC=mpiifort
OPT= $(BASE)
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -ffloat-store -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG= -qmkl -g -traceback -check all -fp-stack-check#-traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OBJ= constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o mod_potentialSolver.o mod_collisions.o mod_Scheme.o mod_simulation.o BoundPlasmaExample.o ../../Libraries/Nitsol/libnitsol.a mod_nonLinSolvers.o
RM = rm -f

all: cleanData clean BoundPlasmaExample

BoundPlasmaExample: $(OBJ)
	$(FC) $(OPT) -o $@ $(OBJ)

BoundPlasmaExample.o: BoundPlasmaExample.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o mod_potentialSolver.o mod_collisions.o mod_Scheme.o mod_simulation.o mod_nonLinSolvers.o
	$(FC) $(OPT) -c  -o $@ $<

mod_simulation.o: ../../SharedModules/mod_simulation.f90 mod_Scheme.o mod_collisions.o mod_potentialSolver.o mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o mod_nonLinSolvers.o
	$(FC) $(OPT) -c  -o $@ $<

mod_Scheme.o: mod_Scheme.f90 mod_collisions.o mod_potentialSolver.o mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o mod_nonLinSolvers.o
	$(FC) $(OPT) -c  -o $@ $<

mod_nonLinSolvers.o: ../../SharedModules/mod_nonLinSolvers.f90 ../../Libraries/Nitsol/libnitsol.a mod_collisions.o mod_potentialSolver.o mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(OPT) -c  -o $@ $<

mod_potentialSolver.o: mod_potentialSolver.f90 mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(OPT) -c  -o $@ $<

mod_domain.o: mod_domain.f90 constants.o mod_BasicFunctions.o
	$(FC) $(OPT) -c  -o $@ $<

mod_collisions.o: ../../SharedModules/mod_collisions.f90 constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(OPT) -c  -o $@ $<

mod_particle.o: ../../SharedModules/mod_particle.f90 constants.o mod_BasicFunctions.o
	$(FC) $(OPT) -c  -o $@ $<

mod_BasicFunctions.o: ../../SharedModules/mod_BasicFunctions.f90 constants.o
	$(FC) $(OPT) -c  -o $@ $<

constants.o: ../../SharedModules/constants.f90
	$(FC) $(OPT) -c  -o $@ $<

../../Libraries/Nitsol/libnitsol.a:
	cd ../../Libraries/Nitsol; $(MAKE) "FFLAGS=$(OPT)" "FC=$(FC)"

cleanNitsol:
	cd ../../Libraries/Nitsol; $(MAKE) clean

clean:
	$(RM) BoundPlasmaExample *.o *.mod

cleanData:
	$(RM) ../Data/*.dat
	$(RM) ../Data/Density/*.dat
	$(RM) ../Data/Phi/*.dat
	$(RM) ../Data/PhaseSpace/*.dat
	$(RM) ../Data/ElectronTemperature/*.dat


