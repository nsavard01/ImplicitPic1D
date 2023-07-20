FC = ifort
CC=mpiifort
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -ffloat-store -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG = -qmkl -g -traceback -check all -fp-stack-check -warn all #-traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= BoundPlasmaExample.o mod_simulation.o mod_nonLinSolvers.o ../../Libraries/Nitsol/libnitsol.a mod_particleMover.o mod_potentialSolver.o mod_Scheme.o mod_domain.o mod_collisions.o mod_particle.o mod_BasicFunctions.o constants.o
RM = rm -f

all: cleanData clean cleanNitsol BoundPlasmaExample

debug: cleanData clean cleanNitsol setDebug BoundPlasmaExample

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -g -traceback -check all -fp-stack-check -warn all)
	

BoundPlasmaExample: $(OBJ)
	$(FC) $(BASE) -o $@ $(OBJ)

BoundPlasmaExample.o: ../../SharedModules/src/BoundPlasmaExample.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o mod_Scheme.o mod_potentialSolver.o mod_collisions.o mod_particleMover.o mod_simulation.o mod_nonLinSolvers.o ../../Libraries/Nitsol/libnitsol.a
	$(FC) $(BASE) -c  -o $@ $<

mod_simulation.o: ../../SharedModules/src/mod_simulation.f90  mod_collisions.o mod_potentialSolver.o mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o mod_Scheme.o mod_particleMover.o mod_nonLinSolvers.o ../../Libraries/Nitsol/libnitsol.a
	$(FC) $(BASE) -c  -o $@ $<

mod_nonLinSolvers.o: ../../SharedModules/src/mod_nonLinSolvers.f90 ../../Libraries/Nitsol/libnitsol.a mod_particleMover.o mod_Scheme.o mod_potentialSolver.o mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o 
	$(FC) $(BASE) -c  -o $@ $<

mod_particleMover.o: mod_particleMover.f90 mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o mod_potentialSolver.o
	$(FC) $(BASE) -c  -o $@ $<

mod_potentialSolver.o: ../../SharedModules/src/mod_potentialSolver.f90 mod_domain.o constants.o mod_BasicFunctions.o
	$(FC) $(BASE) -c  -o $@ $<

mod_Scheme.o: mod_Scheme.f90 mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(BASE) -c  -o $@ $<

mod_domain.o: ../../SharedModules/src/mod_domain.f90 constants.o mod_BasicFunctions.o
	$(FC) $(BASE) -c  -o $@ $<

mod_collisions.o: ../../SharedModules/src/mod_collisions.f90 constants.o mod_BasicFunctions.o mod_particle.o
	$(FC) $(BASE) -c  -o $@ $<

mod_particle.o: ../../SharedModules/src/mod_particle.f90 constants.o mod_BasicFunctions.o
	$(FC) $(BASE) -c  -o $@ $<

mod_BasicFunctions.o: ../../SharedModules/src/mod_BasicFunctions.f90 constants.o
	$(FC) $(BASE) -c  -o $@ $<

constants.o: ../../SharedModules/src/constants.f90
	$(FC) $(BASE) -c  -o $@ $<

../../Libraries/Nitsol/libnitsol.a:
	cd ../../Libraries/Nitsol; $(MAKE) "FFLAGS=$(BASE)" "FC=$(FC)"

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


