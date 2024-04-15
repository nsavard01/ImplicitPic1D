FC = ifort
CC=mpiifort
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG= -traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= MCCExample.o mod_readInputs.o mod_NullCollision.o mod_targetParticle.o mod_particle.o mod_BasicFunctions.o constants.o #    mod_collisions.o   mod_potentialSolver.o 
RM = rm -f

all: cleanData clean MCCExample

debug: cleanData clean setDebug MCCExample

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -qopenmp -g -CB -traceback -check all -fp-stack-check -warn all -heap-arrays)

MCCExample: $(OBJ)
	$(FC) $(BASE) -o $@ $(OBJ)

MCCExample.o: MCCExample.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_targetParticle.o mod_NullCollision.o mod_readInputs.o
	$(FC) $(BASE) -c  -o $@ $<

mod_readInputs.o: mod_readInputs.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_targetParticle.o mod_NullCollision.o
	$(FC) $(BASE) -c  -o $@ $<

mod_NullCollision.o: ../SharedModules/src/mod_NullCollision.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_targetParticle.o
	$(FC) $(BASE) -c  -o $@ $<

mod_targetParticle.o: ../SharedModules/src/mod_targetParticle.f90 constants.o mod_BasicFunctions.o
	$(FC) $(BASE) -c  -o $@ $<

mod_particle.o: ../SharedModules/src/mod_particle.f90 constants.o mod_BasicFunctions.o
	$(FC) $(BASE) -c  -o $@ $<

mod_BasicFunctions.o: ../SharedModules/src/mod_BasicFunctions.f90 constants.o
	$(FC) $(BASE) -c  -o $@ $<

constants.o: ../SharedModules/src/constants.f90
	$(FC) $(BASE) -c  -o $@ $<



clean:
	$(RM) MCCExample *.o *.mod

cleanData:
	$(RM) ../Data/*.dat
	$(RM) ../Data/Density/*.dat
	$(RM) ../Data/Phi/*.dat
	$(RM) ../Data/PhaseSpace/*.dat
	$(RM) ../Data/ElectronTemperature/*.dat


