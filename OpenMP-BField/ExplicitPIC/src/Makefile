FC = ifort
CC=mpiifort
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG= -traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= BoundPlasmaExample.o mod_simulation.o mod_potentialSolver.o mod_particleInjection.o mod_NullCollision.o mod_targetParticle.o mod_particle.o mod_domain.o mod_BasicFunctions.o constants.o #    mod_collisions.o   mod_potentialSolver.o 
RM = rm -f

all: cleanData clean BoundPlasmaExample

debug: cleanData clean setDebug BoundPlasmaExample

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -qopenmp -g -CB -traceback -check all -fp-stack-check -warn all -heap-arrays)

BoundPlasmaExample: $(OBJ)
	$(FC) $(BASE) -o $@ $(OBJ)

BoundPlasmaExample.o: BoundPlasmaExample.f90 constants.o mod_BasicFunctions.o mod_domain.o mod_particle.o mod_targetParticle.o mod_NullCollision.o mod_particleInjection.o mod_potentialSolver.o mod_simulation.o# #mod_collisions.o 
	$(FC) $(BASE) -c  -o $@ $<

mod_simulation.o: mod_simulation.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o mod_particleInjection.o mod_potentialSolver.o
	$(FC) $(BASE) -c  -o $@ $<

# mod_readInputs.o: ../../SharedModules/src/mod_readInputs.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_targetParticle.o mod_NullCollision.o mod_domain.o mod_particleInjection.o mod_potentialSolver.o
# 	$(FC) $(BASE) -c  -o $@ $<

mod_potentialSolver.o: mod_potentialSolver.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o
	$(FC) $(BASE) -c  -o $@ $<

mod_particleInjection.o: ../../SharedModules/src/mod_particleInjection.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o
	$(FC) $(BASE) -c  -o $@ $<

mod_NullCollision.o: ../../SharedModules/src/mod_NullCollision.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_targetParticle.o
	$(FC) $(BASE) -c  -o $@ $<

mod_targetParticle.o: ../../SharedModules/src/mod_targetParticle.f90 constants.o mod_BasicFunctions.o
	$(FC) $(BASE) -c  -o $@ $<

mod_particle.o: ../../SharedModules/src/mod_particle.f90 constants.o mod_BasicFunctions.o mod_domain.o
	$(FC) $(BASE) -c  -o $@ $<

mod_domain.o: mod_domain.f90 constants.o mod_BasicFunctions.o
	$(FC) $(BASE) -c  -o $@ $<

mod_BasicFunctions.o: ../../SharedModules/src/mod_BasicFunctions.f90 constants.o
	$(FC) $(BASE) -c  -o $@ $<

constants.o: ../../SharedModules/src/constants.f90
	$(FC) $(BASE) -c  -o $@ $<



clean:
	$(RM) BoundPlasmaExample *.o *.mod

cleanData:
	$(RM) ../Data/*.dat
	$(RM) ../Data/Density/*.dat
	$(RM) ../Data/Phi/*.dat
	$(RM) ../Data/PhaseSpace/*.dat
	$(RM) ../Data/ElectronTemperature/*.dat


