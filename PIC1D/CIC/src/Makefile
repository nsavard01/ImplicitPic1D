FC = ifort
CC=mpiifort
BASE= -qopenmp -qmkl -O3 -ipo -no-prec-div -m64 -auto -fp-model fast=2 -unroll-aggressive -mtune=skylake -heap-arrays 8192 -xHost  #-xCORE-AVX512 -mt_mpi -qopt-report=1
DEBUG = -qmkl -g -traceback -check all -fp-stack-check -warn all #-traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= main.o mod_simulation.o mod_nonLinSolvers.o mod_nitsol.o mod_particleMover.o mod_potentialSolver.o mod_particleInjection.o mod_Scheme.o mod_NullCollision.o mod_targetParticle.o mod_particle.o mod_domain.o mod_BasicFunctions.o constants.o 
RM = rm -f

all: clean cleanNitsol main

debug: clean cleanNitsol setDebug main

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -qopenmp -g -heap-arrays 8192 -traceback -check all -fp-stack-check -warn all)
	

main: $(OBJ)
	$(FC) $(BASE) -o $@ $(OBJ)

main.o: ../../SharedModules/src/mainImplicit.f90 constants.o mod_BasicFunctions.o mod_domain.o mod_particle.o mod_targetParticle.o mod_NullCollision.o mod_particleInjection.o mod_Scheme.o mod_potentialSolver.o mod_particleMover.o mod_nitsol.o mod_nonLinSolvers.o mod_simulation.o
	$(FC) $(BASE) -c  -o $@ $<

mod_simulation.o: ../../SharedModules/src/mod_simulation.f90 constants.o mod_BasicFunctions.o mod_domain.o mod_particle.o mod_targetParticle.o mod_NullCollision.o mod_Scheme.o mod_particleInjection.o mod_potentialSolver.o mod_particleMover.o mod_nitsol.o mod_nonLinSolvers.o
	$(FC) $(BASE) -c  -o $@ $<

mod_nonLinSolvers.o: ../../SharedModules/src/mod_nonLinSolvers.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o mod_Scheme.o mod_potentialSolver.o mod_particleMover.o mod_nitsol.o
	$(FC) $(BASE) -c  -o $@ $<

mod_particleMover.o: mod_particleMoverPicard.f90 mod_domain.o constants.o mod_BasicFunctions.o mod_particle.o mod_potentialSolver.o
	$(FC) $(BASE) -c  -o $@ $<

mod_potentialSolver.o: mod_potentialSolver.f90 constants.o mod_BasicFunctions.o mod_domain.o
	$(FC) $(BASE) -c  -o $@ $<

mod_particleInjection.o: ../../SharedModules/src/mod_particleInjection.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o
	$(FC) $(BASE) -c  -o $@ $<

mod_Scheme.o: mod_Scheme.f90 constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o
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

mod_nitsol.o: ../../SharedModules/src/mod_nitsol.f90
	$(FC) $(BASE) -c  -o $@ $<

cleanNitsol:
	cd ../../Libraries/Nitsol; $(MAKE) clean

clean:
	$(RM) main *.o *.mod



