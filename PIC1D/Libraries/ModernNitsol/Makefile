
FC = ifort
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG= -traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= BoundPlasmaExample.o mod_simulation.o mod_potentialSolver.o mod_particleInjection.o mod_NullCollision.o mod_targetParticle.o mod_particle.o mod_domain.o mod_BasicFunctions.o constants.o #    mod_collisions.o   mod_potentialSolver.o 
RM = rm -f

all: clean mod_nitsol.o

debug: clean setDebug mod_nitsol.o

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -qopenmp -g -CB -traceback -check all -fp-stack-check -warn all -heap-arrays)

mod_nitsol.o: mod_nitsol.f90
	$(FC) $(BASE) -c  -o $@ $<

clean:
	$(RM) *.o *.mod
