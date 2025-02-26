
COMPILER_C = icpx
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -xHost #-xCORE-AVX512 -mt_mpi 
BASE_C = -O3 -ipo -xHost -fast -qmkl -qopenmp
OBJ= main.o null_collision.o target_particle.o particle.o domain.o pcg_rng.o global_inputs.o Constants.o
RM = rm -f

all: main

debug: setDebug main

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE_C = -qmkl -qopenmp -g -Wall -fexceptions -fsanitize=undefined)

main: $(OBJ)
	$(COMPILER_C) $(BASE_C) -o $@ $(OBJ)

main.o: main.cpp null_collision.o target_particle.o particle.o domain.o pcg_rng.o global_inputs.o Constants.o
	$(COMPILER_C) $(BASE_C) -c -o $@ $<	

null_collision.o: null_collision.cpp particle.o domain.o global_inputs.o pcg_rng.o Constants.o
	$(COMPILER_C) $(BASE_C) -c -o $@ $<	

target_particle.o: target_particle.cpp particle.o domain.o global_inputs.o pcg_rng.o Constants.o
	$(COMPILER_C) $(BASE_C) -c -o $@ $<	

particle.o: particle.cpp domain.o global_inputs.o pcg_rng.o Constants.o
	$(COMPILER_C) $(BASE_C) -c -o $@ $<	

domain.o: domain.cpp global_inputs.o Constants.o
	$(COMPILER_C) $(BASE_C) -c -o $@ $<	

pcg_rng.o: pcg_rng.cpp Constants.o
	$(COMPILER_C) $(BASE_C) -c -o $@ $<	

global_inputs.o: global_inputs.cpp Constants.o
	$(COMPILER_C) $(BASE_C) -c -o $@ $<	

Constants.o: Constants.cpp
	$(COMPILER_C) $(BASE_C) -c -o $@ $<	


clean:
	$(RM) main *.o *.mod



