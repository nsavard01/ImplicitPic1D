FC = ifort
CC=mpiifort
BASE= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -qmkl -mtune=skylake -heap-arrays 8192 -ffloat-store -xHost#-xCORE-AVX512 -mt_mpi 
DEBUG = -qmkl -g -traceback -check all -fp-stack-check -warn all #-traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= main.o mod_MT.o mod_Random.o mod_PCG.o  # mod_collisions.o
RM = rm -f

all: clean mtTest

debug: clean setDebug mtTest

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -qopenmp -g -traceback -check all -fp-stack-check -warn all)
	
mtTest: $(OBJ)
	$(FC) $(BASE) -o $@ $(OBJ)


main.o: mtTest.f90 mod_MT.o mod_Random.o mod_PCG.o
	$(FC) $(BASE) -c  -o $@ $<

mod_Random.o: mod_Random.f90
	$(FC) $(BASE) -c  -o $@ $<

mod_MT.o: mod_MT.f90
	$(FC) $(BASE) -c  -o $@ $<

mod_PCG.o: mod_PCG.f90
	$(FC) $(BASE) -c  -o $@ $<	

clean:
	$(RM) main *.o *.mod






