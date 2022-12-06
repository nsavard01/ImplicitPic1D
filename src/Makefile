FC = ifort
FCFLAGS= -O3 -no-prec-div -m64 -auto -fp-model fast=2 -unroll-aggressive -xHost #-qmkl -qopenmp
OBJS = constants.o mod_BasicFunctions.o mod_particle.o mod_domain.o

.PHONY: all clean
.SUFFIXES: .f90 .o

all: BoundPlasmaExample

BoundPlasmaExample: BoundPlasmaExample.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFlags) $<

%.o: %.mod

clean:
	$(RM) BoundPlasmaExample *.o *.mod
