FC = gfortran
FCFLAGS= -O3 -Wall
OBJS = constants.o mod_particle.o mod_domain.o

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
