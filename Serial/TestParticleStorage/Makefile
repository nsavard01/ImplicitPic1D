FC = ifort
FCFLAGS= -O3 -ipo -no-prec-div -m64 -qopenmp -auto -fp-model fast=2 -unroll-aggressive -xHost -qmkl
OBJS = mod_test.o

.PHONY: all clean
.SUFFIXES: .f90 .o

all: clean main

main: main.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFlags) $<

%.o: %.mod

clean:
	$(RM) main *.o *.mod
