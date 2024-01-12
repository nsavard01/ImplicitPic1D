FC = ifort
FCFLAGS= -qopenmp -qmkl -O3 -ipo -no-prec-div -m64 -auto -fp-model fast=2 -unroll-aggressive  -mtune=skylake -heap-arrays 8192 -xHost

.PHONY: all clean
.SUFFIXES: .f90 .o

all: clean main

main: main.f90
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFlags) $<

%.o: %.mod

clean:
	$(RM) main *.o *.mod *.DAT
