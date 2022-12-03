FC = ifort
FCFLAGS= -qopenmp -O3 -Wall
OBJS = mod_test.o

.PHONY: all clean
.SUFFIXES: .f90 .o

all: main

main: main.f90 $(OBJS)
	$(FC) $(FCFLAGS) $< $(OBJS) -o $@

.f90.o:
	$(FC) -c $(FCFlags) $<

%.o: %.mod

clean:
	$(RM) main *.o *.mod
