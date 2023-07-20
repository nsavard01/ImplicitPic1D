ALL: nitsolTest

# Macro definitions.

FC 		= ifort
FFLAGS 		= -qmkl -g -traceback -check all -fp-stack-check
FLINKER 	= ifort
NITSOL = ../nitsol-master/Nitsol/libnitsol.a


# Default compilation rules.


# Rules to build libraries.

nitsol_lib:
	cd ../nitsol-master/Nitsol; $(MAKE) "FFLAGS=$(FFLAGS)" "FC=$(FC)"

nitsolTest: nitsol_lib nitsolTest.o
	$(FLINKER) $(FFLAGS) -o nitsolTest nitsolTest.o $(NITSOL)

nitsolTest.o: nitsolTest.f90 ../nitsol-master/Nitsol/nitprint.h ../nitsol-master/Nitsol/nitparam.h
	$(FC) $(FFLAGS) -c  -o $@ $<

veryclean:
	make clean
	cd ../nitsol-master/Nitsol; $(MAKE) clean

remake:
	make clean
	make

clean:
	rm -f nitsolTest *.o
