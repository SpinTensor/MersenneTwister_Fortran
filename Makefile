FC = gfortran

FCFLAGS = -O3 -g -fcheck=all -fbacktrace \
	  -fmax-errors=3 -fimplicit-none \
	  -D__HAS_QP -D__HAS_IQP -D__MAX_RANK=15

FCWFLAGS = -Wall -Wextra -Wpedantic -Wsurprising

OUTLIB = ./libMTrandom.a

all: $(OUTLIB)

$(OUTLIB): MTrandom.o \
           kinds.o
	ar rvcs $@ $^

MTrandom.o: MTrandom.F90 kinds.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -c $<

kinds.o: kinds.F90
	$(FC) $(FCFLAGS) $(FCWFLAGS) -c $<

.PHONY: clean distclean

clean:
	rm -f *.o *.mod 

distclean: clean
	rm -f $(OUTLIB)
