include Make.include

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
