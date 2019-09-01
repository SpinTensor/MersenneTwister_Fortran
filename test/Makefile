include ../Make.include

MTINC = -I../
MTLIB = -L../ -lMTrandom

all: test

test: runtest

runtest: runtest_random_integer_restart \
	 runtest_random_uniform_restart \
	 runtest_random_gaussian_restart \
	 runtest_random_integer_isp \
	 runtest_random_integer_idp \
	 runtest_random_integer_iqp \
	 runtest_random_uniform_sp \
	 runtest_random_uniform_dp \
	 runtest_random_uniform_qp \
	 runtest_random_gaussian_sp \
	 runtest_random_gaussian_dp \
	 runtest_random_gaussian_qp

buildtest: test_random_integer_restart.x \
           test_random_uniform_restart.x \
           test_random_gaussian_restart.x \
           test_random_integer_isp.x \
           test_random_integer_idp.x \
           test_random_integer_iqp.x \
           test_random_uniform_sp.x \
           test_random_uniform_dp.x \
           test_random_uniform_qp.x \
           test_random_gaussian_sp.x \
           test_random_gaussian_dp.x \
           test_random_gaussian_qp.x

runtest_random_integer_restart: test_random_integer_restart.x
	@(./test_random_integer_restart.x)
runtest_random_uniform_restart: test_random_uniform_restart.x
	@(./test_random_uniform_restart.x)
runtest_random_gaussian_restart: test_random_gaussian_restart.x
	@(./test_random_gaussian_restart.x)
runtest_random_integer_isp: test_random_integer_isp.x
	@(./test_random_integer_isp.x)
runtest_random_integer_idp: test_random_integer_idp.x
	@(./test_random_integer_idp.x)
runtest_random_integer_iqp: test_random_integer_iqp.x
	@(./test_random_integer_iqp.x)
runtest_random_uniform_sp: test_random_uniform_sp.x
	@(./test_random_uniform_sp.x)
runtest_random_uniform_dp: test_random_uniform_dp.x
	@(./test_random_uniform_dp.x)
runtest_random_uniform_qp: test_random_uniform_qp.x
	@(./test_random_uniform_qp.x)
runtest_random_gaussian_sp: test_random_gaussian_sp.x
	@(./test_random_gaussian_sp.x)
runtest_random_gaussian_dp: test_random_gaussian_dp.x
	@(./test_random_gaussian_dp.x)
runtest_random_gaussian_qp: test_random_gaussian_qp.x
	@(./test_random_gaussian_qp.x)

test_random_integer_restart.x: test_random_integer_restart.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_integer_restart.o: test_random_integer_restart.F90 
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -c $<
test_random_uniform_restart.x: test_random_uniform_restart.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_uniform_restart.o: test_random_uniform_restart.F90 
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -Wno-compare-reals -c $<
test_random_gaussian_restart.x: test_random_gaussian_restart.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_gaussian_restart.o: test_random_gaussian_restart.F90 
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -Wno-compare-reals -c $<
test_random_integer_isp.x: test_random_integer_isp.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_integer_isp.o: test_random_integer_isp.F90 
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -c $<
test_random_integer_idp.x: test_random_integer_idp.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_integer_idp.o: test_random_integer_idp.F90
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -c $<
test_random_integer_iqp.x: test_random_integer_iqp.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_integer_iqp.o: test_random_integer_iqp.F90
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -c $<
test_random_uniform_sp.x: test_random_uniform_sp.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_uniform_sp.o: test_random_uniform_sp.F90
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -Wno-compare-reals -Wno-conversion-extra -c $<
test_random_uniform_dp.x: test_random_uniform_dp.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_uniform_dp.o: test_random_uniform_dp.F90
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -Wno-compare-reals -Wno-conversion-extra -c $<
test_random_uniform_qp.x: test_random_uniform_qp.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_uniform_qp.o: test_random_uniform_qp.F90
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -Wno-compare-reals -Wno-conversion-extra -c $<
test_random_gaussian_sp.x: test_random_gaussian_sp.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_gaussian_sp.o: test_random_gaussian_sp.F90
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -Wno-compare-reals -Wno-conversion-extra -c $<
test_random_gaussian_dp.x: test_random_gaussian_dp.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_gaussian_dp.o: test_random_gaussian_dp.F90
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -Wno-compare-reals -Wno-conversion-extra -c $<
test_random_gaussian_qp.x: test_random_gaussian_qp.o
	$(FC) $(FCFLAGS) $(FCWFLAGS) -o $@ $^ $(MTLIB)
test_random_gaussian_qp.o: test_random_gaussian_qp.F90
	$(FC) $(FCFLAGS) $(FCWFLAGS) $(MTINC) -Wno-compare-reals -Wno-conversion-extra -c $<

.PHONY: test clean distclean

clean:
	rm -f *.o *.mod

distclean: clean
	rm -f *.x
