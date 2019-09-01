# MersenneTwister_Fortran

Description:
   Contains routines to use a non global mersenne twister for 
   uniformly and gaussian distributed pseudo random number generation. 
   The random number states are independen of each other, which allows for 
   arbritrarily many independend streams.

Compilation:
   Adjust the compiler and flags in the Makefile
   type "make" to build
   type "make test" to run tests
   

   There are three adjustment possibilities for your fortran standard.
   Add these to the FCFLAGS in the Makefile
   -D__HAS_QP      : Supports quadruple precision for real numbers 
                     (4 times size of single precision)
   -D__HAS_IQP     : Supports quadruple precision for integer numbers
                     (4 times size of single precision)
   -D__MAX_RANK=10 : Maximum rank of arrays in fortran 
                     (Maximum is  7 with F95)
                     (Maximum is 15 with F08)
                     (Implemented is a maximum of 15)

Usage:
   Import the module:
      USE MT_random
   Define a variable to hold the RNG state for uniform random numbers:
      TYPE(rng_uniform_type) :: rng_stateu
   For gaussian random numbers:
      TYPE(rng_gaussian_type) :: rng_stateg
   For integer random numbers over the complete integer range
      TYPE(rng_integer_type) :: rng_statei
   Initialize the uniform random nuber state:
      CALL init_rng(seed, lowerbound, upperbound, rng_stateu)
   Initialize the gaussian random nuber state:
      CALL init_rng(seed, meanvalue, stddev, rng_stateg)
   Initialize the integer random number state:
      CALL init_rng(seed, rng_statei)
   Get the next random number: 
      random_number might be a scalar or array pointer of any supported dimension
      CALL next_random(rng_state, random_number)
   Store the state in a string for writeout:
      CHARACTER(LEN=default_rng_state_length) :: state_string
      CALL get_random_state(rng_state, state_string)
   Restore a previously extracted state:
      CALL restore_random_state(state_string, rng_state)
   There is no need for deallocating anything

TODO:
   - Add polar method for gaussian random number generation, as it is faster
