!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Description:
!!    Contains routines to use a non global mersenne twister for 
!!    uniformly and gaussian distributed pseudo random number generation. 
!!    The random number states are independen of each other, which allows for 
!!    arbritrarily many independend streams.
!! Usage:
!!    Import the module:
!!       USE MT_random
!!    Define a variable to hold the RNG state for uniform random numbers:
!!       TYPE(rng_uniform_type) :: rng_stateu
!!    For gaussian random numbers:
!!       TYPE(rng_gaussian_type) :: rng_stateg
!!    For integer random numbers over the complete integer range
!!       TYPE(rng_integer_type) :: rng_statei
!!    Initialize the uniform random nuber state:
!!       CALL init_rng(seed, lowerbound, upperbound, rng_stateu)
!!    Initialize the gaussian random nuber state:
!!       CALL init_rng(seed, meanvalue, stddev, rng_stateg)
!!    Initialize the integer random number state:
!!       CALL init_rng(seed, rng_statei)
!!    Get the next random number:
!!    random_number might be a scalar or array pointer of any supported dimension
!!       CALL next_random(rng_state, random_number)
!!    Store the state in a string for writeout:
!!       CHARACTER(LEN=default_rng_state_length) :: state_string
!!       CALL get_random_state(rng_state, state_string)
!!    Restore a previously extracted state:
!!       CALL restore_random_state(state_string, rng_state)
!!    There is no need for deallocating anything
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE MT_random

   USE, INTRINSIC :: ISO_FORTRAN_ENV
#ifdef __HAS_QP
   USE kinds, ONLY : sp, dp, qp, &
#else
   USE kinds, ONLY : sp, dp, &
#endif
#ifdef __HAS_IQP
                     isp, idp, iqp
#else
                     isp, idp
#endif

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: default_rng_state_length
   PUBLIC :: rng_integer_type, &
             rng_gaussian_type, &
             rng_uniform_type
   PUBLIC :: init_rng
   PUBLIC :: next_random
   PUBLIC :: get_random_state
   PUBLIC :: restore_random_state

   REAL(KIND=dp), PARAMETER :: twopi = 6.2831853071795862_dp

   INTEGER(KIND=isp), PARAMETER :: default_rng_state_length = 8192

   INTEGER(KIND=isp), PARAMETER :: seed_def                 = 123456_idp
   INTEGER(KIND=idp), PARAMETER :: nn                       = 312_idp
   INTEGER(KIND=idp), PARAMETER :: mm                       = 156_idp
   INTEGER(KIND=idp), PARAMETER :: matrix_a                 = -5403634167711393303_idp
   INTEGER(KIND=idp), DIMENSION(0:1), PARAMETER :: mag01    = [0_idp, matrix_a]
   INTEGER(KIND=idp), PARAMETER :: um                       = -2147483648_idp ! most significant 33 bits
   INTEGER(KIND=idp), PARAMETER :: lm                       =  2147483647_idp ! least significant 31 bits

   REAL(KIND=dp), PARAMETER :: inv2pow53min1                = 1.0_dp/(2.0_dp**53 - 1.0_dp)
#ifdef __HAS_QP
   REAL(KIND=qp), PARAMETER :: inv2pow53min1qp              = 1.0_qp/(2.0_qp**53 - 1.0_qp)
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Type to hold the state of an integer random number generator
   !! Variables:
   !!    seed: The seed
   !!    mtstate: State of the mersenne twister
   !!    mtidx: index of the current random number of current state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   TYPE rng_integer_type
      INTEGER(KIND=isp) :: seed
      INTEGER(KIND=idp), DIMENSION(nn) :: mtstate
      INTEGER(KIND=idp) :: mtidx
   END TYPE rng_integer_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Type to hold the state of uniform random number generator
   !! Variables:
   !!    seed: The seed
   !!    lower: lower bound of random numbers
   !!    upper: upper bound of random numbers
   !!    mtstate: State of the mersenne twister
   !!    mtidx: index of the current random number of current state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   TYPE rng_uniform_type
      REAL(KIND=dp) :: lower
      REAL(KIND=dp) :: upper
      TYPE(rng_integer_type) :: rng_integer
   END TYPE rng_uniform_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Type to hold the state of gaussian random number generator
   !! Variables:
   !!    mean: mean value of the gaussian
   !!    stddev: standard deviation of the gaussian
   !!    z1: first random number of Box-Muller algorithm
   !!    z2: second random number of Box-Muller algorithm
   !!    generated: Boolean if Box-Muller needs to be reexecuted
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   TYPE rng_gaussian_type
      REAL(KIND=dp) :: mean
      REAL(KIND=dp) :: stddev
      REAL(KIND=dp) :: z1, z2
      LOGICAL :: generated
      TYPE(rng_uniform_type) :: rng_uniform
   END TYPE rng_gaussian_type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a uniform or gaussian random number generator
   !! Variables:
   !!    See individual implementation for description
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTERFACE init_rng
      MODULE PROCEDURE init_rng_integer_isp
      MODULE PROCEDURE init_rng_integer_idp
#ifdef __HAS_IQP
      MODULE PROCEDURE init_rng_integer_iqp
#endif
                       
      MODULE PROCEDURE init_rng_uniform_isp_sp
      MODULE PROCEDURE init_rng_uniform_isp_dp
      MODULE PROCEDURE init_rng_uniform_idp_sp
      MODULE PROCEDURE init_rng_uniform_idp_dp
#ifdef __HAS_QP
      MODULE PROCEDURE init_rng_uniform_isp_qp
      MODULE PROCEDURE init_rng_uniform_idp_qp
#endif
#ifdef __HAS_IQP
      MODULE PROCEDURE init_rng_uniform_iqp_sp
      MODULE PROCEDURE init_rng_uniform_iqp_dp
#endif

      MODULE PROCEDURE init_rng_gaussian_isp_sp
      MODULE PROCEDURE init_rng_gaussian_idp_sp
      MODULE PROCEDURE init_rng_gaussian_isp_dp
      MODULE PROCEDURE init_rng_gaussian_idp_dp
#ifdef __HAS_QP
      MODULE PROCEDURE init_rng_gaussian_isp_qp
      MODULE PROCEDURE init_rng_gaussian_idp_qp
#endif
#ifdef __HAS_IQP
      MODULE PROCEDURE init_rng_gaussian_iqp_sp
      MODULE PROCEDURE init_rng_gaussian_iqp_dp
#ifdef __HAS_QP
      MODULE PROCEDURE init_rng_gaussian_iqp_qp
#endif
#endif
   END INTERFACE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Mersenne twister to produce random numbers for several ranges and types
   !! Variables:
   !!    rng_integer: random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTERFACE mersennetwister
      MODULE PROCEDURE mersennetwister_sp
      MODULE PROCEDURE mersennetwister_dp
#ifdef __HAS_QP
      MODULE PROCEDURE mersennetwister_qp
#endif
      MODULE PROCEDURE mersennetwister_isp
      MODULE PROCEDURE mersennetwister_idp
#ifdef __HAS_IQP
      MODULE PROCEDURE mersennetwister_iqp
#endif
   END INTERFACE

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next random number
   !! Variables:
   !!    rng_uniform/rng_gaussian: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTERFACE next_random
      MODULE PROCEDURE next_random_integer_isp
#if __MAX_RANK > 0
      MODULE PROCEDURE next_random_integer_isp_1d
#if __MAX_RANK > 1
      MODULE PROCEDURE next_random_integer_isp_2d
#if __MAX_RANK > 2
      MODULE PROCEDURE next_random_integer_isp_3d
#if __MAX_RANK > 3
      MODULE PROCEDURE next_random_integer_isp_4d
#if __MAX_RANK > 4
      MODULE PROCEDURE next_random_integer_isp_5d
#if __MAX_RANK > 5
      MODULE PROCEDURE next_random_integer_isp_6d
#if __MAX_RANK > 6
      MODULE PROCEDURE next_random_integer_isp_7d
#if __MAX_RANK > 7
      MODULE PROCEDURE next_random_integer_isp_8d
#if __MAX_RANK > 8
      MODULE PROCEDURE next_random_integer_isp_9d
#if __MAX_RANK > 9
      MODULE PROCEDURE next_random_integer_isp_10d
#if __MAX_RANK > 10
      MODULE PROCEDURE next_random_integer_isp_11d
#if __MAX_RANK > 11
      MODULE PROCEDURE next_random_integer_isp_12d
#if __MAX_RANK > 12
      MODULE PROCEDURE next_random_integer_isp_13d
#if __MAX_RANK > 13
      MODULE PROCEDURE next_random_integer_isp_14d
#if __MAX_RANK > 14
      MODULE PROCEDURE next_random_integer_isp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

      MODULE PROCEDURE next_random_integer_idp
#if __MAX_RANK > 0
      MODULE PROCEDURE next_random_integer_idp_1d
#if __MAX_RANK > 1
      MODULE PROCEDURE next_random_integer_idp_2d
#if __MAX_RANK > 2
      MODULE PROCEDURE next_random_integer_idp_3d
#if __MAX_RANK > 3
      MODULE PROCEDURE next_random_integer_idp_4d
#if __MAX_RANK > 4
      MODULE PROCEDURE next_random_integer_idp_5d
#if __MAX_RANK > 5
      MODULE PROCEDURE next_random_integer_idp_6d
#if __MAX_RANK > 6
      MODULE PROCEDURE next_random_integer_idp_7d
#if __MAX_RANK > 7
      MODULE PROCEDURE next_random_integer_idp_8d
#if __MAX_RANK > 8
      MODULE PROCEDURE next_random_integer_idp_9d
#if __MAX_RANK > 9
      MODULE PROCEDURE next_random_integer_idp_10d
#if __MAX_RANK > 10
      MODULE PROCEDURE next_random_integer_idp_11d
#if __MAX_RANK > 11
      MODULE PROCEDURE next_random_integer_idp_12d
#if __MAX_RANK > 12
      MODULE PROCEDURE next_random_integer_idp_13d
#if __MAX_RANK > 13
      MODULE PROCEDURE next_random_integer_idp_14d
#if __MAX_RANK > 14
      MODULE PROCEDURE next_random_integer_idp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#ifdef __HAS_IQP
      MODULE PROCEDURE next_random_integer_iqp
#if __MAX_RANK > 0
      MODULE PROCEDURE next_random_integer_iqp_1d
#if __MAX_RANK > 1
      MODULE PROCEDURE next_random_integer_iqp_2d
#if __MAX_RANK > 2
      MODULE PROCEDURE next_random_integer_iqp_3d
#if __MAX_RANK > 3
      MODULE PROCEDURE next_random_integer_iqp_4d
#if __MAX_RANK > 4
      MODULE PROCEDURE next_random_integer_iqp_5d
#if __MAX_RANK > 5
      MODULE PROCEDURE next_random_integer_iqp_6d
#if __MAX_RANK > 6
      MODULE PROCEDURE next_random_integer_iqp_7d
#if __MAX_RANK > 7
      MODULE PROCEDURE next_random_integer_iqp_8d
#if __MAX_RANK > 8
      MODULE PROCEDURE next_random_integer_iqp_9d
#if __MAX_RANK > 9
      MODULE PROCEDURE next_random_integer_iqp_10d
#if __MAX_RANK > 10
      MODULE PROCEDURE next_random_integer_iqp_11d
#if __MAX_RANK > 11
      MODULE PROCEDURE next_random_integer_iqp_12d
#if __MAX_RANK > 12
      MODULE PROCEDURE next_random_integer_iqp_13d
#if __MAX_RANK > 13
      MODULE PROCEDURE next_random_integer_iqp_14d
#if __MAX_RANK > 14
      MODULE PROCEDURE next_random_integer_iqp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

      MODULE PROCEDURE next_random_uniform_sp
#if __MAX_RANK > 0
      MODULE PROCEDURE next_random_uniform_sp_1d
#if __MAX_RANK > 1
      MODULE PROCEDURE next_random_uniform_sp_2d
#if __MAX_RANK > 2
      MODULE PROCEDURE next_random_uniform_sp_3d
#if __MAX_RANK > 3
      MODULE PROCEDURE next_random_uniform_sp_4d
#if __MAX_RANK > 4
      MODULE PROCEDURE next_random_uniform_sp_5d
#if __MAX_RANK > 5
      MODULE PROCEDURE next_random_uniform_sp_6d
#if __MAX_RANK > 6
      MODULE PROCEDURE next_random_uniform_sp_7d
#if __MAX_RANK > 7
      MODULE PROCEDURE next_random_uniform_sp_8d
#if __MAX_RANK > 8
      MODULE PROCEDURE next_random_uniform_sp_9d
#if __MAX_RANK > 9
      MODULE PROCEDURE next_random_uniform_sp_10d
#if __MAX_RANK > 10
      MODULE PROCEDURE next_random_uniform_sp_11d
#if __MAX_RANK > 11
      MODULE PROCEDURE next_random_uniform_sp_12d
#if __MAX_RANK > 12
      MODULE PROCEDURE next_random_uniform_sp_13d
#if __MAX_RANK > 13
      MODULE PROCEDURE next_random_uniform_sp_14d
#if __MAX_RANK > 14
      MODULE PROCEDURE next_random_uniform_sp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

      MODULE PROCEDURE next_random_uniform_dp
#if __MAX_RANK > 0
      MODULE PROCEDURE next_random_uniform_dp_1d
#if __MAX_RANK > 1
      MODULE PROCEDURE next_random_uniform_dp_2d
#if __MAX_RANK > 2
      MODULE PROCEDURE next_random_uniform_dp_3d
#if __MAX_RANK > 3
      MODULE PROCEDURE next_random_uniform_dp_4d
#if __MAX_RANK > 4
      MODULE PROCEDURE next_random_uniform_dp_5d
#if __MAX_RANK > 5
      MODULE PROCEDURE next_random_uniform_dp_6d
#if __MAX_RANK > 6
      MODULE PROCEDURE next_random_uniform_dp_7d
#if __MAX_RANK > 7
      MODULE PROCEDURE next_random_uniform_dp_8d
#if __MAX_RANK > 8
      MODULE PROCEDURE next_random_uniform_dp_9d
#if __MAX_RANK > 9
      MODULE PROCEDURE next_random_uniform_dp_10d
#if __MAX_RANK > 10
      MODULE PROCEDURE next_random_uniform_dp_11d
#if __MAX_RANK > 11
      MODULE PROCEDURE next_random_uniform_dp_12d
#if __MAX_RANK > 12
      MODULE PROCEDURE next_random_uniform_dp_13d
#if __MAX_RANK > 13
      MODULE PROCEDURE next_random_uniform_dp_14d
#if __MAX_RANK > 14
      MODULE PROCEDURE next_random_uniform_dp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#ifdef __HAS_QP
      MODULE PROCEDURE next_random_uniform_qp
#if __MAX_RANK > 0
      MODULE PROCEDURE next_random_uniform_qp_1d
#if __MAX_RANK > 1
      MODULE PROCEDURE next_random_uniform_qp_2d
#if __MAX_RANK > 2
      MODULE PROCEDURE next_random_uniform_qp_3d
#if __MAX_RANK > 3
      MODULE PROCEDURE next_random_uniform_qp_4d
#if __MAX_RANK > 4
      MODULE PROCEDURE next_random_uniform_qp_5d
#if __MAX_RANK > 5
      MODULE PROCEDURE next_random_uniform_qp_6d
#if __MAX_RANK > 6
      MODULE PROCEDURE next_random_uniform_qp_7d
#if __MAX_RANK > 7
      MODULE PROCEDURE next_random_uniform_qp_8d
#if __MAX_RANK > 8
      MODULE PROCEDURE next_random_uniform_qp_9d
#if __MAX_RANK > 9
      MODULE PROCEDURE next_random_uniform_qp_10d
#if __MAX_RANK > 10
      MODULE PROCEDURE next_random_uniform_qp_11d
#if __MAX_RANK > 11
      MODULE PROCEDURE next_random_uniform_qp_12d
#if __MAX_RANK > 12
      MODULE PROCEDURE next_random_uniform_qp_13d
#if __MAX_RANK > 13
      MODULE PROCEDURE next_random_uniform_qp_14d
#if __MAX_RANK > 14
      MODULE PROCEDURE next_random_uniform_qp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

      MODULE PROCEDURE next_random_gaussian_sp
#if __MAX_RANK > 0
      MODULE PROCEDURE next_random_gaussian_sp_1d
#if __MAX_RANK > 1
      MODULE PROCEDURE next_random_gaussian_sp_2d
#if __MAX_RANK > 2
      MODULE PROCEDURE next_random_gaussian_sp_3d
#if __MAX_RANK > 3
      MODULE PROCEDURE next_random_gaussian_sp_4d
#if __MAX_RANK > 4
      MODULE PROCEDURE next_random_gaussian_sp_5d
#if __MAX_RANK > 5
      MODULE PROCEDURE next_random_gaussian_sp_6d
#if __MAX_RANK > 6
      MODULE PROCEDURE next_random_gaussian_sp_7d
#if __MAX_RANK > 7
      MODULE PROCEDURE next_random_gaussian_sp_8d
#if __MAX_RANK > 8
      MODULE PROCEDURE next_random_gaussian_sp_9d
#if __MAX_RANK > 9
      MODULE PROCEDURE next_random_gaussian_sp_10d
#if __MAX_RANK > 10
      MODULE PROCEDURE next_random_gaussian_sp_11d
#if __MAX_RANK > 11
      MODULE PROCEDURE next_random_gaussian_sp_12d
#if __MAX_RANK > 12
      MODULE PROCEDURE next_random_gaussian_sp_13d
#if __MAX_RANK > 13
      MODULE PROCEDURE next_random_gaussian_sp_14d
#if __MAX_RANK > 14
      MODULE PROCEDURE next_random_gaussian_sp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

      MODULE PROCEDURE next_random_gaussian_dp
#if __MAX_RANK > 0
      MODULE PROCEDURE next_random_gaussian_dp_1d
#if __MAX_RANK > 1
      MODULE PROCEDURE next_random_gaussian_dp_2d
#if __MAX_RANK > 2
      MODULE PROCEDURE next_random_gaussian_dp_3d
#if __MAX_RANK > 3
      MODULE PROCEDURE next_random_gaussian_dp_4d
#if __MAX_RANK > 4
      MODULE PROCEDURE next_random_gaussian_dp_5d
#if __MAX_RANK > 5
      MODULE PROCEDURE next_random_gaussian_dp_6d
#if __MAX_RANK > 6
      MODULE PROCEDURE next_random_gaussian_dp_7d
#if __MAX_RANK > 7
      MODULE PROCEDURE next_random_gaussian_dp_8d
#if __MAX_RANK > 8
      MODULE PROCEDURE next_random_gaussian_dp_9d
#if __MAX_RANK > 9
      MODULE PROCEDURE next_random_gaussian_dp_10d
#if __MAX_RANK > 10
      MODULE PROCEDURE next_random_gaussian_dp_11d
#if __MAX_RANK > 11
      MODULE PROCEDURE next_random_gaussian_dp_12d
#if __MAX_RANK > 12
      MODULE PROCEDURE next_random_gaussian_dp_13d
#if __MAX_RANK > 13
      MODULE PROCEDURE next_random_gaussian_dp_14d
#if __MAX_RANK > 14
      MODULE PROCEDURE next_random_gaussian_dp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#ifdef __HAS_QP
      MODULE PROCEDURE next_random_gaussian_qp
#if __MAX_RANK > 0
      MODULE PROCEDURE next_random_gaussian_qp_1d
#if __MAX_RANK > 1
      MODULE PROCEDURE next_random_gaussian_qp_2d
#if __MAX_RANK > 2
      MODULE PROCEDURE next_random_gaussian_qp_3d
#if __MAX_RANK > 3
      MODULE PROCEDURE next_random_gaussian_qp_4d
#if __MAX_RANK > 4
      MODULE PROCEDURE next_random_gaussian_qp_5d
#if __MAX_RANK > 5
      MODULE PROCEDURE next_random_gaussian_qp_6d
#if __MAX_RANK > 6
      MODULE PROCEDURE next_random_gaussian_qp_7d
#if __MAX_RANK > 7
      MODULE PROCEDURE next_random_gaussian_qp_8d
#if __MAX_RANK > 8
      MODULE PROCEDURE next_random_gaussian_qp_9d
#if __MAX_RANK > 9
      MODULE PROCEDURE next_random_gaussian_qp_10d
#if __MAX_RANK > 10
      MODULE PROCEDURE next_random_gaussian_qp_11d
#if __MAX_RANK > 11
      MODULE PROCEDURE next_random_gaussian_qp_12d
#if __MAX_RANK > 12
      MODULE PROCEDURE next_random_gaussian_qp_13d
#if __MAX_RANK > 13
      MODULE PROCEDURE next_random_gaussian_qp_14d
#if __MAX_RANK > 14
      MODULE PROCEDURE next_random_gaussian_qp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

   END INTERFACE next_random

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Save rng state to a string of length default_rng_state_length
   !! Variables:
   !!    rng_uniform/rng_gaussian: Random number generator state
   !!    statestr: that contains all the state information
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTERFACE get_random_state
      MODULE PROCEDURE get_random_state_integer, &
                       get_random_state_uniform, &
                       get_random_state_gaussian
   END INTERFACE get_random_state

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Restore rng state from a string of length default_rng_state_length
   !! Variables:
   !!    statestr: that contains all the state information
   !!    rng_uniform/rng_gaussian: Random number generator state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTERFACE restore_random_state
      MODULE PROCEDURE restore_random_state_integer, &
                       restore_random_state_uniform, &
                       restore_random_state_gaussian
   END INTERFACE restore_random_state

CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a integer random number generator
   !! Variables:
   !!    seed: The seed
   !!    rng_integer: integer random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_integer_isp(seed, rng_integer)

      IMPLICIT NONE
   
      INTEGER(KIND=isp), INTENT(IN) :: seed
      TYPE(rng_integer_type), INTENT(OUT) :: rng_integer

      INTEGER(KIND=idp) :: i

      rng_integer%seed = seed
      rng_integer%mtstate(1) = INT(seed, idp)
      DO i = 1, nn-1
         rng_integer%mtstate(i+1) = 6364136223846793005_idp * &
                                    IEOR(rng_integer%mtstate(i), ISHFT(rng_integer%mtstate(i), -62)) + i
      END DO
      rng_integer%mtidx = nn+1

      RETURN

   END SUBROUTINE init_rng_integer_isp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a integer random number generator
   !! Variables:
   !!    seed: The seed
   !!    rng_integer: integer random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_integer_idp(seed, rng_integer)

      IMPLICIT NONE

      INTEGER(KIND=idp), INTENT(IN) :: seed
      TYPE(rng_integer_type), INTENT(OUT) :: rng_integer

      CALL init_rng_integer_isp( &
         INT( &
            MOD(seed, &
               INT(HUGE(1_sp),dp) &
            ) &
         ), rng_integer)

      RETURN

   END SUBROUTINE init_rng_integer_idp

#ifdef __HAS_IQP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a integer random number generator
   !! Variables:
   !!    seed: The seed
   !!    rng_integer: integer random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_integer_iqp(seed, rng_integer)

      IMPLICIT NONE

      INTEGER(KIND=iqp), INTENT(IN) :: seed
      TYPE(rng_integer_type), INTENT(OUT) :: rng_integer

      CALL init_rng_integer_isp( &
         INT( &
            MOD(seed, &
               INT(HUGE(1_sp),iqp) &
            ) &
         ), rng_integer)

      RETURN

   END SUBROUTINE init_rng_integer_iqp
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a uniform random number generator
   !! Variables:
   !!    seed: The seed
   !!    lower: lower bound of random numbers
   !!    upper: upper bound of random numbers
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_uniform_isp_sp(seed, lower, upper, rng_uniform)

      IMPLICIT NONE
   
      INTEGER(KIND=isp), INTENT(IN) :: seed
      REAL(KIND=sp), INTENT(IN) :: lower
      REAL(KIND=sp), INTENT(IN) :: upper
      TYPE(rng_uniform_type), INTENT(OUT) :: rng_uniform

      rng_uniform%lower = REAL(lower, dp)
      rng_uniform%upper = REAL(upper, dp)

      CALL init_rng_integer_isp(seed, rng_uniform%rng_integer)

      RETURN

   END SUBROUTINE init_rng_uniform_isp_sp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a uniform random number generator
   !! Variables:
   !!    seed: The seed
   !!    lower: lower bound of random numbers
   !!    upper: upper bound of random numbers
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_uniform_isp_dp(seed, lower, upper, rng_uniform)

      IMPLICIT NONE
   
      INTEGER(KIND=isp), INTENT(IN) :: seed
      REAL(KIND=dp), INTENT(IN) :: lower
      REAL(KIND=dp), INTENT(IN) :: upper
      TYPE(rng_uniform_type), INTENT(OUT) :: rng_uniform

      rng_uniform%lower = lower
      rng_uniform%upper = upper

      CALL init_rng_integer_isp(seed, rng_uniform%rng_integer)

      RETURN

   END SUBROUTINE init_rng_uniform_isp_dp

#ifdef __HAS_QP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a uniform random number generator
   !! Variables:
   !!    seed: The seed
   !!    lower: lower bound of random numbers
   !!    upper: upper bound of random numbers
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_uniform_isp_qp(seed, lower, upper, rng_uniform)

      IMPLICIT NONE

      INTEGER(KIND=isp), INTENT(IN) :: seed
      REAL(KIND=qp), INTENT(IN) :: lower
      REAL(KIND=qp), INTENT(IN) :: upper
      TYPE(rng_uniform_type), INTENT(OUT) :: rng_uniform

      rng_uniform%lower = REAL(lower, dp)
      rng_uniform%upper = REAL(upper, dp)

      CALL init_rng_integer_isp(seed, rng_uniform%rng_integer)

      RETURN

   END SUBROUTINE init_rng_uniform_isp_qp
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a uniform random number generator
   !! Variables:
   !!    seed: The seed
   !!    lower: lower bound of random numbers
   !!    upper: upper bound of random numbers
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_uniform_idp_sp(seed, lower, upper, rng_uniform)

      IMPLICIT NONE
   
      INTEGER(KIND=idp), INTENT(IN) :: seed
      REAL(KIND=sp), INTENT(IN) :: lower
      REAL(KIND=sp), INTENT(IN) :: upper
      TYPE(rng_uniform_type), INTENT(OUT) :: rng_uniform

      rng_uniform%lower = REAL(lower, dp)
      rng_uniform%upper = REAL(upper, dp)

      CALL init_rng_integer_idp(seed, rng_uniform%rng_integer)

      RETURN

   END SUBROUTINE init_rng_uniform_idp_sp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a uniform random number generator
   !! Variables:
   !!    seed: The seed
   !!    lower: lower bound of random numbers
   !!    upper: upper bound of random numbers
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_uniform_idp_dp(seed, lower, upper, rng_uniform)

      IMPLICIT NONE
   
      INTEGER(KIND=idp), INTENT(IN) :: seed
      REAL(KIND=dp), INTENT(IN) :: lower
      REAL(KIND=dp), INTENT(IN) :: upper
      TYPE(rng_uniform_type), INTENT(OUT) :: rng_uniform

      rng_uniform%lower = lower
      rng_uniform%upper = upper

      CALL init_rng_integer_idp(seed, rng_uniform%rng_integer)

      RETURN

   END SUBROUTINE init_rng_uniform_idp_dp

#ifdef __HAS_QP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a uniform random number generator
   !! Variables:
   !!    seed: The seed
   !!    lower: lower bound of random numbers
   !!    upper: upper bound of random numbers
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_uniform_idp_qp(seed, lower, upper, rng_uniform)

      IMPLICIT NONE

      INTEGER(KIND=idp), INTENT(IN) :: seed
      REAL(KIND=qp), INTENT(IN) :: lower
      REAL(KIND=qp), INTENT(IN) :: upper
      TYPE(rng_uniform_type), INTENT(OUT) :: rng_uniform

      rng_uniform%lower = REAL(lower, dp)
      rng_uniform%upper = REAL(upper, dp)

      CALL init_rng_integer_idp(seed, rng_uniform%rng_integer)

      RETURN

   END SUBROUTINE init_rng_uniform_idp_qp
#endif

#ifdef __HAS_IQP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a uniform random number generator
   !! Variables:
   !!    seed: The seed
   !!    lower: lower bound of random numbers
   !!    upper: upper bound of random numbers
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_uniform_iqp_sp(seed, lower, upper, rng_uniform)

      IMPLICIT NONE
   
      INTEGER(KIND=iqp), INTENT(IN) :: seed
      REAL(KIND=sp), INTENT(IN) :: lower
      REAL(KIND=sp), INTENT(IN) :: upper
      TYPE(rng_uniform_type), INTENT(OUT) :: rng_uniform

      rng_uniform%lower = REAL(lower, dp)
      rng_uniform%upper = REAL(upper, dp)

      CALL init_rng_integer_iqp(seed, rng_uniform%rng_integer)

      RETURN

   END SUBROUTINE init_rng_uniform_iqp_sp
#endif

#ifdef __HAS_IQP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a uniform random number generator
   !! Variables:
   !!    seed: The seed
   !!    lower: lower bound of random numbers
   !!    upper: upper bound of random numbers
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_uniform_iqp_dp(seed, lower, upper, rng_uniform)

      IMPLICIT NONE
   
      INTEGER(KIND=iqp), INTENT(IN) :: seed
      REAL(KIND=dp), INTENT(IN) :: lower
      REAL(KIND=dp), INTENT(IN) :: upper
      TYPE(rng_uniform_type), INTENT(OUT) :: rng_uniform

      rng_uniform%lower = lower
      rng_uniform%upper = upper

      CALL init_rng_integer_iqp(seed, rng_uniform%rng_integer)

      RETURN

   END SUBROUTINE init_rng_uniform_iqp_dp
#endif

#ifdef __HAS_QP
#ifdef __HAS_IQP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a uniform random number generator
   !! Variables:
   !!    seed: The seed
   !!    lower: lower bound of random numbers
   !!    upper: upper bound of random numbers
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_uniform_iqp_qp(seed, lower, upper, rng_uniform)

      IMPLICIT NONE

      INTEGER(KIND=iqp), INTENT(IN) :: seed
      REAL(KIND=qp), INTENT(IN) :: lower
      REAL(KIND=qp), INTENT(IN) :: upper
      TYPE(rng_uniform_type), INTENT(OUT) :: rng_uniform

      rng_uniform%lower = REAL(lower, dp)
      rng_uniform%upper = REAL(upper, dp)

      CALL init_rng_integer_iqp(seed, rng_uniform%rng_integer)

      RETURN

   END SUBROUTINE init_rng_uniform_iqp_qp
#endif
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a gaussian random number generator
   !! Variables:
   !!    seed: The seed
   !!    mean: Mean value of the gaussian
   !!    stddev: standard deviation of the gaussian
   !!    rng_gaussian: gaussian random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_gaussian_isp_sp(seed, mean, stddev, rng_gaussian)

      IMPLICIT NONE

      INTEGER(KIND=isp), INTENT(IN) :: seed
      REAL(KIND=sp), INTENT(IN) :: mean
      REAL(KIND=sp), INTENT(IN) :: stddev
      TYPE(rng_gaussian_type), INTENT(OUT) :: rng_gaussian

      rng_gaussian%mean = REAL(mean, dp)
      rng_gaussian%stddev = REAL(stddev, dp)
      rng_gaussian%z1 = 0.0_dp
      rng_gaussian%z2 = 0.0_dp
      rng_gaussian%generated = .FALSE.
      
      CALL init_rng_uniform_isp_sp(seed, 0.0_sp, 1.0_sp, rng_gaussian%rng_uniform)

      RETURN

   END SUBROUTINE init_rng_gaussian_isp_sp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a gaussian random number generator
   !! Variables:
   !!    seed: The seed
   !!    mean: Mean value of the gaussian
   !!    stddev: standard deviation of the gaussian
   !!    rng_gaussian: gaussian random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_gaussian_isp_dp(seed, mean, stddev, rng_gaussian)

      IMPLICIT NONE

      INTEGER(KIND=isp), INTENT(IN) :: seed
      REAL(KIND=dp), INTENT(IN) :: mean
      REAL(KIND=dp), INTENT(IN) :: stddev
      TYPE(rng_gaussian_type), INTENT(OUT) :: rng_gaussian

      rng_gaussian%mean = mean
      rng_gaussian%stddev = stddev
      rng_gaussian%z1 = 0.0_dp
      rng_gaussian%z2 = 0.0_dp
      rng_gaussian%generated = .FALSE.
      
      CALL init_rng_uniform_isp_dp(seed, 0.0_dp, 1.0_dp, rng_gaussian%rng_uniform)

      RETURN

   END SUBROUTINE init_rng_gaussian_isp_dp

#ifdef __HAS_QP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a gaussian random number generator
   !! Variables:
   !!    seed: The seed
   !!    mean: Mean value of the gaussian
   !!    stddev: standard deviation of the gaussian
   !!    rng_gaussian: gaussian random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_gaussian_isp_qp(seed, mean, stddev, rng_gaussian)

      IMPLICIT NONE

      INTEGER(KIND=isp), INTENT(IN) :: seed
      REAL(KIND=qp), INTENT(IN) :: mean
      REAL(KIND=qp), INTENT(IN) :: stddev
      TYPE(rng_gaussian_type), INTENT(OUT) :: rng_gaussian

      rng_gaussian%mean = REAL(mean, dp)
      rng_gaussian%stddev = REAL(stddev, dp)
      rng_gaussian%z1 = 0.0_dp
      rng_gaussian%z2 = 0.0_dp
      rng_gaussian%generated = .FALSE.
      
      CALL init_rng_uniform_isp_qp(seed, 0.0_qp, 1.0_qp, rng_gaussian%rng_uniform)

      RETURN

   END SUBROUTINE init_rng_gaussian_isp_qp
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a gaussian random number generator
   !! Variables:
   !!    seed: The seed
   !!    mean: Mean value of the gaussian
   !!    stddev: standard deviation of the gaussian
   !!    rng_gaussian: gaussian random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_gaussian_idp_sp(seed, mean, stddev, rng_gaussian)

      IMPLICIT NONE

      INTEGER(KIND=idp), INTENT(IN) :: seed
      REAL(KIND=sp), INTENT(IN) :: mean
      REAL(KIND=sp), INTENT(IN) :: stddev
      TYPE(rng_gaussian_type), INTENT(OUT) :: rng_gaussian

      rng_gaussian%mean = REAL(mean, dp)
      rng_gaussian%stddev = REAL(stddev, dp)
      rng_gaussian%z1 = 0.0_dp
      rng_gaussian%z2 = 0.0_dp
      rng_gaussian%generated = .FALSE.
      
      CALL init_rng_uniform_idp_sp(seed, 0.0_sp, 1.0_sp, rng_gaussian%rng_uniform)

      RETURN

   END SUBROUTINE init_rng_gaussian_idp_sp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a gaussian random number generator
   !! Variables:
   !!    seed: The seed
   !!    mean: Mean value of the gaussian
   !!    stddev: standard deviation of the gaussian
   !!    rng_gaussian: gaussian random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_gaussian_idp_dp(seed, mean, stddev, rng_gaussian)

      IMPLICIT NONE

      INTEGER(KIND=idp), INTENT(IN) :: seed
      REAL(KIND=dp), INTENT(IN) :: mean
      REAL(KIND=dp), INTENT(IN) :: stddev
      TYPE(rng_gaussian_type), INTENT(OUT) :: rng_gaussian

      rng_gaussian%mean = mean
      rng_gaussian%stddev = stddev
      rng_gaussian%z1 = 0.0_dp
      rng_gaussian%z2 = 0.0_dp
      rng_gaussian%generated = .FALSE.
      
      CALL init_rng_uniform_idp_dp(seed, 0.0_dp, 1.0_dp, rng_gaussian%rng_uniform)

      RETURN

   END SUBROUTINE init_rng_gaussian_idp_dp

#ifdef __HAS_QP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a gaussian random number generator
   !! Variables:
   !!    seed: The seed
   !!    mean: Mean value of the gaussian
   !!    stddev: standard deviation of the gaussian
   !!    rng_gaussian: gaussian random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_gaussian_idp_qp(seed, mean, stddev, rng_gaussian)

      IMPLICIT NONE

      INTEGER(KIND=idp), INTENT(IN) :: seed
      REAL(KIND=qp), INTENT(IN) :: mean
      REAL(KIND=qp), INTENT(IN) :: stddev
      TYPE(rng_gaussian_type), INTENT(OUT) :: rng_gaussian

      rng_gaussian%mean = REAL(mean, dp)
      rng_gaussian%stddev = REAL(stddev, dp)
      rng_gaussian%z1 = 0.0_dp
      rng_gaussian%z2 = 0.0_dp
      rng_gaussian%generated = .FALSE.
      
      CALL init_rng_uniform_idp_qp(seed, 0.0_qp, 1.0_qp, rng_gaussian%rng_uniform)

      RETURN

   END SUBROUTINE init_rng_gaussian_idp_qp
#endif

#ifdef __HAS_IQP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a gaussian random number generator
   !! Variables:
   !!    seed: The seed
   !!    mean: Mean value of the gaussian
   !!    stddev: standard deviation of the gaussian
   !!    rng_gaussian: gaussian random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_gaussian_iqp_sp(seed, mean, stddev, rng_gaussian)

      IMPLICIT NONE

      INTEGER(KIND=iqp), INTENT(IN) :: seed
      REAL(KIND=sp), INTENT(IN) :: mean
      REAL(KIND=sp), INTENT(IN) :: stddev
      TYPE(rng_gaussian_type), INTENT(OUT) :: rng_gaussian

      rng_gaussian%mean = REAL(mean, dp)
      rng_gaussian%stddev = REAL(stddev, dp)
      rng_gaussian%z1 = 0.0_dp
      rng_gaussian%z2 = 0.0_dp
      rng_gaussian%generated = .FALSE.
      
      CALL init_rng_uniform_iqp_sp(seed, 0.0_sp, 1.0_sp, rng_gaussian%rng_uniform)

      RETURN

   END SUBROUTINE init_rng_gaussian_iqp_sp
#endif

#ifdef __HAS_IQP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a gaussian random number generator
   !! Variables:
   !!    seed: The seed
   !!    mean: Mean value of the gaussian
   !!    stddev: standard deviation of the gaussian
   !!    rng_gaussian: gaussian random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_gaussian_iqp_dp(seed, mean, stddev, rng_gaussian)

      IMPLICIT NONE

      INTEGER(KIND=iqp), INTENT(IN) :: seed
      REAL(KIND=dp), INTENT(IN) :: mean
      REAL(KIND=dp), INTENT(IN) :: stddev
      TYPE(rng_gaussian_type), INTENT(OUT) :: rng_gaussian

      rng_gaussian%mean = mean
      rng_gaussian%stddev = stddev
      rng_gaussian%z1 = 0.0_dp
      rng_gaussian%z2 = 0.0_dp
      rng_gaussian%generated = .FALSE.
      
      CALL init_rng_uniform_iqp_dp(seed, 0.0_dp, 1.0_dp, rng_gaussian%rng_uniform)

      RETURN

   END SUBROUTINE init_rng_gaussian_iqp_dp
#endif

#ifdef __HAS_QP
#ifdef __HAS_IQP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Initializes a gaussian random number generator
   !! Variables:
   !!    seed: The seed
   !!    mean: Mean value of the gaussian
   !!    stddev: standard deviation of the gaussian
   !!    rng_gaussian: gaussian random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE init_rng_gaussian_iqp_qp(seed, mean, stddev, rng_gaussian)

      IMPLICIT NONE

      INTEGER(KIND=iqp), INTENT(IN) :: seed
      REAL(KIND=qp), INTENT(IN) :: mean
      REAL(KIND=qp), INTENT(IN) :: stddev
      TYPE(rng_gaussian_type), INTENT(OUT) :: rng_gaussian

      rng_gaussian%mean = REAL(mean, dp)
      rng_gaussian%stddev = REAL(stddev, dp)
      rng_gaussian%z1 = 0.0_dp
      rng_gaussian%z2 = 0.0_dp
      rng_gaussian%generated = .FALSE.
      
      CALL init_rng_uniform_iqp_qp(seed, 0.0_qp, 1.0_qp, rng_gaussian%rng_uniform)

      RETURN

   END SUBROUTINE init_rng_gaussian_iqp_qp
#endif
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Mersenne twister to produce random number in [-2^31, 2^31-1]
   !! Variables:
   !!    rng_integer: integer random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE mersennetwister_isp(rng_integer, rng)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(OUT) :: rng

      INTEGER(KIND=idp) :: tmprng

      CALL mersennetwister_idp(rng_integer, tmprng)      
      tmprng = RSHIFT(tmprng, 32)
      rng = INT(tmprng, isp)

      RETURN

   END SUBROUTINE mersennetwister_isp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Mersenne twister to produce random number in [-2^63, 2^63-1]
   !! Variables:
   !!    rng_integer: integer random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE mersennetwister_idp(rng_integer, rng)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(OUT) :: rng

      INTEGER(KIND=idp) :: i

      IF (rng_integer%mtidx >= nn) THEN ! generate nn words at one time


         DO i = 1, nn-mm
            rng = IOR(IAND(rng_integer%mtstate(i),um), IAND(rng_integer%mtstate(i+1), lm))
            rng_integer%mtstate(i) = IEOR(IEOR(rng_integer%mtstate(i+mm), ISHFT(rng, -1)), mag01(IAND(rng, 1_idp)))
         END DO

         DO i = nn-mm+1, nn-1
            rng = IOR(IAND(rng_integer%mtstate(i), um), IAND(rng_integer%mtstate(i+1), lm))
            rng_integer%mtstate(i) = IEOR(IEOR(rng_integer%mtstate(i+mm-nn), ISHFT(rng, -1)), mag01(IAND(rng, 1_idp)))
         END DO

         rng = IOR(IAND(rng_integer%mtstate(nn), um), IAND(rng_integer%mtstate(1), lm))
         rng_integer%mtstate(nn) = IEOR(IEOR(rng_integer%mtstate(mm), ISHFT(rng, -1)), mag01(IAND(rng, 1_idp)))

         rng_integer%mtidx = 0

      END IF

      rng_integer%mtidx = rng_integer%mtidx + 1
      rng = rng_integer%mtstate(rng_integer%mtidx)

      rng = IEOR(rng, IAND(ISHFT(rng,-29), 6148914691236517205_idp))
      rng = IEOR(rng, IAND(ISHFT(rng, 17), 8202884508482404352_idp))
      rng = IEOR(rng, IAND(ISHFT(rng, 37),   -2270628950310912_idp))
      rng = IEOR(rng, ISHFT(rng, -43))

      RETURN

   END SUBROUTINE mersennetwister_idp

#ifdef __HAS_IQP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Mersenne twister to produce random number in [-2^61, 2^61-1]
   !! Variables:
   !!    rng_integer: integer random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE mersennetwister_iqp(rng_integer, rng)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(OUT) :: rng

      INTEGER(KIND=idp) :: tmprng

      CALL mersennetwister_idp(rng_integer, tmprng)      
      rng = INT(tmprng, iqp)

      RETURN

   END SUBROUTINE mersennetwister_iqp
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Mersenne twister to produce random number in [0:1]
   !! Variables:
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE mersennetwister_sp(rng_uniform, rng)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(OUT) :: rng

      INTEGER(KIND=idp) :: x
      
      CALL mersennetwister_idp(rng_uniform%rng_integer, x)
      rng = REAL(REAL(ISHFT(x, -11), KIND=dp) * inv2pow53min1, KIND=sp)

      RETURN

   END SUBROUTINE mersennetwister_sp

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Mersenne twister to produce random number in [0:1]
   !! Variables:
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE mersennetwister_dp(rng_uniform, rng)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(OUT) :: rng

      INTEGER(KIND=idp) :: x

      CALL mersennetwister_idp(rng_uniform%rng_integer, x)
      rng = REAL(ISHFT(x, -11), KIND=dp) * inv2pow53min1

      RETURN

   END SUBROUTINE mersennetwister_dp

#ifdef __HAS_QP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Mersenne twister to produce random number in [0:1]
   !! Variables:
   !!    rng_uniform: uniform random number state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE mersennetwister_qp(rng_uniform, rng)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(OUT) :: rng

      INTEGER(KIND=idp) :: x

      CALL mersennetwister_idp(rng_uniform%rng_integer, x)
      rng = REAL(ISHFT(x, -11), KIND=qp) * inv2pow53min1qp

      RETURN

   END SUBROUTINE mersennetwister_qp
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(OUT) :: randnum
   
      CALL mersennetwister_isp(rng_integer, randnum) 

      RETURN

   END SUBROUTINE next_random_integer_isp

#if __MAX_RANK > 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 1D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_1d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:), POINTER :: randnum
   
      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(1) :: startidx, endidx
      INTEGER :: i1

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)
      
      DO i1 = startidx(1), endidx(1)
         CALL next_random_integer_isp(rng_integer, tmprng)
         randnum(i1) = tmprng
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_1d

#if __MAX_RANK > 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 2D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_2d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:), POINTER :: randnum
   
      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(2) :: startidx, endidx
      INTEGER :: i1, i2

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i2 = startidx(2), endidx(2)
         DO i1 = startidx(1), endidx(1)
            CALL next_random_integer_isp(rng_integer, tmprng)
            randnum(i1,i2) = tmprng
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_2d

#if __MAX_RANK > 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 3D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_3d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: randnum
   
      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(3) :: startidx, endidx
      INTEGER :: i1, i2, i3

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i3 = startidx(3), endidx(3)
         DO i2 = startidx(2), endidx(2)
            DO i1 = startidx(1), endidx(1)
               CALL next_random_integer_isp(rng_integer, tmprng)
               randnum(i1,i2,i3) = tmprng
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_3d

#if __MAX_RANK > 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 4D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_4d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(4) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i4 = startidx(4), endidx(4)
         DO i3 = startidx(3), endidx(3)
            DO i2 = startidx(2), endidx(2)
               DO i1 = startidx(1), endidx(1)
                  CALL next_random_integer_isp(rng_integer, tmprng)
                  randnum(i1,i2,i3,i4) = tmprng
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_4d

#if __MAX_RANK > 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 5D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_5d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(5) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i5 = startidx(5), endidx(5)
         DO i4 = startidx(4), endidx(4)
            DO i3 = startidx(3), endidx(3)
               DO i2 = startidx(2), endidx(2)
                  DO i1 = startidx(1), endidx(1)
                     CALL next_random_integer_isp(rng_integer, tmprng)
                     randnum(i1,i2,i3,i4,i5) = tmprng
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_5d

#if __MAX_RANK > 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 6D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_6d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(6) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i6 = startidx(6), endidx(6)
         DO i5 = startidx(5), endidx(5)
            DO i4 = startidx(4), endidx(4)
               DO i3 = startidx(3), endidx(3)
                  DO i2 = startidx(2), endidx(2)
                     DO i1 = startidx(1), endidx(1)
                        CALL next_random_integer_isp(rng_integer, tmprng)
                        randnum(i1,i2,i3,i4,i5,i6) = tmprng
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_6d

#if __MAX_RANK > 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 7D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_7d(rng_integer, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(7) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i7 = startidx(7), endidx(7)
         DO i6 = startidx(6), endidx(6)
            DO i5 = startidx(5), endidx(5)
               DO i4 = startidx(4), endidx(4)
                  DO i3 = startidx(3), endidx(3)
                     DO i2 = startidx(2), endidx(2)
                        DO i1 = startidx(1), endidx(1)
                           CALL next_random_integer_isp(rng_integer, tmprng)
                           randnum(i1,i2,i3,i4,i5,i6,i7) = tmprng
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_7d

#if __MAX_RANK > 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 8D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_8d(rng_integer, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(8) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i8 = startidx(8), endidx(8)
         DO i7 = startidx(7), endidx(7)
            DO i6 = startidx(6), endidx(6)
               DO i5 = startidx(5), endidx(5)
                  DO i4 = startidx(4), endidx(4)
                     DO i3 = startidx(3), endidx(3)
                        DO i2 = startidx(2), endidx(2)
                           DO i1 = startidx(1), endidx(1)
                              CALL next_random_integer_isp(rng_integer, tmprng)
                              randnum(i1,i2,i3,i4,i5,i6,i7,i8) = tmprng
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_8d

#if __MAX_RANK > 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 9D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_9d(rng_integer, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(9) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i9 = startidx(9), endidx(9)
         DO i8 = startidx(8), endidx(8)
            DO i7 = startidx(7), endidx(7)
               DO i6 = startidx(6), endidx(6)
                  DO i5 = startidx(5), endidx(5)
                     DO i4 = startidx(4), endidx(4)
                        DO i3 = startidx(3), endidx(3)
                           DO i2 = startidx(2), endidx(2)
                              DO i1 = startidx(1), endidx(1)
                                 CALL next_random_integer_isp(rng_integer, tmprng)
                                 randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9) = tmprng
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_9d

#if __MAX_RANK > 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 10D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_10d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(10) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i10 = startidx(10), endidx(10)
         DO i9 = startidx(9), endidx(9)
            DO i8 = startidx(8), endidx(8)
               DO i7 = startidx(7), endidx(7)
                  DO i6 = startidx(6), endidx(6)
                     DO i5 = startidx(5), endidx(5)
                        DO i4 = startidx(4), endidx(4)
                           DO i3 = startidx(3), endidx(3)
                              DO i2 = startidx(2), endidx(2)
                                 DO i1 = startidx(1), endidx(1)
                                    CALL next_random_integer_isp(rng_integer, tmprng)
                                    randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10) = tmprng
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_10d

#if __MAX_RANK > 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 11D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_11d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(11) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i11 = startidx(11), endidx(11)
         DO i10 = startidx(10), endidx(10)
            DO i9 = startidx(9), endidx(9)
               DO i8 = startidx(8), endidx(8)
                  DO i7 = startidx(7), endidx(7)
                     DO i6 = startidx(6), endidx(6)
                        DO i5 = startidx(5), endidx(5)
                           DO i4 = startidx(4), endidx(4)
                              DO i3 = startidx(3), endidx(3)
                                 DO i2 = startidx(2), endidx(2)
                                    DO i1 = startidx(1), endidx(1)
                                       CALL next_random_integer_isp(rng_integer, tmprng)
                                       randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11) = tmprng
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_11d

#if __MAX_RANK > 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 12D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_12d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(12) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i12 = startidx(12), endidx(12)
         DO i11 = startidx(11), endidx(11)
            DO i10 = startidx(10), endidx(10)
               DO i9 = startidx(9), endidx(9)
                  DO i8 = startidx(8), endidx(8)
                     DO i7 = startidx(7), endidx(7)
                        DO i6 = startidx(6), endidx(6)
                           DO i5 = startidx(5), endidx(5)
                              DO i4 = startidx(4), endidx(4)
                                 DO i3 = startidx(3), endidx(3)
                                    DO i2 = startidx(2), endidx(2)
                                       DO i1 = startidx(1), endidx(1)
                                          CALL next_random_integer_isp(rng_integer, tmprng)
                                          randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) = tmprng
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_12d

#if __MAX_RANK > 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 13D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_13d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(13) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i13 = startidx(13), endidx(13)
         DO i12 = startidx(12), endidx(12)
            DO i11 = startidx(11), endidx(11)
               DO i10 = startidx(10), endidx(10)
                  DO i9 = startidx(9), endidx(9)
                     DO i8 = startidx(8), endidx(8)
                        DO i7 = startidx(7), endidx(7)
                           DO i6 = startidx(6), endidx(6)
                              DO i5 = startidx(5), endidx(5)
                                 DO i4 = startidx(4), endidx(4)
                                    DO i3 = startidx(3), endidx(3)
                                       DO i2 = startidx(2), endidx(2)
                                          DO i1 = startidx(1), endidx(1)
                                             CALL next_random_integer_isp(rng_integer, tmprng)
                                             randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13) = tmprng
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_13d

#if __MAX_RANK > 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 14D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_14d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(14) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i14 = startidx(14), endidx(14)
         DO i13 = startidx(13), endidx(13)
            DO i12 = startidx(12), endidx(12)
               DO i11 = startidx(11), endidx(11)
                  DO i10 = startidx(10), endidx(10)
                     DO i9 = startidx(9), endidx(9)
                        DO i8 = startidx(8), endidx(8)
                           DO i7 = startidx(7), endidx(7)
                              DO i6 = startidx(6), endidx(6)
                                 DO i5 = startidx(5), endidx(5)
                                    DO i4 = startidx(4), endidx(4)
                                       DO i3 = startidx(3), endidx(3)
                                          DO i2 = startidx(2), endidx(2)
                                             DO i1 = startidx(1), endidx(1)
                                                CALL next_random_integer_isp(rng_integer, tmprng)
                                                randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) = tmprng
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_14d

#if __MAX_RANK > 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 15D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_isp_15d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=isp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=isp) :: tmprng
      INTEGER, DIMENSION(15) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i15 = startidx(15), endidx(15)
         DO i14 = startidx(14), endidx(14)
            DO i13 = startidx(13), endidx(13)
               DO i12 = startidx(12), endidx(12)
                  DO i11 = startidx(11), endidx(11)
                     DO i10 = startidx(10), endidx(10)
                        DO i9 = startidx(9), endidx(9)
                           DO i8 = startidx(8), endidx(8)
                              DO i7 = startidx(7), endidx(7)
                                 DO i6 = startidx(6), endidx(6)
                                    DO i5 = startidx(5), endidx(5)
                                       DO i4 = startidx(4), endidx(4)
                                          DO i3 = startidx(3), endidx(3)
                                             DO i2 = startidx(2), endidx(2)
                                                DO i1 = startidx(1), endidx(1)
                                                   CALL next_random_integer_isp(rng_integer, tmprng)
                                                   randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) = tmprng
                                                END DO
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_isp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(OUT) :: randnum
   
      CALL mersennetwister(rng_integer, randnum) 

      RETURN

   END SUBROUTINE next_random_integer_idp

#if __MAX_RANK > 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 1D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_1d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:), POINTER :: randnum
   
      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(1) :: startidx, endidx
      INTEGER :: i1

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)
      
      DO i1 = startidx(1), endidx(1)
         CALL next_random_integer_idp(rng_integer, tmprng)
         randnum(i1) = tmprng
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_1d

#if __MAX_RANK > 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 2D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_2d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:), POINTER :: randnum
   
      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(2) :: startidx, endidx
      INTEGER :: i1, i2

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i2 = startidx(2), endidx(2)
         DO i1 = startidx(1), endidx(1)
            CALL next_random_integer_idp(rng_integer, tmprng)
            randnum(i1,i2) = tmprng
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_2d

#if __MAX_RANK > 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 3D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_3d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: randnum
   
      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(3) :: startidx, endidx
      INTEGER :: i1, i2, i3

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i3 = startidx(3), endidx(3)
         DO i2 = startidx(2), endidx(2)
            DO i1 = startidx(1), endidx(1)
               CALL next_random_integer_idp(rng_integer, tmprng)
               randnum(i1,i2,i3) = tmprng
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_3d

#if __MAX_RANK > 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 4D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_4d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(4) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i4 = startidx(4), endidx(4)
         DO i3 = startidx(3), endidx(3)
            DO i2 = startidx(2), endidx(2)
               DO i1 = startidx(1), endidx(1)
                  CALL next_random_integer_idp(rng_integer, tmprng)
                  randnum(i1,i2,i3,i4) = tmprng
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_4d

#if __MAX_RANK > 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 5D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_5d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(5) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i5 = startidx(5), endidx(5)
         DO i4 = startidx(4), endidx(4)
            DO i3 = startidx(3), endidx(3)
               DO i2 = startidx(2), endidx(2)
                  DO i1 = startidx(1), endidx(1)
                     CALL next_random_integer_idp(rng_integer, tmprng)
                     randnum(i1,i2,i3,i4,i5) = tmprng
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_5d

#if __MAX_RANK > 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 6D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_6d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(6) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i6 = startidx(6), endidx(6)
         DO i5 = startidx(5), endidx(5)
            DO i4 = startidx(4), endidx(4)
               DO i3 = startidx(3), endidx(3)
                  DO i2 = startidx(2), endidx(2)
                     DO i1 = startidx(1), endidx(1)
                        CALL next_random_integer_idp(rng_integer, tmprng)
                        randnum(i1,i2,i3,i4,i5,i6) = tmprng
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_6d

#if __MAX_RANK > 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 7D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_7d(rng_integer, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(7) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i7 = startidx(7), endidx(7)
         DO i6 = startidx(6), endidx(6)
            DO i5 = startidx(5), endidx(5)
               DO i4 = startidx(4), endidx(4)
                  DO i3 = startidx(3), endidx(3)
                     DO i2 = startidx(2), endidx(2)
                        DO i1 = startidx(1), endidx(1)
                           CALL next_random_integer_idp(rng_integer, tmprng)
                           randnum(i1,i2,i3,i4,i5,i6,i7) = tmprng
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_7d

#if __MAX_RANK > 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 8D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_8d(rng_integer, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(8) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i8 = startidx(8), endidx(8)
         DO i7 = startidx(7), endidx(7)
            DO i6 = startidx(6), endidx(6)
               DO i5 = startidx(5), endidx(5)
                  DO i4 = startidx(4), endidx(4)
                     DO i3 = startidx(3), endidx(3)
                        DO i2 = startidx(2), endidx(2)
                           DO i1 = startidx(1), endidx(1)
                              CALL next_random_integer_idp(rng_integer, tmprng)
                              randnum(i1,i2,i3,i4,i5,i6,i7,i8) = tmprng
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_8d

#if __MAX_RANK > 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 9D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_9d(rng_integer, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(9) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i9 = startidx(9), endidx(9)
         DO i8 = startidx(8), endidx(8)
            DO i7 = startidx(7), endidx(7)
               DO i6 = startidx(6), endidx(6)
                  DO i5 = startidx(5), endidx(5)
                     DO i4 = startidx(4), endidx(4)
                        DO i3 = startidx(3), endidx(3)
                           DO i2 = startidx(2), endidx(2)
                              DO i1 = startidx(1), endidx(1)
                                 CALL next_random_integer_idp(rng_integer, tmprng)
                                 randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9) = tmprng
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_9d

#if __MAX_RANK > 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 10D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_10d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(10) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i10 = startidx(10), endidx(10)
         DO i9 = startidx(9), endidx(9)
            DO i8 = startidx(8), endidx(8)
               DO i7 = startidx(7), endidx(7)
                  DO i6 = startidx(6), endidx(6)
                     DO i5 = startidx(5), endidx(5)
                        DO i4 = startidx(4), endidx(4)
                           DO i3 = startidx(3), endidx(3)
                              DO i2 = startidx(2), endidx(2)
                                 DO i1 = startidx(1), endidx(1)
                                    CALL next_random_integer_idp(rng_integer, tmprng)
                                    randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10) = tmprng
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_10d

#if __MAX_RANK > 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 11D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_11d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(11) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i11 = startidx(11), endidx(11)
         DO i10 = startidx(10), endidx(10)
            DO i9 = startidx(9), endidx(9)
               DO i8 = startidx(8), endidx(8)
                  DO i7 = startidx(7), endidx(7)
                     DO i6 = startidx(6), endidx(6)
                        DO i5 = startidx(5), endidx(5)
                           DO i4 = startidx(4), endidx(4)
                              DO i3 = startidx(3), endidx(3)
                                 DO i2 = startidx(2), endidx(2)
                                    DO i1 = startidx(1), endidx(1)
                                       CALL next_random_integer_idp(rng_integer, tmprng)
                                       randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11) = tmprng
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_11d

#if __MAX_RANK > 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 12D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_12d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(12) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i12 = startidx(12), endidx(12)
         DO i11 = startidx(11), endidx(11)
            DO i10 = startidx(10), endidx(10)
               DO i9 = startidx(9), endidx(9)
                  DO i8 = startidx(8), endidx(8)
                     DO i7 = startidx(7), endidx(7)
                        DO i6 = startidx(6), endidx(6)
                           DO i5 = startidx(5), endidx(5)
                              DO i4 = startidx(4), endidx(4)
                                 DO i3 = startidx(3), endidx(3)
                                    DO i2 = startidx(2), endidx(2)
                                       DO i1 = startidx(1), endidx(1)
                                          CALL next_random_integer_idp(rng_integer, tmprng)
                                          randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) = tmprng
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_12d

#if __MAX_RANK > 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 13D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_13d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(13) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i13 = startidx(13), endidx(13)
         DO i12 = startidx(12), endidx(12)
            DO i11 = startidx(11), endidx(11)
               DO i10 = startidx(10), endidx(10)
                  DO i9 = startidx(9), endidx(9)
                     DO i8 = startidx(8), endidx(8)
                        DO i7 = startidx(7), endidx(7)
                           DO i6 = startidx(6), endidx(6)
                              DO i5 = startidx(5), endidx(5)
                                 DO i4 = startidx(4), endidx(4)
                                    DO i3 = startidx(3), endidx(3)
                                       DO i2 = startidx(2), endidx(2)
                                          DO i1 = startidx(1), endidx(1)
                                             CALL next_random_integer_idp(rng_integer, tmprng)
                                             randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13) = tmprng
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_13d

#if __MAX_RANK > 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 14D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_14d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(14) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i14 = startidx(14), endidx(14)
         DO i13 = startidx(13), endidx(13)
            DO i12 = startidx(12), endidx(12)
               DO i11 = startidx(11), endidx(11)
                  DO i10 = startidx(10), endidx(10)
                     DO i9 = startidx(9), endidx(9)
                        DO i8 = startidx(8), endidx(8)
                           DO i7 = startidx(7), endidx(7)
                              DO i6 = startidx(6), endidx(6)
                                 DO i5 = startidx(5), endidx(5)
                                    DO i4 = startidx(4), endidx(4)
                                       DO i3 = startidx(3), endidx(3)
                                          DO i2 = startidx(2), endidx(2)
                                             DO i1 = startidx(1), endidx(1)
                                                CALL next_random_integer_idp(rng_integer, tmprng)
                                                randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) = tmprng
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_14d

#if __MAX_RANK > 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 15D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_idp_15d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=idp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=idp) :: tmprng
      INTEGER, DIMENSION(15) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i15 = startidx(15), endidx(15)
         DO i14 = startidx(14), endidx(14)
            DO i13 = startidx(13), endidx(13)
               DO i12 = startidx(12), endidx(12)
                  DO i11 = startidx(11), endidx(11)
                     DO i10 = startidx(10), endidx(10)
                        DO i9 = startidx(9), endidx(9)
                           DO i8 = startidx(8), endidx(8)
                              DO i7 = startidx(7), endidx(7)
                                 DO i6 = startidx(6), endidx(6)
                                    DO i5 = startidx(5), endidx(5)
                                       DO i4 = startidx(4), endidx(4)
                                          DO i3 = startidx(3), endidx(3)
                                             DO i2 = startidx(2), endidx(2)
                                                DO i1 = startidx(1), endidx(1)
                                                   CALL next_random_integer_idp(rng_integer, tmprng)
                                                   randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) = tmprng
                                                END DO
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_idp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#ifdef __HAS_IQP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(OUT) :: randnum

      CALL mersennetwister_iqp(rng_integer, randnum) 

      RETURN

   END SUBROUTINE next_random_integer_iqp

#if __MAX_RANK > 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 1D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_1d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:), POINTER :: randnum
   
      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(1) :: startidx, endidx
      INTEGER :: i1

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)
      
      DO i1 = startidx(1), endidx(1)
         CALL next_random_integer_iqp(rng_integer, tmprng)
         randnum(i1) = tmprng
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_1d

#if __MAX_RANK > 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 2D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_2d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:), POINTER :: randnum
   
      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(2) :: startidx, endidx
      INTEGER :: i1, i2

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i2 = startidx(2), endidx(2)
         DO i1 = startidx(1), endidx(1)
            CALL next_random_integer_iqp(rng_integer, tmprng)
            randnum(i1,i2) = tmprng
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_2d

#if __MAX_RANK > 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 3D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_3d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: randnum
   
      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(3) :: startidx, endidx
      INTEGER :: i1, i2, i3

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i3 = startidx(3), endidx(3)
         DO i2 = startidx(2), endidx(2)
            DO i1 = startidx(1), endidx(1)
               CALL next_random_integer_iqp(rng_integer, tmprng)
               randnum(i1,i2,i3) = tmprng
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_3d

#if __MAX_RANK > 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 4D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_4d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(4) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i4 = startidx(4), endidx(4)
         DO i3 = startidx(3), endidx(3)
            DO i2 = startidx(2), endidx(2)
               DO i1 = startidx(1), endidx(1)
                  CALL next_random_integer_iqp(rng_integer, tmprng)
                  randnum(i1,i2,i3,i4) = tmprng
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_4d

#if __MAX_RANK > 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 5D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_5d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(5) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i5 = startidx(5), endidx(5)
         DO i4 = startidx(4), endidx(4)
            DO i3 = startidx(3), endidx(3)
               DO i2 = startidx(2), endidx(2)
                  DO i1 = startidx(1), endidx(1)
                     CALL next_random_integer_iqp(rng_integer, tmprng)
                     randnum(i1,i2,i3,i4,i5) = tmprng
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_5d

#if __MAX_RANK > 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 6D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_6d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(6) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i6 = startidx(6), endidx(6)
         DO i5 = startidx(5), endidx(5)
            DO i4 = startidx(4), endidx(4)
               DO i3 = startidx(3), endidx(3)
                  DO i2 = startidx(2), endidx(2)
                     DO i1 = startidx(1), endidx(1)
                        CALL next_random_integer_iqp(rng_integer, tmprng)
                        randnum(i1,i2,i3,i4,i5,i6) = tmprng
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_6d

#if __MAX_RANK > 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 7D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_7d(rng_integer, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(7) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i7 = startidx(7), endidx(7)
         DO i6 = startidx(6), endidx(6)
            DO i5 = startidx(5), endidx(5)
               DO i4 = startidx(4), endidx(4)
                  DO i3 = startidx(3), endidx(3)
                     DO i2 = startidx(2), endidx(2)
                        DO i1 = startidx(1), endidx(1)
                           CALL next_random_integer_iqp(rng_integer, tmprng)
                           randnum(i1,i2,i3,i4,i5,i6,i7) = tmprng
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_7d

#if __MAX_RANK > 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 8D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_8d(rng_integer, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(8) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i8 = startidx(8), endidx(8)
         DO i7 = startidx(7), endidx(7)
            DO i6 = startidx(6), endidx(6)
               DO i5 = startidx(5), endidx(5)
                  DO i4 = startidx(4), endidx(4)
                     DO i3 = startidx(3), endidx(3)
                        DO i2 = startidx(2), endidx(2)
                           DO i1 = startidx(1), endidx(1)
                              CALL next_random_integer_iqp(rng_integer, tmprng)
                              randnum(i1,i2,i3,i4,i5,i6,i7,i8) = tmprng
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_8d

#if __MAX_RANK > 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 9D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_9d(rng_integer, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: randnum
   
      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(9) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i9 = startidx(9), endidx(9)
         DO i8 = startidx(8), endidx(8)
            DO i7 = startidx(7), endidx(7)
               DO i6 = startidx(6), endidx(6)
                  DO i5 = startidx(5), endidx(5)
                     DO i4 = startidx(4), endidx(4)
                        DO i3 = startidx(3), endidx(3)
                           DO i2 = startidx(2), endidx(2)
                              DO i1 = startidx(1), endidx(1)
                                 CALL next_random_integer_iqp(rng_integer, tmprng)
                                 randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9) = tmprng
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_9d

#if __MAX_RANK > 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 10D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_10d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(10) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i10 = startidx(10), endidx(10)
         DO i9 = startidx(9), endidx(9)
            DO i8 = startidx(8), endidx(8)
               DO i7 = startidx(7), endidx(7)
                  DO i6 = startidx(6), endidx(6)
                     DO i5 = startidx(5), endidx(5)
                        DO i4 = startidx(4), endidx(4)
                           DO i3 = startidx(3), endidx(3)
                              DO i2 = startidx(2), endidx(2)
                                 DO i1 = startidx(1), endidx(1)
                                    CALL next_random_integer_iqp(rng_integer, tmprng)
                                    randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10) = tmprng
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_10d

#if __MAX_RANK > 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 11D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_11d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(11) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i11 = startidx(11), endidx(11)
         DO i10 = startidx(10), endidx(10)
            DO i9 = startidx(9), endidx(9)
               DO i8 = startidx(8), endidx(8)
                  DO i7 = startidx(7), endidx(7)
                     DO i6 = startidx(6), endidx(6)
                        DO i5 = startidx(5), endidx(5)
                           DO i4 = startidx(4), endidx(4)
                              DO i3 = startidx(3), endidx(3)
                                 DO i2 = startidx(2), endidx(2)
                                    DO i1 = startidx(1), endidx(1)
                                       CALL next_random_integer_iqp(rng_integer, tmprng)
                                       randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11) = tmprng
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_11d

#if __MAX_RANK > 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 12D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_12d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(12) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i12 = startidx(12), endidx(12)
         DO i11 = startidx(11), endidx(11)
            DO i10 = startidx(10), endidx(10)
               DO i9 = startidx(9), endidx(9)
                  DO i8 = startidx(8), endidx(8)
                     DO i7 = startidx(7), endidx(7)
                        DO i6 = startidx(6), endidx(6)
                           DO i5 = startidx(5), endidx(5)
                              DO i4 = startidx(4), endidx(4)
                                 DO i3 = startidx(3), endidx(3)
                                    DO i2 = startidx(2), endidx(2)
                                       DO i1 = startidx(1), endidx(1)
                                          CALL next_random_integer_iqp(rng_integer, tmprng)
                                          randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) = tmprng
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_12d

#if __MAX_RANK > 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 13D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_13d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(13) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i13 = startidx(13), endidx(13)
         DO i12 = startidx(12), endidx(12)
            DO i11 = startidx(11), endidx(11)
               DO i10 = startidx(10), endidx(10)
                  DO i9 = startidx(9), endidx(9)
                     DO i8 = startidx(8), endidx(8)
                        DO i7 = startidx(7), endidx(7)
                           DO i6 = startidx(6), endidx(6)
                              DO i5 = startidx(5), endidx(5)
                                 DO i4 = startidx(4), endidx(4)
                                    DO i3 = startidx(3), endidx(3)
                                       DO i2 = startidx(2), endidx(2)
                                          DO i1 = startidx(1), endidx(1)
                                             CALL next_random_integer_iqp(rng_integer, tmprng)
                                             randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13) = tmprng
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_13d

#if __MAX_RANK > 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 14D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_14d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(14) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i14 = startidx(14), endidx(14)
         DO i13 = startidx(13), endidx(13)
            DO i12 = startidx(12), endidx(12)
               DO i11 = startidx(11), endidx(11)
                  DO i10 = startidx(10), endidx(10)
                     DO i9 = startidx(9), endidx(9)
                        DO i8 = startidx(8), endidx(8)
                           DO i7 = startidx(7), endidx(7)
                              DO i6 = startidx(6), endidx(6)
                                 DO i5 = startidx(5), endidx(5)
                                    DO i4 = startidx(4), endidx(4)
                                       DO i3 = startidx(3), endidx(3)
                                          DO i2 = startidx(2), endidx(2)
                                             DO i1 = startidx(1), endidx(1)
                                                CALL next_random_integer_iqp(rng_integer, tmprng)
                                                randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) = tmprng
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_14d

#if __MAX_RANK > 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next integer random number for an 15D array
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_integer_iqp_15d(rng_integer, randnum)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(INOUT) :: rng_integer
      INTEGER(KIND=iqp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      INTEGER(KIND=iqp) :: tmprng
      INTEGER, DIMENSION(15) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i15 = startidx(15), endidx(15)
         DO i14 = startidx(14), endidx(14)
            DO i13 = startidx(13), endidx(13)
               DO i12 = startidx(12), endidx(12)
                  DO i11 = startidx(11), endidx(11)
                     DO i10 = startidx(10), endidx(10)
                        DO i9 = startidx(9), endidx(9)
                           DO i8 = startidx(8), endidx(8)
                              DO i7 = startidx(7), endidx(7)
                                 DO i6 = startidx(6), endidx(6)
                                    DO i5 = startidx(5), endidx(5)
                                       DO i4 = startidx(4), endidx(4)
                                          DO i3 = startidx(3), endidx(3)
                                             DO i2 = startidx(2), endidx(2)
                                                DO i1 = startidx(1), endidx(1)
                                                   CALL next_random_integer_iqp(rng_integer, tmprng)
                                                   randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) = tmprng
                                                END DO
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_integer_iqp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(OUT) :: randnum
   
      CALL mersennetwister_sp(rng_uniform, randnum) 
      randnum = randnum * REAL(rng_uniform%upper - rng_uniform%lower,sp) + REAL(rng_uniform%lower,sp)

      RETURN

   END SUBROUTINE next_random_uniform_sp

#if __MAX_RANK > 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 1D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_1d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(1) :: startidx, endidx
      INTEGER :: i1

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)
      
      DO i1 = startidx(1), endidx(1)
         CALL next_random_uniform_sp(rng_uniform, tmprng)
         randnum(i1) = tmprng
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_1d

#if __MAX_RANK > 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 2D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_2d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(2) :: startidx, endidx
      INTEGER :: i1, i2

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i2 = startidx(2), endidx(2)
         DO i1 = startidx(1), endidx(1)
            CALL next_random_uniform_sp(rng_uniform, tmprng)
            randnum(i1,i2) = tmprng
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_2d

#if __MAX_RANK > 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 3D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_3d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(3) :: startidx, endidx
      INTEGER :: i1, i2, i3

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i3 = startidx(3), endidx(3)
         DO i2 = startidx(2), endidx(2)
            DO i1 = startidx(1), endidx(1)
               CALL next_random_uniform_sp(rng_uniform, tmprng)
               randnum(i1,i2,i3) = tmprng
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_3d

#if __MAX_RANK > 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 4D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_4d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(4) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i4 = startidx(4), endidx(4)
         DO i3 = startidx(3), endidx(3)
            DO i2 = startidx(2), endidx(2)
               DO i1 = startidx(1), endidx(1)
                  CALL next_random_uniform_sp(rng_uniform, tmprng)
                  randnum(i1,i2,i3,i4) = tmprng
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_4d

#if __MAX_RANK > 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 5D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_5d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(5) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i5 = startidx(5), endidx(5)
         DO i4 = startidx(4), endidx(4)
            DO i3 = startidx(3), endidx(3)
               DO i2 = startidx(2), endidx(2)
                  DO i1 = startidx(1), endidx(1)
                     CALL next_random_uniform_sp(rng_uniform, tmprng)
                     randnum(i1,i2,i3,i4,i5) = tmprng
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_5d

#if __MAX_RANK > 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 6D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_6d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(6) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i6 = startidx(6), endidx(6)
         DO i5 = startidx(5), endidx(5)
            DO i4 = startidx(4), endidx(4)
               DO i3 = startidx(3), endidx(3)
                  DO i2 = startidx(2), endidx(2)
                     DO i1 = startidx(1), endidx(1)
                        CALL next_random_uniform_sp(rng_uniform, tmprng)
                        randnum(i1,i2,i3,i4,i5,i6) = tmprng
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_6d

#if __MAX_RANK > 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 7D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_7d(rng_uniform, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(7) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i7 = startidx(7), endidx(7)
         DO i6 = startidx(6), endidx(6)
            DO i5 = startidx(5), endidx(5)
               DO i4 = startidx(4), endidx(4)
                  DO i3 = startidx(3), endidx(3)
                     DO i2 = startidx(2), endidx(2)
                        DO i1 = startidx(1), endidx(1)
                           CALL next_random_uniform_sp(rng_uniform, tmprng)
                           randnum(i1,i2,i3,i4,i5,i6,i7) = tmprng
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_7d

#if __MAX_RANK > 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 8D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_8d(rng_uniform, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(8) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i8 = startidx(8), endidx(8)
         DO i7 = startidx(7), endidx(7)
            DO i6 = startidx(6), endidx(6)
               DO i5 = startidx(5), endidx(5)
                  DO i4 = startidx(4), endidx(4)
                     DO i3 = startidx(3), endidx(3)
                        DO i2 = startidx(2), endidx(2)
                           DO i1 = startidx(1), endidx(1)
                              CALL next_random_uniform_sp(rng_uniform, tmprng)
                              randnum(i1,i2,i3,i4,i5,i6,i7,i8) = tmprng
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_8d

#if __MAX_RANK > 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 9D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_9d(rng_uniform, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(9) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i9 = startidx(9), endidx(9)
         DO i8 = startidx(8), endidx(8)
            DO i7 = startidx(7), endidx(7)
               DO i6 = startidx(6), endidx(6)
                  DO i5 = startidx(5), endidx(5)
                     DO i4 = startidx(4), endidx(4)
                        DO i3 = startidx(3), endidx(3)
                           DO i2 = startidx(2), endidx(2)
                              DO i1 = startidx(1), endidx(1)
                                 CALL next_random_uniform_sp(rng_uniform, tmprng)
                                 randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9) = tmprng
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_9d

#if __MAX_RANK > 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 10D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_10d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(10) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i10 = startidx(10), endidx(10)
         DO i9 = startidx(9), endidx(9)
            DO i8 = startidx(8), endidx(8)
               DO i7 = startidx(7), endidx(7)
                  DO i6 = startidx(6), endidx(6)
                     DO i5 = startidx(5), endidx(5)
                        DO i4 = startidx(4), endidx(4)
                           DO i3 = startidx(3), endidx(3)
                              DO i2 = startidx(2), endidx(2)
                                 DO i1 = startidx(1), endidx(1)
                                    CALL next_random_uniform_sp(rng_uniform, tmprng)
                                    randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10) = tmprng
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_10d

#if __MAX_RANK > 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 11D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_11d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(11) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i11 = startidx(11), endidx(11)
         DO i10 = startidx(10), endidx(10)
            DO i9 = startidx(9), endidx(9)
               DO i8 = startidx(8), endidx(8)
                  DO i7 = startidx(7), endidx(7)
                     DO i6 = startidx(6), endidx(6)
                        DO i5 = startidx(5), endidx(5)
                           DO i4 = startidx(4), endidx(4)
                              DO i3 = startidx(3), endidx(3)
                                 DO i2 = startidx(2), endidx(2)
                                    DO i1 = startidx(1), endidx(1)
                                       CALL next_random_uniform_sp(rng_uniform, tmprng)
                                       randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11) = tmprng
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_11d

#if __MAX_RANK > 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 12D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_12d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(12) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i12 = startidx(12), endidx(12)
         DO i11 = startidx(11), endidx(11)
            DO i10 = startidx(10), endidx(10)
               DO i9 = startidx(9), endidx(9)
                  DO i8 = startidx(8), endidx(8)
                     DO i7 = startidx(7), endidx(7)
                        DO i6 = startidx(6), endidx(6)
                           DO i5 = startidx(5), endidx(5)
                              DO i4 = startidx(4), endidx(4)
                                 DO i3 = startidx(3), endidx(3)
                                    DO i2 = startidx(2), endidx(2)
                                       DO i1 = startidx(1), endidx(1)
                                          CALL next_random_uniform_sp(rng_uniform, tmprng)
                                          randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) = tmprng
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_12d

#if __MAX_RANK > 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 13D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_13d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(13) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i13 = startidx(13), endidx(13)
         DO i12 = startidx(12), endidx(12)
            DO i11 = startidx(11), endidx(11)
               DO i10 = startidx(10), endidx(10)
                  DO i9 = startidx(9), endidx(9)
                     DO i8 = startidx(8), endidx(8)
                        DO i7 = startidx(7), endidx(7)
                           DO i6 = startidx(6), endidx(6)
                              DO i5 = startidx(5), endidx(5)
                                 DO i4 = startidx(4), endidx(4)
                                    DO i3 = startidx(3), endidx(3)
                                       DO i2 = startidx(2), endidx(2)
                                          DO i1 = startidx(1), endidx(1)
                                             CALL next_random_uniform_sp(rng_uniform, tmprng)
                                             randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13) = tmprng
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_13d

#if __MAX_RANK > 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 14D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_14d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(14) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i14 = startidx(14), endidx(14)
         DO i13 = startidx(13), endidx(13)
            DO i12 = startidx(12), endidx(12)
               DO i11 = startidx(11), endidx(11)
                  DO i10 = startidx(10), endidx(10)
                     DO i9 = startidx(9), endidx(9)
                        DO i8 = startidx(8), endidx(8)
                           DO i7 = startidx(7), endidx(7)
                              DO i6 = startidx(6), endidx(6)
                                 DO i5 = startidx(5), endidx(5)
                                    DO i4 = startidx(4), endidx(4)
                                       DO i3 = startidx(3), endidx(3)
                                          DO i2 = startidx(2), endidx(2)
                                             DO i1 = startidx(1), endidx(1)
                                                CALL next_random_uniform_sp(rng_uniform, tmprng)
                                                randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) = tmprng
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_14d

#if __MAX_RANK > 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 15D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_sp_15d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(15) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i15 = startidx(15), endidx(15)
         DO i14 = startidx(14), endidx(14)
            DO i13 = startidx(13), endidx(13)
               DO i12 = startidx(12), endidx(12)
                  DO i11 = startidx(11), endidx(11)
                     DO i10 = startidx(10), endidx(10)
                        DO i9 = startidx(9), endidx(9)
                           DO i8 = startidx(8), endidx(8)
                              DO i7 = startidx(7), endidx(7)
                                 DO i6 = startidx(6), endidx(6)
                                    DO i5 = startidx(5), endidx(5)
                                       DO i4 = startidx(4), endidx(4)
                                          DO i3 = startidx(3), endidx(3)
                                             DO i2 = startidx(2), endidx(2)
                                                DO i1 = startidx(1), endidx(1)
                                                   CALL next_random_uniform_sp(rng_uniform, tmprng)
                                                   randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) = tmprng
                                                END DO
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_sp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(OUT) :: randnum
   
      CALL mersennetwister(rng_uniform, randnum) 
      randnum = randnum * &
                (rng_uniform%upper - rng_uniform%lower) + &
                rng_uniform%lower

      RETURN

   END SUBROUTINE next_random_uniform_dp

#if __MAX_RANK > 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 1D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_1d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(1) :: startidx, endidx
      INTEGER :: i1

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)
      
      DO i1 = startidx(1), endidx(1)
         CALL next_random_uniform_dp(rng_uniform, tmprng)
         randnum(i1) = tmprng
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_1d

#if __MAX_RANK > 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 2D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_2d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(2) :: startidx, endidx
      INTEGER :: i1, i2

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i2 = startidx(2), endidx(2)
         DO i1 = startidx(1), endidx(1)
            CALL next_random_uniform_dp(rng_uniform, tmprng)
            randnum(i1,i2) = tmprng
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_2d

#if __MAX_RANK > 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 3D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_3d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(3) :: startidx, endidx
      INTEGER :: i1, i2, i3

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i3 = startidx(3), endidx(3)
         DO i2 = startidx(2), endidx(2)
            DO i1 = startidx(1), endidx(1)
               CALL next_random_uniform_dp(rng_uniform, tmprng)
               randnum(i1,i2,i3) = tmprng
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_3d

#if __MAX_RANK > 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 4D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_4d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(4) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i4 = startidx(4), endidx(4)
         DO i3 = startidx(3), endidx(3)
            DO i2 = startidx(2), endidx(2)
               DO i1 = startidx(1), endidx(1)
                  CALL next_random_uniform_dp(rng_uniform, tmprng)
                  randnum(i1,i2,i3,i4) = tmprng
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_4d

#if __MAX_RANK > 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 5D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_5d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(5) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i5 = startidx(5), endidx(5)
         DO i4 = startidx(4), endidx(4)
            DO i3 = startidx(3), endidx(3)
               DO i2 = startidx(2), endidx(2)
                  DO i1 = startidx(1), endidx(1)
                     CALL next_random_uniform_dp(rng_uniform, tmprng)
                     randnum(i1,i2,i3,i4,i5) = tmprng
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_5d

#if __MAX_RANK > 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 6D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_6d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(6) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i6 = startidx(6), endidx(6)
         DO i5 = startidx(5), endidx(5)
            DO i4 = startidx(4), endidx(4)
               DO i3 = startidx(3), endidx(3)
                  DO i2 = startidx(2), endidx(2)
                     DO i1 = startidx(1), endidx(1)
                        CALL next_random_uniform_dp(rng_uniform, tmprng)
                        randnum(i1,i2,i3,i4,i5,i6) = tmprng
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_6d

#if __MAX_RANK > 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 7D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_7d(rng_uniform, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(7) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i7 = startidx(7), endidx(7)
         DO i6 = startidx(6), endidx(6)
            DO i5 = startidx(5), endidx(5)
               DO i4 = startidx(4), endidx(4)
                  DO i3 = startidx(3), endidx(3)
                     DO i2 = startidx(2), endidx(2)
                        DO i1 = startidx(1), endidx(1)
                           CALL next_random_uniform_dp(rng_uniform, tmprng)
                           randnum(i1,i2,i3,i4,i5,i6,i7) = tmprng
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_7d

#if __MAX_RANK > 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 8D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_8d(rng_uniform, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(8) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i8 = startidx(8), endidx(8)
         DO i7 = startidx(7), endidx(7)
            DO i6 = startidx(6), endidx(6)
               DO i5 = startidx(5), endidx(5)
                  DO i4 = startidx(4), endidx(4)
                     DO i3 = startidx(3), endidx(3)
                        DO i2 = startidx(2), endidx(2)
                           DO i1 = startidx(1), endidx(1)
                              CALL next_random_uniform_dp(rng_uniform, tmprng)
                              randnum(i1,i2,i3,i4,i5,i6,i7,i8) = tmprng
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_8d

#if __MAX_RANK > 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 9D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_9d(rng_uniform, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(9) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i9 = startidx(9), endidx(9)
         DO i8 = startidx(8), endidx(8)
            DO i7 = startidx(7), endidx(7)
               DO i6 = startidx(6), endidx(6)
                  DO i5 = startidx(5), endidx(5)
                     DO i4 = startidx(4), endidx(4)
                        DO i3 = startidx(3), endidx(3)
                           DO i2 = startidx(2), endidx(2)
                              DO i1 = startidx(1), endidx(1)
                                 CALL next_random_uniform_dp(rng_uniform, tmprng)
                                 randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9) = tmprng
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_9d

#if __MAX_RANK > 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 10D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_10d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(10) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i10 = startidx(10), endidx(10)
         DO i9 = startidx(9), endidx(9)
            DO i8 = startidx(8), endidx(8)
               DO i7 = startidx(7), endidx(7)
                  DO i6 = startidx(6), endidx(6)
                     DO i5 = startidx(5), endidx(5)
                        DO i4 = startidx(4), endidx(4)
                           DO i3 = startidx(3), endidx(3)
                              DO i2 = startidx(2), endidx(2)
                                 DO i1 = startidx(1), endidx(1)
                                    CALL next_random_uniform_dp(rng_uniform, tmprng)
                                    randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10) = tmprng
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_10d

#if __MAX_RANK > 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 11D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_11d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(11) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i11 = startidx(11), endidx(11)
         DO i10 = startidx(10), endidx(10)
            DO i9 = startidx(9), endidx(9)
               DO i8 = startidx(8), endidx(8)
                  DO i7 = startidx(7), endidx(7)
                     DO i6 = startidx(6), endidx(6)
                        DO i5 = startidx(5), endidx(5)
                           DO i4 = startidx(4), endidx(4)
                              DO i3 = startidx(3), endidx(3)
                                 DO i2 = startidx(2), endidx(2)
                                    DO i1 = startidx(1), endidx(1)
                                       CALL next_random_uniform_dp(rng_uniform, tmprng)
                                       randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11) = tmprng
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_11d

#if __MAX_RANK > 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 12D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_12d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(12) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i12 = startidx(12), endidx(12)
         DO i11 = startidx(11), endidx(11)
            DO i10 = startidx(10), endidx(10)
               DO i9 = startidx(9), endidx(9)
                  DO i8 = startidx(8), endidx(8)
                     DO i7 = startidx(7), endidx(7)
                        DO i6 = startidx(6), endidx(6)
                           DO i5 = startidx(5), endidx(5)
                              DO i4 = startidx(4), endidx(4)
                                 DO i3 = startidx(3), endidx(3)
                                    DO i2 = startidx(2), endidx(2)
                                       DO i1 = startidx(1), endidx(1)
                                          CALL next_random_uniform_dp(rng_uniform, tmprng)
                                          randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) = tmprng
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_12d

#if __MAX_RANK > 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 13D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_13d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(13) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i13 = startidx(13), endidx(13)
         DO i12 = startidx(12), endidx(12)
            DO i11 = startidx(11), endidx(11)
               DO i10 = startidx(10), endidx(10)
                  DO i9 = startidx(9), endidx(9)
                     DO i8 = startidx(8), endidx(8)
                        DO i7 = startidx(7), endidx(7)
                           DO i6 = startidx(6), endidx(6)
                              DO i5 = startidx(5), endidx(5)
                                 DO i4 = startidx(4), endidx(4)
                                    DO i3 = startidx(3), endidx(3)
                                       DO i2 = startidx(2), endidx(2)
                                          DO i1 = startidx(1), endidx(1)
                                             CALL next_random_uniform_dp(rng_uniform, tmprng)
                                             randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13) = tmprng
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_13d

#if __MAX_RANK > 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 14D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_14d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(14) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i14 = startidx(14), endidx(14)
         DO i13 = startidx(13), endidx(13)
            DO i12 = startidx(12), endidx(12)
               DO i11 = startidx(11), endidx(11)
                  DO i10 = startidx(10), endidx(10)
                     DO i9 = startidx(9), endidx(9)
                        DO i8 = startidx(8), endidx(8)
                           DO i7 = startidx(7), endidx(7)
                              DO i6 = startidx(6), endidx(6)
                                 DO i5 = startidx(5), endidx(5)
                                    DO i4 = startidx(4), endidx(4)
                                       DO i3 = startidx(3), endidx(3)
                                          DO i2 = startidx(2), endidx(2)
                                             DO i1 = startidx(1), endidx(1)
                                                CALL next_random_uniform_dp(rng_uniform, tmprng)
                                                randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) = tmprng
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_14d

#if __MAX_RANK > 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 15D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_dp_15d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(15) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i15 = startidx(15), endidx(15)
         DO i14 = startidx(14), endidx(14)
            DO i13 = startidx(13), endidx(13)
               DO i12 = startidx(12), endidx(12)
                  DO i11 = startidx(11), endidx(11)
                     DO i10 = startidx(10), endidx(10)
                        DO i9 = startidx(9), endidx(9)
                           DO i8 = startidx(8), endidx(8)
                              DO i7 = startidx(7), endidx(7)
                                 DO i6 = startidx(6), endidx(6)
                                    DO i5 = startidx(5), endidx(5)
                                       DO i4 = startidx(4), endidx(4)
                                          DO i3 = startidx(3), endidx(3)
                                             DO i2 = startidx(2), endidx(2)
                                                DO i1 = startidx(1), endidx(1)
                                                   CALL next_random_uniform_dp(rng_uniform, tmprng)
                                                   randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) = tmprng
                                                END DO
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_dp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#ifdef __HAS_QP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(OUT) :: randnum
   
      CALL mersennetwister(rng_uniform, randnum) 
      randnum = randnum * &
                REAL(rng_uniform%upper - rng_uniform%lower, qp) + &
                REAL(rng_uniform%lower,qp)

      RETURN

   END SUBROUTINE next_random_uniform_qp

#if __MAX_RANK > 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 1D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_1d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(1) :: startidx, endidx
      INTEGER :: i1

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)
      
      DO i1 = startidx(1), endidx(1)
         CALL next_random_uniform_qp(rng_uniform, tmprng)
         randnum(i1) = tmprng
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_1d

#if __MAX_RANK > 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 2D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_2d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(2) :: startidx, endidx
      INTEGER :: i1, i2

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i2 = startidx(2), endidx(2)
         DO i1 = startidx(1), endidx(1)
            CALL next_random_uniform_qp(rng_uniform, tmprng)
            randnum(i1,i2) = tmprng
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_2d

#if __MAX_RANK > 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 3D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_3d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(3) :: startidx, endidx
      INTEGER :: i1, i2, i3

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i3 = startidx(3), endidx(3)
         DO i2 = startidx(2), endidx(2)
            DO i1 = startidx(1), endidx(1)
               CALL next_random_uniform_qp(rng_uniform, tmprng)
               randnum(i1,i2,i3) = tmprng
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_3d

#if __MAX_RANK > 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 4D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_4d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(4) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i4 = startidx(4), endidx(4)
         DO i3 = startidx(3), endidx(3)
            DO i2 = startidx(2), endidx(2)
               DO i1 = startidx(1), endidx(1)
                  CALL next_random_uniform_qp(rng_uniform, tmprng)
                  randnum(i1,i2,i3,i4) = tmprng
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_4d

#if __MAX_RANK > 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 5D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_5d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(5) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i5 = startidx(5), endidx(5)
         DO i4 = startidx(4), endidx(4)
            DO i3 = startidx(3), endidx(3)
               DO i2 = startidx(2), endidx(2)
                  DO i1 = startidx(1), endidx(1)
                     CALL next_random_uniform_qp(rng_uniform, tmprng)
                     randnum(i1,i2,i3,i4,i5) = tmprng
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_5d

#if __MAX_RANK > 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 6D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_6d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(6) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i6 = startidx(6), endidx(6)
         DO i5 = startidx(5), endidx(5)
            DO i4 = startidx(4), endidx(4)
               DO i3 = startidx(3), endidx(3)
                  DO i2 = startidx(2), endidx(2)
                     DO i1 = startidx(1), endidx(1)
                        CALL next_random_uniform_qp(rng_uniform, tmprng)
                        randnum(i1,i2,i3,i4,i5,i6) = tmprng
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_6d

#if __MAX_RANK > 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniformly distributed random number for an 7D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_7d(rng_uniform, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(7) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i7 = startidx(7), endidx(7)
         DO i6 = startidx(6), endidx(6)
            DO i5 = startidx(5), endidx(5)
               DO i4 = startidx(4), endidx(4)
                  DO i3 = startidx(3), endidx(3)
                     DO i2 = startidx(2), endidx(2)
                        DO i1 = startidx(1), endidx(1)
                           CALL next_random_uniform_qp(rng_uniform, tmprng)
                           randnum(i1,i2,i3,i4,i5,i6,i7) = tmprng
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_7d

#if __MAX_RANK > 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 8D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_8d(rng_uniform, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(8) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i8 = startidx(8), endidx(8)
         DO i7 = startidx(7), endidx(7)
            DO i6 = startidx(6), endidx(6)
               DO i5 = startidx(5), endidx(5)
                  DO i4 = startidx(4), endidx(4)
                     DO i3 = startidx(3), endidx(3)
                        DO i2 = startidx(2), endidx(2)
                           DO i1 = startidx(1), endidx(1)
                              CALL next_random_uniform_qp(rng_uniform, tmprng)
                              randnum(i1,i2,i3,i4,i5,i6,i7,i8) = tmprng
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_8d

#if __MAX_RANK > 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 9D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_9d(rng_uniform, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(9) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i9 = startidx(9), endidx(9)
         DO i8 = startidx(8), endidx(8)
            DO i7 = startidx(7), endidx(7)
               DO i6 = startidx(6), endidx(6)
                  DO i5 = startidx(5), endidx(5)
                     DO i4 = startidx(4), endidx(4)
                        DO i3 = startidx(3), endidx(3)
                           DO i2 = startidx(2), endidx(2)
                              DO i1 = startidx(1), endidx(1)
                                 CALL next_random_uniform_qp(rng_uniform, tmprng)
                                 randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9) = tmprng
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_9d

#if __MAX_RANK > 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 10D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_10d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(10) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i10 = startidx(10), endidx(10)
         DO i9 = startidx(9), endidx(9)
            DO i8 = startidx(8), endidx(8)
               DO i7 = startidx(7), endidx(7)
                  DO i6 = startidx(6), endidx(6)
                     DO i5 = startidx(5), endidx(5)
                        DO i4 = startidx(4), endidx(4)
                           DO i3 = startidx(3), endidx(3)
                              DO i2 = startidx(2), endidx(2)
                                 DO i1 = startidx(1), endidx(1)
                                    CALL next_random_uniform_qp(rng_uniform, tmprng)
                                    randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10) = tmprng
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_10d

#if __MAX_RANK > 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 11D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_11d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(11) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i11 = startidx(11), endidx(11)
         DO i10 = startidx(10), endidx(10)
            DO i9 = startidx(9), endidx(9)
               DO i8 = startidx(8), endidx(8)
                  DO i7 = startidx(7), endidx(7)
                     DO i6 = startidx(6), endidx(6)
                        DO i5 = startidx(5), endidx(5)
                           DO i4 = startidx(4), endidx(4)
                              DO i3 = startidx(3), endidx(3)
                                 DO i2 = startidx(2), endidx(2)
                                    DO i1 = startidx(1), endidx(1)
                                       CALL next_random_uniform_qp(rng_uniform, tmprng)
                                       randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11) = tmprng
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_11d

#if __MAX_RANK > 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 12D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_12d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(12) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i12 = startidx(12), endidx(12)
         DO i11 = startidx(11), endidx(11)
            DO i10 = startidx(10), endidx(10)
               DO i9 = startidx(9), endidx(9)
                  DO i8 = startidx(8), endidx(8)
                     DO i7 = startidx(7), endidx(7)
                        DO i6 = startidx(6), endidx(6)
                           DO i5 = startidx(5), endidx(5)
                              DO i4 = startidx(4), endidx(4)
                                 DO i3 = startidx(3), endidx(3)
                                    DO i2 = startidx(2), endidx(2)
                                       DO i1 = startidx(1), endidx(1)
                                          CALL next_random_uniform_qp(rng_uniform, tmprng)
                                          randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) = tmprng
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_12d

#if __MAX_RANK > 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 13D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_13d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(13) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i13 = startidx(13), endidx(13)
         DO i12 = startidx(12), endidx(12)
            DO i11 = startidx(11), endidx(11)
               DO i10 = startidx(10), endidx(10)
                  DO i9 = startidx(9), endidx(9)
                     DO i8 = startidx(8), endidx(8)
                        DO i7 = startidx(7), endidx(7)
                           DO i6 = startidx(6), endidx(6)
                              DO i5 = startidx(5), endidx(5)
                                 DO i4 = startidx(4), endidx(4)
                                    DO i3 = startidx(3), endidx(3)
                                       DO i2 = startidx(2), endidx(2)
                                          DO i1 = startidx(1), endidx(1)
                                             CALL next_random_uniform_qp(rng_uniform, tmprng)
                                             randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13) = tmprng
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_13d

#if __MAX_RANK > 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 14D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_14d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(14) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i14 = startidx(14), endidx(14)
         DO i13 = startidx(13), endidx(13)
            DO i12 = startidx(12), endidx(12)
               DO i11 = startidx(11), endidx(11)
                  DO i10 = startidx(10), endidx(10)
                     DO i9 = startidx(9), endidx(9)
                        DO i8 = startidx(8), endidx(8)
                           DO i7 = startidx(7), endidx(7)
                              DO i6 = startidx(6), endidx(6)
                                 DO i5 = startidx(5), endidx(5)
                                    DO i4 = startidx(4), endidx(4)
                                       DO i3 = startidx(3), endidx(3)
                                          DO i2 = startidx(2), endidx(2)
                                             DO i1 = startidx(1), endidx(1)
                                                CALL next_random_uniform_qp(rng_uniform, tmprng)
                                                randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) = tmprng
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_14d

#if __MAX_RANK > 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next uniform random number for an 15D array
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_uniform_qp_15d(rng_uniform, randnum)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(INOUT) :: rng_uniform
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(15) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i15 = startidx(15), endidx(15)
         DO i14 = startidx(14), endidx(14)
            DO i13 = startidx(13), endidx(13)
               DO i12 = startidx(12), endidx(12)
                  DO i11 = startidx(11), endidx(11)
                     DO i10 = startidx(10), endidx(10)
                        DO i9 = startidx(9), endidx(9)
                           DO i8 = startidx(8), endidx(8)
                              DO i7 = startidx(7), endidx(7)
                                 DO i6 = startidx(6), endidx(6)
                                    DO i5 = startidx(5), endidx(5)
                                       DO i4 = startidx(4), endidx(4)
                                          DO i3 = startidx(3), endidx(3)
                                             DO i2 = startidx(2), endidx(2)
                                                DO i1 = startidx(1), endidx(1)
                                                   CALL next_random_uniform_qp(rng_uniform, tmprng)
                                                   randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) = tmprng
                                                END DO
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_uniform_qp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number by Box-Muller
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp(rng_gaussian, randnum)
      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(OUT) :: randnum

      REAL(KIND=dp) :: tmprnd

      CALL next_random_gaussian_dp(rng_gaussian, tmprnd)

      randnum = REAL(tmprnd, sp)

      RETURN
   END SUBROUTINE next_random_gaussian_sp

#if __MAX_RANK > 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 1D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_1d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(1) :: startidx, endidx
      INTEGER :: i1

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)
      
      DO i1 = startidx(1), endidx(1)
         CALL next_random_gaussian_sp(rng_gaussian, tmprng)
         randnum(i1) = tmprng
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_1d

#if __MAX_RANK > 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 2D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_2d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(2) :: startidx, endidx
      INTEGER :: i1, i2

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i2 = startidx(2), endidx(2)
         DO i1 = startidx(1), endidx(1)
            CALL next_random_gaussian_sp(rng_gaussian, tmprng)
            randnum(i1,i2) = tmprng
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_2d

#if __MAX_RANK > 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 3D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_3d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(3) :: startidx, endidx
      INTEGER :: i1, i2, i3

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i3 = startidx(3), endidx(3)
         DO i2 = startidx(2), endidx(2)
            DO i1 = startidx(1), endidx(1)
               CALL next_random_gaussian_sp(rng_gaussian, tmprng)
               randnum(i1,i2,i3) = tmprng
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_3d

#if __MAX_RANK > 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 4D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_4d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(4) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i4 = startidx(4), endidx(4)
         DO i3 = startidx(3), endidx(3)
            DO i2 = startidx(2), endidx(2)
               DO i1 = startidx(1), endidx(1)
                  CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                  randnum(i1,i2,i3,i4) = tmprng
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_4d

#if __MAX_RANK > 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 5D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_5d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(5) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i5 = startidx(5), endidx(5)
         DO i4 = startidx(4), endidx(4)
            DO i3 = startidx(3), endidx(3)
               DO i2 = startidx(2), endidx(2)
                  DO i1 = startidx(1), endidx(1)
                     CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                     randnum(i1,i2,i3,i4,i5) = tmprng
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_5d

#if __MAX_RANK > 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 6D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_6d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(6) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i6 = startidx(6), endidx(6)
         DO i5 = startidx(5), endidx(5)
            DO i4 = startidx(4), endidx(4)
               DO i3 = startidx(3), endidx(3)
                  DO i2 = startidx(2), endidx(2)
                     DO i1 = startidx(1), endidx(1)
                        CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                        randnum(i1,i2,i3,i4,i5,i6) = tmprng
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_6d

#if __MAX_RANK > 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 7D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_7d(rng_gaussian, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(7) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i7 = startidx(7), endidx(7)
         DO i6 = startidx(6), endidx(6)
            DO i5 = startidx(5), endidx(5)
               DO i4 = startidx(4), endidx(4)
                  DO i3 = startidx(3), endidx(3)
                     DO i2 = startidx(2), endidx(2)
                        DO i1 = startidx(1), endidx(1)
                           CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                           randnum(i1,i2,i3,i4,i5,i6,i7) = tmprng
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_7d

#if __MAX_RANK > 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 8D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_8d(rng_gaussian, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(8) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i8 = startidx(8), endidx(8)
         DO i7 = startidx(7), endidx(7)
            DO i6 = startidx(6), endidx(6)
               DO i5 = startidx(5), endidx(5)
                  DO i4 = startidx(4), endidx(4)
                     DO i3 = startidx(3), endidx(3)
                        DO i2 = startidx(2), endidx(2)
                           DO i1 = startidx(1), endidx(1)
                              CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                              randnum(i1,i2,i3,i4,i5,i6,i7,i8) = tmprng
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_8d

#if __MAX_RANK > 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 9D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_9d(rng_gaussian, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(9) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i9 = startidx(9), endidx(9)
         DO i8 = startidx(8), endidx(8)
            DO i7 = startidx(7), endidx(7)
               DO i6 = startidx(6), endidx(6)
                  DO i5 = startidx(5), endidx(5)
                     DO i4 = startidx(4), endidx(4)
                        DO i3 = startidx(3), endidx(3)
                           DO i2 = startidx(2), endidx(2)
                              DO i1 = startidx(1), endidx(1)
                                 CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                                 randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9) = tmprng
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_9d

#if __MAX_RANK > 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 10D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_10d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(10) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i10 = startidx(10), endidx(10)
         DO i9 = startidx(9), endidx(9)
            DO i8 = startidx(8), endidx(8)
               DO i7 = startidx(7), endidx(7)
                  DO i6 = startidx(6), endidx(6)
                     DO i5 = startidx(5), endidx(5)
                        DO i4 = startidx(4), endidx(4)
                           DO i3 = startidx(3), endidx(3)
                              DO i2 = startidx(2), endidx(2)
                                 DO i1 = startidx(1), endidx(1)
                                    CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                                    randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10) = tmprng
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_10d

#if __MAX_RANK > 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 11D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_11d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(11) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i11 = startidx(11), endidx(11)
         DO i10 = startidx(10), endidx(10)
            DO i9 = startidx(9), endidx(9)
               DO i8 = startidx(8), endidx(8)
                  DO i7 = startidx(7), endidx(7)
                     DO i6 = startidx(6), endidx(6)
                        DO i5 = startidx(5), endidx(5)
                           DO i4 = startidx(4), endidx(4)
                              DO i3 = startidx(3), endidx(3)
                                 DO i2 = startidx(2), endidx(2)
                                    DO i1 = startidx(1), endidx(1)
                                       CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                                       randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11) = tmprng
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_11d

#if __MAX_RANK > 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 12D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_12d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(12) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i12 = startidx(12), endidx(12)
         DO i11 = startidx(11), endidx(11)
            DO i10 = startidx(10), endidx(10)
               DO i9 = startidx(9), endidx(9)
                  DO i8 = startidx(8), endidx(8)
                     DO i7 = startidx(7), endidx(7)
                        DO i6 = startidx(6), endidx(6)
                           DO i5 = startidx(5), endidx(5)
                              DO i4 = startidx(4), endidx(4)
                                 DO i3 = startidx(3), endidx(3)
                                    DO i2 = startidx(2), endidx(2)
                                       DO i1 = startidx(1), endidx(1)
                                          CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                                          randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) = tmprng
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_12d

#if __MAX_RANK > 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 13D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_13d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(13) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i13 = startidx(13), endidx(13)
         DO i12 = startidx(12), endidx(12)
            DO i11 = startidx(11), endidx(11)
               DO i10 = startidx(10), endidx(10)
                  DO i9 = startidx(9), endidx(9)
                     DO i8 = startidx(8), endidx(8)
                        DO i7 = startidx(7), endidx(7)
                           DO i6 = startidx(6), endidx(6)
                              DO i5 = startidx(5), endidx(5)
                                 DO i4 = startidx(4), endidx(4)
                                    DO i3 = startidx(3), endidx(3)
                                       DO i2 = startidx(2), endidx(2)
                                          DO i1 = startidx(1), endidx(1)
                                             CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                                             randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13) = tmprng
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_13d

#if __MAX_RANK > 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 14D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_14d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(14) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i14 = startidx(14), endidx(14)
         DO i13 = startidx(13), endidx(13)
            DO i12 = startidx(12), endidx(12)
               DO i11 = startidx(11), endidx(11)
                  DO i10 = startidx(10), endidx(10)
                     DO i9 = startidx(9), endidx(9)
                        DO i8 = startidx(8), endidx(8)
                           DO i7 = startidx(7), endidx(7)
                              DO i6 = startidx(6), endidx(6)
                                 DO i5 = startidx(5), endidx(5)
                                    DO i4 = startidx(4), endidx(4)
                                       DO i3 = startidx(3), endidx(3)
                                          DO i2 = startidx(2), endidx(2)
                                             DO i1 = startidx(1), endidx(1)
                                                CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                                                randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) = tmprng
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_14d

#if __MAX_RANK > 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 15D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_sp_15d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=sp) :: tmprng
      INTEGER, DIMENSION(15) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i15 = startidx(15), endidx(15)
         DO i14 = startidx(14), endidx(14)
            DO i13 = startidx(13), endidx(13)
               DO i12 = startidx(12), endidx(12)
                  DO i11 = startidx(11), endidx(11)
                     DO i10 = startidx(10), endidx(10)
                        DO i9 = startidx(9), endidx(9)
                           DO i8 = startidx(8), endidx(8)
                              DO i7 = startidx(7), endidx(7)
                                 DO i6 = startidx(6), endidx(6)
                                    DO i5 = startidx(5), endidx(5)
                                       DO i4 = startidx(4), endidx(4)
                                          DO i3 = startidx(3), endidx(3)
                                             DO i2 = startidx(2), endidx(2)
                                                DO i1 = startidx(1), endidx(1)
                                                   CALL next_random_gaussian_sp(rng_gaussian, tmprng)
                                                   randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) = tmprng
                                                END DO
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_sp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number by Box-Muller
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp(rng_gaussian, randnum)
      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(OUT) :: randnum

      REAL(KIND=dp) :: u1, u2

      IF (rng_gaussian%generated) THEN
         rng_gaussian%generated = .FALSE.
         randnum = rng_gaussian%z2 * rng_gaussian%stddev + rng_gaussian%mean
         RETURN
      END IF

      CALL next_random_uniform_dp(rng_gaussian%rng_uniform, u1)
      CALL next_random_uniform_dp(rng_gaussian%rng_uniform, u2)

      DO WHILE (u1 <= EPSILON(1.0_dp))
         CALL next_random_uniform_dp(rng_gaussian%rng_uniform, u1)
         CALL next_random_uniform_dp(rng_gaussian%rng_uniform, u2)
      END DO

      rng_gaussian%z1 = SQRT(-2.0_dp*LOG(u1)) * COS(twopi*u2)
      rng_gaussian%z2 = SQRT(-2.0_dp*LOG(u1)) * SIN(twopi*u2)
      rng_gaussian%generated = .TRUE.

      randnum = rng_gaussian%z1 * rng_gaussian%stddev + rng_gaussian%mean

      RETURN
   END SUBROUTINE next_random_gaussian_dp

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number by polar Box-Muller
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_polar(rng_gaussian, randnum)
      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(OUT) :: randnum

      REAL(KIND=dp) :: u1, u2
      REAL(KIND=dp) :: q, p

      IF (rng_gaussian%generated) THEN
         rng_gaussian%generated = .FALSE.
         randnum = rng_gaussian%z2 * rng_gaussian%stddev + rng_gaussian%mean
         RETURN
      END IF

      ! need to uniformly random number such that u1 and u2 are from [-1,1]
      ! and u1**2+u2**2 < 1 but larger than 0
      DO 
         CALL next_random_uniform_dp(rng_gaussian%rng_uniform, u1)
         CALL next_random_uniform_dp(rng_gaussian%rng_uniform, u2)

         q = u1**2 + u2**2
         IF (q < EPSILON(1.0_dp)) CYCLE
         IF (q < 1.0_dp) EXIT
      END DO

      p = SQRT(-2.0_dp * LOG(q)/q)

      rng_gaussian%z1 = u1 * p
      rng_gaussian%z2 = u2 * p                                                                          
      rng_gaussian%generated = .TRUE.

      randnum = rng_gaussian%z1 * rng_gaussian%stddev + rng_gaussian%mean

      RETURN
   END SUBROUTINE next_random_gaussian_dp_polar

#if __MAX_RANK > 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 1D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_1d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(1) :: startidx, endidx
      INTEGER :: i1

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)
      
      DO i1 = startidx(1), endidx(1)
         CALL next_random_gaussian_dp(rng_gaussian, tmprng)
         randnum(i1) = tmprng
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_1d

#if __MAX_RANK > 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 2D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_2d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(2) :: startidx, endidx
      INTEGER :: i1, i2

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i2 = startidx(2), endidx(2)
         DO i1 = startidx(1), endidx(1)
            CALL next_random_gaussian_dp(rng_gaussian, tmprng)
            randnum(i1,i2) = tmprng
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_2d

#if __MAX_RANK > 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 3D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_3d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(3) :: startidx, endidx
      INTEGER :: i1, i2, i3

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i3 = startidx(3), endidx(3)
         DO i2 = startidx(2), endidx(2)
            DO i1 = startidx(1), endidx(1)
               CALL next_random_gaussian_dp(rng_gaussian, tmprng)
               randnum(i1,i2,i3) = tmprng
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_3d

#if __MAX_RANK > 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 4D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_4d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(4) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i4 = startidx(4), endidx(4)
         DO i3 = startidx(3), endidx(3)
            DO i2 = startidx(2), endidx(2)
               DO i1 = startidx(1), endidx(1)
                  CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                  randnum(i1,i2,i3,i4) = tmprng
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_4d

#if __MAX_RANK > 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 5D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_5d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(5) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i5 = startidx(5), endidx(5)
         DO i4 = startidx(4), endidx(4)
            DO i3 = startidx(3), endidx(3)
               DO i2 = startidx(2), endidx(2)
                  DO i1 = startidx(1), endidx(1)
                     CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                     randnum(i1,i2,i3,i4,i5) = tmprng
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_5d

#if __MAX_RANK > 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 6D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_6d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(6) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i6 = startidx(6), endidx(6)
         DO i5 = startidx(5), endidx(5)
            DO i4 = startidx(4), endidx(4)
               DO i3 = startidx(3), endidx(3)
                  DO i2 = startidx(2), endidx(2)
                     DO i1 = startidx(1), endidx(1)
                        CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                        randnum(i1,i2,i3,i4,i5,i6) = tmprng
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_6d

#if __MAX_RANK > 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 7D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_7d(rng_gaussian, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(7) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i7 = startidx(7), endidx(7)
         DO i6 = startidx(6), endidx(6)
            DO i5 = startidx(5), endidx(5)
               DO i4 = startidx(4), endidx(4)
                  DO i3 = startidx(3), endidx(3)
                     DO i2 = startidx(2), endidx(2)
                        DO i1 = startidx(1), endidx(1)
                           CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                           randnum(i1,i2,i3,i4,i5,i6,i7) = tmprng
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_7d

#if __MAX_RANK > 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 8D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_8d(rng_gaussian, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(8) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i8 = startidx(8), endidx(8)
         DO i7 = startidx(7), endidx(7)
            DO i6 = startidx(6), endidx(6)
               DO i5 = startidx(5), endidx(5)
                  DO i4 = startidx(4), endidx(4)
                     DO i3 = startidx(3), endidx(3)
                        DO i2 = startidx(2), endidx(2)
                           DO i1 = startidx(1), endidx(1)
                              CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                              randnum(i1,i2,i3,i4,i5,i6,i7,i8) = tmprng
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_8d

#if __MAX_RANK > 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 9D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_9d(rng_gaussian, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(9) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i9 = startidx(9), endidx(9)
         DO i8 = startidx(8), endidx(8)
            DO i7 = startidx(7), endidx(7)
               DO i6 = startidx(6), endidx(6)
                  DO i5 = startidx(5), endidx(5)
                     DO i4 = startidx(4), endidx(4)
                        DO i3 = startidx(3), endidx(3)
                           DO i2 = startidx(2), endidx(2)
                              DO i1 = startidx(1), endidx(1)
                                 CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                                 randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9) = tmprng
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_9d

#if __MAX_RANK > 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 10D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_10d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(10) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i10 = startidx(10), endidx(10)
         DO i9 = startidx(9), endidx(9)
            DO i8 = startidx(8), endidx(8)
               DO i7 = startidx(7), endidx(7)
                  DO i6 = startidx(6), endidx(6)
                     DO i5 = startidx(5), endidx(5)
                        DO i4 = startidx(4), endidx(4)
                           DO i3 = startidx(3), endidx(3)
                              DO i2 = startidx(2), endidx(2)
                                 DO i1 = startidx(1), endidx(1)
                                    CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                                    randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10) = tmprng
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_10d

#if __MAX_RANK > 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 11D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_11d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(11) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i11 = startidx(11), endidx(11)
         DO i10 = startidx(10), endidx(10)
            DO i9 = startidx(9), endidx(9)
               DO i8 = startidx(8), endidx(8)
                  DO i7 = startidx(7), endidx(7)
                     DO i6 = startidx(6), endidx(6)
                        DO i5 = startidx(5), endidx(5)
                           DO i4 = startidx(4), endidx(4)
                              DO i3 = startidx(3), endidx(3)
                                 DO i2 = startidx(2), endidx(2)
                                    DO i1 = startidx(1), endidx(1)
                                       CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                                       randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11) = tmprng
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_11d

#if __MAX_RANK > 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 12D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_12d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(12) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i12 = startidx(12), endidx(12)
         DO i11 = startidx(11), endidx(11)
            DO i10 = startidx(10), endidx(10)
               DO i9 = startidx(9), endidx(9)
                  DO i8 = startidx(8), endidx(8)
                     DO i7 = startidx(7), endidx(7)
                        DO i6 = startidx(6), endidx(6)
                           DO i5 = startidx(5), endidx(5)
                              DO i4 = startidx(4), endidx(4)
                                 DO i3 = startidx(3), endidx(3)
                                    DO i2 = startidx(2), endidx(2)
                                       DO i1 = startidx(1), endidx(1)
                                          CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                                          randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) = tmprng
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_12d

#if __MAX_RANK > 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 13D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_13d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(13) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i13 = startidx(13), endidx(13)
         DO i12 = startidx(12), endidx(12)
            DO i11 = startidx(11), endidx(11)
               DO i10 = startidx(10), endidx(10)
                  DO i9 = startidx(9), endidx(9)
                     DO i8 = startidx(8), endidx(8)
                        DO i7 = startidx(7), endidx(7)
                           DO i6 = startidx(6), endidx(6)
                              DO i5 = startidx(5), endidx(5)
                                 DO i4 = startidx(4), endidx(4)
                                    DO i3 = startidx(3), endidx(3)
                                       DO i2 = startidx(2), endidx(2)
                                          DO i1 = startidx(1), endidx(1)
                                             CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                                             randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13) = tmprng
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_13d

#if __MAX_RANK > 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 14D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_14d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(14) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i14 = startidx(14), endidx(14)
         DO i13 = startidx(13), endidx(13)
            DO i12 = startidx(12), endidx(12)
               DO i11 = startidx(11), endidx(11)
                  DO i10 = startidx(10), endidx(10)
                     DO i9 = startidx(9), endidx(9)
                        DO i8 = startidx(8), endidx(8)
                           DO i7 = startidx(7), endidx(7)
                              DO i6 = startidx(6), endidx(6)
                                 DO i5 = startidx(5), endidx(5)
                                    DO i4 = startidx(4), endidx(4)
                                       DO i3 = startidx(3), endidx(3)
                                          DO i2 = startidx(2), endidx(2)
                                             DO i1 = startidx(1), endidx(1)
                                                CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                                                randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) = tmprng
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_14d

#if __MAX_RANK > 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 15D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_dp_15d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=dp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=dp) :: tmprng
      INTEGER, DIMENSION(15) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i15 = startidx(15), endidx(15)
         DO i14 = startidx(14), endidx(14)
            DO i13 = startidx(13), endidx(13)
               DO i12 = startidx(12), endidx(12)
                  DO i11 = startidx(11), endidx(11)
                     DO i10 = startidx(10), endidx(10)
                        DO i9 = startidx(9), endidx(9)
                           DO i8 = startidx(8), endidx(8)
                              DO i7 = startidx(7), endidx(7)
                                 DO i6 = startidx(6), endidx(6)
                                    DO i5 = startidx(5), endidx(5)
                                       DO i4 = startidx(4), endidx(4)
                                          DO i3 = startidx(3), endidx(3)
                                             DO i2 = startidx(2), endidx(2)
                                                DO i1 = startidx(1), endidx(1)
                                                   CALL next_random_gaussian_dp(rng_gaussian, tmprng)
                                                   randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) = tmprng
                                                END DO
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_dp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#ifdef __HAS_QP
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number by Box-Muller
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp(rng_gaussian, randnum)
      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(OUT) :: randnum

      REAL(KIND=dp) :: tmprnd

      CALL next_random_gaussian_dp(rng_gaussian, tmprnd)

      randnum = REAL(tmprnd, qp)

      RETURN
   END SUBROUTINE next_random_gaussian_qp

#if __MAX_RANK > 0
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 1D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_1d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(1) :: startidx, endidx
      INTEGER :: i1

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)
      
      DO i1 = startidx(1), endidx(1)
         CALL next_random_gaussian_qp(rng_gaussian, tmprng)
         randnum(i1) = tmprng
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_1d

#if __MAX_RANK > 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 2D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_2d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(2) :: startidx, endidx
      INTEGER :: i1, i2

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i2 = startidx(2), endidx(2)
         DO i1 = startidx(1), endidx(1)
            CALL next_random_gaussian_qp(rng_gaussian, tmprng)
            randnum(i1,i2) = tmprng
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_2d

#if __MAX_RANK > 2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 3D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_3d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(3) :: startidx, endidx
      INTEGER :: i1, i2, i3

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i3 = startidx(3), endidx(3)
         DO i2 = startidx(2), endidx(2)
            DO i1 = startidx(1), endidx(1)
               CALL next_random_gaussian_qp(rng_gaussian, tmprng)
               randnum(i1,i2,i3) = tmprng
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_3d

#if __MAX_RANK > 3
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 4D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_4d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(4) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i4 = startidx(4), endidx(4)
         DO i3 = startidx(3), endidx(3)
            DO i2 = startidx(2), endidx(2)
               DO i1 = startidx(1), endidx(1)
                  CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                  randnum(i1,i2,i3,i4) = tmprng
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_4d

#if __MAX_RANK > 4
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 5D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_5d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(5) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i5 = startidx(5), endidx(5)
         DO i4 = startidx(4), endidx(4)
            DO i3 = startidx(3), endidx(3)
               DO i2 = startidx(2), endidx(2)
                  DO i1 = startidx(1), endidx(1)
                     CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                     randnum(i1,i2,i3,i4,i5) = tmprng
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_5d

#if __MAX_RANK > 5
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 6D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_6d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(6) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i6 = startidx(6), endidx(6)
         DO i5 = startidx(5), endidx(5)
            DO i4 = startidx(4), endidx(4)
               DO i3 = startidx(3), endidx(3)
                  DO i2 = startidx(2), endidx(2)
                     DO i1 = startidx(1), endidx(1)
                        CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                        randnum(i1,i2,i3,i4,i5,i6) = tmprng
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_6d

#if __MAX_RANK > 6
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 7D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_7d(rng_gaussian, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(7) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i7 = startidx(7), endidx(7)
         DO i6 = startidx(6), endidx(6)
            DO i5 = startidx(5), endidx(5)
               DO i4 = startidx(4), endidx(4)
                  DO i3 = startidx(3), endidx(3)
                     DO i2 = startidx(2), endidx(2)
                        DO i1 = startidx(1), endidx(1)
                           CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                           randnum(i1,i2,i3,i4,i5,i6,i7) = tmprng
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_7d

#if __MAX_RANK > 7
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 8D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_8d(rng_gaussian, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(8) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i8 = startidx(8), endidx(8)
         DO i7 = startidx(7), endidx(7)
            DO i6 = startidx(6), endidx(6)
               DO i5 = startidx(5), endidx(5)
                  DO i4 = startidx(4), endidx(4)
                     DO i3 = startidx(3), endidx(3)
                        DO i2 = startidx(2), endidx(2)
                           DO i1 = startidx(1), endidx(1)
                              CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                              randnum(i1,i2,i3,i4,i5,i6,i7,i8) = tmprng
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_8d

#if __MAX_RANK > 8
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 9D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_9d(rng_gaussian, randnum)
   
      IMPLICIT NONE
   
      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: randnum
   
      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(9) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i9 = startidx(9), endidx(9)
         DO i8 = startidx(8), endidx(8)
            DO i7 = startidx(7), endidx(7)
               DO i6 = startidx(6), endidx(6)
                  DO i5 = startidx(5), endidx(5)
                     DO i4 = startidx(4), endidx(4)
                        DO i3 = startidx(3), endidx(3)
                           DO i2 = startidx(2), endidx(2)
                              DO i1 = startidx(1), endidx(1)
                                 CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                                 randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9) = tmprng
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_9d

#if __MAX_RANK > 9
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 10D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_10d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(10) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i10 = startidx(10), endidx(10)
         DO i9 = startidx(9), endidx(9)
            DO i8 = startidx(8), endidx(8)
               DO i7 = startidx(7), endidx(7)
                  DO i6 = startidx(6), endidx(6)
                     DO i5 = startidx(5), endidx(5)
                        DO i4 = startidx(4), endidx(4)
                           DO i3 = startidx(3), endidx(3)
                              DO i2 = startidx(2), endidx(2)
                                 DO i1 = startidx(1), endidx(1)
                                    CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                                    randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10) = tmprng
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_10d

#if __MAX_RANK > 10
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 11D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_11d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(11) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i11 = startidx(11), endidx(11)
         DO i10 = startidx(10), endidx(10)
            DO i9 = startidx(9), endidx(9)
               DO i8 = startidx(8), endidx(8)
                  DO i7 = startidx(7), endidx(7)
                     DO i6 = startidx(6), endidx(6)
                        DO i5 = startidx(5), endidx(5)
                           DO i4 = startidx(4), endidx(4)
                              DO i3 = startidx(3), endidx(3)
                                 DO i2 = startidx(2), endidx(2)
                                    DO i1 = startidx(1), endidx(1)
                                       CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                                       randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11) = tmprng
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_11d

#if __MAX_RANK > 11
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 12D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_12d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(12) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i12 = startidx(12), endidx(12)
         DO i11 = startidx(11), endidx(11)
            DO i10 = startidx(10), endidx(10)
               DO i9 = startidx(9), endidx(9)
                  DO i8 = startidx(8), endidx(8)
                     DO i7 = startidx(7), endidx(7)
                        DO i6 = startidx(6), endidx(6)
                           DO i5 = startidx(5), endidx(5)
                              DO i4 = startidx(4), endidx(4)
                                 DO i3 = startidx(3), endidx(3)
                                    DO i2 = startidx(2), endidx(2)
                                       DO i1 = startidx(1), endidx(1)
                                          CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                                          randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) = tmprng
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_12d

#if __MAX_RANK > 12
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 13D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_13d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(13) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i13 = startidx(13), endidx(13)
         DO i12 = startidx(12), endidx(12)
            DO i11 = startidx(11), endidx(11)
               DO i10 = startidx(10), endidx(10)
                  DO i9 = startidx(9), endidx(9)
                     DO i8 = startidx(8), endidx(8)
                        DO i7 = startidx(7), endidx(7)
                           DO i6 = startidx(6), endidx(6)
                              DO i5 = startidx(5), endidx(5)
                                 DO i4 = startidx(4), endidx(4)
                                    DO i3 = startidx(3), endidx(3)
                                       DO i2 = startidx(2), endidx(2)
                                          DO i1 = startidx(1), endidx(1)
                                             CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                                             randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13) = tmprng
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_13d

#if __MAX_RANK > 13
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 14D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_14d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(14) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i14 = startidx(14), endidx(14)
         DO i13 = startidx(13), endidx(13)
            DO i12 = startidx(12), endidx(12)
               DO i11 = startidx(11), endidx(11)
                  DO i10 = startidx(10), endidx(10)
                     DO i9 = startidx(9), endidx(9)
                        DO i8 = startidx(8), endidx(8)
                           DO i7 = startidx(7), endidx(7)
                              DO i6 = startidx(6), endidx(6)
                                 DO i5 = startidx(5), endidx(5)
                                    DO i4 = startidx(4), endidx(4)
                                       DO i3 = startidx(3), endidx(3)
                                          DO i2 = startidx(2), endidx(2)
                                             DO i1 = startidx(1), endidx(1)
                                                CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                                                randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) = tmprng
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_14d

#if __MAX_RANK > 14
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Generates the next gaussian distributed random number for an 15D array
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    randnum: next random number or array
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE next_random_gaussian_qp_15d(rng_gaussian, randnum)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(INOUT) :: rng_gaussian
      REAL(KIND=qp), INTENT(INOUT), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: randnum

      REAL(KIND=qp) :: tmprng
      INTEGER, DIMENSION(15) :: startidx, endidx
      INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15

      startidx = LBOUND(randnum)
      endidx =   UBOUND(randnum)

      DO i15 = startidx(15), endidx(15)
         DO i14 = startidx(14), endidx(14)
            DO i13 = startidx(13), endidx(13)
               DO i12 = startidx(12), endidx(12)
                  DO i11 = startidx(11), endidx(11)
                     DO i10 = startidx(10), endidx(10)
                        DO i9 = startidx(9), endidx(9)
                           DO i8 = startidx(8), endidx(8)
                              DO i7 = startidx(7), endidx(7)
                                 DO i6 = startidx(6), endidx(6)
                                    DO i5 = startidx(5), endidx(5)
                                       DO i4 = startidx(4), endidx(4)
                                          DO i3 = startidx(3), endidx(3)
                                             DO i2 = startidx(2), endidx(2)
                                                DO i1 = startidx(1), endidx(1)
                                                   CALL next_random_gaussian_qp(rng_gaussian, tmprng)
                                                   randnum(i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) = tmprng
                                                END DO
                                             END DO
                                          END DO
                                       END DO
                                    END DO
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN

   END SUBROUTINE next_random_gaussian_qp_15d
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Save rng state to a string of length default_rng_state_length
   !! Variables:
   !!    rng_integer: Random number generator state
   !!    statestr: that contains all the state information
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE get_random_state_integer(rng_integer, statestr)

      IMPLICIT NONE

      TYPE(rng_integer_type), INTENT(IN) :: rng_integer
      CHARACTER(LEN=default_rng_state_length), INTENT(OUT) :: statestr

      WRITE(UNIT=statestr, FMT=*) rng_integer

      RETURN

   END SUBROUTINE get_random_state_integer

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Save rng state to a string of length default_rng_state_length
   !! Variables:
   !!    rng_uniform: Random number generator state
   !!    statestr: that contains all the state information
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE get_random_state_uniform(rng_uniform, statestr)

      IMPLICIT NONE

      TYPE(rng_uniform_type), INTENT(IN) :: rng_uniform
      CHARACTER(LEN=default_rng_state_length), INTENT(OUT) :: statestr
      
      WRITE(UNIT=statestr, FMT=*) rng_uniform

      RETURN

   END SUBROUTINE get_random_state_uniform

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Save rng state to a string of length default_rng_state_length
   !! Variables:
   !!    rng_gaussian: Random number generator state
   !!    statestr: that contains all the state information
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE get_random_state_gaussian(rng_gaussian, statestr)

      IMPLICIT NONE

      TYPE(rng_gaussian_type), INTENT(IN) :: rng_gaussian
      CHARACTER(LEN=default_rng_state_length), INTENT(OUT) :: statestr

      WRITE(UNIT=statestr, FMT=*) rng_gaussian

      RETURN

   END SUBROUTINE get_random_state_gaussian

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Restore rng state from a string of length default_rng_state_length
   !! Variables:
   !!    statestr: that contains all the state information
   !!    rng_integer: Random number generator state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE restore_random_state_integer(statestr, rng_integer)

      IMPLICIT NONE
      
      CHARACTER(LEN=default_rng_state_length), INTENT(IN) :: statestr
      TYPE(rng_integer_type), INTENT(OUT) :: rng_integer

      READ(UNIT=statestr, FMT=*) rng_integer

      RETURN

   END SUBROUTINE restore_random_state_integer

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Restore rng state from a string of length default_rng_state_length
   !! Variables:
   !!    statestr: that contains all the state information
   !!    rng_uniform: Random number generator state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE restore_random_state_uniform(statestr, rng_uniform)

      IMPLICIT NONE
      
      CHARACTER(LEN=default_rng_state_length), INTENT(IN) :: statestr
      TYPE(rng_uniform_type), INTENT(OUT) :: rng_uniform

      READ(UNIT=statestr, FMT=*) rng_uniform

      RETURN

   END SUBROUTINE restore_random_state_uniform

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Restore rng state from a string of length default_rng_state_length
   !! Variables:
   !!    statestr: that contains all the state information
   !!    rng_gaussian: Random number generator state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE restore_random_state_gaussian(statestr, rng_gaussian)

      IMPLICIT NONE
      
      CHARACTER(LEN=default_rng_state_length), INTENT(IN) :: statestr
      TYPE(rng_gaussian_type), INTENT(OUT) :: rng_gaussian

      READ(UNIT=statestr, FMT=*) rng_gaussian

      RETURN

   END SUBROUTINE restore_random_state_gaussian

END MODULE MT_random
