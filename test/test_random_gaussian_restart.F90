PROGRAM test_random_gaussian_restart

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   USE kinds, ONLY : dp
   USE MTrandom, ONLY : rng_gaussian_type, &
                      default_rng_state_length, &
                      init_rng, &
                      next_random, &
                      get_random_state, &
                      restore_random_state

   IMPLICIT NONE

   CHARACTER(LEN=80) :: cmdname
   
   TYPE(rng_gaussian_type) :: rng_state_cont
   TYPE(rng_gaussian_type) :: rng_state_restr

   CHARACTER(LEN=default_rng_state_length) :: state_string

   REAL(KIND=dp) :: random1, random2
   INTEGER :: i

   LOGICAL :: passed

   CALL GET_COMMAND_ARGUMENT(0, cmdname)

   CALL init_rng(1337, 3.0_dp, 2.0_dp,rng_state_cont)

   ! advance random state
   DO i = 1, 100000
      CALL next_random(rng_state_cont, random1)
   END DO

   !transfer random state via restarting
   CALL get_random_state(rng_state_cont, state_string)
   CALL restore_random_state(state_string, rng_state_restr)

   ! propagate both states independently from each other
   passed = .TRUE.
   DO i = 1, 100000
      CALL next_random(rng_state_cont, random1)
      CALL next_random(rng_state_restr, random2)
      passed = (random1 == random2).AND.(passed)
   END DO

   IF (passed) THEN
      WRITE(UNIT=OUTPUT_UNIT, FMT="(A70,2X,A8)") ADJUSTL(cmdname), "[PASSED]"
   ELSE
      WRITE(UNIT=OUTPUT_UNIT, FMT="(A70,2X,A8)") ADJUSTL(cmdname), "  [FAIL]"
   END IF
   
END PROGRAM test_random_gaussian_restart
