PROGRAM test_random_uniform_dp

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   USE kinds, ONLY : dp
   USE MTrandom, ONLY : rng_uniform_type, &
                      init_rng, &
                      next_random

   IMPLICIT NONE

   CHARACTER(LEN=80) :: cmdname
   
   TYPE(rng_uniform_type) :: rng_state

   REAL(KIND=dp), PARAMETER :: eps = 1e-14_dp

   REAL(KIND=dp) :: random
#if __MAX_RANK > 0
   REAL(KIND=dp), DIMENSION(:), POINTER :: random1d
#if __MAX_RANK > 1
   REAL(KIND=dp), DIMENSION(:,:), POINTER :: random2d
#if __MAX_RANK > 2
   REAL(KIND=dp), DIMENSION(:,:,:), POINTER :: random3d
#if __MAX_RANK > 3
   REAL(KIND=dp), DIMENSION(:,:,:,:), POINTER :: random4d
#if __MAX_RANK > 4
   REAL(KIND=dp), DIMENSION(:,:,:,:,:), POINTER :: random5d
#if __MAX_RANK > 5
   REAL(KIND=dp), DIMENSION(:,:,:,:,:,:), POINTER :: random6d
#if __MAX_RANK > 6
   REAL(KIND=dp), DIMENSION(:,:,:,:,:,:,:), POINTER :: random7d
#if __MAX_RANK > 7
   REAL(KIND=dp), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: random8d
#if __MAX_RANK > 8
   REAL(KIND=dp), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: random9d
#if __MAX_RANK > 9
   REAL(KIND=dp), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: random10d
#if __MAX_RANK > 10
   REAL(KIND=dp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: random11d
#if __MAX_RANK > 11
   REAL(KIND=dp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random12d
#if __MAX_RANK > 12
   REAL(KIND=dp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random13d
#if __MAX_RANK > 13
   REAL(KIND=dp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random14d
#if __MAX_RANK > 14
   REAL(KIND=dp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random15d
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

   LOGICAL :: passed

   CALL GET_COMMAND_ARGUMENT(0, cmdname)

   CALL init_rng(1337, 42.0, 137.0, rng_state)

   passed = .TRUE.

   ! 0D array
   IF (passed) THEN
      random = 0.0_dp
      CALL next_random(rng_state, random)
      passed = ABS(random - 108.50245417806168_dp) < eps
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp of module random"
      END IF
   END IF
   
#if __MAX_RANK > 0
   ! 1D array
   IF (passed) THEN
      NULLIFY(random1d)
      ALLOCATE(random1d(3))
      random1d(:) = 0.0_dp
      CALL next_random(rng_state, random1d)
      passed = (.NOT.ANY(random1d == 0.0_dp)) &
               .AND. &
               (ABS(random1d(3) - 109.71054484167301_dp) < eps)
      DEALLOCATE(random1d)
      NULLIFY(random1d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_1d of module random"
      END IF
   END IF

#if __MAX_RANK > 1
   ! 2D array
   IF (passed) THEN
      NULLIFY(random2d)
      ALLOCATE(random2d(3,3))
      random2d(:,:) = 0.0_dp
      CALL next_random(rng_state, random2d)
      passed = (.NOT.ANY(random2d == 0.0_dp)) &
               .AND. &
               (ABS(random2d(3,3) - 123.92918250253740_dp) < eps)
      DEALLOCATE(random2d)
      NULLIFY(random2d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_2d of module random"
      END IF
   END IF

#if __MAX_RANK > 2
   ! 3D array
   IF (passed) THEN
      NULLIFY(random3d)
      ALLOCATE(random3d(3,3,3))
      random3d(:,:,:) = 0.0_dp
      CALL next_random(rng_state, random3d)
      passed = (.NOT.ANY(random3d == 0.0_dp)) &
               .AND. &
               (ABS(random3d(3,3,3) - 69.721358114589833_dp) < eps)
      DEALLOCATE(random3d)
      NULLIFY(random3d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_3d of module random"
      END IF
   END IF

#if __MAX_RANK > 3
   ! 4D array
   IF (passed) THEN
      NULLIFY(random4d)
      ALLOCATE(random4d(3,3,3,3))
      random4d(:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random4d)
      passed = (.NOT.ANY(random4d == 0.0_dp)) &
               .AND. &
               (ABS(random4d(3,3,3,3) - 89.184162119968647_dp) < eps)
      DEALLOCATE(random4d)
      NULLIFY(random4d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_4d of module random"
      END IF
   END IF

#if __MAX_RANK > 4
   ! 5D array
   IF (passed) THEN
      NULLIFY(random5d)
      ALLOCATE(random5d(3,3,3,3,3))
      random5d(:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random5d)
      passed = (.NOT.ANY(random5d == 0.0_dp)) &
               .AND. &
               (ABS(random5d(3,3,3,3,3) - 129.69589270435972_dp) < eps)
      DEALLOCATE(random5d)
      NULLIFY(random5d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_5d of module random"
      END IF
   END IF

#if __MAX_RANK > 5
   ! 6D array
   IF (passed) THEN
      NULLIFY(random6d)
      ALLOCATE(random6d(3,3,3,3,3,3))
      random6d(:,:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random6d)
      passed = (.NOT.ANY(random6d == 0.0_dp)) &
               .AND. &
               (ABS(random6d(3,3,3,3,3,3) - 82.877968155095743_dp) < eps)
      DEALLOCATE(random6d)
      NULLIFY(random6d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_6d of module random"
      END IF
   END IF

#if __MAX_RANK > 6
   ! 7D array
   IF (passed) THEN
      NULLIFY(random7d)
      ALLOCATE(random7d(3,3,3,3,3,3,3))
      random7d(:,:,:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random7d)
      passed = (.NOT.ANY(random7d == 0.0_dp)) &
               .AND. &
               (ABS(random7d(3,3,3,3,3,3,3) - 43.183187588782133_dp) < eps)
      DEALLOCATE(random7d)
      NULLIFY(random7d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_7d of module random"
      END IF
   END IF

#if __MAX_RANK > 7
   ! 8D array
   IF (passed) THEN
      NULLIFY(random8d)
      ALLOCATE(random8d(3,3,3,3,3,3,3,3))
      random8d(:,:,:,:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random8d)
      passed = (.NOT.ANY(random8d == 0.0_dp)) &
               .AND. &
               (ABS(random8d(3,3,3,3,3,3,3,3) - 85.701156252421569_dp) < eps)
      DEALLOCATE(random8d)
      NULLIFY(random8d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_8d of module random"
      END IF
   END IF

#if __MAX_RANK > 8
   ! 9D array
   IF (passed) THEN
      NULLIFY(random9d)
      ALLOCATE(random9d(3,3,3,3,3,3,3,3,3))
      random9d(:,:,:,:,:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random9d)
      passed = (.NOT.ANY(random9d == 0.0_dp)) &
               .AND. &
               (ABS(random9d(3,3,3,3,3,3,3,3,3) - 126.32548842872599_dp) < eps)
      DEALLOCATE(random9d)
      NULLIFY(random9d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_9d of module random"
      END IF
   END IF

#if __MAX_RANK > 9
   ! 10D array
   IF (passed) THEN
      NULLIFY(random10d)
      ALLOCATE(random10d(3,3,3,3,3,3,3,3,3,3))
      random10d(:,:,:,:,:,:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random10d)
      passed = (.NOT.ANY(random10d == 0.0_dp)) &
               .AND. &
               (ABS(random10d(3,3,3,3,3,3,3,3,3,3) - 59.419562013761563_dp) < eps)
      DEALLOCATE(random10d)
      NULLIFY(random10d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_10d of module random"
      END IF
   END IF

#if __MAX_RANK > 10
   ! 11D array
   IF (passed) THEN
      NULLIFY(random11d)
      ALLOCATE(random11d(3,3,3,3,3,3,3,3,3,3,3))
      random11d(:,:,:,:,:,:,:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random11d)
      passed = (.NOT.ANY(random11d == 0.0_dp)) &
               .AND. &
               (ABS(random11d(3,3,3,3,3,3,3,3,3,3,3) - 125.23791093297741_dp) < eps)
      DEALLOCATE(random11d)
      NULLIFY(random11d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_11d of module random"
      END IF
   END IF

#if __MAX_RANK > 11
   ! 12D array
   IF (passed) THEN
      NULLIFY(random12d)
      ALLOCATE(random12d(3,3,3,3,3,3,3,3,3,3,3,3))
      random12d(:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random12d)
      passed = (.NOT.ANY(random12d == 0.0_dp)) &
               .AND. &
               (ABS(random12d(3,3,3,3,3,3,3,3,3,3,3,3) - 86.473265710980058_dp) < eps)
      DEALLOCATE(random12d)
      NULLIFY(random12d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_12d of module random"
      END IF
   END IF

#if __MAX_RANK > 12
   ! 13D array
   IF (passed) THEN
      NULLIFY(random13d)
      ALLOCATE(random13d(3,3,3,3,3,3,3,3,3,3,3,3,3))
      random13d(:,:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random13d)
      passed = (.NOT.ANY(random13d == 0.0_dp)) &
               .AND. &
               (ABS(random13d(3,3,3,3,3,3,3,3,3,3,3,3,3) - 132.55110440598307_dp) < eps)
      DEALLOCATE(random13d)
      NULLIFY(random13d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_13d of module random"
      END IF
   END IF

#if __MAX_RANK > 13
   ! 14D array
   IF (passed) THEN
      NULLIFY(random14d)
      ALLOCATE(random14d(3,3,3,3,3,3,3,3,3,3,3,3,3,3))
      random14d(:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random14d)
      passed = (.NOT.ANY(random14d == 0.0_dp)) &
               .AND. &
               (ABS(random14d(3,3,3,3,3,3,3,3,3,3,3,3,3,3) - 115.99875463285895_dp) < eps)
      DEALLOCATE(random14d)
      NULLIFY(random14d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_14d of module random"
      END IF
   END IF

#if __MAX_RANK > 14
   ! 15D array
   IF (passed) THEN
      NULLIFY(random15d)
      ALLOCATE(random15d(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3))
      random15d(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_dp
      CALL next_random(rng_state, random15d)
      passed = (.NOT.ANY(random15d == 0.0_dp)) &
               .AND. &
               (ABS(random15d(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3) - 48.278816525412616_dp) < eps)
      DEALLOCATE(random15d)
      NULLIFY(random15d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_dp_15d of module random"
      END IF
   END IF
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

   IF (passed) THEN
      WRITE(UNIT=OUTPUT_UNIT, FMT="(A70,2X,A8)") ADJUSTL(cmdname), "[PASSED]"
   ELSE
      WRITE(UNIT=OUTPUT_UNIT, FMT="(A70,2X,A8)") ADJUSTL(cmdname), "  [FAIL]"
   END IF
   
END PROGRAM test_random_uniform_dp
