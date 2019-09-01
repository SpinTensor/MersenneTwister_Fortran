PROGRAM test_random_gaussian_sp

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   USE kinds, ONLY : sp
   USE MTrandom, ONLY : rng_gaussian_type, &
                      init_rng, &
                      next_random

   IMPLICIT NONE

   CHARACTER(LEN=80) :: cmdname
   
   TYPE(rng_gaussian_type) :: rng_state

   REAL(KIND=sp), PARAMETER :: eps = 1.0e-7_sp

   REAL(KIND=sp) :: random
#if __MAX_RANK > 0
   REAL(KIND=sp), DIMENSION(:), POINTER :: random1d
#if __MAX_RANK > 1
   REAL(KIND=sp), DIMENSION(:,:), POINTER :: random2d
#if __MAX_RANK > 2
   REAL(KIND=sp), DIMENSION(:,:,:), POINTER :: random3d
#if __MAX_RANK > 3
   REAL(KIND=sp), DIMENSION(:,:,:,:), POINTER :: random4d
#if __MAX_RANK > 4
   REAL(KIND=sp), DIMENSION(:,:,:,:,:), POINTER :: random5d
#if __MAX_RANK > 5
   REAL(KIND=sp), DIMENSION(:,:,:,:,:,:), POINTER :: random6d
#if __MAX_RANK > 6
   REAL(KIND=sp), DIMENSION(:,:,:,:,:,:,:), POINTER :: random7d
#if __MAX_RANK > 7
   REAL(KIND=sp), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: random8d
#if __MAX_RANK > 8
   REAL(KIND=sp), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: random9d
#if __MAX_RANK > 9
   REAL(KIND=sp), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: random10d
#if __MAX_RANK > 10
   REAL(KIND=sp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: random11d
#if __MAX_RANK > 11
   REAL(KIND=sp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random12d
#if __MAX_RANK > 12
   REAL(KIND=sp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random13d
#if __MAX_RANK > 13
   REAL(KIND=sp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random14d
#if __MAX_RANK > 14
   REAL(KIND=sp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random15d
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
      random = 0.0_sp
      CALL next_random(rng_state, random)
      passed = ABS(random + 62.6903648_sp) < eps
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp of module random"
      END IF
   END IF
   
#if __MAX_RANK > 0
   ! 1D array
   IF (passed) THEN
      NULLIFY(random1d)
      ALLOCATE(random1d(3))
      random1d(:) = 0.0_sp
      CALL next_random(rng_state, random1d)
      passed = (.NOT.ANY(random1d == 0.0_sp)) &
               .AND. &
               (ABS(random1d(3) + 73.7698517_sp) < eps)
      DEALLOCATE(random1d)
      NULLIFY(random1d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_1d of module random"
      END IF
   END IF

#if __MAX_RANK > 1
   ! 2D array
   IF (passed) THEN
      NULLIFY(random2d)
      ALLOCATE(random2d(3,3))
      random2d(:,:) = 0.0_sp
      CALL next_random(rng_state, random2d)
      passed = (.NOT.ANY(random2d == 0.0_sp)) &
               .AND. &
               (ABS(random2d(3,3) - 42.1143456_sp) < eps)
      DEALLOCATE(random2d)
      NULLIFY(random2d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_2d of module random"
      END IF
   END IF

#if __MAX_RANK > 2
   ! 3D array
   IF (passed) THEN
      NULLIFY(random3d)
      ALLOCATE(random3d(3,3,3))
      random3d(:,:,:) = 0.0_sp
      CALL next_random(rng_state, random3d)
      passed = (.NOT.ANY(random3d == 0.0_sp)) &
               .AND. &
               (ABS(random3d(3,3,3) - 110.235695_sp) < eps)
      DEALLOCATE(random3d)
      NULLIFY(random3d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_3d of module random"
      END IF
   END IF

#if __MAX_RANK > 3
   ! 4D array
   IF (passed) THEN
      NULLIFY(random4d)
      ALLOCATE(random4d(3,3,3,3))
      random4d(:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random4d)
      passed = (.NOT.ANY(random4d == 0.0_sp)) &
               .AND. &
               (ABS(random4d(3,3,3,3) + 115.012611_sp) < eps)
      DEALLOCATE(random4d)
      NULLIFY(random4d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_4d of module random"
      END IF
   END IF

#if __MAX_RANK > 4
   ! 5D array
   IF (passed) THEN
      NULLIFY(random5d)
      ALLOCATE(random5d(3,3,3,3,3))
      random5d(:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random5d)
      passed = (.NOT.ANY(random5d == 0.0_sp)) &
               .AND. &
               (ABS(random5d(3,3,3,3,3) + 14.4771872_sp) < eps)
      DEALLOCATE(random5d)
      NULLIFY(random5d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_5d of module random"
      END IF
   END IF

#if __MAX_RANK > 5
   ! 6D array
   IF (passed) THEN
      NULLIFY(random6d)
      ALLOCATE(random6d(3,3,3,3,3,3))
      random6d(:,:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random6d)
      passed = (.NOT.ANY(random6d == 0.0_sp)) &
               .AND. &
               (ABS(random6d(3,3,3,3,3,3) + 127.068863_sp) < eps)
      DEALLOCATE(random6d)
      NULLIFY(random6d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_6d of module random"
      END IF
   END IF

#if __MAX_RANK > 6
   ! 7D array
   IF (passed) THEN
      NULLIFY(random7d)
      ALLOCATE(random7d(3,3,3,3,3,3,3))
      random7d(:,:,:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random7d)
      passed = (.NOT.ANY(random7d == 0.0_sp)) &
               .AND. &
               (ABS(random7d(3,3,3,3,3,3,3) - 69.0604858_sp) < eps)
      DEALLOCATE(random7d)
      NULLIFY(random7d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_7d of module random"
      END IF
   END IF

#if __MAX_RANK > 7
   ! 8D array
   IF (passed) THEN
      NULLIFY(random8d)
      ALLOCATE(random8d(3,3,3,3,3,3,3,3))
      random8d(:,:,:,:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random8d)
      passed = (.NOT.ANY(random8d == 0.0_sp)) &
               .AND. &
               (ABS(random8d(3,3,3,3,3,3,3,3) - 212.672592_sp) < eps)
      DEALLOCATE(random8d)
      NULLIFY(random8d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_8d of module random"
      END IF
   END IF

#if __MAX_RANK > 8
   ! 9D array
   IF (passed) THEN
      NULLIFY(random9d)
      ALLOCATE(random9d(3,3,3,3,3,3,3,3,3))
      random9d(:,:,:,:,:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random9d)
      passed = (.NOT.ANY(random9d == 0.0_sp)) &
               .AND. &
               (ABS(random9d(3,3,3,3,3,3,3,3,3) + 138.763245_sp) < eps)
      DEALLOCATE(random9d)
      NULLIFY(random9d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_9d of module random"
      END IF
   END IF

#if __MAX_RANK > 9
   ! 10D array
   IF (passed) THEN
      NULLIFY(random10d)
      ALLOCATE(random10d(3,3,3,3,3,3,3,3,3,3))
      random10d(:,:,:,:,:,:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random10d)
      passed = (.NOT.ANY(random10d == 0.0_sp)) &
               .AND. &
               (ABS(random10d(3,3,3,3,3,3,3,3,3,3) + 206.054062_sp) < eps)
      DEALLOCATE(random10d)
      NULLIFY(random10d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_10d of module random"
      END IF
   END IF

#if __MAX_RANK > 10
   ! 11D array
   IF (passed) THEN
      NULLIFY(random11d)
      ALLOCATE(random11d(3,3,3,3,3,3,3,3,3,3,3))
      random11d(:,:,:,:,:,:,:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random11d)
      passed = (.NOT.ANY(random11d == 0.0_sp)) &
               .AND. &
               (ABS(random11d(3,3,3,3,3,3,3,3,3,3,3) + 87.7085953_sp) < eps)
      DEALLOCATE(random11d)
      NULLIFY(random11d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_11d of module random"
      END IF
   END IF

#if __MAX_RANK > 11
   ! 12D array
   IF (passed) THEN
      NULLIFY(random12d)
      ALLOCATE(random12d(3,3,3,3,3,3,3,3,3,3,3,3))
      random12d(:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random12d)
      passed = (.NOT.ANY(random12d == 0.0_sp)) &
               .AND. &
               (ABS(random12d(3,3,3,3,3,3,3,3,3,3,3,3) + 113.784981_sp) < eps)
      DEALLOCATE(random12d)
      NULLIFY(random12d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_12d of module random"
      END IF
   END IF

#if __MAX_RANK > 12
   ! 13D array
   IF (passed) THEN
      NULLIFY(random13d)
      ALLOCATE(random13d(3,3,3,3,3,3,3,3,3,3,3,3,3))
      random13d(:,:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random13d)
      passed = (.NOT.ANY(random13d == 0.0_sp)) &
               .AND. &
               (ABS(random13d(3,3,3,3,3,3,3,3,3,3,3,3,3) + 21.9586010_sp) < eps)
      DEALLOCATE(random13d)
      NULLIFY(random13d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_13d of module random"
      END IF
   END IF

#if __MAX_RANK > 13
   ! 14D array
   IF (passed) THEN
      NULLIFY(random14d)
      ALLOCATE(random14d(3,3,3,3,3,3,3,3,3,3,3,3,3,3))
      random14d(:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random14d)
      passed = (.NOT.ANY(random14d == 0.0_sp)) &
               .AND. &
               (ABS(random14d(3,3,3,3,3,3,3,3,3,3,3,3,3,3) - 22.9003124_sp) < eps)
      DEALLOCATE(random14d)
      NULLIFY(random14d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_14d of module random"
      END IF
   END IF

#if __MAX_RANK > 14
   ! 15D array
   IF (passed) THEN
      NULLIFY(random15d)
      ALLOCATE(random15d(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3))
      random15d(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_sp
      CALL next_random(rng_state, random15d)
      passed = (.NOT.ANY(random15d == 0.0_sp)) &
               .AND. &
               (ABS(random15d(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3) - 127.097763_sp) < eps)
      DEALLOCATE(random15d)
      NULLIFY(random15d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_gaussian_sp_15d of module random"
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
   
END PROGRAM test_random_gaussian_sp
