PROGRAM test_random_integer_isp

   USE, INTRINSIC :: ISO_FORTRAN_ENV
   USE kinds, ONLY : isp
   USE MTrandom, ONLY : rng_integer_type, &
                      init_rng, &
                      next_random

   IMPLICIT NONE

   CHARACTER(LEN=80) :: cmdname
   
   TYPE(rng_integer_type) :: rng_state

   INTEGER(KIND=isp) :: random
#if __MAX_RANK > 0
   INTEGER(KIND=isp), DIMENSION(:), POINTER :: random1d
#if __MAX_RANK > 1
   INTEGER(KIND=isp), DIMENSION(:,:), POINTER :: random2d
#if __MAX_RANK > 2
   INTEGER(KIND=isp), DIMENSION(:,:,:), POINTER :: random3d
#if __MAX_RANK > 3
   INTEGER(KIND=isp), DIMENSION(:,:,:,:), POINTER :: random4d
#if __MAX_RANK > 4
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:), POINTER :: random5d
#if __MAX_RANK > 5
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:,:), POINTER :: random6d
#if __MAX_RANK > 6
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:,:,:), POINTER :: random7d
#if __MAX_RANK > 7
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: random8d
#if __MAX_RANK > 8
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: random9d
#if __MAX_RANK > 9
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: random10d
#if __MAX_RANK > 10
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: random11d
#if __MAX_RANK > 11
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random12d
#if __MAX_RANK > 12
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random13d
#if __MAX_RANK > 13
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random14d
#if __MAX_RANK > 14
   INTEGER(KIND=isp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random15d
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

   CALL init_rng(1337, rng_state)

   passed = .TRUE.

   ! 0D array
   IF (passed) THEN
      random = 0_isp
      CALL next_random(rng_state, random)
      passed = random == -1288379235_isp
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp of module random"
      END IF
   END IF
   
#if __MAX_RANK > 0
   ! 1D array
   IF (passed) THEN
      NULLIFY(random1d)
      ALLOCATE(random1d(3))
      random1d(:) = 0_isp
      CALL next_random(rng_state, random1d)
      passed = (.NOT.ANY(random1d == 0_isp)) &
               .AND. &
               (random1d(3) == -1233761237_isp)
      DEALLOCATE(random1d)
      NULLIFY(random1d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_1d of module random"
      END IF
   END IF

#if __MAX_RANK > 1
   ! 2D array
   IF (passed) THEN
      NULLIFY(random2d)
      ALLOCATE(random2d(3,3))
      random2d(:,:) = 0_isp
      CALL next_random(rng_state, random2d)
      passed = (.NOT.ANY(random2d == 0_isp)) &
               .AND. &
               (random2d(3,3) == -590934039_isp)
      DEALLOCATE(random2d)
      NULLIFY(random2d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_2d of module random"
      END IF
   END IF

#if __MAX_RANK > 2
   ! 3D array
   IF (passed) THEN
      NULLIFY(random3d)
      ALLOCATE(random3d(3,3,3))
      random3d(:,:,:) = 0_isp
      CALL next_random(rng_state, random3d)
      passed = (.NOT.ANY(random3d == 0_isp)) &
               .AND. &
               (random3d(3,3,3) == 1253287647_isp)
      DEALLOCATE(random3d)
      NULLIFY(random3d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_3d of module random"
      END IF
   END IF

#if __MAX_RANK > 3
   ! 4D array
   IF (passed) THEN
      NULLIFY(random4d)
      ALLOCATE(random4d(3,3,3,3))
      random4d(:,:,:,:) = 0_isp
      CALL next_random(rng_state, random4d)
      passed = (.NOT.ANY(random4d == 0_isp)) &
               .AND. &
               (random4d(3,3,3,3) == 2133204559_isp)
      DEALLOCATE(random4d)
      NULLIFY(random4d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_4d of module random"
      END IF
   END IF

#if __MAX_RANK > 4
   ! 5D array
   IF (passed) THEN
      NULLIFY(random5d)
      ALLOCATE(random5d(3,3,3,3,3))
      random5d(:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random5d)
      passed = (.NOT.ANY(random5d == 0_isp)) &
               .AND. &
               (random5d(3,3,3,3,3) == -330220021_isp)
      DEALLOCATE(random5d)
      NULLIFY(random5d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_5d of module random"
      END IF
   END IF

#if __MAX_RANK > 5
   ! 6D array
   IF (passed) THEN
      NULLIFY(random6d)
      ALLOCATE(random6d(3,3,3,3,3,3))
      random6d(:,:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random6d)
      passed = (.NOT.ANY(random6d == 0_isp)) &
               .AND. &
               (random6d(3,3,3,3,3,3) == 1848100382_isp)
      DEALLOCATE(random6d)
      NULLIFY(random6d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_6d of module random"
      END IF
   END IF

#if __MAX_RANK > 6
   ! 7D array
   IF (passed) THEN
      NULLIFY(random7d)
      ALLOCATE(random7d(3,3,3,3,3,3,3))
      random7d(:,:,:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random7d)
      passed = (.NOT.ANY(random7d == 0_isp)) &
               .AND. &
               (random7d(3,3,3,3,3,3,3) == 53492126_isp)
      DEALLOCATE(random7d)
      NULLIFY(random7d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_7d of module random"
      END IF
   END IF

#if __MAX_RANK > 7
   ! 8D array
   IF (passed) THEN
      NULLIFY(random8d)
      ALLOCATE(random8d(3,3,3,3,3,3,3,3))
      random8d(:,:,:,:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random8d)
      passed = (.NOT.ANY(random8d == 0_isp)) &
               .AND. &
               (random8d(3,3,3,3,3,3,3,3) == 1975737230_isp)
      DEALLOCATE(random8d)
      NULLIFY(random8d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_8d of module random"
      END IF
   END IF

#if __MAX_RANK > 8
   ! 9D array
   IF (passed) THEN
      NULLIFY(random9d)
      ALLOCATE(random9d(3,3,3,3,3,3,3,3,3))
      random9d(:,:,:,:,:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random9d)
      passed = (.NOT.ANY(random9d == 0_isp)) &
               .AND. &
               (random9d(3,3,3,3,3,3,3,3,3) == -482596612_isp)
      DEALLOCATE(random9d)
      NULLIFY(random9d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_9d of module random"
      END IF
   END IF

#if __MAX_RANK > 9
   ! 10D array
   IF (passed) THEN
      NULLIFY(random10d)
      ALLOCATE(random10d(3,3,3,3,3,3,3,3,3,3))
      random10d(:,:,:,:,:,:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random10d)
      passed = (.NOT.ANY(random10d == 0_isp)) &
               .AND. &
               (random10d(3,3,3,3,3,3,3,3,3,3) == 787541570_isp)
      DEALLOCATE(random10d)
      NULLIFY(random10d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_10d of module random"
      END IF
   END IF

#if __MAX_RANK > 10
   ! 11D array
   IF (passed) THEN
      NULLIFY(random11d)
      ALLOCATE(random11d(3,3,3,3,3,3,3,3,3,3,3))
      random11d(:,:,:,:,:,:,:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random11d)
      passed = (.NOT.ANY(random11d == 0_isp)) &
               .AND. &
               (random11d(3,3,3,3,3,3,3,3,3,3,3) == -531766189_isp)
      DEALLOCATE(random11d)
      NULLIFY(random11d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_11d of module random"
      END IF
   END IF

#if __MAX_RANK > 11
   ! 12D array
   IF (passed) THEN
      NULLIFY(random12d)
      ALLOCATE(random12d(3,3,3,3,3,3,3,3,3,3,3,3))
      random12d(:,:,:,:,:,:,:,:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random12d)
      passed = (.NOT.ANY(random12d == 0_isp)) &
               .AND. &
               (random12d(3,3,3,3,3,3,3,3,3,3,3,3) == 2010644439_isp)
      DEALLOCATE(random12d)
      NULLIFY(random12d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_12d of module random"
      END IF
   END IF

#if __MAX_RANK > 12
   ! 13D array
   IF (passed) THEN
      NULLIFY(random13d)
      ALLOCATE(random13d(3,3,3,3,3,3,3,3,3,3,3,3,3))
      random13d(:,:,:,:,:,:,:,:,:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random13d)
      passed = (.NOT.ANY(random13d == 0_isp)) &
               .AND. &
               (random13d(3,3,3,3,3,3,3,3,3,3,3,3,3) == -201135380_isp)
      DEALLOCATE(random13d)
      NULLIFY(random13d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_13d of module random"
      END IF
   END IF

#if __MAX_RANK > 13
   ! 14D array
   IF (passed) THEN
      NULLIFY(random14d)
      ALLOCATE(random14d(3,3,3,3,3,3,3,3,3,3,3,3,3,3))
      random14d(:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random14d)
      passed = (.NOT.ANY(random14d == 0_isp)) &
               .AND. &
               (random14d(3,3,3,3,3,3,3,3,3,3,3,3,3,3) == -949470127_isp)
      DEALLOCATE(random14d)
      NULLIFY(random14d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_14d of module random"
      END IF
   END IF

#if __MAX_RANK > 14
   ! 15D array
   IF (passed) THEN
      NULLIFY(random15d)
      ALLOCATE(random15d(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3))
      random15d(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0_isp
      CALL next_random(rng_state, random15d)
      passed = (.NOT.ANY(random15d == 0_isp)) &
               .AND. &
               (random15d(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3) == 283866438_isp)
      DEALLOCATE(random15d)
      NULLIFY(random15d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_isp_15d of module random"
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
   
END PROGRAM test_random_integer_isp
