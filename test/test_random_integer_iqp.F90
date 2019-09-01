PROGRAM test_random_integer_iqp

#ifdef __HAS_IQP
   USE, INTRINSIC :: ISO_FORTRAN_ENV
   USE kinds, ONLY : iqp
   USE MTrandom, ONLY : rng_integer_type, &
                      init_rng, &
                      next_random

   IMPLICIT NONE

   CHARACTER(LEN=80) :: cmdname
   
   TYPE(rng_integer_type) :: rng_state

   INTEGER(KIND=iqp) :: random
#if __MAX_RANK > 0
   INTEGER(KIND=iqp), DIMENSION(:), POINTER :: random1d
#if __MAX_RANK > 1
   INTEGER(KIND=iqp), DIMENSION(:,:), POINTER :: random2d
#if __MAX_RANK > 2
   INTEGER(KIND=iqp), DIMENSION(:,:,:), POINTER :: random3d
#if __MAX_RANK > 3
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:), POINTER :: random4d
#if __MAX_RANK > 4
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:), POINTER :: random5d
#if __MAX_RANK > 5
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:,:), POINTER :: random6d
#if __MAX_RANK > 6
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:,:,:), POINTER :: random7d
#if __MAX_RANK > 7
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: random8d
#if __MAX_RANK > 8
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: random9d
#if __MAX_RANK > 9
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: random10d
#if __MAX_RANK > 10
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: random11d
#if __MAX_RANK > 11
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random12d
#if __MAX_RANK > 12
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random13d
#if __MAX_RANK > 13
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random14d
#if __MAX_RANK > 14
   INTEGER(KIND=iqp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random15d
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
      random = 0_iqp
      CALL next_random(rng_state, random)
      passed = random == -5533546679011654786_iqp
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp of module random"
      END IF
   END IF
   
#if __MAX_RANK > 0
   ! 1D array
   IF (passed) THEN
      NULLIFY(random1d)
      ALLOCATE(random1d(3))
      random1d(:) = 0_iqp
      CALL next_random(rng_state, random1d)
      passed = (.NOT.ANY(random1d == 0_iqp)) &
               .AND. &
               (random1d(3) == -5298964160175064061_iqp)
      DEALLOCATE(random1d)
      NULLIFY(random1d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_1d of module random"
      END IF
   END IF

#if __MAX_RANK > 1
   ! 2D array
   IF (passed) THEN
      NULLIFY(random2d)
      ALLOCATE(random2d(3,3))
      random2d(:,:) = 0_iqp
      CALL next_random(rng_state, random2d)
      passed = (.NOT.ANY(random2d == 0_iqp)) &
               .AND. &
               (random2d(3,3) == -2538042370630080369_iqp)
      DEALLOCATE(random2d)
      NULLIFY(random2d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_2d of module random"
      END IF
   END IF

#if __MAX_RANK > 2
   ! 3D array
   IF (passed) THEN
      NULLIFY(random3d)
      ALLOCATE(random3d(3,3,3))
      random3d(:,:,:) = 0_iqp
      CALL next_random(rng_state, random3d)
      passed = (.NOT.ANY(random3d == 0_iqp)) &
               .AND. &
               (random3d(3,3,3) == 5382829458057791205_iqp)
      DEALLOCATE(random3d)
      NULLIFY(random3d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_3d of module random"
      END IF
   END IF

#if __MAX_RANK > 3
   ! 4D array
   IF (passed) THEN
      NULLIFY(random4d)
      ALLOCATE(random4d(3,3,3,3))
      random4d(:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random4d)
      passed = (.NOT.ANY(random4d == 0_iqp)) &
               .AND. &
               (random4d(3,3,3,3) == 9162043820626129831_iqp)
      DEALLOCATE(random4d)
      NULLIFY(random4d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_4d of module random"
      END IF
   END IF

#if __MAX_RANK > 4
   ! 5D array
   IF (passed) THEN
      NULLIFY(random5d)
      ALLOCATE(random5d(3,3,3,3,3))
      random5d(:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random5d)
      passed = (.NOT.ANY(random5d == 0_iqp)) &
               .AND. &
               (random5d(3,3,3,3,3) == -1418284189153594140_iqp)
      DEALLOCATE(random5d)
      NULLIFY(random5d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_5d of module random"
      END IF
   END IF

#if __MAX_RANK > 5
   ! 6D array
   IF (passed) THEN
      NULLIFY(random6d)
      ALLOCATE(random6d(3,3,3,3,3,3))
      random6d(:,:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random6d)
      passed = (.NOT.ANY(random6d == 0_iqp)) &
               .AND. &
               (random6d(3,3,3,3,3,3) == 7937530703266315534_iqp)
      DEALLOCATE(random6d)
      NULLIFY(random6d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_6d of module random"
      END IF
   END IF

#if __MAX_RANK > 6
   ! 7D array
   IF (passed) THEN
      NULLIFY(random7d)
      ALLOCATE(random7d(3,3,3,3,3,3,3))
      random7d(:,:,:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random7d)
      passed = (.NOT.ANY(random7d == 0_iqp)) &
               .AND. &
               (random7d(3,3,3,3,3,3,3) == 229746933067933294_iqp)
      DEALLOCATE(random7d)
      NULLIFY(random7d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_7d of module random"
      END IF
   END IF

#if __MAX_RANK > 7
   ! 8D array
   IF (passed) THEN
      NULLIFY(random8d)
      ALLOCATE(random8d(3,3,3,3,3,3,3,3))
      random8d(:,:,:,:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random8d)
      passed = (.NOT.ANY(random8d == 0_iqp)) &
               .AND. &
               (random8d(3,3,3,3,3,3,3,3) == 8485726790669606640_iqp)
      DEALLOCATE(random8d)
      NULLIFY(random8d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_8d of module random"
      END IF
   END IF

#if __MAX_RANK > 8
   ! 9D array
   IF (passed) THEN
      NULLIFY(random9d)
      ALLOCATE(random9d(3,3,3,3,3,3,3,3,3))
      random9d(:,:,:,:,:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random9d)
      passed = (.NOT.ANY(random9d == 0_iqp)) &
               .AND. &
               (random9d(3,3,3,3,3,3,3,3,3) == -2072736663864665851_iqp)
      DEALLOCATE(random9d)
      NULLIFY(random9d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_9d of module random"
      END IF
   END IF

#if __MAX_RANK > 9
   ! 10D array
   IF (passed) THEN
      NULLIFY(random10d)
      ALLOCATE(random10d(3,3,3,3,3,3,3,3,3,3))
      random10d(:,:,:,:,:,:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random10d)
      passed = (.NOT.ANY(random10d == 0_iqp)) &
               .AND. &
               (random10d(3,3,3,3,3,3,3,3,3,3) == 3382465287831286048_iqp)
      DEALLOCATE(random10d)
      NULLIFY(random10d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_10d of module random"
      END IF
   END IF

#if __MAX_RANK > 10
   ! 11D array
   IF (passed) THEN
      NULLIFY(random11d)
      ALLOCATE(random11d(3,3,3,3,3,3,3,3,3,3,3))
      random11d(:,:,:,:,:,:,:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random11d)
      passed = (.NOT.ANY(random11d == 0_iqp)) &
               .AND. &
               (random11d(3,3,3,3,3,3,3,3,3,3,3) == -2283918387279399934_iqp)
      DEALLOCATE(random11d)
      NULLIFY(random11d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_11d of module random"
      END IF
   END IF

#if __MAX_RANK > 11
   ! 12D array
   IF (passed) THEN
      NULLIFY(random12d)
      ALLOCATE(random12d(3,3,3,3,3,3,3,3,3,3,3,3))
      random12d(:,:,:,:,:,:,:,:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random12d)
      passed = (.NOT.ANY(random12d == 0_iqp)) &
               .AND. &
               (random12d(3,3,3,3,3,3,3,3,3,3,3,3) == 8635652112552965933_iqp)
      DEALLOCATE(random12d)
      NULLIFY(random12d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_12d of module random"
      END IF
   END IF

#if __MAX_RANK > 12
   ! 13D array
   IF (passed) THEN
      NULLIFY(random13d)
      ALLOCATE(random13d(3,3,3,3,3,3,3,3,3,3,3,3,3))
      random13d(:,:,:,:,:,:,:,:,:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random13d)
      passed = (.NOT.ANY(random13d == 0_iqp)) &
               .AND. &
               (random13d(3,3,3,3,3,3,3,3,3,3,3,3,3) == -863869878247205681_iqp)
      DEALLOCATE(random13d)
      NULLIFY(random13d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_13d of module random"
      END IF
   END IF

#if __MAX_RANK > 13
   ! 14D array
   IF (passed) THEN
      NULLIFY(random14d)
      ALLOCATE(random14d(3,3,3,3,3,3,3,3,3,3,3,3,3,3))
      random14d(:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random14d)
      passed = (.NOT.ANY(random14d == 0_iqp)) &
               .AND. &
               (random14d(3,3,3,3,3,3,3,3,3,3,3,3,3,3) == -4077943142282418122_iqp)
      DEALLOCATE(random14d)
      NULLIFY(random14d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_14d of module random"
      END IF
   END IF

#if __MAX_RANK > 14
   ! 15D array
   IF (passed) THEN
      NULLIFY(random15d)
      ALLOCATE(random15d(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3))
      random15d(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0_iqp
      CALL next_random(rng_state, random15d)
      passed = (.NOT.ANY(random15d == 0_iqp)) &
               .AND. &
               (random15d(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3) == 1219197068737524773_iqp)
      DEALLOCATE(random15d)
      NULLIFY(random15d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_integer_iqp_15d of module random"
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

#else

   USE, INTRINSIC :: ISO_FORTRAN_ENV

   IMPLICIT NONE

   CHARACTER(LEN=80) :: cmdname

   CALL GET_COMMAND_ARGUMENT(0, cmdname)

   WRITE(UNIT=OUTPUT_UNIT, FMT="(A70,2X,A8)") ADJUSTL(cmdname), "[NO IQP]"

#endif
END PROGRAM test_random_integer_iqp
