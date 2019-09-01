PROGRAM test_random_uniform_qp

#ifdef __HAS_QP
   USE, INTRINSIC :: ISO_FORTRAN_ENV
   USE kinds, ONLY : qp
   USE MTrandom, ONLY : rng_uniform_type, &
                      init_rng, &
                      next_random

   IMPLICIT NONE

   CHARACTER(LEN=80) :: cmdname
   
   TYPE(rng_uniform_type) :: rng_state

   REAL(KIND=qp), PARAMETER :: eps = 1.0e-14_qp

   REAL(KIND=qp) :: random
#if __MAX_RANK > 0
   REAL(KIND=qp), DIMENSION(:), POINTER :: random1d
#if __MAX_RANK > 1
   REAL(KIND=qp), DIMENSION(:,:), POINTER :: random2d
#if __MAX_RANK > 2
   REAL(KIND=qp), DIMENSION(:,:,:), POINTER :: random3d
#if __MAX_RANK > 3
   REAL(KIND=qp), DIMENSION(:,:,:,:), POINTER :: random4d
#if __MAX_RANK > 4
   REAL(KIND=qp), DIMENSION(:,:,:,:,:), POINTER :: random5d
#if __MAX_RANK > 5
   REAL(KIND=qp), DIMENSION(:,:,:,:,:,:), POINTER :: random6d
#if __MAX_RANK > 6
   REAL(KIND=qp), DIMENSION(:,:,:,:,:,:,:), POINTER :: random7d
#if __MAX_RANK > 7
   REAL(KIND=qp), DIMENSION(:,:,:,:,:,:,:,:), POINTER :: random8d
#if __MAX_RANK > 8
   REAL(KIND=qp), DIMENSION(:,:,:,:,:,:,:,:,:), POINTER :: random9d
#if __MAX_RANK > 9
   REAL(KIND=qp), DIMENSION(:,:,:,:,:,:,:,:,:,:), POINTER :: random10d
#if __MAX_RANK > 10
   REAL(KIND=qp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:), POINTER :: random11d
#if __MAX_RANK > 11
   REAL(KIND=qp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random12d
#if __MAX_RANK > 12
   REAL(KIND=qp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random13d
#if __MAX_RANK > 13
   REAL(KIND=qp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random14d
#if __MAX_RANK > 14
   REAL(KIND=qp), DIMENSION(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), POINTER :: random15d
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
      random = 0.0_qp
      CALL next_random(rng_state, random)
      passed = ABS(random - 108.502454178061676970662905235185061_qp) < eps
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp of module random"
      END IF
   END IF
   
#if __MAX_RANK > 0
   ! 1D array
   IF (passed) THEN
      NULLIFY(random1d)
      ALLOCATE(random1d(3))
      random1d(:) = 0.0_qp
      CALL next_random(rng_state, random1d)
      passed = (.NOT.ANY(random1d == 0.0_qp)) &
               .AND. &
               (ABS(random1d(3) - 109.710544841673005213083507322414128_qp) < eps)
      DEALLOCATE(random1d)
      NULLIFY(random1d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_1d of module random"
      END IF
   END IF

#if __MAX_RANK > 1
   ! 2D array
   IF (passed) THEN
      NULLIFY(random2d)
      ALLOCATE(random2d(3,3))
      random2d(:,:) = 0.0_qp
      CALL next_random(rng_state, random2d)
      passed = (.NOT.ANY(random2d == 0.0_qp)) &
               .AND. &
               (ABS(random2d(3,3) - 123.929182502537391755612905812332843_qp) < eps)
      DEALLOCATE(random2d)
      NULLIFY(random2d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_2d of module random"
      END IF
   END IF

#if __MAX_RANK > 2
   ! 3D array
   IF (passed) THEN
      NULLIFY(random3d)
      ALLOCATE(random3d(3,3,3))
      random3d(:,:,:) = 0.0_qp
      CALL next_random(rng_state, random3d)
      passed = (.NOT.ANY(random3d == 0.0_qp)) &
               .AND. &
               (ABS(random3d(3,3,3) - 69.7213581145898255410473647974784552_qp) < eps)
      DEALLOCATE(random3d)
      NULLIFY(random3d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_3d of module random"
      END IF
   END IF

#if __MAX_RANK > 3
   ! 4D array
   IF (passed) THEN
      NULLIFY(random4d)
      ALLOCATE(random4d(3,3,3,3))
      random4d(:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random4d)
      passed = (.NOT.ANY(random4d == 0.0_qp)) &
               .AND. &
               (ABS(random4d(3,3,3,3) - 89.1841621199686464187730416917794683_qp) < eps)
      DEALLOCATE(random4d)
      NULLIFY(random4d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_4d of module random"
      END IF
   END IF

#if __MAX_RANK > 4
   ! 5D array
   IF (passed) THEN
      NULLIFY(random5d)
      ALLOCATE(random5d(3,3,3,3,3))
      random5d(:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random5d)
      passed = (.NOT.ANY(random5d == 0.0_qp)) &
               .AND. &
               (ABS(random5d(3,3,3,3,3) - 129.695892704359701152571371847133229_qp) < eps)
      DEALLOCATE(random5d)
      NULLIFY(random5d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_5d of module random"
      END IF
   END IF

#if __MAX_RANK > 5
   ! 6D array
   IF (passed) THEN
      NULLIFY(random6d)
      ALLOCATE(random6d(3,3,3,3,3,3))
      random6d(:,:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random6d)
      passed = (.NOT.ANY(random6d == 0.0_qp)) &
               .AND. &
               (ABS(random6d(3,3,3,3,3,3) - 82.8779681550957303268720893880383938_qp) < eps)
      DEALLOCATE(random6d)
      NULLIFY(random6d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_6d of module random"
      END IF
   END IF

#if __MAX_RANK > 6
   ! 7D array
   IF (passed) THEN
      NULLIFY(random7d)
      ALLOCATE(random7d(3,3,3,3,3,3,3))
      random7d(:,:,:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random7d)
      passed = (.NOT.ANY(random7d == 0.0_qp)) &
               .AND. &
               (ABS(random7d(3,3,3,3,3,3,3) - 43.1831875887821331686730706758285013_qp) < eps)
      DEALLOCATE(random7d)
      NULLIFY(random7d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_7d of module random"
      END IF
   END IF

#if __MAX_RANK > 7
   ! 8D array
   IF (passed) THEN
      NULLIFY(random8d)
      ALLOCATE(random8d(3,3,3,3,3,3,3,3))
      random8d(:,:,:,:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random8d)
      passed = (.NOT.ANY(random8d == 0.0_qp)) &
               .AND. &
               (ABS(random8d(3,3,3,3,3,3,3,3) - 85.7011562524215611770413150629145627_qp) < eps)
      DEALLOCATE(random8d)
      NULLIFY(random8d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_8d of module random"
      END IF
   END IF

#if __MAX_RANK > 8
   ! 9D array
   IF (passed) THEN
      NULLIFY(random9d)
      ALLOCATE(random9d(3,3,3,3,3,3,3,3,3))
      random9d(:,:,:,:,:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random9d)
      passed = (.NOT.ANY(random9d == 0.0_qp)) &
               .AND. &
               (ABS(random9d(3,3,3,3,3,3,3,3,3) - 126.325488428725977532783678461522926_qp) < eps)
      DEALLOCATE(random9d)
      NULLIFY(random9d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_9d of module random"
      END IF
   END IF

#if __MAX_RANK > 9
   ! 10D array
   IF (passed) THEN
      NULLIFY(random10d)
      ALLOCATE(random10d(3,3,3,3,3,3,3,3,3,3))
      random10d(:,:,:,:,:,:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random10d)
      passed = (.NOT.ANY(random10d == 0.0_qp)) &
               .AND. &
               (ABS(random10d(3,3,3,3,3,3,3,3,3,3) - 59.4195620137615655908593940994284033_qp) < eps)
      DEALLOCATE(random10d)
      NULLIFY(random10d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_10d of module random"
      END IF
   END IF

#if __MAX_RANK > 10
   ! 11D array
   IF (passed) THEN
      NULLIFY(random11d)
      ALLOCATE(random11d(3,3,3,3,3,3,3,3,3,3,3))
      random11d(:,:,:,:,:,:,:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random11d)
      passed = (.NOT.ANY(random11d == 0.0_qp)) &
               .AND. &
               (ABS(random11d(3,3,3,3,3,3,3,3,3,3,3) - 125.237910932977407312040144857490907_qp) < eps)
      DEALLOCATE(random11d)
      NULLIFY(random11d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_11d of module random"
      END IF
   END IF

#if __MAX_RANK > 11
   ! 12D array
   IF (passed) THEN
      NULLIFY(random12d)
      ALLOCATE(random12d(3,3,3,3,3,3,3,3,3,3,3,3))
      random12d(:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random12d)
      passed = (.NOT.ANY(random12d == 0.0_qp)) &
               .AND. &
               (ABS(random12d(3,3,3,3,3,3,3,3,3,3,3,3) - 86.4732657109800627008040175813866233_qp) < eps)
      DEALLOCATE(random12d)
      NULLIFY(random12d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_12d of module random"
      END IF
   END IF

#if __MAX_RANK > 12
   ! 13D array
   IF (passed) THEN
      NULLIFY(random13d)
      ALLOCATE(random13d(3,3,3,3,3,3,3,3,3,3,3,3,3))
      random13d(:,:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random13d)
      passed = (.NOT.ANY(random13d == 0.0_qp)) &
               .AND. &
               (ABS(random13d(3,3,3,3,3,3,3,3,3,3,3,3,3) - 132.551104405983062822905128140045546_qp) < eps)
      DEALLOCATE(random13d)
      NULLIFY(random13d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_13d of module random"
      END IF
   END IF

#if __MAX_RANK > 13
   ! 14D array
   IF (passed) THEN
      NULLIFY(random14d)
      ALLOCATE(random14d(3,3,3,3,3,3,3,3,3,3,3,3,3,3))
      random14d(:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random14d)
      passed = (.NOT.ANY(random14d == 0.0_qp)) &
               .AND. &
               (ABS(random14d(3,3,3,3,3,3,3,3,3,3,3,3,3,3) - 115.998754632858937285572701643025760_qp) < eps)
      DEALLOCATE(random14d)
      NULLIFY(random14d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_14d of module random"
      END IF
   END IF

#if __MAX_RANK > 14
   ! 15D array
   IF (passed) THEN
      NULLIFY(random15d)
      ALLOCATE(random15d(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3))
      random15d(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) = 0.0_qp
      CALL next_random(rng_state, random15d)
      passed = (.NOT.ANY(random15d == 0.0_qp)) &
               .AND. &
               (ABS(random15d(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3) - 48.2788165254126156210397497373216386_qp) < eps)
      DEALLOCATE(random15d)
      NULLIFY(random15d)
      IF (.NOT.passed) THEN
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Fail in routine next_random_uniform_qp_15d of module random"
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

   WRITE(UNIT=OUTPUT_UNIT, FMT="(A70,2X,A8)") ADJUSTL(cmdname), " [NO QP]"

#endif
END PROGRAM test_random_uniform_qp
