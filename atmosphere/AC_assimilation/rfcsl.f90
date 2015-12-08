! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation

!LL   SUBROUTINE RFCSL
!LL   RECURSIVE FILTER - CONSTANT COEFF - FOR A SCALAR LTD.AREA FIELD
!LL
!LL   THIS VERSION FOR A VARIABLE RESOLUTION GRID AND DLAT AND DLONG
!LL   ARE VECTORS OF THE GRID SPACING
!LL
!*L   ARGUMENTS:
MODULE rfcsl_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE RFCSL(F,N_ROWS,ROW_LENGTH,M_GRID,F_BOUNDARY,           &
     & COSLAT,DLAT,DLONG,SCALE,NPASS                                    &
     &  )

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
!$    USE omp_lib 
      IMPLICIT NONE
!---- INPUT
      INTEGER N_ROWS      ! NUMBER OF ROWS IN GRID
      INTEGER ROW_LENGTH  ! NUMBER OF POINTS IN EACH ROW
      INTEGER M_GRID      ! 0=USE F_BOUNDARY, -1=NO SPECIAL BOUNDARY
      REAL F_BOUNDARY     ! F ASSUMED OUTSIDE LTD.AREA
      REAL COSLAT(N_ROWS) ! COSINE(LATITUDE)
      REAL DLAT(N_ROWS)             ! ROW SPACING (RADIANS)
      REAL DLONG(ROW_LENGTH)        ! POINT SPACING (RADIANS)
      REAL SCALE                    ! FILTER SCALE (RADIANS)
      INTEGER NPASS       ! NUMBER OF PASSES IN EACH DIRECTION
!             NPASS=2 WILL GIVE EQUIVALENT OF SOAR FUNCTION (USUAL)
!             NPASS=LARGE WILL GIVE EQUIVALENT OF GAUSSIAN FUNCTION
!---- IN-OUT
      REAL F(ROW_LENGTH,N_ROWS)     ! FIELD TO BE FILTERED
!
!---- WORK SPACE
      REAL A_NS(N_ROWS)                    ! FILTER COEFFICIENTS N-S-N
      REAL A_WE(ROW_LENGTH,N_ROWS)         ! W-E-W FILTER COEFFICIENTS
      REAL E     ! INTERMEDIATE VALUE IN FORMULA FOR FILTER COEFFICIENT
      REAL Z     ! WORK VARIABLE FOR CALCULATING COEFFS
      INTEGER  JPASS,JROW,JPT     ! LOOP COUNTERS FOR PASS,ROW,PT IN ROW
      INTEGER  JJ, OMP_BLOCK      ! FOR OMP LOOP BLOCKING 

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!*----------------------------------------------------------------------
!           PROGRAMMING NOTES
!     THE ORDER OF NESTING OF DOUBLE LOOPS ON JPT & JROW IS OPTIONAL.
!     SUCH LOOPS ARE INDICATED BY !)) ON THE ENDDO STATEMENT.
!     THE LOOP IN THE DIRECTION OF FILTERING IS RECURSIVE, AND SO CANNOT
!     BE VECTORIZED; HERE IT IS MADE THE OUTER LOOP.
!     A SCALAR COMPUTER MAY BE MORE EFFICIENT WITH A RECURSIVE INNER
!-----------------------------------------------------------------------
!***        BRIEF DESCRIPTION OF METHOD
!     THE INPUT FIELD IS F
!     FIRST PASS IN ONE DIRECTION GIVES G (OVERWRITING ARRAY F):
!      G(I)=A(I)*G(I-1)+B(I)*F(I)  :  I=1,N
!     PASS IN OPPOSITE DIRECTION GIVES H (OVERWRITING ARRAY F):
!      H(I)=A(I)*H(I+1)+B(I)*G(I)  :  I=N,1
!     FOR FILTER ALONG ROW (CALLED E-W-E IN COMMENTS) THE GRID IS
!      REGULAR, SO B=1-A, AND A IS PRECALCULATED.
!     FOR FILTER OF LIMITED AREA THE GRID IS
!      ASSUMED REGULAR, SO B=1-A, AND A IS PRECALCULATED.
!      BECAUSE OF THIS APPROX, THE GRID SHOULD NOT GO NEAR POLE.
!-----------------------------------------------------------------------
!*L   NO EXTERNAL SUBROUTINE CALLS
!-----------------------------------------------------------------------
!L*** PRECALCULATE COEFFICIENTS
!
!---- CALCULATE COEFFS FOR E->W SWEEP (ALSO USED FOR W->E)
!       USING 6.5.13 & 6.5.9
      IF (lhook) CALL dr_hook('RFCSL',zhook_in,zhook_handle)

!$OMP  PARALLEL DEFAULT(NONE) SHARED(row_length, n_rows, coslat, dlong,  &
!$OMP& a_we, a_ns, dlat, npass, scale) PRIVATE(z, e, jrow, jpt)

!$OMP DO SCHEDULE(STATIC)
      DO JPT=1,ROW_LENGTH
        DO JROW=1,N_ROWS
          Z=NPASS*(DLONG(JPT)*COSLAT(JROW))**2*0.25
          E=Z/SCALE          **2
          A_WE(JPT,JROW)= 1.+E-SQRT(E*(E+2.))
        ENDDO ! JROW
      ENDDO   ! JPT
!$OMP END DO

!
!---- CALCULATE COEFFS FOR N->S SWEEP (ALSO USED FOR S->N)
!       USING 6.5.13 & 6.5.9

!$OMP DO SCHEDULE(STATIC)
      DO JROW=1,N_ROWS
        Z=NPASS*DLAT(JROW)**2*0.25
        E=Z/SCALE          **2
        A_NS(JROW)      = 1.+E-SQRT(E*(E+2.))
      ENDDO ! JROW
!$OMP END DO

!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(jrow, jpt, jpass, omp_block, jj)

!-----------------------------------------------------------------------
!L
!L*** START LOOP OVER PASSES
      DO JPASS=1,NPASS     ! ===========================================
!L
!L***  FILTER W->E              -----------------------------------
!L



       IF(M_GRID == 0)THEN ! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
!      DO W BC
!      THE BOUNDARY CONDITIONS DEPEND ON THE NUMBER OF PREVIOUS PASSES
!      IN THE OPPOSITE DIRECTION.
!      IF F_BOUNDARY=0 THE FORMULAE FOR 0,1,2 PREVIOUS PASSES ARE:-
!        G1=(1-A)F1                          (6.5.26)
!        G1=(1/(1+A))F1                      (6.5.27)
!        G1=(1/(1+A)(1-A**2))(F1-A**3*F2)    (6.5.28)
!      IF F_BOUNDARY /= 0 IT MUST BE SUBTRACTED BEFORE USING FORMULAE.
       IF     (JPASS == 1)THEN     !     NO PREVIOUS E->W PASS

!$OMP DO SCHEDULE(STATIC)
         DO      JROW=1,N_ROWS
          F(1,JROW)=F_BOUNDARY+(F(1,JROW)-F_BOUNDARY)                   &
     &    *(1.-A_WE(1,JROW))
         ENDDO ! JROW
!$OMP END DO

       ELSE IF(JPASS == 2)THEN     !    ONE PREVIOUS E->W PASS

!$OMP DO SCHEDULE(STATIC)
         DO      JROW=1,N_ROWS
          F(1,JROW)=F_BOUNDARY+(F(1,JROW)-F_BOUNDARY)                   &
     &                  /(1.+A_WE(1,JROW))
         ENDDO ! JROW
!$OMP END DO
       ELSE IF(JPASS >= 3)THEN     !    TWO PREVIOUS E->W PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2

!$OMP DO SCHEDULE(STATIC)
         DO      JROW=1,N_ROWS
          F(1,JROW)=F_BOUNDARY+(F(1,JROW)-F_BOUNDARY-                   &
     &         A_WE(1,JROW)**3*(F(2,JROW)-F_BOUNDARY))                  &
     &     /((1.+A_WE(1,JROW))*(1.-A_WE(1,JROW)**2))
         ENDDO ! JROW
!$OMP END DO

       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER W->E PASS      EQN 6.5.1

       OMP_BLOCK = N_ROWS
!      ENSURE AS LARGE A STRIDE AS POSSIBLE FOR EACH OMP THREAD 
!$     OMP_BLOCK = CEILING(REAL(N_ROWS)/OMP_GET_NUM_THREADS())

!$OMP DO SCHEDULE(STATIC)
       DO        JJ=1, N_ROWS, OMP_BLOCK
       DO        JPT=2,ROW_LENGTH
       DO        JROW=JJ,MIN(JJ+OMP_BLOCK-1,N_ROWS)
         F(JPT,JROW)=F(JPT,JROW)+                                       &
     &       (F(JPT-1,JROW)-F(JPT,JROW))*A_WE(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
       ENDDO !)) JJ
!$OMP END DO

!L
!L***  FILTER E->W              -----------------------------------
!L
       IF(M_GRID == 0)THEN ! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
!      DO E BC
       IF     (JPASS == 1)THEN     !    ONE PREVIOUS W->E PASS

!$OMP DO SCHEDULE(STATIC)
         DO      JROW=1,N_ROWS
          F(ROW_LENGTH,JROW)=F_BOUNDARY+(F(ROW_LENGTH,JROW)-F_BOUNDARY) &
     &                           /(1.+A_WE(ROW_LENGTH,JROW))
         ENDDO ! JROW
!$OMP END DO
       ELSE IF(JPASS >= 2)THEN     !    TWO PREVIOUS W->E PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
!$OMP DO SCHEDULE(STATIC)
         DO      JROW=1,N_ROWS
          F(ROW_LENGTH,JROW)=F_BOUNDARY+(F(ROW_LENGTH,JROW)-F_BOUNDARY- &
     &     A_WE(ROW_LENGTH,JROW)**3*(F(ROW_LENGTH-1,JROW)-F_BOUNDARY))  &
     &     /((1.+A_WE(ROW_LENGTH,JROW))*(1.-A_WE(ROW_LENGTH,JROW)**2))
         ENDDO ! JROW
!$OMP END DO 

       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER E->W PASS

       OMP_BLOCK = N_ROWS
!      ENSURE AS LARGE A STRIDE AS POSSIBLE FOR EACH OMP THREAD 
!$     OMP_BLOCK = CEILING(REAL(N_ROWS)/OMP_GET_NUM_THREADS())

!$OMP DO SCHEDULE(STATIC)
       DO        JJ=1,N_ROWS,OMP_BLOCK
       DO        JPT=ROW_LENGTH-1,1,-1
       DO        JROW=JJ,MIN(JJ+OMP_BLOCK-1, N_ROWS)
         F(JPT,JROW)=F(JPT,JROW)+                                       &
     &       (F(JPT+1,JROW)-F(JPT,JROW))*A_WE(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
       ENDDO !)) JJ
!$OMP END DO

!L
!L***   FILTER N->S             -----------------------------------
!L
       IF(M_GRID == 0)THEN ! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
!      DO N BC
       IF     (JPASS == 1)THEN     !     NO PREVIOUS S->N PASS

!$OMP DO SCHEDULE(STATIC)
         DO      JPT=1,ROW_LENGTH
          F(JPT,1)=F_BOUNDARY+(F(JPT,1)-F_BOUNDARY)                     &
     &     *(1.-A_NS(1))
         ENDDO ! JPT
!$OMP END DO

       ELSE IF(JPASS == 2)THEN     !    ONE PREVIOUS S->N PASS

!$OMP DO SCHEDULE(STATIC)
         DO      JPT=1,ROW_LENGTH
          F(JPT,1)=F_BOUNDARY+(F(JPT,1)-F_BOUNDARY)                     &
     &                 /(1.+A_NS(1))
         ENDDO ! JPT
!$OMP END DO
       ELSE IF(JPASS >= 3)THEN     !    TWO PREVIOUS S->N PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
!$OMP DO SCHEDULE(STATIC)
         DO      JPT=1,ROW_LENGTH
          F(JPT,1)=F_BOUNDARY+(F(JPT,1)-F_BOUNDARY-                     &
     &     A_NS(1)    **3*(F(JPT,2)-F_BOUNDARY))                        &
     &     /((1.+A_NS(1))*(1.-A_NS(1)**2))
         ENDDO ! JPT
!$OMP END DO
       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!

       OMP_BLOCK = ROW_LENGTH
!      ENSURE AS LARGE A STRIDE AS POSSIBLE FOR EACH OMP THREAD 
!$     OMP_BLOCK = CEILING(REAL(ROW_LENGTH)/OMP_GET_NUM_THREADS())

!L         FILTER N->S PASS
!$OMP DO SCHEDULE(STATIC)
       DO JJ=1, ROW_LENGTH, OMP_BLOCK
       DO        JROW=2,N_ROWS
       DO        JPT=JJ,MIN(JJ+OMP_BLOCK-1,ROW_LENGTH)
         F(JPT,JROW)=F(JPT,JROW)+                                       &
     &       (F(JPT,JROW-1)-F(JPT,JROW))*A_NS(JROW)
       ENDDO !)) JPT
       ENDDO !)) JROW
       ENDDO !)) JJ
!$OMP END DO

!L
!L***   FILTER S->N             -----------------------------------
!L
       IF(M_GRID == 0)THEN ! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
!      DO S BC
       IF     (JPASS == 1)THEN     !    ONE PREVIOUS N->S PASS
!        BC FOR AFTER A SINGLE SWEEP IS USED FOR ALL LATER SWEEPS
!$OMP DO SCHEDULE(STATIC)
         DO      JPT=1,ROW_LENGTH
          F(JPT,N_ROWS)=F_BOUNDARY+(F(JPT,N_ROWS)-F_BOUNDARY)           &
     &                      /(1.+A_NS(N_ROWS)    )
         ENDDO ! JPT
!$OMP END DO
       ELSE IF(JPASS >= 2)THEN     !    TWO PREVIOUS S->N PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
!$OMP DO SCHEDULE(STATIC)
         DO      JPT=1,ROW_LENGTH
          F(JPT,N_ROWS)=F_BOUNDARY+(F(JPT,N_ROWS)-F_BOUNDARY-           &
     &     A_NS(N_ROWS)    **3*(F(JPT,N_ROWS-1)-F_BOUNDARY))            &
     &     /((1.+A_NS(N_ROWS)    )*(1.-A_NS(N_ROWS)    **2))
         ENDDO ! JPT
!$OMP END DO
       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER N->S PASS

       OMP_BLOCK = ROW_LENGTH
!      ENSURE AS LARGE A STRIDE AS POSSIBLE FOR EACH OMP THREAD 
!$     OMP_BLOCK = CEILING(REAL(ROW_LENGTH)/OMP_GET_NUM_THREADS())

!$OMP DO SCHEDULE(STATIC)
       DO        JJ=1,ROW_LENGTH,OMP_BLOCK
       DO        JROW=N_ROWS-1,1,-1
       DO        JPT=JJ,MIN(JJ+OMP_BLOCK-1,ROW_LENGTH)
         F(JPT,JROW)=F(JPT,JROW)+                                       &
     &       (F(JPT,JROW+1)-F(JPT,JROW))*A_NS(JROW)
       ENDDO !)) JPT
       ENDDO !)) JROW
       ENDDO !)) JJ
!$OMP END DO


!
      ENDDO ! JPASS ====================================================

!$OMP END PARALLEL
!
      IF (lhook) CALL dr_hook('RFCSL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE RFCSL

END MODULE rfcsl_mod
