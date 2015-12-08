! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   SUBROUTINE RFCSLR
!LL   RECURSIVE FILTER - CONSTANT COEFF - FOR A SCALAR LTD.AREA FIELD
!LL
!LL   THIS VERSION IS FOR A REGULAR GRID AND DLAT AND DLONG ARE 
!LL   SCALARS. 
!LL
!LL   DOCUMENTATION:
!LL   THIS IS EXACTLY LIKE RFVSL (SEE BELOW),
!LL   EXCEPT THAT SCALE IS A SINGLE VALUE, APPLICABLE EVERYWHERE.
!
!*L   ARGUMENTS:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation

MODULE rfcslr_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE RFCSLR(F,N_ROWS,ROW_LENGTH,M_GRID,F_BOUNDARY,          &
     & COSLAT,DLAT,DLONG,SCALE,NPASS                                    &
     &  )
!FPP$ NOCONCUR R

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!---- INPUT
      INTEGER N_ROWS      ! NUMBER OF ROWS IN GRID
      INTEGER ROW_LENGTH  ! NUMBER OF POINTS IN EACH ROW
      INTEGER M_GRID      ! 0=USE F_BOUNDARY, -1=NO SPECIAL BOUNDARY
      REAL F_BOUNDARY     ! F ASSUMED OUTSIDE LTD.AREA
      REAL COSLAT(N_ROWS) ! COSINE(LATITUDE)
      REAL DLAT                     ! ROW SPACING (RADIANS)
      REAL DLONG                    ! POINT SPACING (RADIANS)
      REAL SCALE                    ! FILTER SCALE (RADIANS)
      INTEGER NPASS       ! NUMBER OF PASSES IN EACH DIRECTION
!             NPASS=2 WILL GIVE EQUIVALENT OF SOAR FUNCTION (USUAL)
!             NPASS=LARGE WILL GIVE EQUIVALENT OF GAUSSIAN FUNCTION
!---- IN-OUT
      REAL F(ROW_LENGTH,N_ROWS)     ! FIELD TO BE FILTERED
!
!---- WORK SPACE
      REAL A_NS                            ! FILTER COEFFICIENTS N-S-N
      REAL A_WE(           N_ROWS)         ! W-E-W FILTER COEFFICIENTS
      REAL E     ! INTERMEDIATE VALUE IN FORMULA FOR FILTER COEFFICIENT
      REAL Z     ! WORK VARIABLE FOR CALCULATING COEFFS
      INTEGER  JPASS,JROW,JPT     ! LOOP COUNTERS FOR PASS,ROW,PT IN ROW

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
      IF (lhook) CALL dr_hook('RFCSLR',zhook_in,zhook_handle)
      DO      JROW=1,N_ROWS
        Z=NPASS*(DLONG*COSLAT(JROW))**2*0.25
          E=Z/SCALE          **2
          A_WE(    JROW)= 1.+E-SQRT(E*(E+2.))
      ENDDO ! JROW
!
!---- CALCULATE COEFFS FOR N->S SWEEP (ALSO USED FOR S->N)
!       USING 6.5.13 & 6.5.9
      Z=NPASS*DLAT**2*0.25
          E=Z/SCALE          **2
          A_NS          = 1.+E-SQRT(E*(E+2.))

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
         DO      JROW=1,N_ROWS
          F(1,JROW)=F_BOUNDARY+(F(1,JROW)-F_BOUNDARY)                   &
     &    *(1.-A_WE(  JROW))
         ENDDO ! JROW
       ELSE IF(JPASS == 2)THEN     !    ONE PREVIOUS E->W PASS
         DO      JROW=1,N_ROWS
          F(1,JROW)=F_BOUNDARY+(F(1,JROW)-F_BOUNDARY)                   &
     &                  /(1.+A_WE(  JROW))
         ENDDO ! JROW
       ELSE IF(JPASS >= 3)THEN     !    TWO PREVIOUS E->W PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
         DO      JROW=1,N_ROWS
          F(1,JROW)=F_BOUNDARY+(F(1,JROW)-F_BOUNDARY-                   &
     &         A_WE(  JROW)**3*(F(2,JROW)-F_BOUNDARY))                  &
     &     /((1.+A_WE(  JROW))*(1.-A_WE(  JROW)**2))
         ENDDO ! JROW
       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER W->E PASS      EQN 6.5.1
       DO        JPT=2,ROW_LENGTH
       DO        JROW=1,N_ROWS
         F(JPT,JROW)=F(JPT,JROW)+                                       &
     &       (F(JPT-1,JROW)-F(JPT,JROW))*A_WE(    JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
!L
!L***  FILTER E->W              -----------------------------------
!L
       IF(M_GRID == 0)THEN ! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
!      DO E BC
       IF     (JPASS == 1)THEN     !    ONE PREVIOUS W->E PASS
         DO      JROW=1,N_ROWS
          F(ROW_LENGTH,JROW)=F_BOUNDARY+(F(ROW_LENGTH,JROW)-F_BOUNDARY) &
     &                           /(1.+A_WE(           JROW))
         ENDDO ! JROW
       ELSE IF(JPASS >= 2)THEN     !    TWO PREVIOUS W->E PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
         DO      JROW=1,N_ROWS
          F(ROW_LENGTH,JROW)=F_BOUNDARY+(F(ROW_LENGTH,JROW)-F_BOUNDARY- &
     &     A_WE(           JROW)**3*(F(ROW_LENGTH-1,JROW)-F_BOUNDARY))  &
     &     /((1.+A_WE(           JROW))*(1.-A_WE(           JROW)**2))
         ENDDO ! JROW
       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER E->W PASS
       DO        JPT=ROW_LENGTH-1,1,-1
       DO        JROW=1,N_ROWS
         F(JPT,JROW)=F(JPT,JROW)+                                       &
     &       (F(JPT+1,JROW)-F(JPT,JROW))*A_WE(    JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
!L
!L***   FILTER N->S             -----------------------------------
!L
       IF(M_GRID == 0)THEN ! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
!      DO N BC
       IF     (JPASS == 1)THEN     !     NO PREVIOUS S->N PASS
         DO      JPT=1,ROW_LENGTH
          F(JPT,1)=F_BOUNDARY+(F(JPT,1)-F_BOUNDARY)                     &
     &     *(1.-A_NS       )
         ENDDO ! JPT
       ELSE IF(JPASS == 2)THEN     !    ONE PREVIOUS S->N PASS
         DO      JPT=1,ROW_LENGTH
          F(JPT,1)=F_BOUNDARY+(F(JPT,1)-F_BOUNDARY)                     &
     &                 /(1.+A_NS       )
         ENDDO ! JPT
       ELSE IF(JPASS >= 3)THEN     !    TWO PREVIOUS S->N PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
         DO      JPT=1,ROW_LENGTH
          F(JPT,1)=F_BOUNDARY+(F(JPT,1)-F_BOUNDARY-                     &
     &     A_NS       **3*(F(JPT,2)-F_BOUNDARY))                        &
     &     /((1.+A_NS       )*(1.-A_NS       **2))
         ENDDO ! JPT
       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER N->S PASS
       DO        JROW=2,N_ROWS
       DO        JPT=1,ROW_LENGTH
         F(JPT,JROW)=F(JPT,JROW)+                                       &
     &       (F(JPT,JROW-1)-F(JPT,JROW))*A_NS
       ENDDO !)) JPT
       ENDDO !)) JROW
!L
!L***   FILTER S->N             -----------------------------------
!L
       IF(M_GRID == 0)THEN ! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
!      DO S BC
       IF     (JPASS == 1)THEN     !    ONE PREVIOUS N->S PASS
!        BC FOR AFTER A SINGLE SWEEP IS USED FOR ALL LATER SWEEPS
         DO      JPT=1,ROW_LENGTH
          F(JPT,N_ROWS)=F_BOUNDARY+(F(JPT,N_ROWS)-F_BOUNDARY)           &
     &                      /(1.+A_NS            )
         ENDDO ! JPT
       ELSE IF(JPASS >= 2)THEN     !    TWO PREVIOUS S->N PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
         DO      JPT=1,ROW_LENGTH
          F(JPT,N_ROWS)=F_BOUNDARY+(F(JPT,N_ROWS)-F_BOUNDARY-           &
     &     A_NS            **3*(F(JPT,N_ROWS-1)-F_BOUNDARY))            &
     &     /((1.+A_NS            )*(1.-A_NS            **2))
         ENDDO ! JPT
       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER N->S PASS
       DO        JROW=N_ROWS-1,1,-1
       DO        JPT=1,ROW_LENGTH
         F(JPT,JROW)=F(JPT,JROW)+                                       &
     &       (F(JPT,JROW+1)-F(JPT,JROW))*A_NS
       ENDDO !)) JPT
       ENDDO !)) JROW
!
      ENDDO ! JPASS ====================================================
!
      IF (lhook) CALL dr_hook('RFCSLR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE RFCSLR

END MODULE rfcslr_mod
