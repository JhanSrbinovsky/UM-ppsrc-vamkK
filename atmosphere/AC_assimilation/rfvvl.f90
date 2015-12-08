! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL    The U.M. version of this deck was made from the PP module PPRFILT
!LL    To keep the two codes equivalent,  changes should be made to
!LL    PPRFILT, and the U.M. deck re-created from it.
!LL    Instructions for doing this are in PPRFILT.
!LL
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation

!LL   SUBROUTINE RFVVL
!LL   RECURSIVE FILTER - VARIABLE COEFF - FOR A VECTOR LTD.AREA FIELD

MODULE rfvvl_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE RFVVL(U,V,N_ROWS,ROW_LENGTH,M_GRID,                    &
     & COSLAT,FIRSTROW,DLAT,DLONG,SCALE,NPASS                           &
     & )

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!---- INPUT
      INTEGER N_ROWS      ! NUMBER OF ROWS IN GRID
      INTEGER ROW_LENGTH  ! NUMBER OF POINTS IN EACH ROW
      INTEGER M_GRID      ! 0=ASSUME 0 OUTSIDE. -1=NO SPECIAL BOUNDARY
      REAL COSLAT(N_ROWS) ! COSINE(LATITUDE)
!      EXCEPT FOR ROWS AT POLE, WHICH HAVE COSLAT(1)=COSLAT(2)/8
!      & COSLAT(NNS)=COSLAT(NNS-1)/8 IF M_GRID=1.
!      TO ALLOW FOR THE NON-ZERO AREA OF GRID TRIANGLE SEGMENTS AT POLE
      REAL FIRSTROW                 ! CO-LATITUDE OF FIRST ROW (RADIANS)
      REAL DLAT                     ! ROW SPACING (RADIANS)
      REAL DLONG                    ! POINT SPACING (RADIANS)
      REAL SCALE(ROW_LENGTH,N_ROWS) ! FILTER SCALE (RADIANS)
      INTEGER NPASS       ! NUMBER OF PASSES IN EACH DIRECTION
!             NPASS=2 WILL GIVE EQUIVALENT OF SOAR FUNCTION (USUAL)
!             NPASS=LARGE WILL GIVE EQUIVALENT OF GAUSSIAN FUNCTION
!---- IN-OUT
      REAL U(ROW_LENGTH,N_ROWS)     ! U FIELD TO BE FILTERED
      REAL V(ROW_LENGTH,N_ROWS)     ! V FIELD TO BE FILTERED
!
!---- WORK SPACE
!
      REAL A_NS(ROW_LENGTH,N_ROWS)         ! FILTER COEFFICIENTS N-S-N
      REAL A_WE(ROW_LENGTH,N_ROWS)         ! W-E-W FILTER COEFFICIENTS
      REAL ZU(N_ROWS)   ! RECURSIVELY CALCULATED U VALUE IN W-E-W SWEEPS
      REAL ZV(N_ROWS)   ! RECURSIVELY CALCULATED V VALUE IN W-E-W SWEEPS
      REAL C_ROT(N_ROWS)  ! COS(ROTATION) : ROTATION=SIN(LAT)*DLONG
      REAL S_ROT(N_ROWS)  ! SIN(ROTATION) : ROTATION=SIN(LAT)*DLONG

      REAL ZZU,ZZV     ! TEMP U,V IN W-E-W FILTER
      REAL E     ! INTERMEDIATE VALUE IN FORMULA FOR FILTER COEFFICIENT
      REAL Z     ! WORK VARIABLE FOR CALCULATING COEFFS
      INTEGER  JPASS,JROW,JPT     ! LOOP COUNTERS FOR PASS,ROW,PT IN ROW

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!*----------------------------------------------------------------------
!           PROGRAMMING NOTES
!     THE ORDER OF NESTING OF DOUBLE LOOPS ON JPT & JROW IS OPTIONAL.
!     SUCH LOOPS ARE INDICATED BY !)) ON THE ENDDO STATEMENT.
!     THE LOOP IN THE DIRECTION OF FILTERING IS RECURSIVE, AND SO CANNOT
!     BE VECTORIZED; HERE IT IS MADE THE OUTER LOOP.
!     A SCALAR COMPUTER MAY BE MORE EFFICIENT WITH A RECURSIVE INNER
!-----------------------------------------------------------------------
!***        BRIEF DESCRIPTION OF METHOD
!     THE INPUT FIELD IS VECTOR F=(U,V)
!     FIRST PASS IN ONE DIRECTION GIVES VECTOR G (OVERWRITING ARRAY F):
!      G(I)=A(I)*G(I-1)+B(I)*F(I)  :  I=1,N
!     PASS IN OPPOSITE DIRECTION GIVES VECTOR H (OVERWRITING ARRAY F):
!      H(I)=A(I)*H(I+1)+B(I)*G(I)  :  I=N,1
!     FOR FILTER ALONG ROW (CALLED E-W-E IN COMMENTS) THE GRID IS
!      REGULAR, SO B=1-A, AND A IS PRECALCULATED.
!     FOR FILTER ALONG COLUMNS (CALLED N-S-N IN COMMENTS) THE GRID IS
!      ASSUMED REGULAR, SO B=1-A, AND A IS PRECALCULATED.
!      BECAUSE OF THIS APPROX, THE GRID SHOULD NOT GO NEAR POLE.
!-----------------------------------------------------------------------
!*L   NO EXTERNAL SUBROUTINE CALLS
!-----------------------------------------------------------------------
!L*** PRECALCULATE COEFFICIENTS
      IF (lhook) CALL dr_hook('RFVVL',zhook_in,zhook_handle)
      Z=FIRSTROW       ! N.B. THIS IS CO-LATITUDE
      DO JROW=1,N_ROWS ! CALCULATE ROTATION COEFFICIENTS
        E=COS(Z)*DLONG
        C_ROT(JROW)=COS(E)
        S_ROT(JROW)=SIN(E)
        Z=Z+DLAT
      ENDDO ! JROW
!
!---- CALCULATE COEFFS FOR E->W SWEEP (ALSO USED FOR W->E)
      DO      JROW=1,N_ROWS
        Z=NPASS*(DLONG*COSLAT(JROW))**2*0.25
        DO      JPT=1,ROW_LENGTH
!      STRICTLY, SCALE SHOULD BE AVERAGE JPT & JPT+(-)1 FOR W->E (E->W)
!       NOT DOING SO MAKES THE CODE SIMPLER.
!       SINCE SCALE IS SMOOTH, THE DIFFERENCE IS NEGLIGABLE.
          E=Z/SCALE(JPT,JROW)**2
          A_WE(JPT,JROW)= 1.+E-SQRT(E*(E+2.))
        ENDDO ! JPT
      ENDDO ! JROW
!
!---- CALCULATE COEFFS FOR N->S SWEEP (ALSO USED FOR S->N)
      Z=NPASS*DLAT**2*0.25
      DO      JROW=1,N_ROWS
        DO      JPT=1,ROW_LENGTH
!      STRICTLY, SCALE SHOULD BE AVERAGE JROW & JROW+(-)1
!       NOT DOING SO MAKES THE CODE SIMPLER.
!       SINCE SCALE IS SMOOTH, THE DIFFERENCE IS NEGLIGABLE.
          E=Z/SCALE(JPT,JROW)**2
          A_NS(JPT,JROW)= 1.+E-SQRT(E*(E+2.))
        ENDDO ! JPT
      ENDDO ! JROW
!-----------------------------------------------------------------------
!L
!L*** START LOOP OVER PASSES
      DO JPASS=1,NPASS     ! ===========================================
!L
!L***  FILTER W->E              -----------------------------------
!L
!
       IF(M_GRID == 0)THEN ! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
!      DO W BC
!      THE BOUNDARY CONDITIONS DEPEND ON THE NUMBER OF PREVIOUS PASSES
!      IN THE OPPOSITE DIRECTION.
!      THE FORMULAE FOR 0,1,2 PREVIOUS PASSES ARE:-
!        G1=(1-A)F1
!        G1=(1/(1+A))F1
!        G1=(1/(1+A)(1-A**2))(F1-A**3*F2)
!      SEE METO11 WORKING PAPER 91 BY R J PURSER (1987) FOR DETAILS
       IF     (JPASS == 1)THEN     !     NO PREVIOUS E->W PASS
         DO      JROW=1,N_ROWS
          U(1,JROW)=U(1,JROW)*(1.-A_WE(1,JROW))
          V(1,JROW)=V(1,JROW)*(1.-A_WE(1,JROW))
          ZU(JROW)=U(1,JROW)
          ZV(JROW)=V(1,JROW)
         ENDDO ! JROW
       ELSE IF(JPASS == 2)THEN     !    ONE PREVIOUS E->W PASS
         DO      JROW=1,N_ROWS
          U(1,JROW)=U(1,JROW)/(1.+A_WE(1,JROW))
          V(1,JROW)=V(1,JROW)/(1.+A_WE(1,JROW))
          ZU(JROW)=U(1,JROW)
          ZV(JROW)=V(1,JROW)
         ENDDO ! JROW
       ELSE IF(JPASS >= 3)THEN     !    TWO PREVIOUS E->W PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
         DO      JROW=1,N_ROWS     !    EQN 6.5.28
          U(1,JROW)=(U(1,JROW)-A_WE(1,JROW)**3*U(2,JROW))               &
     &     /((1.+A_WE(1,JROW))*(1.-A_WE(1,JROW)**2))
          V(1,JROW)=(V(1,JROW)-A_WE(1,JROW)**3*V(2,JROW))               &
     &     /((1.+A_WE(1,JROW))*(1.-A_WE(1,JROW)**2))
          ZU(JROW)=U(1,JROW)
          ZV(JROW)=V(1,JROW)
         ENDDO ! JROW
       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER W->E PASS
       DO        JPT=2,ROW_LENGTH
       DO        JROW=1,N_ROWS
         U(JPT,JROW)=U(JPT,JROW)+A_WE(JPT,JROW)*                        &
     &    (C_ROT(JROW)*ZU(JROW)+S_ROT(JROW)*ZV(JROW)-U(JPT,JROW))
         V(JPT,JROW)=V(JPT,JROW)+A_WE(JPT,JROW)*                        &
     &    (C_ROT(JROW)*ZV(JROW)-S_ROT(JROW)*ZU(JROW)-V(JPT,JROW))
         ZU(JROW)=U(JPT,JROW)
         ZV(JROW)=V(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
!L
!L***  FILTER E->W              -----------------------------------
!L
       IF(M_GRID == 0)THEN ! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
!      DO E BC
       IF     (JPASS == 1)THEN     !    ONE PREVIOUS W->E PASS
         DO      JROW=1,N_ROWS
          U(ROW_LENGTH,JROW)=U(ROW_LENGTH,JROW)                         &
     &                           /(1.+A_WE(ROW_LENGTH,JROW))
          V(ROW_LENGTH,JROW)=V(ROW_LENGTH,JROW)                         &
     &                           /(1.+A_WE(ROW_LENGTH,JROW))
          ZU(JROW)=U(ROW_LENGTH,JROW)
          ZV(JROW)=V(ROW_LENGTH,JROW)
         ENDDO ! JROW
       ELSE IF(JPASS >= 2)THEN     !    TWO PREVIOUS W->E PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
         DO      JROW=1,N_ROWS     !    EQN 6.5.28
          U(ROW_LENGTH,JROW)=(U(ROW_LENGTH,JROW)-                       &
     &     A_WE(ROW_LENGTH,JROW)**3*U(ROW_LENGTH-1,JROW))               &
     &     /((1.+A_WE(ROW_LENGTH,JROW))*(1.-A_WE(ROW_LENGTH,JROW)**2))
          V(ROW_LENGTH,JROW)=(V(ROW_LENGTH,JROW)-                       &
     &     A_WE(ROW_LENGTH,JROW)**3*V(ROW_LENGTH-1,JROW))               &
     &     /((1.+A_WE(ROW_LENGTH,JROW))*(1.-A_WE(ROW_LENGTH,JROW)**2))
          ZU(JROW)=U(ROW_LENGTH,JROW)
          ZV(JROW)=V(ROW_LENGTH,JROW)
         ENDDO ! JROW
       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER E->W PASS
       DO        JPT=ROW_LENGTH-1,1,-1
       DO        JROW=1,N_ROWS
         U(JPT,JROW)=U(JPT,JROW)+A_WE(JPT,JROW)*                        &
     &    (C_ROT(JROW)*ZU(JROW)-S_ROT(JROW)*ZV(JROW)-U(JPT,JROW))
         V(JPT,JROW)=V(JPT,JROW)+A_WE(JPT,JROW)*                        &
     &    (C_ROT(JROW)*ZV(JROW)+S_ROT(JROW)*ZU(JROW)-V(JPT,JROW))
         ZU(JROW)=U(JPT,JROW)
         ZV(JROW)=V(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
!L
!L***   FILTER N->S             -----------------------------------
!L
       IF(M_GRID == 0)THEN ! BCS SPECIFIED TO ALLOW FOR CONSTANT OUTSIDE
!      DO N BC
       IF     (JPASS == 1)THEN     !     NO PREVIOUS S->N PASS
         DO      JPT=1,ROW_LENGTH
          U(JPT,1)=U(JPT,1)*(1.-A_NS(JPT,1))
          V(JPT,1)=V(JPT,1)*(1.-A_NS(JPT,1))
         ENDDO ! JPT
       ELSE IF(JPASS == 2)THEN     !    ONE PREVIOUS S->N PASS
         DO      JPT=1,ROW_LENGTH
          U(JPT,1)=U(JPT,1)/(1.+A_NS(JPT,1))
          V(JPT,1)=V(JPT,1)/(1.+A_NS(JPT,1))
         ENDDO ! JPT
       ELSE IF(JPASS >= 3)THEN     !    TWO PREVIOUS S->N PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
         DO      JPT=1,ROW_LENGTH  !    EQN 6.5.28
          U(JPT,1)=(U(JPT,1)-A_NS(JPT,1)**3*U(JPT,2))                   &
     &     /((1.+A_NS(JPT,1))*(1.-A_NS(JPT,1)**2))
          V(JPT,1)=(V(JPT,1)-A_NS(JPT,1)**3*V(JPT,2))                   &
     &     /((1.+A_NS(JPT,1))*(1.-A_NS(JPT,1)**2))
         ENDDO ! JPT
       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER N->S PASS
       DO        JROW=2,N_ROWS
       DO        JPT=1,ROW_LENGTH
         U(JPT,JROW)=U(JPT,JROW)+                                       &
     &              (U(JPT,JROW-1)-U(JPT,JROW))*A_NS(JPT,JROW)
         V(JPT,JROW)=V(JPT,JROW)+                                       &
     &              (V(JPT,JROW-1)-V(JPT,JROW))*A_NS(JPT,JROW)
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
          U(JPT,N_ROWS)=U(JPT,N_ROWS)/(1.+A_NS(JPT,N_ROWS))
          V(JPT,N_ROWS)=V(JPT,N_ROWS)/(1.+A_NS(JPT,N_ROWS))
         ENDDO ! JPT
       ELSE IF(JPASS >= 2)THEN     !    TWO PREVIOUS S->N PASSES
!                                       BC FOR 2 IS ALSO USED FOR >2
         DO      JPT=1,ROW_LENGTH  !    EQN 6.5.28
          U(JPT,N_ROWS)=                                                &
     &      (U(JPT,N_ROWS)-A_NS(JPT,N_ROWS)**3*U(JPT,N_ROWS-1))         &
     &     /((1.+A_NS(JPT,N_ROWS))*(1.-A_NS(JPT,N_ROWS)**2))
          V(JPT,N_ROWS)=                                                &
     &      (V(JPT,N_ROWS)-A_NS(JPT,N_ROWS)**3*V(JPT,N_ROWS-1))         &
     &     /((1.+A_NS(JPT,N_ROWS))*(1.-A_NS(JPT,N_ROWS)**2))
         ENDDO ! JPT
       ENDIF ! JPASS
!      ELSE  ! M_GRID /= 0  NO SPECIAL TREATMENT OF BOUNDARY
       ENDIF ! M_GRID
!
!L         FILTER N->S PASS
       DO        JROW=N_ROWS-1,1,-1
       DO        JPT=1,ROW_LENGTH
         U(JPT,JROW)=U(JPT,JROW)+                                       &
     &       (U(JPT,JROW+1)-U(JPT,JROW))*A_NS(JPT,JROW)
         V(JPT,JROW)=V(JPT,JROW)+                                       &
     &       (V(JPT,JROW+1)-V(JPT,JROW))*A_NS(JPT,JROW)
       ENDDO !)) JPT
       ENDDO !)) JROW
!
      ENDDO ! JPASS ====================================================
!
      IF (lhook) CALL dr_hook('RFVVL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE RFVVL
END MODULE rfvvl_mod
