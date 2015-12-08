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

!LL   SUBROUTINE RFVVG
!LL   RECURSIVE FILTER - VARIABLE COEFF - FOR A VECTOR GLOBAL FIELD
!LL   DOCUMENTATION:
!
MODULE rfvvg_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE RFVVG(U,V,N_ROWS,ROW_LENGTH,M_GRID,                    &
     & COSLAT,FIRSTROW,DLAT,DLONG,SCALE,NPASS                           &
     & )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!---- INPUT
      INTEGER N_ROWS      ! NUMBER OF ROWS IN GRID
      INTEGER ROW_LENGTH  ! NUMBER OF POINTS IN EACH ROW
      INTEGER M_GRID      ! 1=(PTS AT POLE), 2=(STAG'D FROM POLE)
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
      REAL A(ROW_LENGTH,N_ROWS)         ! FILTER COEFFICIENTS N-S-N
      REAL B(ROW_LENGTH,N_ROWS)         ! FILTER COEFFICIENTS N-S-N
      REAL A_WE(ROW_LENGTH,N_ROWS)         ! W-E-W FILTER COEFFICIENTS
      REAL A_PRODUCT(N_ROWS)     ! PRODUCT OF A_WE ROUND LATITUDE CIRCLE
      REAL ZU(N_ROWS)   ! RECURSIVELY CALCULATED U VALUE IN W-E-W SWEEPS
      REAL ZV(N_ROWS)   ! RECURSIVELY CALCULATED V VALUE IN W-E-W SWEEPS
      REAL C_ROT(N_ROWS)  ! COS(ROTATION) : ROTATION=SIN(LAT)*DLONG
      REAL S_ROT(N_ROWS)  ! SIN(ROTATION) : ROTATION=SIN(LAT)*DLONG
      REAL ONES(ROW_LENGTH) ! RECURSIVE RESULT OF FILTERING ONES   N_S_N

      REAL SCALE_N     ! FILTER SCALE AT N POLE
      REAL SCALE_S     ! FILTER SCALE AT S POLE
      REAL E     ! INTERMEDIATE VALUE IN FORMULA FOR FILTER COEFFICIENT
      REAL Z     ! WORK VARIABLE FOR CALCULATING COEFFS
      REAL AREG     ! A COEFF IF THE GRIDBOXES DID NOT CHANGE SIZE N-S
      REAL U_POLE,V_POLE     ! FILTERED VALUE OF FIELD AT POLE
      REAL ZZU,ZZV     ! TEMP U,V IN W-E-W FILTER
      REAL RW1,RW2     ! RATIOS OF WEIGHTS BETWEEN ADJACENT ROWS
      REAL COTLAT     ! COTAN(GRID LATITUDE)
      REAL RSCALE     ! 1/SCALE FOR EXTRA SWEEP E-W-E
      INTEGER IRF,IRG,IRH,IRL ! LIMITS TO ROWS DONE IN W-E-W SWEEPS
      INTEGER  JPASS,JROW,JPT     ! LOOP COUNTERS FOR PASS,ROW,PT IN ROW

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!*----------------------------------------------------------------------
!           PROGRAMMING NOTES
!     THE ORDER OF NESTING OF DOUBLE LOOPS ON JPT & JROW IS OPTIONAL.
!     SUCH LOOPS ARE INDICATED BY !)) ON THE ENDDO STATEMENT.
!     THE LOOP IN THE DIRECTION OF FILTERING IS RECURSIVE, AND SO CANNOT
!     BE VECTORIZED; HERE IT IS MADE THE OUTER LOOP, AND ZU IS A
!     TEMPORARY VECTOR, SO THE INNER LOOP CAN BE VECTORIZED.
!     A SCALAR COMPUTER MAY BE MORE EFFICIENT WITH A RECURSIVE INNER
!     LOOP, WITH ZU IN A TEMPORARY SCALAR.
!-----------------------------------------------------------------------
!***        BRIEF DESCRIPTION OF METHOD
!     THIS FILTER IS THE SAME AS THAT OF A SCALAR FIELD (SEE RFVSG)
!     EXCEPT FOR A ROTATION COSLAT*DLONG FROM PT TO PT ALONG EACH ROW.
!
!     THE INPUT FIELD IS A VECTOR F=(U,V)
!     FIRST PASS IN ONE DIRECTION GIVES VECTOR G (OVERWRITING ARRAY F):
!      G(I)=A(I)*G(I-1)+B(I)*F(I)  :  I=1,N
!     PASS IN OPPOSITE DIRECTION GIVES VECTOR H (OVERWRITING ARRAY F):
!      H(I)=C(I)*H(I+1)+D(I)*G(I)  :  I=N,1
!     FOR FILTER ALONG ROW (CALLED E-W-E IN COMMENTS) THE GRID IS
!      REGULAR, SO D=B=1-C=1-A, AND A IS PRECALCULATED.
!     FOR FILTER ALONG COLUMNS (CALLED N-S-N IN COMMENTS) THE GRID BOX
!      SIZE CHANGES, SO A B C D ARE CALCULATED TO ALLOW FOR THIS.
!      A B C D ARE ALPHA BETA GAMA DELTA IN DOCUMENTATION
!      C IS STORED IN A, D IS STORED IN B.
!-----------------------------------------------------------------------
!*L   NO EXTERNAL SUBROUTINE CALLS
!-----------------------------------------------------------------------
!L*** PRECALCULATE COEFFICIENTS
      IF (lhook) CALL dr_hook('RFVVG',zhook_in,zhook_handle)
      Z=FIRSTROW       ! N.B. THIS IS CO-LATITUDE
      DO JROW=1,N_ROWS ! CALCULATE ROTATION COEFFICIENTS R(COSLAT DLONG)
        E=COS(Z)*DLONG ! FOR ROTATION MATRIX 6.3.8
        C_ROT(JROW)=COS(E)
        S_ROT(JROW)=SIN(E)
        Z=Z+DLAT
      ENDDO ! JROW

      SCALE_N=0.
      SCALE_S=0.
      DO      JPT=1,ROW_LENGTH
        SCALE_N=SCALE_N+SCALE(JPT,1  )
        SCALE_S=SCALE_S+SCALE(JPT,N_ROWS)
      ENDDO ! JPT
      SCALE_N=SCALE_N/ROW_LENGTH
      SCALE_S=SCALE_S/ROW_LENGTH
!
      IF(M_GRID == 1)THEN ! ROWS 1&N_ROWS ARE POLES: NO W->E FILTER
        IRF=2
        IRL=N_ROWS-1
      ELSE                ! FIRST & LAST ROWS ARE 1/2 G-L FROM POLES
        IRF=1
        IRL=N_ROWS
      ENDIF ! M_GRID
!
!---- CALCULATE COEFFS FOR E->W SWEEP (ALSO USED FOR W->E)
!               EQNS. 6.5.13 & 6.5.9
      DO      JROW=IRF,IRL
        Z=NPASS*(DLONG*COSLAT(JROW))**2*0.25
        DO      JPT=1,ROW_LENGTH
!      STRICTLY, SCALE SHOULD BE AVERAGE JPT & JPT+(-)1 FOR W->E (E->W)
!       NOT DOING SO MAKES THE CODE SIMPLER.
!       SINCE SCALE IS SMOOTH, THE DIFFERENCE IS NEGLIGABLE.
          E=Z/SCALE(JPT,JROW)**2
          A_WE(JPT,JROW)= 1.+E-SQRT(E*(E+2.))
        ENDDO ! JPT
        A_PRODUCT(JROW)=A_WE(1 ,JROW)
      ENDDO ! JROW
!              FOR USE IN 6.5.23
      DO        JPT=2,ROW_LENGTH
      DO        JROW=IRF,IRL
        A_PRODUCT(JROW)=A_PRODUCT(JROW)*A_WE(JPT,JROW)
      ENDDO !)) JROW
      ENDDO !)) JPT
!-----------------------------------------------------------------------
!L
!L*** START LOOP OVER PASSES
      DO JPASS=1,NPASS     ! ===========================================
!L
!L***  FILTER W->E              -----------------------------------
!L
       IF(M_GRID == 1)THEN ! REPLACE ROWS 1 & N_ROWS BY AVERAGE
         U_POLE=0.
         V_POLE=0.
         DO      JPT=1,ROW_LENGTH
           ZZU    =C_ROT(1)*U_POLE     +S_ROT(1)*V_POLE
           ZZV    =C_ROT(1)*V_POLE     -S_ROT(1)*U_POLE
           U_POLE=ZZU   +U(JPT,1  )
           V_POLE=ZZV   +V(JPT,1  )
         ENDDO ! JPT
         U_POLE=U_POLE/ROW_LENGTH
         V_POLE=V_POLE/ROW_LENGTH
         DO      JPT=ROW_LENGTH,1,-1
           U(JPT,1  )=U_POLE
           V(JPT,1  )=V_POLE
           ZZU    =C_ROT(1)*U_POLE     -S_ROT(1)*V_POLE
           V_POLE =C_ROT(1)*V_POLE     +S_ROT(1)*U_POLE
           U_POLE =ZZU
         ENDDO ! JPT

         U_POLE=0.
         V_POLE=0.
         DO      JPT=1,ROW_LENGTH
           ZZU    =C_ROT(N_ROWS)*U_POLE +S_ROT(N_ROWS)*V_POLE
           ZZV    =C_ROT(N_ROWS)*V_POLE -S_ROT(N_ROWS)*U_POLE
           U_POLE=ZZU   +U(JPT,N_ROWS )
           V_POLE=ZZV   +V(JPT,N_ROWS )
         ENDDO ! JPT
         U_POLE=U_POLE/ROW_LENGTH
         V_POLE=V_POLE/ROW_LENGTH
         DO      JPT=ROW_LENGTH,1,-1
           U(JPT,N_ROWS )=U_POLE
           V(JPT,N_ROWS )=V_POLE
           ZZU    =C_ROT(N_ROWS)*U_POLE -S_ROT(N_ROWS)*V_POLE
           V_POLE =C_ROT(N_ROWS)*V_POLE +S_ROT(N_ROWS)*U_POLE
           U_POLE =ZZU
         ENDDO ! JPT
       ENDIF ! M_GRID
!
!      SET UP W BC FOR WRAP-ROUND BY FILTERING ACROSS GRID
       DO   JROW=IRF,IRL
         ZU(JROW)=(1.-A_WE(1 ,JROW))*U(1 ,JROW)
         ZV(JROW)=(1.-A_WE(1 ,JROW))*V(1 ,JROW)
       ENDDO ! JROW
       DO        JPT=2,ROW_LENGTH
       DO        JROW=IRF,IRL
           ZZU    =C_ROT(JROW)*ZU(JROW)+S_ROT(JROW)*ZV(JROW)
           ZZV    =C_ROT(JROW)*ZV(JROW)-S_ROT(JROW)*ZU(JROW)
         ZU(JROW)=U(JPT,JROW)+(ZZU     -U(JPT,JROW))*A_WE(JPT,JROW)
         ZV(JROW)=V(JPT,JROW)+(ZZV     -V(JPT,JROW))*A_WE(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
!      SOLVE IMPLICIT EQ. 6.5.23 FOR WRAP-ROUND OF LATITUDE CIRCLE
       DO      JROW=IRF,IRL
         ZU(JROW)=ZU(JROW)/(1.-A_PRODUCT(JROW))
         ZV(JROW)=ZV(JROW)/(1.-A_PRODUCT(JROW))
       ENDDO ! JROW
!
!L         FILTER W->E PASS   (5.6.1A)
       DO        JPT=1,ROW_LENGTH
       DO        JROW=IRF,IRL
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
!      SET UP E BC FOR WRAP-ROUND BY FILTERING ACROSS GRID
       DO      JROW=IRF,IRL
         ZU(JROW)=(1.-A_WE(ROW_LENGTH,JROW))*U(ROW_LENGTH,JROW)
         ZV(JROW)=(1.-A_WE(ROW_LENGTH,JROW))*V(ROW_LENGTH,JROW)
       ENDDO ! JROW
       DO        JPT=ROW_LENGTH-1,1,-1
       DO        JROW=IRF,IRL
           ZZU    =C_ROT(JROW)*ZU(JROW)-S_ROT(JROW)*ZV(JROW)
           ZZV    =C_ROT(JROW)*ZV(JROW)+S_ROT(JROW)*ZU(JROW)
          ZU(JROW)=U(JPT,JROW)+(ZZU     -U(JPT,JROW))*A_WE(JPT,JROW)
          ZV(JROW)=V(JPT,JROW)+(ZZV     -V(JPT,JROW))*A_WE(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
!      SOLVE IMPLICIT EQ. 6.5.23 FOR WRAP-ROUND OF LATITUDE CIRCLE
       DO      JROW=IRF,IRL
         ZU(JROW)=ZU(JROW)/(1.-A_PRODUCT(JROW))
         ZV(JROW)=ZV(JROW)/(1.-A_PRODUCT(JROW))
       ENDDO ! JROW
!
!          FILTER E->W PASS   (5.6.1B)
       DO        JPT=ROW_LENGTH,1,-1
       DO        JROW=IRF,IRL
         U(JPT,JROW)=U(JPT,JROW)+A_WE(JPT,JROW)*                        &
     &    (C_ROT(JROW)*ZU(JROW)-S_ROT(JROW)*ZV(JROW)-U(JPT,JROW))
         V(JPT,JROW)=V(JPT,JROW)+A_WE(JPT,JROW)*                        &
     &    (C_ROT(JROW)*ZV(JROW)+S_ROT(JROW)*ZU(JROW)-V(JPT,JROW))
         ZU(JROW)=U(JPT,JROW)
         ZV(JROW)=V(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
!
!L     N-S-N & S-N-S SWEEPS ARE ALTERNATED TO EQUALIZE FILTER AT POLES
       IF(MOD(JPASS,2) == 1)THEN ! SWEEP N->S->N
!L***   CALCULATE COEFFICIENTS FOR A N->S SWEEP
!         B FOR S.POLE ROW
        DO JPT=1,ROW_LENGTH
          B(JPT,N_ROWS)= 1.
        ENDDO ! JPT
!L      CALCULATE COEFFS USING A RECURSIVE RELATIONSHIP TO CONSERVE
!       THE (AREA WEIGHTED) OUTPUT FROM EACH GRIDPOINT.
!       A FROM 6.5.22, B FROM 6.5.18
        DO      JROW=N_ROWS,2,-1
          RW1=COSLAT(JROW)/COSLAT(JROW-1)
          RW2=MAX(RW1,1./RW1)
          DO      JPT=1,ROW_LENGTH
           E=NPASS*DLAT**2/(SCALE(JPT,JROW-1)+SCALE(JPT,JROW))**2
           AREG= 1.+E-SQRT(E*(E+2.))
           A(JPT,JROW)=RW2*AREG*(B(JPT,JROW)/(1.-AREG))
           B(JPT,JROW-1)=1./(1.+RW1*A(JPT,JROW)/B(JPT,JROW))
          ENDDO ! JPT
        ENDDO ! JROW
!L
!L***   FILTER N->S             -----------------------------------
!L
!L      AT THE SAME TIME CALCULATE COEFFICIENTS FOR A S->N SWEEP
!L      SPECIAL TREATMENT FOR ROW 1
        DO      JPT=1,ROW_LENGTH
          U(JPT,1)=U(JPT,1)*B(JPT,1)
          V(JPT,1)=V(JPT,1)*B(JPT,1)
          ONES(JPT)=B(JPT,1)
          B(JPT,1)=1.             ! D FOR NEXT (S->N) SWEEP 6.5.19
          A(JPT,1)=1.-ONES(JPT)   ! C FOR NEXT (S->N) SWEEP 6.5.20
        ENDDO ! JPT
!L      LOOP OVER ROWS 2,N_ROWS
        DO        JROW=2,N_ROWS
        DO        JPT=1,ROW_LENGTH
           U(JPT,JROW)=U(JPT,JROW)*B(JPT,JROW)+U(JPT,JROW-1)*A(JPT,JROW)
           V(JPT,JROW)=V(JPT,JROW)*B(JPT,JROW)+V(JPT,JROW-1)*A(JPT,JROW)
           ONES(JPT)=            B(JPT,JROW)+ONES(JPT)  *A(JPT,JROW)
!          B(S->N) IS CALCULATED TO CONSERVE EACH INPUT GRID VALUE
           B(JPT,JROW)=1./ (1.+(COSLAT(JROW-1)*A(JPT,JROW-1))/          &
     &                           (COSLAT(JROW)*B(JPT,JROW-1)))
!          A(S->N) IS CALCULATED SUCH THAT A FIELD OF ONES IS UNCHANGED
           A(JPT,JROW)=1.-ONES(JPT)*B(JPT,JROW)
        ENDDO !)) JPT
        ENDDO !)) JROW
!L
!L***   FILTER S->N             -----------------------------------
!L
!L      SPECIAL TREATMENT FOR ROW N_ROWS
        DO      JPT=1,ROW_LENGTH
          U(JPT,N_ROWS)=U(JPT,N_ROWS)/ONES(JPT)
          V(JPT,N_ROWS)=V(JPT,N_ROWS)/ONES(JPT)
        ENDDO ! JPT
        DO        JROW=N_ROWS-1,1,-1
        DO        JPT=1,ROW_LENGTH
           U(JPT,JROW)=U(JPT,JROW)*B(JPT,JROW)+U(JPT,JROW+1)*A(JPT,JROW)
           V(JPT,JROW)=V(JPT,JROW)*B(JPT,JROW)+V(JPT,JROW+1)*A(JPT,JROW)
        ENDDO !)) JPT
        ENDDO !)) JROW
!
!L     N-S-N & S-N-S SWEEPS ARE ALTERNATED TO EQUALIZE FILTER AT POLES
       ELSE  ! SWEEP S->N->S ON EVEN SWEEPS TO EQUALIZE FILTER OF POLES
!
!L***   CALCULATE COEFFICIENTS FOR A S->N SWEEP
        DO      JPT=1,ROW_LENGTH
          B(JPT,1  )= 1.
        ENDDO ! JPT
!L      CALCULATE COEFFS USING A RECURSIVE RELATIONSHIP TO CONSERVE
!       THE (AREA WEIGHTED) OUTPUT FROM EACH GRIDPOINT.
        DO      JROW=1,N_ROWS-1
          RW1=COSLAT(JROW)/COSLAT(JROW+1)
          RW2=MAX(RW1,1./RW1)
          DO      JPT=1,ROW_LENGTH
           E=NPASS*DLAT**2/(SCALE(JPT,JROW+1)+SCALE(JPT,JROW))**2
           AREG= 1.+E-SQRT(E*(E+2.))
           A(JPT,JROW)=RW2*AREG*(B(JPT,JROW)/(1.-AREG))
           B(JPT,JROW+1)=1./(1.+RW1*A(JPT,JROW)/B(JPT,JROW))
          ENDDO ! JPT
        ENDDO ! JROW
!L
!L***   FILTER S->N             -----------------------------------
!L
!L      AT THE SAME TIME CALCULATE COEFFICIENTS FOR A N->S SWEEP
!L      SPECIAL TREATMENT FOR ROW N_ROWS
        DO      JPT=1,ROW_LENGTH
          U(JPT,N_ROWS)=U(JPT,N_ROWS)*B(JPT,N_ROWS)
          V(JPT,N_ROWS)=V(JPT,N_ROWS)*B(JPT,N_ROWS)
          ONES(JPT)=B(JPT,N_ROWS)
          B(JPT,N_ROWS)=1.
          A(JPT,N_ROWS)=1.-ONES(JPT)
        ENDDO ! JPT
!L      LOOP OVER ROWS N_ROWS-1,1
        DO        JROW=N_ROWS-1,1,-1
        DO        JPT=1,ROW_LENGTH
           U(JPT,JROW)=U(JPT,JROW)*B(JPT,JROW)+U(JPT,JROW+1)*A(JPT,JROW)
           V(JPT,JROW)=V(JPT,JROW)*B(JPT,JROW)+V(JPT,JROW+1)*A(JPT,JROW)
           ONES(JPT)=            B(JPT,JROW)+ONES(JPT)  *A(JPT,JROW)
           B(JPT,JROW)=1./ (1.+(COSLAT(JROW+1)*A(JPT,JROW+1))/          &
     &                           (COSLAT(JROW)*B(JPT,JROW+1)))
           A(JPT,JROW)=1.-ONES(JPT)*B(JPT,JROW)
        ENDDO !)) JPT
        ENDDO !)) JROW
!L
!L***   FILTER N->S             -----------------------------------
!L
!L      SPECIAL TREATMENT FOR ROW 1
        DO      JPT=1,ROW_LENGTH
          U(JPT,1  )=U(JPT,1  )/ONES(JPT)
          V(JPT,1  )=V(JPT,1  )/ONES(JPT)
        ENDDO ! JPT
        DO        JROW=2,N_ROWS
        DO        JPT=1,ROW_LENGTH
           U(JPT,JROW)=U(JPT,JROW)*B(JPT,JROW)+U(JPT,JROW-1)*A(JPT,JROW)
           V(JPT,JROW)=V(JPT,JROW)*B(JPT,JROW)+V(JPT,JROW-1)*A(JPT,JROW)
        ENDDO !)) JPT
        ENDDO !)) JROW
!L     N-S-N & S-N-S SWEEPS ARE ALTERNATED TO EQUALIZE FILTER AT POLES
       ENDIF ! END OF (ODD) N->S->N OR (EVEN) S->N->S SWEEP
!
      ENDDO ! JPASS ====================================================
!
!L*** EXTRA W-E-W PASS NEAR POLES
!L***  FILTER W->E              -----------------------------------
       IF(M_GRID == 1)THEN !REPLACE  1&N_ROWS BY AVERAGE
         U_POLE=0.
         V_POLE=0.
         DO      JPT=1,ROW_LENGTH
           ZZU    =C_ROT(1)*U_POLE     +S_ROT(1)*V_POLE
           ZZV    =C_ROT(1)*V_POLE     -S_ROT(1)*U_POLE
           U_POLE=ZZU   +U(JPT,1  )
           V_POLE=ZZV   +V(JPT,1  )
         ENDDO ! JPT
         U_POLE=U_POLE/ROW_LENGTH
         V_POLE=V_POLE/ROW_LENGTH
         DO      JPT=ROW_LENGTH,1,-1
           U(JPT,1  )=U_POLE
           V(JPT,1  )=V_POLE
           ZZU    =C_ROT(1)*U_POLE     -S_ROT(1)*V_POLE
           V_POLE =C_ROT(1)*V_POLE     +S_ROT(1)*U_POLE
           U_POLE =ZZU
         ENDDO ! JPT

         U_POLE=0.
         V_POLE=0.
         DO      JPT=1,ROW_LENGTH
           ZZU    =C_ROT(N_ROWS)*U_POLE +S_ROT(N_ROWS)*V_POLE
           ZZV    =C_ROT(N_ROWS)*V_POLE -S_ROT(N_ROWS)*U_POLE
           U_POLE=ZZU   +U(JPT,N_ROWS )
           V_POLE=ZZV   +V(JPT,N_ROWS )
         ENDDO ! JPT
         U_POLE=U_POLE/ROW_LENGTH
         V_POLE=V_POLE/ROW_LENGTH
         DO      JPT=ROW_LENGTH,1,-1
           U(JPT,N_ROWS )=U_POLE
           V(JPT,N_ROWS )=V_POLE
           ZZU    =C_ROT(N_ROWS)*U_POLE -S_ROT(N_ROWS)*V_POLE
           V_POLE =C_ROT(N_ROWS)*V_POLE +S_ROT(N_ROWS)*U_POLE
           U_POLE =ZZU
         ENDDO ! JPT
       ENDIF ! M_GRID
       IRG=N_ROWS/2
       IRH=IRG+1
       DO JROW=IRF,IRL
!       CALCULATE COEFF FOR THIS EXTRA SWEEP USING 6.5.24 & 6.5.25
        COTLAT=COSLAT(JROW)/MAX(.000001,1.-COSLAT(JROW)**2)
        IF(JROW <  N_ROWS/2)THEN  ! N.POLE
          RSCALE=MAX(COTLAT/SCALE_N,1.)/SCALE_N
        ELSE ! S.POLE
          RSCALE=MAX(COTLAT/SCALE_S,1.)/SCALE_S
        ENDIF ! N.OR S.POLE
!       A_WE IS ASSUMED CONSTANT ALONE ROW FOR THIS EXTRA SWEEP
!       SO ONLY THE FIRST POINT IN ROW IS USED.
        A_WE(1,JROW)=MAX(0.,1.-2.000*DLONG*COSLAT(JROW)*RSCALE)
        IF(A_WE(1,JROW) >  0.)THEN ! NEAR TO POLE
          A_PRODUCT(JROW)=A_WE(1,JROW)**ROW_LENGTH
!         SET UP W BC FOR WRAP-ROUND BY FILTERING ACROSS GRID
          ZU(JROW)=(1.-A_WE(1,JROW))*U(1 ,JROW)
          ZV(JROW)=(1.-A_WE(1,JROW))*V(1 ,JROW)
        ELSE ! DETERMINE LIMITS WHERE FILTER IS DOING NOTHING
          IRG=MIN(IRG,JROW-1)
          IRH=MAX(IRH,JROW+1)
        ENDIF ! NEAR TO POLE
       ENDDO ! JROW
       DO        JPT=2,ROW_LENGTH
       DO        JROW=IRF,IRG
           ZZU    =C_ROT(JROW)*ZU(JROW)+S_ROT(JROW)*ZV(JROW)
           ZZV    =C_ROT(JROW)*ZV(JROW)-S_ROT(JROW)*ZU(JROW)
         ZU(JROW)=U(JPT,JROW)+(ZZU     -U(JPT,JROW))*A_WE(1  ,JROW)
         ZV(JROW)=V(JPT,JROW)+(ZZV     -V(JPT,JROW))*A_WE(1  ,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
       DO        JPT=2,ROW_LENGTH
       DO        JROW=IRH,IRL
           ZZU    =C_ROT(JROW)*ZU(JROW)+S_ROT(JROW)*ZV(JROW)
           ZZV    =C_ROT(JROW)*ZV(JROW)-S_ROT(JROW)*ZU(JROW)
         ZU(JROW)=U(JPT,JROW)+(ZZU     -U(JPT,JROW))*A_WE(1  ,JROW)
         ZV(JROW)=V(JPT,JROW)+(ZZV     -V(JPT,JROW))*A_WE(1  ,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
!      SOLVE IMPLICIT EQUATION FOR WRAP-ROUND FILTER OF LATITUDE CIRCLE
       DO      JROW=IRF,IRG
         ZU(JROW)=ZU(JROW)/(1.-A_PRODUCT(JROW))
         ZV(JROW)=ZV(JROW)/(1.-A_PRODUCT(JROW))
       ENDDO ! JROW
       DO      JROW=IRH,IRL
         ZU(JROW)=ZU(JROW)/(1.-A_PRODUCT(JROW))
         ZV(JROW)=ZV(JROW)/(1.-A_PRODUCT(JROW))
       ENDDO ! JROW
!L         FILTER W->E PASS
       DO        JPT=1,ROW_LENGTH
       DO        JROW=IRF,IRG
         U(JPT,JROW)=U(JPT,JROW)+A_WE(1  ,JROW)*                        &
     &    (C_ROT(JROW)*ZU(JROW)+S_ROT(JROW)*ZV(JROW)-U(JPT,JROW))
         V(JPT,JROW)=V(JPT,JROW)+A_WE(1  ,JROW)*                        &
     &    (C_ROT(JROW)*ZV(JROW)-S_ROT(JROW)*ZU(JROW)-V(JPT,JROW))
         ZU(JROW)=U(JPT,JROW)
         ZV(JROW)=V(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
       DO        JPT=1,ROW_LENGTH
       DO        JROW=IRH,IRL
         U(JPT,JROW)=U(JPT,JROW)+A_WE(1  ,JROW)*                        &
     &    (C_ROT(JROW)*ZU(JROW)+S_ROT(JROW)*ZV(JROW)-U(JPT,JROW))
         V(JPT,JROW)=V(JPT,JROW)+A_WE(1  ,JROW)*                        &
     &    (C_ROT(JROW)*ZV(JROW)-S_ROT(JROW)*ZU(JROW)-V(JPT,JROW))
         ZU(JROW)=U(JPT,JROW)
         ZV(JROW)=V(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
!L***  FILTER E->W              -----------------------------------
!      SET UP E BC FOR WRAP-ROUND BY FILTERING ACROSS GRID
       DO      JROW=IRF,IRG
         ZU(JROW)=(1.-A_WE(1,JROW))*U(ROW_LENGTH,JROW)
         ZV(JROW)=(1.-A_WE(1,JROW))*V(ROW_LENGTH,JROW)
       ENDDO ! JROW
       DO      JROW=IRH,IRL
         ZU(JROW)=(1.-A_WE(1,JROW))*U(ROW_LENGTH,JROW)
         ZV(JROW)=(1.-A_WE(1,JROW))*V(ROW_LENGTH,JROW)
       ENDDO ! JROW
       DO        JPT=ROW_LENGTH-1,1,-1
       DO        JROW=IRF,IRG
           ZZU    =C_ROT(JROW)*ZU(JROW)-S_ROT(JROW)*ZV(JROW)
           ZZV    =C_ROT(JROW)*ZV(JROW)+S_ROT(JROW)*ZU(JROW)
          ZU(JROW)=U(JPT,JROW)+(ZZU     -U(JPT,JROW))*A_WE(1  ,JROW)
          ZV(JROW)=V(JPT,JROW)+(ZZV     -V(JPT,JROW))*A_WE(1  ,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
       DO        JPT=ROW_LENGTH-1,1,-1
       DO        JROW=IRH,IRL
           ZZU    =C_ROT(JROW)*ZU(JROW)-S_ROT(JROW)*ZV(JROW)
           ZZV    =C_ROT(JROW)*ZV(JROW)+S_ROT(JROW)*ZU(JROW)
          ZU(JROW)=U(JPT,JROW)+(ZZU     -U(JPT,JROW))*A_WE(1  ,JROW)
          ZV(JROW)=V(JPT,JROW)+(ZZV     -V(JPT,JROW))*A_WE(1  ,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
!      SOLVE IMPLICIT EQUATION FOR WRAP-ROUND FILTER OF LATITUDE CIRCLE
       DO      JROW=IRF,IRG
         ZU(JROW)=ZU(JROW)/(1.-A_PRODUCT(JROW))
         ZV(JROW)=ZV(JROW)/(1.-A_PRODUCT(JROW))
       ENDDO ! JROW
       DO      JROW=IRH,IRL
         ZU(JROW)=ZU(JROW)/(1.-A_PRODUCT(JROW))
         ZV(JROW)=ZV(JROW)/(1.-A_PRODUCT(JROW))
       ENDDO ! JROW
!L         FILTER E->W PASS
       DO        JPT=ROW_LENGTH,1,-1
       DO        JROW=IRF,IRG
         U(JPT,JROW)=U(JPT,JROW)+A_WE(1  ,JROW)*                        &
     &    (C_ROT(JROW)*ZU(JROW)-S_ROT(JROW)*ZV(JROW)-U(JPT,JROW))
         V(JPT,JROW)=V(JPT,JROW)+A_WE(1  ,JROW)*                        &
     &    (C_ROT(JROW)*ZV(JROW)+S_ROT(JROW)*ZU(JROW)-V(JPT,JROW))
         ZU(JROW)=U(JPT,JROW)
         ZV(JROW)=V(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
       DO        JPT=ROW_LENGTH,1,-1
       DO        JROW=IRH,IRL
         U(JPT,JROW)=U(JPT,JROW)+A_WE(1  ,JROW)*                        &
     &    (C_ROT(JROW)*ZU(JROW)-S_ROT(JROW)*ZV(JROW)-U(JPT,JROW))
         V(JPT,JROW)=V(JPT,JROW)+A_WE(1  ,JROW)*                        &
     &    (C_ROT(JROW)*ZV(JROW)+S_ROT(JROW)*ZU(JROW)-V(JPT,JROW))
         ZU(JROW)=U(JPT,JROW)
         ZV(JROW)=V(JPT,JROW)
       ENDDO !)) JROW
       ENDDO !)) JPT
! *** END OF EXTRA W-E-W PASS NEAR POLES
       IF (lhook) CALL dr_hook('RFVVG',zhook_out,zhook_handle)
       RETURN
      END SUBROUTINE RFVVG

END MODULE rfvvg_mod
