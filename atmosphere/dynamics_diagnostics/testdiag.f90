! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate simple test diagnostics based on a simple analytic formula
!
! Subroutine Interface:
      SUBROUTINE Testdiag(                                              &
     &  p_field,v_field,rows,n_rows,row_length                          &
     & ,EW_SPACE,NS_SPACE,FIRST_LAT,FIRST_LONG,PHI_POLE,LAMBDA_POLE     &
     & ,ELF                                                             &
     & ,PRESS_LEVELS_LIST,NO_PRESS_LEVELS                               &
     & ,MODEL_LEVELS_LIST,NO_MODEL_LEVELS,FORECAST_HRS                  &
     & ,DIAG1,DIAG2,DIAG3,DIAG4                                         &
     & ,qdia1,qdia2,qdia3,qdia4)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE eqtoll_mod, ONLY: eqtoll
      IMPLICIT NONE
!
! Description:

! Calculate simple test diagnostics based on a simple analytic formula:
!
!    VALUE=A*(LATITUDE+90.)+B*LONGITUDE+C*LEVEL+D*FORECAST_HRS
!    where A=1.0, B=1.0E2, C=1.0E3, D=1.0E4
!    and (LAT,LONG) are in degrees, actual position (rotated for LAM),
!    LEVEL is either the model level (real number) or
!                        pressure level (mb),
!    FORECAST_HRS in T+hours after assimilation time.
!    Theses diagnostics are to be used for checking output procedures
!    for various post-processing routes.
!    Four diagnostics are supported:
!    1. single-level FIELD (LEVEL=0.) at V points of Arakawa C grid
!    2. single-level FIELD (LEVEL=0.) at P points of Arakawa C grid
!    3. multi -level FIELD (LEVEL=press level) at P points of C grid
!    4. multi -level FIELD (LEVEL=model level) at P points of C grid
!
! Method:
!      Calculate lat,long offsets for MPP
!   1. Calculate FIRST DIAGNOSTIC  (V GRID SINGLE LEVEL)
!   2. Calculate ACTUAL LATITUDES, LONGITUDES FOR P FIELDS (DIAG 2-4)
!   3. Calculate SECOND DIAGNOSTIC (P GRID SINGLE LEVEL)
!   4. Calculate THIRD  DIAGNOSTIC (P GRID PRESSURE LEVELS)
!   5. Calculate FOURTH DIAGNOSTIC (P GRID MODEL    LEVELS)
!
!    DOCUMENTATION:  UM Doc Paper D7
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     &  p_field                                                         &
                         !IN   Horizontal field size p points
     &, v_field                                                         &
                         !IN   Horizontal field size v points
     &, rows                                                            &
                         !IN   No. of rows for p field
     &, n_rows                                                          &
                         !IN   No. of rows for v field
     &, row_length                                                      &
                         !IN   No. of points per row
     &, NO_MODEL_LEVELS                                                 &
                         !IN   model levels for output
     &, NO_PRESS_LEVELS                                                 
                         !IN   press levels for output

! UM6.5 - MODEL_ANALYSIS_HRS chnged to real -
!                   requires FORECAST_HRS changed to REAL also
      REAL FORECAST_HRS     !IN   FORECAST HOURS T+0, etc

      LOGICAL                                                           &
     &  ELF                                                             &
                         !IN  TRUE IF MODEL IS LAM WITH ROTATED GRID
     & ,qdia1                                                           &
                         !IN  STASHflag for DIAG1
     & ,qdia2                                                           &
                         !IN  STASHflag for DIAG2
     & ,qdia3                                                           &
                         !IN  STASHflag for DIAG3
     & ,qdia4            !IN  STASHflag for DIAG4

      REAL                                                              &
     &  EW_SPACE                                                        &
                         !IN  DELTA LONGITUDE (DEGREES)
     &, NS_SPACE                                                        &
                         !IN  DELTA  LATITUDE (DEGREES)
     &, FIRST_LAT                                                       &
                         !IN  latitude  of first p row (degrees)
     &, FIRST_LONG                                                      &
                         !IN  longitude of first p col (degrees)
     &, PHI_POLE                                                        &
                         !IN  latitude  of the pseudo pole
     &, LAMBDA_POLE      !IN  longitude of the pseudo pole

!   Array  arguments with intent(in):
      REAL                                                              &
     &  MODEL_LEVELS_LIST(NO_MODEL_LEVELS)                              &
                                           !IN LEVELS list (for DIAG3)
     &, PRESS_LEVELS_LIST(NO_PRESS_LEVELS) !IN LEVELS list (for DIAG4)

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL                                                              &
     &  DIAG1(v_field)                                                  &
                                        !OUT DIAGNOSTIC 1
     &, DIAG2(p_field)                                                  &
                                        !OUT DIAGNOSTIC 2
     &, DIAG3(p_field,NO_PRESS_LEVELS)                                  &
                                        !OUT DIAGNOSTIC 3
     &, DIAG4(p_field,NO_MODEL_LEVELS)  !OUT DIAGNOSTIC 4


! Local parameters:
      REAL                                                              &
     &  A,B,C,D  ! COEFFICIENTS FOR CALCULATING VALUES OF FIELD
!
      PARAMETER(A=1.0,B=1.0E2,C=1.0E3,D=1.0E4)

! Local scalars:
      INTEGER                                                           &
     &  I,J,K                                                           &
                     !  LOOP COUNTERS
     & ,L                                                               &
                     !  LOOP INDEX
     & ,OFFSETX                                                         &
                     !  INDEX OFFSETs (for calculating MPP lat,longs)
     & ,OFFSETY

! Local dynamic arrays:
      REAL                                                              &
     &  LATITUDE(p_field)                                               &
                          ! latitude  in degrees
     &, LONGITUDE(p_field)                                              &
                          ! longitude in degrees
     &, LAT(p_field)                                                    &
                          ! latitude  in degrees on equatorial grid
     &, LONG(p_field)     ! longitude in degrees on equatorial grid

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!- End of header


!
!  Calculate lat,long offsets for MPP
!

       IF (lhook) CALL dr_hook('TESTDIAG',zhook_in,zhook_handle)
      OFFSETX=datastart(1)- Offx - 1
      OFFSETY=datastart(2)- Offy - 1

!-----------------------------------------------------------------------
!   1. CALCULATE FIRST DIAGNOSTIC  (V GRID SINGLE LEVEL)
!-----------------------------------------------------------------------
      IF(qdia1) THEN
!
!   1a. FIND EQUATORIAL LATITUDES,LONGITUDES
!
         DO J=1,n_rows
           DO I=1,row_length
             L= I + (J-1)*row_length
             LAT (L)=  FIRST_LAT  - NS_SPACE*(J+OFFSETY-0.5)
             LONG(L)=  FIRST_LONG + EW_SPACE*(I+OFFSETX-0.5)
           ENDDO
         ENDDO
!
!   1b. CONVERT TO ACTUAL LATITUDE,LONGITUDE IF ELF GRID
!
         IF(ELF) THEN
           CALL EQTOLL(LAT,LONG,LATITUDE,LONGITUDE,PHI_POLE,LAMBDA_POLE,&
     &                 v_field)
         ELSE
           DO I=1,v_field
             LATITUDE(I) =LAT(I)
             LONGITUDE(I)=LONG(I)
           ENDDO
         ENDIF
!
!   1c. CALCULATE VALUE FROM ANALYTIC FUNCTION
!

         DO I=1,v_field
           DIAG1(I)=A*(LATITUDE(I)+90.0) + B*LONGITUDE(I) +             &
     &              D*FORECAST_HRS
         ENDDO

      ENDIF               ! END OF qdia1 TEST

!-----------------------------------------------------------------------
!   2. CALCULATE ACTUAL LATITUDES, LONGITUDES FOR P FIELDS (DIAG 2-4)
!-----------------------------------------------------------------------
      IF(qdia2.OR.qdia3.OR.qdia4) THEN
!
!   2a. FIND EQUATORIAL LATITUDES,LONGITUDES
!
         DO J=1,rows
           DO I=1,row_length
             L= I + (J-1)*row_length
             LAT (L)=  FIRST_LAT  - NS_SPACE*(J+OFFSETY-1)
             LONG(L)=  FIRST_LONG + EW_SPACE*(I+OFFSETX-1)
           ENDDO
         ENDDO
!
!   2b. CONVERT TO ACTUAL LATITUDE,LONGITUDE IF ELF GRID
!
         IF(ELF) THEN
           CALL EQTOLL(LAT,LONG,LATITUDE,LONGITUDE,PHI_POLE,LAMBDA_POLE,&
     &                 p_field)
         ELSE
           DO I=1,p_field
             LATITUDE(I) =LAT(I)
             LONGITUDE(I)=LONG(I)
           ENDDO
         ENDIF
      ENDIF                      ! END OF qdia2-4 TEST
!-----------------------------------------------------------------------
!   3. CALCULATE SECOND DIAGNOSTIC (P GRID SINGLE LEVEL)
!-----------------------------------------------------------------------
      IF(qdia2) THEN

         DO I=1,p_field
           DIAG2(I)=A*(LATITUDE(I)+90.0) + B*LONGITUDE(I) +             &
     &              D*FORECAST_HRS
         ENDDO

      ENDIF               ! END OF qdia2 TEST
!-----------------------------------------------------------------------
!   4. CALCULATE THIRD  DIAGNOSTIC (P GRID PRESSURE LEVELS)
!-----------------------------------------------------------------------
      IF(qdia3) THEN

         DO K=1,NO_PRESS_LEVELS
           DO I=1,p_field
             DIAG3(I,K)=A*(LATITUDE(I)+90.0)   + B*LONGITUDE(I) +       &
     &                  C*PRESS_LEVELS_LIST(K) + D*FORECAST_HRS
           ENDDO
         ENDDO

      ENDIF               ! END OF qdia3 TEST
!-----------------------------------------------------------------------
!   5. CALCULATE FOURTH DIAGNOSTIC (P GRID MODEL    LEVELS)
!-----------------------------------------------------------------------
      IF(qdia4) THEN

         DO K=1,NO_MODEL_LEVELS
           DO I=1,p_field
             DIAG4(I,K)=A*(LATITUDE(I)+90.0)   + B*LONGITUDE(I) +       &
     &                  C*MODEL_LEVELS_LIST(K) + D*FORECAST_HRS
           ENDDO
         ENDDO

      ENDIF               ! END OF qdia4 TEST

      IF (lhook) CALL dr_hook('TESTDIAG',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Testdiag
!=======================================================================
