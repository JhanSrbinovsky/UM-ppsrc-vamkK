! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to calculate Clear Air Turbulence diagnostics.



!=======================================================================

SUBROUTINE DuttonCAT( UStdLev,     &  ! in
                      VStdLev,     &  ! in
                      dUdX,        &  ! in
                      dUdY,        &  ! in
                      dVdX,        &  ! in
                      dVdY,        &  ! in
                      dWdZ,        &  ! in
                      CATPred,     &  ! inout
                      ErrorStatus )   ! inout

! Description:
!   Routine to convert values of horizontal and vertical wind shear into
!   a CAT predictor between 0 and 7.5
!
! Method:
!   Calculate HWS and VWS from derivatives of wind components.  Use
!   Dutton's empirical formula to combine HWS and VWS (see Dutton, 1980:
!   Probability forecasts of clear air turbulence based on numerical
!   output: Met Mag v109 pp293-310).  Interpolate to values between 0
!   and 7.5.
!
! IMPORTANT NOTE:
!   The interpolation knots and Dutton's formula are empirically
!   derived, so adjustments are necessary with significant model changes
!   to maintain consistent output on AvSigWX charts.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: UStdLev       ! U & V wind comps on
TYPE(PP_Field_type), INTENT(IN) :: VStdLev       !   std level
TYPE(PP_Field_type), INTENT(IN) :: dUdX, dUdY    ! Horiz. derivs of u
TYPE(PP_Field_type), INTENT(IN) :: dVdX, dVdY    ! Horiz. derivs of v
TYPE(PP_Field_type), INTENT(IN) :: dWdZ          ! Vert. deriv of wspeed

TYPE(PP_Field_type), INTENT(INOUT) :: CATPred    ! CAT diagnostic
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DuttonCAT"
INTEGER,PARAMETER :: knots=3     ! No knots used to interp. from Dutton

! Local Variables:
INTEGER :: i, j                  ! Loop counters
INTEGER :: kk, kref, itest       ! For finding relevant knots for interp
REAL :: alat                     ! Latitude
REAL :: alpha                    ! Interpolation coefficient
REAL :: vmag                     ! Wind speed
REAL :: CATProb = 7.5            ! Final CAT predictor
REAL :: empire                   ! Dutton's empirical CAT predictor
REAL :: ShearH                   ! Horizontal wind shear (HWS)
REAL :: ShearV                   ! Vertical wind shear   (VWS)

REAL :: eval(3)   = (/5.0, 7.5 , 65.0/) ! Values of Dutton's
                                        ! empirical indicator.
REAL :: catval(3) = (/0.0, 1.75,  7.5/) ! Corresponding values of
                                        ! final CAT predictor

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

DO j = 1,CATPred % Hdr % NumRows
!--------------------------------------- calculate latitude
  alat = UStdLev % Hdr % ZerothLat + &
         UStdLev % Hdr % LatInt * FLOAT(j)
  DO i = 1,CATPred % Hdr % NumCols

    IF ( (dUdX % RData(i,j) == dUdX % Hdr % BMDI) .OR.  &
         (dUdY % RData(i,j) == dUdY % Hdr % BMDI) .OR.  &
         (dVdX % RData(i,j) == dVdX % Hdr % BMDI) .OR.  &
         (dVdY % RData(i,j) == dVdY % Hdr % BMDI) ) THEN
      CATPred % RData(i,j) = CATPred % Hdr % BMDI
    ELSE
      !  Requires ShearH in units of (m/s)/(100km),ShearV in units of
      !  (m/s)/km.  SI units are (m/s)/m, so need to take
      !  100000.*ShearH and 1000.*ShearV.
      ShearV = 1000. * dWdZ % RData(i,j)
      vmag  = UStdLev % RData(i,j)**2 + VStdLev % RData(i,j)**2
      IF ( vmag > 0.0 ) THEN
        ShearH = (100000./vmag) * &
         (  UStdLev%RData(i,j) * VStdLev%RData(i,j) * dUdX%RData(i,j) &
          - UStdLev%RData(i,j) * UStdLev%RData(i,j) * dUdY%RData(i,j) &
          + VStdLev%RData(i,j) * VStdLev%RData(i,j) * dVdX%RData(i,j) &
          - UStdLev%RData(i,j) * VStdLev%RData(i,j) * dVdY%RData(i,j) )
      ELSE
        ShearH = 0.0
      END IF
      ! vector form for shearH implies a change of sign if in S.Hemi
      IF ( alat < 0.0 ) THEN
        ShearH = -ShearH
      END IF

      empire = 10.5+1.5*(ShearH+ShearV)
      CATProb = catval(knots)
      itest = 0
! kref must be initialised here for vectorisation
      kref=0

      IF ( empire < eval(1) ) THEN
        CATProb = catval(1)
        itest = 2
      END IF
!      DO kk = 2,knots
!        IF ( (empire < eval(kk)) .AND.  &
!             (itest == 0) ) THEN
!          kref  = kk
!          itest = 1
!        END IF
!      END DO
! the loop is manually unrolled for vectorization, knowing that knots=3
!
      IF (( empire < eval(2)).and.( itest==0)) THEN
         kref=2
         itest=1
      END IF

      IF (( empire < eval(3)).and.( itest==0)) THEN
         kref=3
         itest=1
      END IF

      IF ( itest == 1 ) THEN
        alpha   = (empire-eval(kref-1))/(eval(kref)-eval(kref-1))
        CATProb = alpha*catval(kref) + (1.0-alpha)*catval(kref-1)
      END IF

      CATPred % RData(i,j) = CATProb
    END IF
  END DO
END DO

9999 CONTINUE

END SUBROUTINE DuttonCAT

!=======================================================================


!=======================================================================


