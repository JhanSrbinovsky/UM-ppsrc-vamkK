! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate thermal advection diagnostic

SUBROUTINE ThermAdv( UField,       &  ! in
                     VField,       &  ! in
                     TField,       &  ! in
                     AdvField,     &  ! inout
                     ErrorStatus )    ! inout

! Description:
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE conversions_mod, ONLY: pi_over_180

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE FldCodes_Mod, ONLY:   &
  ST_ThAdv,  MO8_ThAdv,  PP_ThAdv

USE earth_constants_mod, ONLY: earth_radius

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: UField      ! U Component of wind
TYPE(PP_Field_type), INTENT(IN) :: VField      ! V Component of wind
TYPE(PP_Field_type), INTENT(IN) :: TField      ! Temperature

TYPE(PP_Field_type), INTENT(INOUT) :: AdvField ! Thermal Advection field
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "ThermAdv"
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Local Variables:
INTEGER :: NumCols, NumRows
INTEGER :: i, ip1, im1, j
REAL :: HdelX, HdelY, alat
REAL,ALLOCATABLE :: UdTdX(:,:)  ! u*longitudinal derivative of t
REAL,ALLOCATABLE :: VdTdY(:,:)  ! v*latitudinal derivative of t

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( TField % Hdr % LBHEM /= 0 ) THEN  ! Not global field - error
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, &
                "This routine is for Global fields only" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF
! Must be regular, un-rotated lat/lon grid only
IF ( TField % Hdr % LBCODE /= 1 ) THEN
  ErrorStatus = TField % Hdr % LBCODE

  CALL EReport( RoutineName, ErrorStatus,                        &
                "Cannot calculate thermal advection - " // &
                "incorrect grid type" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED( AdvField % RData ) ) THEN
  DEALLOCATE( AdvField % RData )
END IF
AdvField % Hdr = TField % Hdr
AdvField % Hdr % PPcode  =  PP_ThAdv
AdvField % Hdr % MO8Type = MO8_ThAdv
AdvField % Hdr % STCode  =  ST_ThAdv
AdvField % Hdr % RLevel  = -1.0       ! This will change at UM5.3
AdvField % Hdr % BMDI    = RMDI
ALLOCATE( AdvField % RData(AdvField % Hdr % NumCols, &
                           AdvField % Hdr % NumRows) )

NumCols = AdvField % Hdr % NumCols
NumRows = AdvField % Hdr % NumRows

ALLOCATE( UdTdX(UField % Hdr % NumCols, &
                UField % Hdr % NumRows) )
ALLOCATE( VdTdY(UField % Hdr % NumCols, &
                UField % Hdr % NumRows) )

! calculate latitude, u*dtdx and v*dtdy on UV grid (T grid larger)
HdelY = 0.5/( Earth_Radius * Pi_Over_180 * TField % Hdr % LatInt )
DO j = 1,UField % Hdr % NumRows
  alat  = TField % Hdr % ZerothLat + &
          TField % Hdr % LatInt * (FLOAT(j)+0.5)
  HdelX = 0.5/( Earth_Radius * Pi_Over_180 * TField % Hdr % LonInt * &
                   COS( alat * Pi_Over_180) )
  DO i = 1,UField % Hdr % NumCols
    ip1 = i+1
    IF ( i >= UField % Hdr % NumCols ) THEN
      ip1 = ip1 - UField % Hdr % NumCols   ! wrap around
    END IF
    UdTdX(i,j) = UField % RData(i,j) * HdelX *                       &
            ( (TField % RData(ip1,j  ) + TField % RData(ip1,j+1)) -  &
              (TField % RData(i  ,j  ) + TField % RData(i  ,j+1)) )
    VdTdY(i,j) = VField % RData(i,j) * HdelY *                       &
            ( (TField % RData(i  ,j+1) + TField % RData(ip1,j+1)) -  &
              (TField % RData(i  ,j  ) + TField % RData(ip1,j  )) )
  END DO
END DO

!------------------------------------------------------------------
! Calculation of thermal advection using UM4.5 formulae
! Thermal advection=-grad(T).V
! Average values of dtdx and dtdy are mutiplied by values of u and v
! respectively at the 4 U/V points surrounding the T point of interest.
! THADV is set to the average of the 4 values
!------------------------------------------------------------

DO j = 2,NumRows-1
  im1 = NumCols
  DO i = 1, NumCols

    AdvField % RData(i,j) = -0.25 * ( UdTdX(i  ,j) + UdTdX(i  ,j-1) +  &
                                      UdTdX(im1,j) + UdTdX(im1,j-1) +  &
                                      VdTdY(i  ,j) + VdTdY(i  ,j-1) +  &
                                      VdTdY(im1,j) + VdTdY(im1,j-1) )
    im1 = i
  END DO
END DO

IF ( TField % Hdr % LBHEM == 0 ) THEN
  ! average 1st and last rows of UV grid to copy into N & S poles
  AdvField % RData(1:NumCols,      1) =                                &
                                - (  SUM(UdTdX(1:NumCols,1))           &
                                   + SUM(VdTdY(1:NumCols,1)) )         &
                                         / FLOAT(TField % Hdr % NumCols)
  AdvField % RData(1:NumCols,NumRows) =                                &
                                - (  SUM(UdTdX(1:NumCols,NumRows-1))   &
                                   + SUM(VdTdY(1:NumCols,NumRows-1)) ) &
                                         / FLOAT(TField % Hdr % NumCols)
ELSE
  ! Limited area - set N & S rows to missing data
  AdvField % RData(1:NumCols,      1) = AdvField % Hdr % BMDI
  AdvField % RData(1:NumCols,NumRows) = AdvField % Hdr % BMDI
END IF
!IF ( TField % Hdr % LBHEM == 3 ) THEN
  ! No longitude wrap - set E & W cols to missing data
! See original code in DIA_THADV (UM) for afternative method that
! uses extrapolation
!  AdvField % RData(      1,1:NumRows) = AdvField % Hdr % BMDI
!  AdvField % RData(NumCols,1:NumRows) = AdvField % Hdr % BMDI
!END IF

DEALLOCATE( UdTdX )
DEALLOCATE( VdTdY )

9999 CONTINUE

END SUBROUTINE ThermAdv

