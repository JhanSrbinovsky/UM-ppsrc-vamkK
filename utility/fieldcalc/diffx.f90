! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to calculate horizontal derivatives.


SUBROUTINE DiffX( FField,       &  ! in
                  dFdX,         &  ! inout
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

USE earth_constants_mod, ONLY: earth_radius

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: FField   ! Input variable, F

TYPE(PP_Field_type), INTENT(INOUT) :: dFdX  ! Lon derivative of F
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DiffX"
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
INTEGER :: j
INTEGER :: NumCols
REAL :: HdelX, alat      ! temporary scalar

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Must have a regular, non-rotated lat/lon grid
IF ( FField % Hdr % LBCODE /= 1 ) THEN
  ErrorStatus = FField % Hdr % LBCODE

  CALL EReport( RoutineName, ErrorStatus,                        &
                "Cannot calculate longitudinal derivative - " // &
                "incorrect grid type" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED( dFdX % RData ) ) THEN
  DEALLOCATE( dFdX % RData )
END IF
dFdX % Hdr = FField % Hdr
ALLOCATE( dFdX % RData(dFdX % Hdr % NumCols, &
                       dFdX % Hdr % NumRows) )
dFdX % Hdr % BMDI = RMDI
dFdX % Hdr % STCode = IMDI

NumCols = FField % Hdr % NumCols
DO j = 1,FField % Hdr % NumRows
  alat  = FField % Hdr % ZerothLat + &
          FField % Hdr % LatInt * FLOAT(j)
  HdelX = 0.5/( Earth_Radius * Pi_Over_180 * FField % Hdr % LonInt * &
                          COS(Pi_Over_180*alat))
  dFdX % RData(1,j)          = HdelX * ( FField % RData(2          ,j) &
                                       - FField % RData(NumCols    ,j) )
  dFdX % RData(2:NumCols-1,j)= HdelX * ( FField % RData(3:NumCols  ,j) &
                                       - FField % RData(1:NumCols-2,j) )
  dFdX % RData(NumCols,j)    = HdelX * ( FField % RData(1          ,j) &
                                       - FField % RData(NumCols-1  ,j) )
END DO

IF ( FField % Hdr % LBHEM == 3 ) THEN
  WRITE(6,*) "DiffX Warning: No wrap-around - ", &
             "setting E & W boundaries to missing"
  dFdX % RData(1,      :) = RMDI
  dFdX % Rdata(NumCols,:) = RMDI
END IF

9999 CONTINUE

END SUBROUTINE DiffX

!=======================================================================


!=======================================================================


