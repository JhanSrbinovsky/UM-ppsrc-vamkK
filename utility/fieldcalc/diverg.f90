! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate divergence diagnostic

SUBROUTINE Diverg( UField,       &  ! in
                   VField,       &  ! in
                   DivField,     &  ! inout
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

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE FldCodes_Mod, ONLY:   &
  ST_Diverg, MO8_Diverg, PP_Diverg
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: UField       ! U Component of wind
TYPE(PP_Field_type), INTENT(IN) :: VField       ! V Component of wind

TYPE(PP_Field_type), INTENT(INOUT) :: DivField  ! Divergence field
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Diverg"
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
TYPE(PP_Field_type) :: DuDx      ! longitudinal derivative of u
TYPE(PP_Field_type) :: DvDy      ! latitudinal derivative of v

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

Nullify ( DuDx % RData )
Nullify ( DvDy % RData )

IF ( ( UField % Hdr % NumCols /= VField % Hdr % NumCols )  .OR.   &
     ( UField % Hdr % NumRows /= VField % Hdr % NumRows ) )     THEN
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, &
               "U&V Fields must have the same dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED( DivField % RData ) ) THEN
  DEALLOCATE( DivField % RData )
END IF
DivField % Hdr = UField % Hdr
DivField % Hdr % PPCode  =  PP_Diverg
DivField % Hdr % MO8Type = MO8_Diverg
DivField % Hdr % STCode  =  ST_Diverg
DivField % Hdr % BMDI    = RMDI
ALLOCATE( DivField % RData(DivField % Hdr % NumCols, &
                           DivField % Hdr % NumRows) )

! DEPENDS ON: diffx
CALL DiffX   ( UField, DuDx, ErrorStatus )  ! lat derivative
! DEPENDS ON: diffcosy
CALL DiffCOSY( VField, DvDy, ErrorStatus )  ! lon derivative
IF ( ErrorStatus /= StatusOK ) THEN
  GO TO 9999
END IF

WHERE( (DuDx % RData /= DuDx % Hdr % BMDI) .AND. &
       (DvDy % RData /= DvDy % Hdr % BMDI) )
  DivField % RData = DuDx % RData + DvDy % RData
ELSEWHERE
  DivField % RData = RMDI
END WHERE

DEALLOCATE( DuDx % RData )
DEALLOCATE( DvDy % RData )

9999 CONTINUE

END SUBROUTINE Diverg

