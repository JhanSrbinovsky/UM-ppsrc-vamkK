! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generic routines for manipulating pp-fields within Fieldcalc

!=======================================================================
! Vector direction - NOT wind direction
SUBROUTINE VecDir( Field1,      &  ! in
                   Field2,      &  ! in
                   VDField,     &  ! inout
                   ErrorStatus )   ! inout

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

USE conversions_mod, ONLY: pi

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN)    :: Field1
TYPE(PP_Field_type), INTENT(IN)    :: Field2
TYPE(PP_Field_type), INTENT(INOUT) :: VDField

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "VecDir"
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

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ( Field1 % Hdr % NumCols /= Field2 % Hdr % NumCols )  .OR.   &
     ( Field1 % Hdr % NumRows /= Field2 % Hdr % NumRows ) )     THEN
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus,   &
               "VecDir : Fields must have the same dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED(VDField % RData) ) THEN
  DEALLOCATE( VDField % RData )
END IF
VDField % Hdr = Field1 % Hdr
ALLOCATE( VDField % RData(VDField % Hdr % NumCols, &
                          VDField % Hdr % NumRows) )
VDField % Hdr % BMDI = RMDI
VDField % Hdr % STCode = IMDI

WHERE( (Field1 % RData == 0.0) .AND. (Field2 % RData == 0.0) )
  VDField % RData = 0.0
ELSEWHERE ( (Field1 % RData /= Field1 % Hdr % BMDI) .AND.  &
            (Field2 % RData /= Field2 % Hdr % BMDI) )
  VDField % RData = ATAN2( Field1 % RData, Field2 % RData )*(180.0/Pi)
ELSEWHERE
  VDField % RData = RMDI
END WHERE

9999 CONTINUE

END SUBROUTINE VecDir

