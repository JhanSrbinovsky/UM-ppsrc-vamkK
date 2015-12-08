! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generic routines for manipulating pp-fields within Fieldcalc


!=======================================================================

!=======================================================================

!=======================================================================

!=======================================================================
SUBROUTINE VecMag( Field1,      &  ! in
                   Field2,      &  ! in
                   VMField,     &  ! inout
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

TYPE(PP_Field_type), INTENT(INOUT) :: VMField
INTEGER, INTENT(INOUT) :: ErrorStatus
INTEGER :: i,j

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "VecMag"
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

IF ( ( Field1 % Hdr % NumCols /= Field2 % Hdr % NumCols )  .OR. &
     ( Field1 % Hdr % NumRows /= Field2 % Hdr % NumRows ) ) THEN
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, &
               "Vector fields must have the same dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED( VMField % RData ) ) THEN
  DEALLOCATE( VMField % RData )
END IF
VMField % Hdr = Field1 % Hdr
VMField % Hdr % STCode  = IMDI
VMField % Hdr % BMDI    = RMDI
ALLOCATE( VMField % RData(VMField % Hdr % NumCols, &
                          VMField % Hdr % NumRows) )

DO j = 1, VMField % Hdr % NumRows
  DO i = 1, VMField % Hdr % NumCols
    IF (Field1 % RData(i,j) /= Field1 % Hdr % BMDI ) THEN
      IF (Field2 % RData(i,j) /= Field2 % Hdr % BMDI ) THEN
        VMField % RData(i,j) = SQRT((Field1 % RData(i,j))**2 +  &
                                    (Field2 % RData(i,j))**2 )
      ELSE
        VMField % RData(i,j) = RMDI
      END IF
    ELSE
      VMField % RData(i,j) = RMDI
    END IF
  END DO
END DO


9999 CONTINUE

END SUBROUTINE VecMag

!=======================================================================
! Vector direction - NOT wind direction

