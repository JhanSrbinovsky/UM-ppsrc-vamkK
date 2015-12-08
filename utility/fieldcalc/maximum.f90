! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE Maximum( Field1,       &  ! in
                    Field2,       &  ! in
                    MaxField,     &  ! inout
                    ErrorStatus )    ! inout

! Description:
!   Finds the maximum of two fields at each point and outputs a
!   field holding the maximum values.
!
! Method:
!   1) Allocates the maximum field structure
!   2) Where there are data values in at least one field, it takes
!      the maximum of the two input fields using the Fortran
!      MAX function.
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
TYPE(PP_Field_type), INTENT(IN)    :: Field1      ! Fields to find
TYPE(PP_Field_type), INTENT(IN)    :: Field2      !    maximum of
TYPE(PP_Field_type), INTENT(INOUT) :: MaxField    ! Max of two fields
INTEGER, INTENT(INOUT)             :: ErrorStatus
INTEGER :: i,j

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Maximum"
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

IF ( (Field1 % Hdr % NumCols /= Field2 % Hdr % NumCols) .OR. &
     (Field1 % Hdr % NumRows /= Field2 % Hdr % NumRows) ) THEN
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, &
                "Cannot find maximum of fields of different dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! If there is data in MaxField, and it is neither of the two input
! fields, get rid of it
IF ( ASSOCIATED( MaxField % RData ) .AND. &
    (MaxField % ArrayPos /= Field1 % ArrayPos) .AND.  &
    (MaxField % ArrayPos /= Field2 % ArrayPos) ) THEN
  DEALLOCATE( MaxField % RData )
  NULLIFY( MaxField % RData )
END IF

IF ( .NOT. ASSOCIATED( MaxField % RData ) ) THEN
  ALLOCATE( MaxField % RData(Field1 % Hdr % NumCols, &
                             Field1 % Hdr % NumRows) )
  MaxField % Hdr = Field1 % Hdr
END IF

MaxField % Hdr % BMDI = RMDI

!Find the maximum of the two fields
DO j = 1, Field1 % Hdr % NumRows
  DO i = 1, Field1 % Hdr % NumCols

    IF (Field1 % Rdata(i,j) /= Field1 % Hdr % BMDI .AND. &
        Field2 % Rdata(i,j) /= Field2 % Hdr % BMDI ) THEN
      MaxField % RData(i,j) = MAX( Field1 % RData(i,j), Field2 % Rdata(i,j) )
    ELSE
      MaxField % Rdata(i,j) = RMDI
    END IF

  END DO
END DO


9999 CONTINUE

END SUBROUTINE Maximum
