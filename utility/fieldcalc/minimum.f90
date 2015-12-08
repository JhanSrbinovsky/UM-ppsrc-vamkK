! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE Minimum( Field1,       &  ! in
                    Field2,       &  ! in
                    MinField,     &  ! inout
                    ErrorStatus )    ! inout

! Description:
!   Finds the minimum of two fields at each point and outputs a
!   field holding the minimum values.
!
! Method:
!   1) Allocates the minimum field structure
!   2) Where there are data values in at least one field, it takes
!      the minimum of the two input fields using the Fortran
!      MIN function.
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

USE Ereport_mod, ONLY : Ereport
IMPLICIT NONE

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN)    :: Field1      ! Fields to find
TYPE(PP_Field_type), INTENT(IN)    :: Field2      !    minimum of
TYPE(PP_Field_type), INTENT(INOUT) :: MinField    ! Max of two fields
INTEGER, INTENT(INOUT)             :: ErrorStatus
INTEGER :: i,j

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "minimum"
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
                "Cannot find minimum of fields of different dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! If there is data in MinField, and it is neither of the two input
! fields, get rid of it
IF ( ASSOCIATED( MinField % RData ) .AND. &
    (MinField % ArrayPos /= Field1 % ArrayPos) .AND.  &
    (MinField % ArrayPos /= Field2 % ArrayPos) ) THEN
  DEALLOCATE( MinField % RData )
  NULLIFY( MinField % RData )
END IF

IF ( .NOT. ASSOCIATED( MinField % RData ) ) THEN
  ALLOCATE( MinField % RData(Field1 % Hdr % NumCols, &
                             Field1 % Hdr % NumRows) )
  MinField % Hdr = Field1 % Hdr
END IF

MinField % Hdr % BMDI = RMDI

!Find the minimum of the two fields
DO j = 1, Field1 % Hdr % NumRows
  DO i = 1, Field1 % Hdr % NumCols

    IF (Field1 % Rdata(i,j) /= Field1 % Hdr % BMDI .AND. &
        Field2 % Rdata(i,j) /= Field2 % Hdr % BMDI ) THEN
      MinField % RData(i,j) = MIN( Field1 % RData(i,j), Field2 % Rdata(i,j) )
    ELSE
      MinField % Rdata(i,j) = RMDI
    END IF

  END DO
END DO


9999 CONTINUE

END SUBROUTINE minimum
