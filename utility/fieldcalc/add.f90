! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generic routines for manipulating pp-fields within Fieldcalc

!=======================================================================
SUBROUTINE Add( Factor,       &  ! in
                Field,        &  ! in
                AddField,     &  ! inout
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
  StatusOK
USE FldCodes_mod, ONLY:   &
  UMSectnNo
IMPLICIT None

! Subroutine Arguments:
REAL, INTENT(IN) :: Factor
TYPE(PP_Field_type), INTENT(IN) :: Field

TYPE(PP_Field_type), INTENT(INOUT) :: AddField ! Location of result
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Add"
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

INTEGER :: i, j

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! If there is data in MaxField, and it is neither of the two input
! fields, get rid of it
IF ( ASSOCIATED( AddField % RData ) .AND. &
    (AddField % ArrayPos /= Field % ArrayPos) ) THEN
  DEALLOCATE( AddField % RData )
  NULLIFY( AddField % RData )
END IF

IF ( .NOT. ASSOCIATED( AddField % RData ) ) THEN
  ALLOCATE( AddField % RData(Field % Hdr % NumCols, &
                             Field % Hdr % NumRows) )
  AddField % Hdr = Field % Hdr
END IF

AddField % Hdr % BMDI = RMDI
IF ( INT( Field % Hdr % STCode / 1000 ) /= UMSectnNo ) THEN
  AddField % Hdr % STCode = IMDI
ENDIF

!Add factor to field
DO j = 1, Field % Hdr % NumRows
  DO i = 1, Field % Hdr % NumCols

    IF (Field % Rdata(i,j) /= Field % Hdr % BMDI) THEN
      AddField % RData(i,j) = Factor + Field % RData(i,j)
    ELSE
      AddField % Rdata(i,j) = RMDI
    END IF

  END DO
END DO

9999 CONTINUE

END SUBROUTINE Add

