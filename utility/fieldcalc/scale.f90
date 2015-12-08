! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generic routines for manipulating pp-fields within Fieldcalc

!=======================================================================
SUBROUTINE Scale( Factor,       &  ! in
                  Field,        &  ! in
                  ScField,      &  ! inout
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

TYPE(PP_Field_type), INTENT(INOUT) :: ScField ! Location of result
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
INTEGER :: i, j
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Scale"
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

IF ( ASSOCIATED( ScField % RData )            .AND.  &
    (ScField % ArrayPos /= Field % ArrayPos) ) THEN
  DEALLOCATE( ScField % RData )
  NULLIFY( ScField % RData )
END IF
! IF SumField is now empty, allocate some space
IF ( .NOT. ASSOCIATED( ScField % RData ) ) THEN
  ALLOCATE( ScField % RData(Field % Hdr % NumCols, &
                            Field % Hdr % NumRows) )
  ScField % Hdr = Field % Hdr
END IF

ScField % Hdr % BMDI = RMDI
IF ( INT( Field % Hdr % STCode / 1000 ) /= UMSectnNo ) THEN
  ScField % Hdr % STCode = IMDI
ENDIF

DO j = 1, Field % Hdr % NumRows
  DO i = 1, Field % Hdr % NumCols

    IF (Field % Rdata(i,j) /= Field % Hdr % BMDI) THEN
      ScField % RData(i,j) = Factor * Field % Rdata(i,j)
    ELSE
      ScField % Rdata(i,j) = RMDI
    END IF

  END DO
END DO

9999 CONTINUE

END SUBROUTINE Scale

