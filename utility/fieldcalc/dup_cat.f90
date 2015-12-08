! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE Dup_CAT( Field1,       &  ! in
                    Field2,       &  ! inout
                    ErrorStatus )    ! inout

! Description:
!   Duplicates an input CAT field and modifies the parameters
!   in the header to label it as a maximum rather than mean field.
!   (Although the fields will be identical on output from FieldCalc,
!   they will be converted into mean and maximum fields when they
!   are converted to thinned GRIB fields by software on COSMOS.)
!
!   Note: this subroutine was written for use within the WAFC
!         turbulence action (WAFC_turb) only.
!
! Method:
!    1) Set up header for duplicate field (held in Field2) using values
!       in the header of the field to be duplicated (held in Field1).
!    2) Re-set the fieldcodes in the header for the duplicate field to
!       those for maximum CAT
!    3) Copy the data in Field1 to Field2
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
USE FldCodes_mod, ONLY:   &
  ST_MxCATPt, MO8_MxCATPt,&
  PP_MxCATPt, VC_MxCATPt
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning

IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN)    :: Field1     ! input CAT field
TYPE(PP_Field_type), INTENT(INOUT) :: Field2     ! output duplicate CAT fie
INTEGER, INTENT(INOUT)             :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Dup_CAT"

! Local Variables:
INTEGER :: i,j               !Loop counters

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

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF


!-----------------------------------------------------------------------
! Allocate duplicate field
IF ( ASSOCIATED( Field2 % RData )  ) THEN
  DEALLOCATE( Field2 % RData )
END IF
Field2 % Hdr = Field1 % Hdr
Field2 % Hdr % PPCode  =  PP_MxCATPt
Field2 % Hdr % MO8Type = MO8_MxCATPt
Field2 % Hdr % STCode  =  ST_MxCATPt
Field2 % Hdr % LBVC    =  VC_MxCATPt
Field2 % Hdr % BMDI    = RMDI

ALLOCATE( Field2 % RData(Field2 % Hdr % NumCols, &
                         Field2 % Hdr % NumRows) )


!----------------------------------------------------------------------
! Copy Field1 to Field2
DO j = 1, Field1 % Hdr % NumRows
  DO i = 1, Field1 % Hdr % NumCols
     Field2 % RData(i,j) = Field1 % RData(i,j)
  END DO
END DO

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Dup_CAT
