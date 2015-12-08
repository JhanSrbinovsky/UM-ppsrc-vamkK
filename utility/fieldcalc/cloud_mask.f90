! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate mask field from input field given a specified range.

SUBROUTINE Cloud_Mask( Range,                &  ! in
                       NumLevs,              &  ! in
                       IFields,              &  ! in
                       MFields,              &  ! inout
                       ErrorStatus )            ! inout

! Description:
!   Outputs a mask where cloud fraction is greater than 0. This comprises
!   a set of fields corresponding to the input cloud fields. In these
!   fields every gridbox with a value within the given range has a 1
!   assigned to it. Other gridboxes have 0's assigned to them.
!
! Method:
!   The output fields are initialised with the same headers as the
!   input fields. The stash codes are then adjusted so that the output
!   fields are not confused with input fields.
!   The output fields are then set to 1 where the corresponding input
!   value lies within the specified range, or is equal to the maximum
!   value, and are set to 0 elsewhere. (i.e. if there is any cloud
!   in the gridbox the mask field is set to 1.)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod
USE Err_Mod, ONLY:        &
  StatusOK
IMPLICIT None

! Subroutine Arguments:
REAL,    INTENT(IN) :: Range(2)                        ! Masking range
INTEGER, INTENT(IN) :: NumLevs                         ! Number of levels
TYPE(PP_Field_type), INTENT(IN)    :: IFields(NumLevs) ! Input fields
TYPE(PP_Field_type), INTENT(INOUT) :: MFields(NumLevs) ! Mask fields
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Cloud_Mask"
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
INTEGER :: i, j, k

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF


!Allocate output mask fields
DO k = 1, NumLevs

  IF ( ASSOCIATED( MFields(k) % RData ) ) THEN
    DEALLOCATE( MFields(k) % RData )
    NULLIFY( MFields(k) % RData )
  END IF

  IF ( .NOT. ASSOCIATED( MFields(k) % RData ) ) THEN
    ALLOCATE( MFields(k) % RData(IFields(k) % Hdr % NumCols,     &
                                 IFields(k) % Hdr % NumRows) )
  END IF

END DO


MFields % Hdr = IFields % Hdr
MFields % Hdr % BMDI    = RMDI

WHERE ( (IFields % Hdr % STCode < 50000) )
  MFields % Hdr % STCode = IFields % Hdr % STCode + 50000
END WHERE


!Fill in mask fields
DO k = 1, NumLevs
  DO j = 1, IFields(k) % Hdr % NumRows
    DO i = 1, IFields(k) % Hdr % NumCols

      IF (IFields(k) % RData(i,j) == IFields(k) % Hdr % BMDI) THEN

        MFields(k) % RData(i,j) = MFields(k) % Hdr % BMDI

      ELSE IF ( (IFields(k) % RData(i,j) /= IFields(k) % Hdr % BMDI) &
                 .AND. (IFields(k) % RData(i,j) > Range(1))  &
                 .AND. (IFields(k) % RData(i,j) <= Range(2))) THEN

        MFields(k) % RData(i,j) = 1.0

      ELSE

        MFields(k) % RData(i,j) = 0.0

      END IF

    END DO
  END DO

END DO


9999 CONTINUE


! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Cloud_Mask
