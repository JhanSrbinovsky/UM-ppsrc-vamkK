! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE Temperature_Mask( TempRange,    &  ! in
                     NumLevs,              &  ! in
                     TFields,              &  ! in
                     MFields,              &  ! inout
                     ErrorStatus )            ! inout

! Description:
!   Outputs a mask for a given temperature range. This comprises a set
!   of fields corresponding to the input temperature fields. In these
!   fields every gridbox with a temperature within the given range has
!   a 1 assigned to it. Other gridboxes have 0's assigned to them.
!
! Method:
!   The output fields are initialised with the same headers as the
!   input temperature fields. The stash codes are then adjusted and the
!   so that the output fields are not confused with temperature fields.
!   The output fields are then set to 1 where the corresponding input
!   temperature lies between a supplied minumum and maxium temperature,
!   which should be in Kelvin, and are set to 0 elsewhere.
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
IMPLICIT None

! Subroutine Arguments:
REAL,    INTENT(IN) :: TempRange(2)     ! Temperature range in Kelvin
INTEGER, INTENT(IN) :: NumLevs
TYPE(PP_Field_type), INTENT(IN)    :: TFields(NumLevs) ! Temp fields
TYPE(PP_Field_type), INTENT(INOUT) :: MFields(NumLevs) ! Mask fields
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Temperature_Mask"
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

DO k = 1, NumLevs

  IF ( ASSOCIATED( MFields(k) % RData ) ) THEN
    DEALLOCATE( MFields(k) % RData )
    NULLIFY( MFields(k) % RData )
  END IF

  IF ( .NOT. ASSOCIATED( MFields(k) % RData ) ) THEN
    ALLOCATE( MFields(k) % RData(TFields(k) % Hdr % NumCols,     &
                                 TFields(k) % Hdr % NumRows) )
  END IF

END DO

MFields % Hdr = TFields % Hdr
MFields % Hdr % BMDI    = RMDI

WHERE ( (TFields % Hdr % STCode < 50000) )
  MFields % Hdr % STCode = TFields % Hdr % STCode + 50000
END WHERE

DO k = 1, NumLevs
  DO j = 1, TFields(k) % Hdr % NumRows
    DO i = 1, TFields(k) % Hdr % NumCols

      IF (TFields(k) % RData(i,j) == TFields(k) % Hdr % BMDI) THEN

        MFields(k) % RData(i,j) = MFields(k) % Hdr % BMDI

      ELSE IF ( (TFields(k) % RData(i,j) /= TFields(k) % Hdr % BMDI) &
                 .AND. (TFields(k) % RData(i,j) > TempRange(1))  &
                 .AND. (TFields(k) % RData(i,j) < TempRange(2))) THEN

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

END SUBROUTINE Temperature_Mask
