! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate pressure layer fields from model layer fields
! Used in the calculation of icing potential diagnostics.

SUBROUTINE LIOnPLyr( NumLevs ,      &  ! in
                     NumLyrs,       &  ! in
                     LyrMO8L,       &  ! in
                     LyrLwrB,       &  ! in
                     LyrUprB,       &  ! in
                     LevFields ,    &  ! in
                     PFields ,      &  ! in
                     LyrFields,     &  ! inout
                     ErrorStatus)      ! inout

! Description:
!   Calculates a vector of pressure level fields from a vector of model
!   level fields.
!
! Method:
!   Two vectors of pressures are supplied.  These contain the lower and
!   upper boundaries of a set of layers of the atmosphere.  For each of
!   these layers the corresponding data on the model level fields is
!   found using pressures on model levels supplied.  To each gridpoint
!   in the fields comprising the output, the mean of the corresponding
!   data is assigned.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:

! Number of levels for input fields
INTEGER, INTENT(IN) :: NumLevs

! Number of output (std level) levels
INTEGER, INTENT(IN) :: NumLyrs

! MO8 level code of output layers
INTEGER, INTENT(IN) :: LyrMO8L(NumLyrs)

! Lower boundary of output layers
REAL,    INTENT(IN) :: LyrLwrB(NumLyrs)

! Upper boundary of output layers
REAL,    INTENT(IN) :: LyrUprB(NumLyrs)

! Model level fields of parameter to be determined at std levels
TYPE(PP_Field_type), INTENT(IN)    :: LevFields(NumLevs)

! Model level Pressure
TYPE(PP_Field_type), INTENT(IN)    :: PFields(NumLevs)

! Output fields on std levels
TYPE(PP_Field_type), INTENT(INOUT) :: LyrFields(NumLyrs)

! Error return code
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "LIOnPLyr"
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
INTEGER :: i, j, k              ! Loop counters
INTEGER :: n                    ! Number of model level fields lying
                                !   within atmospheric layer/slab
REAL    :: RData(NumLevs)       ! Holds data part of LevFields(k)
REAL    :: Press(NumLevs)       ! Holds pressure on model levels in hPa
LOGICAL :: NotMDI(NumLevs)      ! Mask showing where LevFields and corresponding
                                !   pressure fields do not have missing data
LOGICAL :: Mask(NumLevs)        ! Holds mask showing which model level pressure
                                !   fields are within the slab of atmosphere

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )


IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF (NumLevs < NumLyrs) THEN

  CALL EReport( RoutineName, ErrorStatus,                           &
                "Too few model level fields supplied")
  ErrorStatus = StatusWarning
  GO TO 9999
END IF


!-----------------------------------------------------------------------
! Check input model level fields and pressure fields are on the same grid
DO k = 1, NumLevs

  IF ( (LevFields(k) % Hdr % NumCols   /=                            &
        PFields(k) % Hdr % NumCols) .OR.                             &
       (LevFields(k) % Hdr % NumRows   /=                            &
        PFields(k) % Hdr % NumRows) .OR.                             &
       (LevFields(k) % Hdr % ZerothLon /=                            &
        PFields(k) % Hdr % ZerothLon) .OR.                           &
       (LevFields(k) % Hdr % ZerothLat /=                            &
        PFields(k) % Hdr % ZerothLat) ) THEN

    CALL EReport( RoutineName, ErrorStatus,                          &
                "Model level fields supplied are on different grids")
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF

END DO


!-----------------------------------------------------------------------
! Allocate output fields
DO k = 1, NumLyrs

  IF ( ASSOCIATED( LyrFields(k) % RData ) ) THEN
    DEALLOCATE( LyrFields(k) % RData )
    NULLIFY( LyrFields(k) % RData )
  END IF

  IF ( .NOT. ASSOCIATED( LyrFields(k) % RData ) ) THEN
    ALLOCATE( LyrFields(k) % RData(LevFields(1) % Hdr % NumCols,     &
                                   LevFields(1) % Hdr % NumRows) )
  END IF

END DO

LyrFields(1:NumLyrs) % Hdr = LevFields(1:NumLyrs) % Hdr
LyrFields(1:NumLyrs) % Hdr % MO8Level    = LyrMO8L
LyrFields(1:NumLyrs) % Hdr % RLevel      = REAL(LyrMO8L)
LyrFields(1:NumLyrs) % Hdr % RefLevel    = LyrLwrB

LyrFields % Hdr % LBVC   = 8
LyrFields % Hdr % LBRVC  = 8

LyrFields % Hdr % BULev  = 0.0
LyrFields % Hdr % BHULev = 0.0
LyrFields % Hdr % BHLev  = 0.0
LyrFields % Hdr % BHRLev = 0.0



!-----------------------------------------------------------------------
! Convert model level fields to pressure levels by determining the
! number of model pressure levels which are in the surrounding slab
! of atmosphere and calculating the mean of the model level fields
! in that slab of atmosphere

!First loop over gridpoints
DO j = 1, LyrFields(1) % Hdr % NumRows
  DO i = 1, LyrFields(1) % Hdr % NumCols

    DO k = 1, NumLevs
      RData(k) = LevFields(k) % RData(i,j)
      Press(k) = 0.01*PFields(k) % RData(i,j)
    END DO

    NotMDI = (RData /= LevFields(1) % Hdr % BMDI) .AND.            &
             (Press /= 0.01*PFields(1) % Hdr % BMDI)

   !Now loop over desired atmosphere layers
    DO k = 1, NumLyrs

      Mask = NotMDI .AND.                                          &
             (Press <= LyrLwrB(k)) .AND. (Press >= LyrUprB(k))
      n = COUNT( Mask )

      IF (n > 0) THEN
        LyrFields(k) % RData(i,j) = SUM( RData, Mask = Mask )
        LyrFields(k) % RData(i,j) =                                &
                       LyrFields(k) % RData(i,j) / REAL( n )
      ELSE

        LyrFields(k) % RData(i,j) = LyrFields(k) % Hdr % BMDI

      END IF

    END DO

  END DO
END DO


9999 CONTINUE


! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE LIOnPLyr
