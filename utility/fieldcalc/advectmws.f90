! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to calculate Mountain Wave Predictor diagnostic
! Note WindDir is wind direction in vector sense, not meteorological


!=======================================================================

SUBROUTINE AdvectMWS( NPts,         &  ! in
                      CutOff,       &  ! in
                      MinWSpeed,    &  ! in
                      PtPosn,       &  ! in
                      MaxStress,    &  ! in
                      WindSpd,      &  ! in
                      WindDir,      &  ! in
                      AdvPtPosn,    &  ! inout
                      MWTPred,      &  ! inout
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

USE conversions_mod, ONLY: pi_over_180

USE IO_mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal
IMPLICIT None

!     Description:
!       The input arrays wdir and PtPosn are used to calculate
!       the output arrays PtPosn1 and mwt_pt1 which contain the 2 and 1
!       dimension locations of the next point downwind.
!
!     Method:
!       See Documentation...

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NPts        ! The size of the arrays
REAL,    INTENT(IN) :: CutOff      ! Threshold value for this call
REAL,    INTENT(IN) :: MinWSpeed   ! WSpeed needed for advection
INTEGER, INTENT(IN) :: PtPosn(2,NPts)
TYPE(PP_Field_type), INTENT(IN) :: MaxStress  ! Maximum stress magnitude
TYPE(PP_Field_type), INTENT(IN) :: WindSpd    ! Wind speed field
TYPE(PP_Field_type), INTENT(IN) :: WindDir    ! Wind direction field

INTEGER, INTENT(INOUT) :: AdvPtPosn(2,NPts)   ! Lat Lon array of points
TYPE(PP_Field_type), INTENT(INOUT) :: MWTPred
INTEGER, INTENT(INOUT) :: ErrorStatus         ! Return error status

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "AdvectStress"

! Local Variables
INTEGER :: NumCols, NumRows
INTEGER :: i                 ! DO loop variable
INTEGER :: row,col           !
INTEGER :: new_row,new_col   !  Temps for storing row and col indicators
INTEGER :: first_row,first_col
REAL :: latitude             ! Latitude in degrees
REAL :: cosang               ! COS(latitude*deg2r)
REAL :: alpha1               ! The first comparison angle
REAL :: alpha2               ! The second comparison
REAL :: mod_deltalat         ! The absolute value of lat_interval
REAL :: wdir
REAL :: wspd

! End of header --------------------------------------------------------

IF( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

NumCols = MWTPred % Hdr % NumCols
NumRows = MWTPred % Hdr % NumRows

!     Check that there are some points. If so perform the calculation
!     for each point.
IF ( (NPts <= 0) .OR. (NumCols <= 0) .OR. (NumRows <= 0) ) THEN
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

DO i = 1, NPts
  first_col = PtPosn(1,i)
  first_row = PtPosn(2,i)
  IF ( (first_col == 0) .OR. (first_row == 0) ) THEN
    wspd = 0.0
  ELSE IF( (first_col*first_row > 0) .AND. &
           (first_col <= NumCols)    .AND. &
           (first_row <= NumRows) ) THEN
    wspd = WindSpd % RData(first_col,first_row)
  ELSE
    ErrorStatus = StatusWarning
    WRITE(6,*) "Counter not in range ", first_col,first_row
    EXIT
  END IF

  col = AdvPtPosn(1,i)
  row = AdvPtPosn(2,i)
  IF ( (col == 0) .OR. (row == 0) ) THEN
    wdir = 0.0
  ELSE IF( (col*row > 0)    .AND. &
           (col <= NumCols) .AND. &
           (row <= NumRows) ) THEN
    wdir = WindDir % RData(col,row)
  ELSE
    ErrorStatus = StatusWarning
    WRITE(6,*) "Counter not in range ", col,row
    EXIT
  END IF

  ! Calculate the wind modification if needed. Point 5 in Appendix.
  ! If the wind speed is small, the arrays are not used (reset to zero)
  IF ( (wspd < MinWSpeed) ) THEN
    new_col = 0
    new_row = 0
  ELSE

    IF ( wdir < 0.0 ) THEN  ! wdir is the vector (not met.) direction
      wdir = wdir + 360.0
    END IF

    ! Determine the comparison angles alpha1 and alpha2
    ! Angles alpha1 and alpha2 equal 27.0 and 63.0 for a square
    ! grid on the equator.
    latitude = MWTPred % Hdr % ZerothLat +  &
               MWTPred % Hdr % LatInt * row
    cosang = COS(latitude*Pi_Over_180)

    ! Determine the angles alpha1 and alpha2
    mod_deltalat = ABS( MWTPred % Hdr % LatInt )
    alpha1 = ATAN( MWTPred % Hdr % LonInt * cosang/(2*mod_deltalat)) &
                                                          / Pi_Over_180
    alpha2 = ATAN(2* MWTPred % Hdr % LonInt *cosang/mod_deltalat)    &
                                                          / Pi_Over_180

    ! Determine the adjacent point in longitude direction
    new_col = col                      ! Longitude the same
    IF ( (wdir > alpha1) .AND. &
         (wdir < (180.0 - alpha1)) ) THEN  ! NW, W or SW
      new_col = col + 1                ! Longitude increased
      IF ( col == NumCols ) THEN
        new_col = 1                    ! Wrap around ***NOT MPP***
      END IF

    ELSE IF ( (wdir > (180.0 + alpha1)) .AND. &
              (wdir < (360.0 - alpha1)) ) THEN  ! NE, E or SE
      new_col = col - 1                ! Longitude decreased
      IF ( col == 1 ) THEN
        new_col = NumCols              ! Wrap around  ***NOT MPP***
      END IF

    END IF

    ! Determine the adjacent point in latitude direction
    new_row = row                      !  Latitude the same
    IF ( (wdir > (360.0 - alpha2)) .OR.  &
         (wdir < alpha2) ) THEN        !NW, N or NE
      IF ( row > 1 ) THEN
        new_row = row - 1              ! Latitude decreased
      END IF

    ELSE IF ( (wdir > (180.0 - alpha2)) .AND. &
              (wdir < (180.0 + alpha2)) ) THEN  !SW, S or SE
      IF ( row < NumRows ) THEN
        new_row = row + 1              ! Latitude increased
      END IF

    END IF
  END IF

  ! Reassign value of mwt_pred if point has been set up.
  ! Zero points imply do nothing
  IF ( (new_col > 0) .AND. (first_col > 0) .AND. &
       (new_row > 0) .AND. (first_row > 0) ) THEN
    IF ( MWTPred % RData(new_col,new_row) < CutOff ) THEN
      MWTPred % RData(new_col,new_row) = &
          MaxStress % RData(first_col,first_row)
    END IF
  END IF

  AdvPtPosn(1:2,i) = (/new_col, new_row/)

END DO

9999 CONTINUE

END SUBROUTINE AdvectMWS

