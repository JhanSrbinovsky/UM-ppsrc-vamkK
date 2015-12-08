! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to calculate Mountain Wave Predictor diagnostic
! Note WindDir is wind direction in vector sense, not meteorological

SUBROUTINE MtnStress( NumLevs,      &  ! in
                      WindU200,     &  ! in
                      WindV200,     &  ! in
                      UStress,      &  ! inout
                      VStress,      &  ! inout
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

USE IO_mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
USE FldCodes_Mod, ONLY:   &
  ST_MWT,   MO8_MWT,      &
  PP_MWT,   VC_Turb,      &
  LV_Special
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
TYPE(PP_Field_type), INTENT(IN) :: WindU200  ! wind U-component at 200mb
TYPE(PP_Field_type), INTENT(IN) :: WindV200  ! wind V-component at 200mb
TYPE(PP_Field_type), INTENT(IN) :: UStress(NumLevs)
TYPE(PP_Field_type), INTENT(IN) :: VStress(NumLevs)

TYPE(PP_Field_type), INTENT(INOUT) :: MWTPred
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "MtnStress"
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
INTEGER, PARAMETER :: NumContours = 3
REAL,    PARAMETER :: Contours(NumContours) =  &
                          (/0.007,0.065,0.25/) ! The contour thresholds

! Local Variables:
INTEGER :: i, j, k        ! DO loop variables
INTEGER :: NumRows        ! The number of latitudes = header(18)
INTEGER :: NumCols        ! The number of longitudes = header(19)
INTEGER :: NumAbvCutOff   ! Number of points above Contours
REAL :: CutOff            ! One of the Contours values
REAL :: MinWSpeed         ! Wind threshold for advection

INTEGER, ALLOCATABLE :: PtPosn(:,:)    ! Lat Lon array of
                                       ! stress > cutoff
INTEGER, ALLOCATABLE :: AdvPtPosn(:,:) ! advected LL array of
                                       ! stress > cutoff

TYPE(PP_Field_type) :: WindF200    ! wind speed at 200mb
TYPE(PP_Field_type) :: WindD200    ! wind direction at 200mb
TYPE(PP_Field_type) :: MaxStress   ! Maximum of gravity wave stress
TYPE(PP_Field_type) :: Stress      ! Magnitude of gravity wave stress

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( MWTPred % RData ) ) THEN
  DEALLOCATE( MWTPred % RData )
END IF
MWTPred % Hdr = WindU200 % Hdr
MWTPred % Hdr % PPCode   =  PP_MWT
MWTPred % Hdr % LBVC     =  VC_Turb
MWTPred % Hdr % MO8Type  = MO8_MWT
MWTPred % Hdr % MO8Level =  LV_Special
MWTPred % Hdr % STCode   =  ST_MWT
MWTPred % Hdr % RLevel   = 0.0
MWTPred % Hdr % RefLevel = 0.0
MWTPred % Hdr % BHLEV    = 0.0
MWTPred % Hdr % BHRLEV   = 0.0
MWTPred % Hdr % BULEV    = 0.0
MWTPred % Hdr % BHULEV   = 0.0
MWTPred % Hdr % BMDI     = RMDI
ALLOCATE( MWTPred % RData(MWTPred % Hdr % NumCols, &
                          MWTPred % Hdr % NumRows) )

! Initialise the local variables
NumCols = MWTPred % Hdr % NumCols
NumRows = MWTPred % Hdr % NumRows
NULLIFY( Stress % RData )
NULLIFY( MaxStress % RData )
NULLIFY( WindF200 % RData )
NULLIFY( WindD200 % RData )

! DEPENDS ON: vecmag
CALL VecMag( WindU200, WindV200, WindF200, ErrorStatus )
! DEPENDS ON: vecdir
CALL VecDir( WindU200, WindV200, WindD200, ErrorStatus )

! Determine the maximum magnitude of the gravity wave stress
! DEPENDS ON: vecmag
CALL VecMag( UStress(1), VStress(1), MaxStress, ErrorStatus )
DO k = 2,NumLevs
! DEPENDS ON: vecmag
  CALL VecMag( UStress(k), VStress(k), Stress, ErrorStatus )
  WHERE ( Stress % RData > MaxStress % RData ) &
          MaxStress % RData = Stress % RData
END DO

DEALLOCATE( Stress % RData )
NULLIFY( Stress % RData )
ALLOCATE( PtPosn(2,NumCols*NumRows) )

MWTPred % RData(:,:) = MaxStress % RData(:,:)

DO k = 1, NumContours     ! loop over the contour threshold values

  !  Determine the arrays of values above the thresholds.
  CutOff = Contours(k)
  NumAbvCutOff = 0       ! Initialise to 0
  PtPosn(:,:) = 0
  DO j = 1, NumRows
    DO i = 1, NumCols
      IF ( MaxStress % RData(i,j) > CutOff ) THEN
        ! Increase number of points by 1 and store in PtPosn
        NumAbvCutOff = NumAbvCutOff + 1
        PtPosn( 1:2, NumAbvCutOff ) = (/i,j/)
      END IF
    END DO
  END DO

  ! If NumAbvCutOff > 0, the 'new' mwt points need calculating.
  IF ( NumAbvCutOff > 0 ) THEN

    ALLOCATE( AdvPtPosn(2,NumAbvCutOff) )
    AdvPtPosn(:,:) = PtPosn(:,1:NumAbvCutOff)

    ! Determine the array of next points. Point 2 in Appendix.
    MinWSpeed = 0.0
! DEPENDS ON: advectmws
    CALL AdvectMWS( NumAbvCutOff, CutOff,            &
                    MinWSpeed,                       &
                    PtPosn    (1:2,1:NumAbvCutOff),  &
                    MaxStress,                       &
                    WindF200,     WindD200,          &
                    AdvPtPosn (1:2,1:NumAbvCutOff),  &
                    MWTPred,                         &
                    ErrorStatus )

    !*** need to swap halos of AdvPtPosn if MPP***
    ! Determine the second array of next points. Point 4 in Appendix.
    MinWSpeed = 30.0  ! This is the wind speed required to advect a
                      ! distance of 2 grid squares.
! DEPENDS ON: advectmws
    CALL AdvectMWS( NumAbvCutOff, CutOff,            &
                    MinWSpeed,                       &
                    PtPosn    (1:2,1:NumAbvCutOff),  &
                    MaxStress,                       &
                    WindF200,     WindD200,          &
                    AdvPtPosn (1:2,1:NumAbvCutOff),  &
                    MWTPred,                         &
                    ErrorStatus )

    ! If an extra point is required downwind for extra resolution
    ! or strong flow, calculate an extra point. Point 6 in Appendix.
    ! Possibly do this by calling AdvectMWS from within a loop.
    ! Calculate a MinWindSpd threshold from the grid resolution and pass
    ! it into the subroutine.

    DEALLOCATE(AdvPtPosn)
  END IF
END DO

DEALLOCATE( WindF200 % RData )
DEALLOCATE( WindD200 % RData )
DEALLOCATE( MaxStress % RData )
DEALLOCATE( PtPosn )

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE MtnStress

!=======================================================================


