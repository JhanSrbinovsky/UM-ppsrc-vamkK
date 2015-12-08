! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate Max Wind value and height

SUBROUTINE MaxWind( NumLevs,      &  ! in
                    UFields,      &  ! in
                    VFields,      &  ! in
                    PFields,      &  ! in
                    MaxWindU,     &  ! inout
                    MaxWindV,     &  ! inout
                    MaxWindP,     &  ! inout
                    ErrorStatus )    ! inout

! Description:
!
! Method:
!   Old version uses U&V on rho levels, and pressure on theta levels to
!   interpolate between.  To save space we will interpolate between rho
!   levels.  This may mean we use a different number of levels, and we
!   may want to widen the 100-700 hPa range.
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
USE FldCodes_Mod, ONLY:         &
  ST_MaxWU,  MO8_MaxWU,  PP_U,  &
  ST_MaxWV,  MO8_MaxWV,  PP_V,  &
  ST_MaxWP,  MO8_MaxWP,  PP_P,  &
  VC_Upper,   VC_Lower,         &
  VC_MaxWind, LV_Special,       &
  ST_MWBase,  ST_MWTop,         &
  MO8_MxWBase, MO8_MxWTop

USE PrintStatus_mod
IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
TYPE(PP_Field_type), INTENT(IN) :: UFields(NumLevs) ! U-wind on B-grid
TYPE(PP_Field_type), INTENT(IN) :: VFields(NumLevs) ! V-wind on B-grid
! Pressure on B-grid U/V points on rho levels:
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs)

TYPE(PP_Field_type), INTENT(INOUT) :: MaxWindU  ! U comp of max wind
TYPE(PP_Field_type), INTENT(INOUT) :: MaxWindV  ! V comp of max wind
TYPE(PP_Field_type), INTENT(INOUT) :: MaxWindP  ! P at max wind level

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "MaxWind"

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
INTEGER, PARAMETER :: kmax=5            ! total levels needed for spline
INTEGER, PARAMETER :: khalf=(kmax+1)/2  ! number of levels each side
INTEGER, PARAMETER :: ninc=16           ! number of increments used

! Local variables:
INTEGER :: i, j, k
INTEGER :: LevBtm, LevTop
INTEGER :: Level_Below, Level_Above
INTEGER :: MaxWLev(PFields(1)%Hdr%NumCols, PFields(1)%Hdr%NumRows)
REAL :: U, V, WS
REAL :: Min_Diff
REAL :: MaxW(PFields(1)%Hdr%NumCols, PFields(1)%Hdr%NumRows)
TYPE(PP_Field_Type) :: WindSpd

! Variables involving spline interpolation

! Local Variables for the vectorised version
REAL    :: Uinc(2*ninc)
REAL    :: Vinc(2*ninc)
REAL    :: Pinc(2*ninc) ! u,v and p at increments
REAL    :: Winc(2*ninc)
INTEGER :: MaxWinc(1)



! End of Header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( MaxWindU % RData ) ) THEN
  DEALLOCATE( MaxWindU % RData )
END IF
IF ( ASSOCIATED( MaxWindV % RData ) ) THEN
  DEALLOCATE( MaxWindV % RData )
END IF
IF ( ASSOCIATED( MaxWindP % RData ) ) THEN
  DEALLOCATE( MaxWindP % RData )
END IF

MaxWindU % Hdr = UFields(1) % Hdr
MaxWindU % Hdr % LBVC     = VC_MaxWind
MaxWindU % Hdr % MO8Level = LV_Special
MaxWindU % Hdr % BULEV    = 0.0
MaxWindU % Hdr % BHULEV   = 0.0
MaxWindU % Hdr % RLevel   = 0.0
MaxWindU % Hdr % RefLevel = 0.0
MaxWindU % Hdr % BHLEV    = 0.0
MaxWindU % Hdr % BHRLEV   = 0.0
MaxWindU % Hdr % BMDI     = RMDI
MaxWindV % Hdr = MaxWindU % Hdr
MaxWindP % Hdr = MaxWindU % Hdr
MaxWindU % Hdr % PPCode   =  PP_U
MaxWindU % Hdr % MO8Type  = MO8_MaxWU
MaxWindU % Hdr % STCode   =  ST_MaxWU
MaxWindV % Hdr % PPCode   =  PP_V
MaxWindV % Hdr % MO8Type  = MO8_MaxWV
MaxWindV % Hdr % STCode   =  ST_MaxWV
MaxWindP % Hdr % PPCode   =  PP_P
MaxWindP % Hdr % MO8Type  = MO8_MaxWP
MaxWindP % Hdr % STcode   =  ST_MaxWP
ALLOCATE( MaxWindU % RData(MaxWindU % Hdr % NumCols, &
                           MaxWindU % Hdr % NumRows) )
ALLOCATE( MaxWindV % RData(MaxWindV % Hdr % NumCols, &
                           MaxWindV % Hdr % NumRows) )
ALLOCATE( MaxWindP % RData(MaxWindP % Hdr % NumCols, &
                           MaxWindP % Hdr % NumRows) )

LevBtm = khalf
LevTop = NumLevs+1 - khalf

NULLIFY( WindSpd % RData )

! Calculate 1st guess of max wind and level
MaxWLev(:,:) = 0
MaxW   (:,:) = RMDI   ! Large and -ve
DO k = LevBtm, LevTop
! DEPENDS ON: vecmag
  CALL VecMag( UFields(k), VFields(k), WindSpd, ErrorStatus )
  ! If pressure is in the correct range and greater than previous:


  DO j = 1, PFields(1) % Hdr % NumRows
    DO i = 1, PFields(1) % Hdr % NumCols
      IF (WindSpd % RData(i,j) /= WindSpd % Hdr % BMDI ) THEN
        IF (PFields(k+1) % RData(i,j) >= 10000.0) THEN
          IF (PFields(k-1) % RData(i,j) <  70000.0) THEN
            IF (WindSpd % RData(i,j) >= MaxW(i,j)) THEN
              MaxWLev(i,j) = k
              MaxW(i,j) = WindSpd % RData(i,j)
            END IF
          END IF
        END IF
      END IF
    END DO
  END DO


END DO

IF ( ErrorStatus /= StatusOK ) THEN
  GO TO 9999
END IF

DO i = 1,MaxWindU % Hdr % NumCols
  DO j = 1,MaxWindU % Hdr % NumRows
    IF ( MaxWLev(i,j) == 0 ) THEN
      ! If U or V is missing MaxWLev will contain 0
      MaxWindU % RData(i,j) = RMDI
      MaxWindV % RData(i,j) = RMDI
      MaxWindP % RData(i,j) = RMDI
    ELSE

! DEPENDS ON: maxwindspline
      Call MaxWindSpline (NumLevs, i, j,                         &
                          UFields, VFields, PFields, MaxWLev,    &
                          Uinc, Vinc, Pinc)


     ! Look along increments for new values for max wind and level
     ! This can be speed squared to remove need to use costly sqrt
      Winc(:) = Uinc(:)**2 + Vinc(:)**2
      MaxWinc = MAXLOC( Winc )
      MaxWindU % RData(i,j) = Uinc(MaxWinc(1))
      MaxWindV % RData(i,j) = Vinc(MaxWinc(1))
      MaxWindP % RData(i,j) = Pinc(MaxWinc(1))
    END IF

  END DO
END DO

DEALLOCATE ( WindSpd % RData )

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE MaxWind

