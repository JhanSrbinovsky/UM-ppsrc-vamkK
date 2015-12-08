! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to calculate weather diagnostics


!=======================================================================

SUBROUTINE WXCode( DynRainRate,  &  ! in
                   ConRainRate,  &  ! in
                   DynSnowRate,  &  ! in
                   ConSnowRate,  &  ! in
                   SurfTemp,     &  ! in
                   DynRainAcc,   &  ! in
                   ConRainAcc,   &  ! in
                   DynSnowAcc,   &  ! in
                   ConSnowAcc,   &  ! in
                   FogProb,      &  ! in
                   Visib,        &  ! in
                   CldCover,     &  ! in
                   WXField,      &  ! inout
                   ErrorStatus )    ! inout

! Description:
!   Calculate present weather code.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE conversions_mod, ONLY: zerodegc

USE IO_Mod, ONLY:         &
  PP_Field_type,          &
  PP_Header_type
USE Err_Mod, ONLY:        &
  StatusOK
USE FldCodes_Mod, ONLY:   &
  ST_WXCode, MO8_WXCode, PP_WXCode
IMPLICIT None

! Subroutine arguments:
TYPE(PP_Field_type), INTENT(IN) :: DynRainRate  ! Dynamic Rain Rate
TYPE(PP_Field_type), INTENT(IN) :: ConRainRate  ! Convective Rain Rate
TYPE(PP_Field_type), INTENT(IN) :: DynSnowRate  ! Dynamic Snow Rate
TYPE(PP_Field_type), INTENT(IN) :: ConSnowRate  ! Convective Snow Rate
TYPE(PP_Field_type), INTENT(IN) :: SurfTemp     ! Surface Temperature
TYPE(PP_Field_type), INTENT(IN) :: DynRainAcc   ! Dyn Rain Accumulation
TYPE(PP_Field_type), INTENT(IN) :: ConRainAcc   ! Conv Rain Accumulation
TYPE(PP_Field_type), INTENT(IN) :: DynSnowAcc   ! Dyn Snow Accumulation
TYPE(PP_Field_type), INTENT(IN) :: ConSnowAcc   ! Conv Snow Accumulation
TYPE(PP_Field_type), INTENT(IN) :: FogProb      ! Fog Probability
TYPE(PP_Field_type), INTENT(IN) :: Visib        ! Visibility
TYPE(PP_Field_type), INTENT(IN) :: CldCover     ! Cloud Cover (fraction)

TYPE(PP_Field_type), INTENT(INOUT) :: WXField   ! Present weather code
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "WXCode"

REAL, PARAMETER :: Hr = 3600.

! Local variables:
INTEGER :: i, j
REAL :: DynRate (SurfTemp%Hdr%NumCols, SurfTemp%Hdr%NumRows)
REAL :: ConRate (SurfTemp%Hdr%NumCols, SurfTemp%Hdr%NumRows)
REAL :: DynAcc  (SurfTemp%Hdr%NumCols, SurfTemp%Hdr%NumRows)
LOGICAL :: LRain(SurfTemp%Hdr%NumCols, SurfTemp%Hdr%NumRows)
LOGICAL :: LSnow(SurfTemp%Hdr%NumCols, SurfTemp%Hdr%NumRows)
LOGICAL :: LConv(SurfTemp%Hdr%NumCols, SurfTemp%Hdr%NumRows)
LOGICAL :: LDyn (SurfTemp%Hdr%NumCols, SurfTemp%Hdr%NumRows)

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Initialise output array
IF ( ASSOCIATED( WXField % RData ) ) THEN
  DEALLOCATE( WXField % RData )
END IF
WXField % Hdr = DynRainRate % Hdr
WXField % Hdr % PPcode   =  PP_WXCode
WXField % Hdr % MO8Type  = MO8_WXCode
WXField % Hdr % MO8Level = 8888
WXField % Hdr % STCode   =  ST_WXCode
ALLOCATE( WXField % RData(DynRainRate % Hdr % NumCols, &
                          DynRainRate % Hdr % NumRows) )
WXField % RData(:,:) = 0.0

! Set Past Weather Codes
WHERE ( DynRainAcc % RData > 0.015 )
  WXField % RData = 20.0          ! Past Drizzle
END WHERE
WHERE ( DynRainAcc % RData > 0.030 )
  WXField % RData = 21.0          ! Past Rain
END WHERE
WHERE ( DynSnowAcc % RData > 0.015 )
  WXField % RData = 22.0          ! Past Snow
END WHERE
WHERE ( ConRainAcc % RData > 0.015 )
  WXField % RData = 25.0          ! Past Rain Shower
END WHERE
WHERE ( ConSnowAcc % RData > 0.015 )
  WXField % RData = 26.0          ! Past Snow Shower
END WHERE
WHERE ( (ConSnowAcc % RData > 0.015) .AND. &
        (SurfTemp % RData > 2.0+ZeroDegC) )
  WXField % RData = 25.0          ! Past Rain Shower
END WHERE

! Set Visibility Codes
WHERE ( (Visib % RData < 5.0) .AND. &
        (WXField % RData < 10.0) )
  WXField % RData = 10.0          ! Mist
END WHERE
WHERE ( FogProb % RData > 0.3 )
  WXField % RData = 45.0          ! Fog
END WHERE
WHERE ( (FogProb % RData > 0.3) .AND. &
        (CldCover % RData <= 0.75 ) )
  WXField % RData = 44.0          ! Fog - Sky Visible
END WHERE
WHERE ( (WXField % RData == 45.0) .AND. &
        (SurfTemp % RData < ZeroDegC) )
  WXField % RData = 49.0          ! Freezing Fog
END WHERE
WHERE ( (WXField % RData == 44.0) .AND. &
        (SurfTemp % RData < ZeroDegC) )
  WXField % RData = 48.0
END WHERE

! Set Precipitation Codes
! Rates are given in mm/hr, but data is mm/s, so multiply by 3600
! Temps are given in Celsius, but data is in Kelvin so add ZeroDegC
DynRate = (DynRainRate % RData + DynSnowRate % RData) * Hr
ConRate = (ConRainRate % RData + ConSnowRate % RData) * Hr
DynAcc  = DynRainAcc % RData  + DynSnowAcc % RData

LRain = ( (DynRainRate % RData > 0.03/Hr) .OR.  &
          (ConRainRate % RData > 0.10/Hr) )
LSnow   = ( (DynSnowRate % RData > 0.03/Hr) .OR.  &
           ((ConSnowRate % RData > 0.10/Hr) .AND.  &
            (SurfTemp % RData <= 2.0+ZeroDegC)) )  ! Isn't 2C a bit low?
LConv   = ( (ConRate > 0.1) .OR. &
            ((LRain .OR. LSnow) .AND. (CldCover % RData <= 0.75)) )
LDyn    = ( DynRate > 0.03 )

DO j = 1, WXField % Hdr % NumRows
  DO i = 1, WXField % Hdr % NumCols

    IF ( LConv(i,j) ) THEN
      IF ( LDyn(i,j) ) THEN                ! Mixed Convective/Dynamic
        IF ( DynRate(i,j) < 0.5 ) THEN
          IF ( LRain(i,j) ) THEN
            WXField % RData(i,j)   = 60.0  ! Slt Intermittent Rain
            IF ( SurfTemp % RData(i,j) < ZeroDegC ) THEN
              WXField % RData(i,j) = 66.0  ! Slt Freezing Rain
            END IF
            IF ( LSnow(i,j) ) THEN
              WXField % RData(i,j) = 68.0  ! Slt Rain & Snow
            END IF
          ELSE
            WXField % RData(i,j)   = 70.0  ! Slt Intermittent Snow
          END IF
        ELSE IF ( DynRate(i,j) < 4.0 ) THEN
          IF ( LRain(i,j) ) THEN
            WXField % RData(i,j)   = 62.0  ! Mod Intermittent Rain
            IF ( SurfTemp % RData(i,j) < ZeroDegC ) THEN
              WXField % RData(i,j) = 67.0  ! Mod Freezing Rain
            END IF
            IF ( LSnow(i,j) ) THEN
              WXField % RData(i,j) = 69.0  ! Mod/Heavy Rain & Snow
            END IF
          ELSE
            WXField % RData(i,j)   = 72.0  ! Mod Intermittent Snow
          END IF
        ELSE
          IF ( LRain(i,j) ) THEN
            WXField % RData(i,j)   = 64.0  ! Heavy Intermittent Rain
            IF ( SurfTemp % RData(i,j) < ZeroDegC ) THEN
              WXField % RData(i,j) = 67.0  ! Mod Freezing Rain
            END IF
            IF ( LSnow(i,j) ) THEN
              WXField % RData(i,j) = 69.0  ! Mod/Heavy Rain & Snow
            END IF
          ELSE
            WXField % RData(i,j)   = 74.0  ! Heavy Intermittent Snow
          END IF
        END IF
      ELSE                                 ! Convective Only
        IF ( ConRate(i,j) >= 2.0 ) THEN
          IF ( LRain(i,j) ) THEN
            WXField % RData(i,j)   = 81.0  ! Mod/Heavy Rain Shower
            IF ( LSnow(i,j) ) THEN
              WXField % RData(i,j) = 84.0  ! Mod/Heavy Rain &
                                           ! Snow Shower
            END IF
          ELSE
            WXField % RData(i,j)   = 86.0  ! Mod/Heavy Snow Shower
          END IF
        ELSE
          IF ( LRain(i,j) ) THEN
            WXField % RData(i,j)   = 80.0  ! Rain Shower
            IF ( LSnow(i,j) ) THEN
              WXField % RData(i,j) = 83.0  ! Rain & Snow Shower
            END IF
          ELSE
            WXField % RData(i,j)   = 85.0  ! Snow Shower
          END IF
        END IF
      END IF
    ELSE IF ( LDyn(i,j) ) THEN             ! Dynamic Only
      IF ( DynRate(i,j) >= 4.0 ) THEN
        IF ( LRain(i,j) ) THEN
          WXField % RData(i,j)     = 65.0  ! Heavy Rain
          IF ( DynAcc(i,j) < 0.05 ) THEN
            WXField % RData(i,j)   = 64.0  ! Heavy Intermittent Rain
          END IF
          IF ( SurfTemp % RData(i,j) < ZeroDegC ) THEN
            WXField % RData(i,j)   = 67.0  ! Mod Freezing Rain
          END IF
          IF ( LSnow(i,j) ) THEN
            WXField % RData(i,j)   = 69.0  ! Heavy Rain & Snow
          END IF
        ELSE
          WXField % RData(i,j)     = 75.0  ! Heavy Snow
          IF ( DynAcc(i,j) < 0.05 ) THEN
            WXField % RData(i,j)   = 74.0  ! Heavy Intermittent Snow
          END IF
        END IF
      ELSE IF ( DynRate(i,j) > 0.5 ) THEN
        IF ( LRain(i,j) ) THEN
          WXField % RData(i,j)     = 63.0  ! Mod Rain
          IF ( DynAcc(i,j) < 0.05 ) THEN
            WXField % RData(i,j)   = 62.0  ! Mod Intermittent Rain
          END IF
          IF ( SurfTemp % RData(i,j) < ZeroDegC ) THEN
            WXField % RData(i,j)   = 67.0  ! Mod Freezing Rain
          END IF
          IF ( LSnow(i,j) ) THEN
            WXField % RData(i,j)   = 69.0  ! Mod Rain & Snow
          END IF
        ELSE
          WXField % RData(i,j)     = 73.0  ! Mod Snow
          IF ( DynAcc(i,j) < 0.05 ) THEN
            WXField % RData(i,j)   = 72.0  ! Mod Intermittent Snow
          END IF
        END IF
      ELSE IF ( DynRate(i,j) < 0.1 ) THEN
        WXField % RData(i,j) = 51.0        ! Drizzle
        IF ( SurfTemp % RData(i,j) < ZeroDegC ) THEN
          WXField % RData(i,j)     = 56.0  ! Slt Freezing Drizzle
        END IF
        IF ( LSnow(i,j) ) THEN
          WXField % RData(i,j)     = 77.0  ! Snow Grains
        END IF
      ELSE
        IF ( LRain(i,j) ) THEN
          WXField % RData(i,j)     = 61.0  ! Light Rain
          IF ( DynAcc(i,j) < 0.05 ) THEN
            WXField % RData(i,j)   = 60.0  ! Light Intermittent Rain
          END IF
          IF ( SurfTemp % RData(i,j) < ZeroDegC ) THEN
            WXField % RData(i,j)   = 66.0  ! Slt. Freezing Rain
          END IF
          IF ( LSnow(i,j) ) THEN
            WXField % RData(i,j)   = 68.0  ! Light Rain & Snow
          END IF
        ELSE
          WXField % RData(i,j)     = 71.0  ! Light Snow
          IF ( DynAcc(i,j) < 0.05 ) THEN
            WXField % RData(i,j)   = 70.0  ! Light Intermittent Snow
          END IF
        END IF
      END IF
    END IF

  END DO
END DO

WHERE ( (ConRate > 15.0) .OR. (DynRate > 15.0) )
  WXField % RData = 95.0          ! Thunderstorm with rain/snow
END WHERE
WHERE ( (ConRate > 20.0) .OR. (DynRate > 20.0) )
  WXField % RData = 96.0          ! Thunderstorm with hail
END WHERE

9999 CONTINUE

END SUBROUTINE WXCode

