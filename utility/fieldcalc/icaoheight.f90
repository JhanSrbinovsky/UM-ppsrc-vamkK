! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to convert pressure fields to ICAO height

SUBROUTINE ICAOHeight( PField,       &  ! in
                       IHField,      &  ! inout
                       ErrorStatus )    ! inout

! Description:
!   Convert pressure (Pa) to height (kft) using ICAO standard atmosphere
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE atmos_constants_mod, ONLY: r

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
USE FldCodes_Mod, ONLY:              &
  ST_CClBP,  ST_CClBI,  MO8_CClBI,   &
  ST_CClTP,  ST_CClTI,  MO8_CClTI,   &
  ST_LCClBP, ST_LCClBI, MO8_LCClBI,  &
  ST_LCClTP, ST_LCClTI, MO8_LCClTI,  &
  ST_MaxWP,  ST_MaxWI,  MO8_MaxWI,   &
  ST_MWBase, MO8_MxWBase,            &
  ST_MWTop,  MO8_MxWTop,             &
  ST_Iso70P, ST_Iso70I, MO8_Iso70I,  &
  ST_Iso20P, ST_Iso20I, MO8_Iso20I,  &
  ST_FreezP, ST_FreezI, MO8_FreezI,  &
  ST_TropP,  ST_TropI,  MO8_TropI,   &
  ST_TropP_GRIB2,                    &
  ST_TropI_GRIB2, MO8_TropI_GRIB2,   &
  PP_ICAOHt,                         &
  ST_P_CBB,   MO8_P_CBB,  PP_P_CBB,  &
  ST_P_CBT,   MO8_P_CBT,  PP_P_CBT,  &
  ST_I_CBB,   MO8_I_CBB,  PP_I_CBB,  &
  ST_I_CBT,   MO8_I_CBT,  PP_I_CBT,  &
  ST_P_ECBB,  MO8_P_ECBB, PP_P_ECBB, &
  ST_P_ECBT,  MO8_P_ECBT, PP_P_ECBT, &
  ST_I_ECBB,  MO8_I_ECBB, PP_I_ECBB, &
  ST_I_ECBT,  MO8_I_ECBT, PP_I_ECBT

USE earth_constants_mod, ONLY: g

IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: PField       !P field for conversion

TYPE(PP_Field_type), INTENT(INOUT) :: IHField   !ICAO height in kft
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "ICAOHeight"
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
REAL, PARAMETER :: FT2M = 0.3048
REAL, PARAMETER :: G_over_R = G / R
REAL, PARAMETER :: mtokft = 1/(FT2M * 1000.0)
REAL, PARAMETER :: Lapse_RateL = 6.5E-03  ! For levels below 11,000 gpm
REAL, PARAMETER :: Lapse_RateU = -1.0E-03 ! For levels above 11,000 gpm
REAL, PARAMETER :: Press_Bot = 101325.    ! ICAO std: surface pressure
REAL, PARAMETER :: Press_Mid = 22632.     !      pressure @ 11,000 gpm
REAL, PARAMETER :: Press_Top = 5474.87    !      pressure @ 20,000 gpm
REAL, PARAMETER :: Temp_Bot = 288.15      ! Surface temperature
REAL, PARAMETER :: Temp_Top = 216.65      ! Temperature of isotherm
REAL, PARAMETER :: Gpm1 = 11000.0  ! Ht limit (gpm) for std lower
                                   ! lapse rate
REAL, PARAMETER :: Gpm2 = 20000.0  ! Ht (gpm) of top of isothermal layer
REAL, PARAMETER :: ZP1 = Lapse_RateL/G_over_R ! Exponents used for
REAL, PARAMETER :: ZP2 = Lapse_RateU/G_over_R ! calculation

! Local Variables:
INTEGER :: i,j      ! loop counters
REAL :: Pressure    ! Local pressure

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( IHField % RData ) ) THEN
  DEALLOCATE( IHField % RData )
END IF
IHField % Hdr = PField % Hdr
IHField % Hdr % PPCode = PP_ICAOHt
IF      ( PField % Hdr % STCode == ST_CClBP  ) THEN
  ! Convective Cloud Base Pressure
  IHField % Hdr % MO8Type = MO8_CClBI
  IHField % Hdr % STCode  =  ST_CClBI
ELSE IF ( PField % Hdr % STCode == ST_CClTP  ) THEN
  ! Convective Cloud Top Pressure
  IHField % Hdr % MO8Type = MO8_CClTI
  IHField % Hdr % STCode  =  ST_CClTI
ELSE IF ( PField % Hdr % STCode == ST_LCClBP ) THEN
  ! Lowest Convective Cloud Base Pressure
  IHField % Hdr % MO8Type = MO8_LCClBI
  IHField % Hdr % STCode  =  ST_LCClBI
ELSE IF ( PField % Hdr % STCode == ST_LCClTP ) THEN
  ! Lowest Convective Cloud Top Pressure
  IHField % Hdr % MO8Type = MO8_LCClTI
  IHField % Hdr % STCode  =  ST_LCClTI
ELSE IF ( PField % Hdr % STCode == ST_P_CBB  ) THEN
  ! Pressure at Cb Base
  IHField % Hdr % PPCode  = PP_I_CBB
  IHField % Hdr % MO8Type = MO8_I_CBB
  IHField % Hdr % STCode  =  ST_I_CBB
ELSE IF ( PField % Hdr % STCode == ST_P_CBT  ) THEN
  ! Pressure at Cb Top
  IHField % Hdr % PPCode  = PP_I_CBT
  IHField % Hdr % MO8Type = MO8_I_CBT
  IHField % Hdr % STCode  =  ST_I_CBT
ELSE IF ( PField % Hdr % STCode == ST_P_ECBB  ) THEN
  ! Pressure at Embedded Cb Base
  IHField % Hdr % PPCode  = PP_I_ECBB
  IHField % Hdr % MO8Type = MO8_I_ECBB
  IHField % Hdr % STCode  =  ST_I_ECBB
ELSE IF ( PField % Hdr % STCode == ST_P_ECBT  ) THEN
  ! Pressure at Embedded Cb Base
  IHField % Hdr % PPCode  = PP_I_ECBT
  IHField % Hdr % MO8Type = MO8_I_ECBT
  IHField % Hdr % STCode  =  ST_I_ECBT
ELSE IF ( PField % Hdr % STCode == ST_MaxWP  ) THEN
  ! MaxWind Pressure
  IHField % Hdr % MO8Type = MO8_MaxWI
  IHField % Hdr % STCode  =  ST_MaxWI
ELSE IF ( PField % Hdr % STCode == ST_MWBase  ) THEN
  ! MaxWind Base Pressure
  IHField % Hdr % MO8Type = MO8_MxWBase
  IHField % Hdr % STCode  =  ST_MWBase
ELSE IF ( PField % Hdr % STCode == ST_MWTop  ) THEN
  ! MaxWind Top Pressure
  IHField % Hdr % MO8Type = MO8_MxWTop
  IHField % Hdr % STCode  =  ST_MWTop
ELSE IF ( PField % Hdr % STCode == ST_Iso70P ) THEN
  ! -70C Isotherm Pressure
  IHField % Hdr % MO8Type = MO8_Iso70I
  IHField % Hdr % STCode  =  ST_Iso70I
ELSE IF ( PField % Hdr % STCode == ST_Iso20P ) THEN
  ! -20C Isotherm Pressure
  IHField % Hdr % MO8Type = MO8_Iso20I
  IHField % Hdr % STCode  =  ST_Iso20I
ELSE IF ( PField % Hdr % STCode == ST_FreezP ) THEN
  ! Freezing Level Pressure
  IHField % Hdr % MO8Type = MO8_FreezI
  IHField % Hdr % STCode  =  ST_FreezI
ELSE IF ( PField % Hdr % STCode == ST_TropP  ) THEN
  ! Tropopause Pressure
  IHField % Hdr % MO8Type = MO8_TropI
  IHField % Hdr % STCode  =  ST_TropI
ELSE IF ( PField % Hdr % STCode == ST_TropP_GRIB2 ) THEN
  ! Tropopause Pressure for GRIB2
  IHField % Hdr % MO8Type = MO8_TropI_GRIB2
  IHField % Hdr % STCode  =  ST_TropI_GRIB2
ELSE
  ! Unrecognised
  IHField % Hdr % MO8Type = IMDI
  IHField % Hdr % STCode  = IMDI
END IF
IHField % Hdr % BMDI   = RMDI
ALLOCATE( IHField % RData(IHField % Hdr % NumCols, &
                          IHField % Hdr % NumRows) )

DO j = 1,IHField % Hdr % NumRows
  DO i = 1,IHField % Hdr % NumCols
    pressure = PField % RData(i,j)
    IF ( (pressure <= 1000) .AND. (pressure >= 0.) ) THEN
      pressure = 1000.
    END IF
    IF ( pressure > Press_Bot) THEN
      pressure = Press_Bot
    END IF

    IF (pressure == PField % Hdr % BMDI) THEN
      IHField % RData(i,j) = RMDI
    ELSE IF (pressure > Press_Mid) THEN ! Hts up to 11,000 GPM
      pressure = pressure/Press_Bot
      pressure = 1.0 - pressure**ZP1
      IHField % RData(i,j) = pressure*Temp_Bot/Lapse_RateL

    ELSE IF (pressure > Press_Top) THEN ! Hts between 11,000
                                        !     and     20,000 GPM
      pressure = pressure/Press_Mid
      pressure = -ALOG(pressure)
      IHField % RData(i,j) = Gpm1 + pressure*Temp_Top/G_over_R

    ELSE                                ! Hts above 20,000 GPM
      pressure = pressure/Press_Top
      pressure = 1.0 - pressure**ZP2
      IHField % RData(i,j) = Gpm2 + pressure*Temp_Top/Lapse_RateU

    END IF

  ENDDO
ENDDO

WHERE( IHField % RData /= RMDI )
  IHField % RData = IHField % RData * MtoKft
END WHERE

9999 CONTINUE

END SUBROUTINE ICAOHeight

