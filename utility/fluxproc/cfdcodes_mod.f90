! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small Execs
MODULE cfdcodes_mod

IMPLICIT NONE
!----------------------------------------------------------------------
! comdeck: CFDCODES
! Purpose: list of parameters for stash and field codes.
!----------------------------------------------------------------------

! fields are listed by order of their numbers in both lists

      ! Stash codes

      ! Surface wind stress
      INTEGER,PARAMETER::StCWindStressU = 03219 !U (west -> east)
      INTEGER,PARAMETER::StCWindStressV = 03220 !V (south -> north)

      INTEGER,PARAMETER::StCWindMixEng  = 03224 !Wind mixing energy

      ! Net down short wave radiation
      INTEGER,PARAMETER::StCSW          = 01203

      ! Net down short wave radiation (band 1)
      INTEGER,PARAMETER::StCSW1         = 01204

      ! Net down long wave radiation
      INTEGER,PARAMETER::StCLongWave    = 02203

      ! Surface sensible heat flux from sea
      INTEGER,PARAMETER::StCSensibleHeat= 03228

      INTEGER,PARAMETER::StCSublim      = 03231 !Sublimation (in)
      INTEGER,PARAMETER::StCTopmelt     = 03235 !Topmelt (in)
      INTEGER,PARAMETER::StCBotmelt     = 03201 !Botmelt (in)
      INTEGER,PARAMETER::StCEvaporation = 03232 !Evaporation from sea
      INTEGER,PARAMETER::StCDrain       = 04201 !Large Scale Rain Amount
      INTEGER,PARAMETER::StCConvrain    = 05201 !Convective Rain Amount
      INTEGER,PARAMETER::StCDsnow       = 04202 !Large Scale Snow Amount
      INTEGER,PARAMETER::StCConvsnow    = 05202 !Convective Snow Amount
      INTEGER,PARAMETER::StCSST         = 00024 !Sea Surface Temp
      INTEGER,PARAMETER::StCSSS         = 00181 !Sea Surface Salinity
      INTEGER,PARAMETER::StCHICE        = 00032 !Ice Depth
      INTEGER,PARAMETER::StCAICE        = 00031 !Ice fraction
      INTEGER,PARAMETER::StCSSP         = 16222 !Sea Surface Pressure
      INTEGER, PARAMETER :: StCWindSpeedU   = 03225
      INTEGER, PARAMETER :: StCWindSpeedV   = 03226
      INTEGER,PARAMETER::OutStCTAUX  = 00150 ! windstress (x)
      INTEGER,PARAMETER::OutStCTAUY  = 00151 ! windstress (y)
      INTEGER,PARAMETER::OutStCWME   = 00152 ! wind mixing energy
      INTEGER,PARAMETER::OutStCSOL   = 00161 ! Solar Radiation
      INTEGER,PARAMETER::OutStCHTN   = 00162 ! net non penetrating heat
      INTEGER,PARAMETER::OutStCPLE   = 00165 ! Precip less evaporation
      INTEGER,PARAMETER::OutStCSNO   = 00171 ! Snowfall Rate
      INTEGER,PARAMETER::OutStCSUB   = 00172 ! Sublimation Rate
      INTEGER,PARAMETER::OutStCTOP   = 00190 ! Topmelt
      INTEGER,PARAMETER::OutStCBOT   = 00191 ! Bottom Melt
      INTEGER,PARAMETER::OutStCSST   = 00180 ! Sea Surface Temp
      INTEGER,PARAMETER::OutStCSSS   = 00181 ! Sea Surface Salinity
      INTEGER,PARAMETER::OutStCHICE  = 00183 ! Ice Depth
      INTEGER,PARAMETER::OutStCSSP   = 16222 ! Sea Surface Pressure
      INTEGER, PARAMETER :: OutStCWSPX      = 03225
      INTEGER, PARAMETER :: OutStCWSPY      = 03226

      ! FF Codes

      INTEGER,PARAMETER::FFTAUX = 254 ! U (west -> east)   surface
      INTEGER,PARAMETER::FFTAUY = 255 ! V (south -> north) wind stress
      INTEGER,PARAMETER::FFWME  = 256 ! Wind mixing energy
      INTEGER,PARAMETER::FFSOL  = 257 ! Net pen solar radiation
      INTEGER,PARAMETER::FFHTN  = 258 ! Net non penetrating heat
      INTEGER,PARAMETER::FFPLE  = 259 ! Precip less evap
      INTEGER,PARAMETER::FFSNO  = 260 ! Snowfall Rate
      INTEGER,PARAMETER::FFSUB  = 261 ! Sublimation Rate
      INTEGER,PARAMETER::FFTOP  = 265 ! Topmelt
      INTEGER,PARAMETER::FFBOT  = 266 ! Bottom Melt
      INTEGER,PARAMETER::FFSST  = 262 ! Sea Surface Temp
      INTEGER,PARAMETER::FFSSS  = 263 ! Sea Surface Salinity
      INTEGER,PARAMETER::FFHICE = 264 ! Ice Depth
      INTEGER,PARAMETER::FFSSP  = 012 ! Sea Surface Pressure
      INTEGER,PARAMETER::FFWSPX = 005 ! U (west -> east)   wind speed
      INTEGER,PARAMETER::FFWSPY = 006 ! V (south -> north) wind speed

      ! PP Codes

      INTEGER,PARAMETER::PPTAUX = 721 ! U (west -> east)   surface
      INTEGER,PARAMETER::PPTAUY = 722 ! V (south -> north) wind stress
      INTEGER,PARAMETER::PPWME  = 627 ! Wind mixing energy
      INTEGER,PARAMETER::PPSOL  = 625 ! Net pen solar radiation
      INTEGER,PARAMETER::PPHTN  = 626 ! Net non penetrating heat
      INTEGER,PARAMETER::PPPLE  = 629 ! Precip less evap
      INTEGER,PARAMETER::PPSNO  = 623 ! Snowfall Rate
      INTEGER,PARAMETER::PPSUB  = 624 ! Sublimation Rate
      INTEGER,PARAMETER::PPTOP  = 681 ! Topmelt
      INTEGER,PARAMETER::PPBOT  = 682 ! Bottom Melt
      INTEGER,PARAMETER::PPSST  = 650 ! Sea Surface Temp
      INTEGER,PARAMETER::PPSSS  = 649 ! Sea Surface Salinity
      INTEGER,PARAMETER::PPHICE = 675 ! Ice Depth
      INTEGER,PARAMETER::PPSSP  = 008 ! Sea surface pressure
      INTEGER,PARAMETER::PPWSPX = 048 ! U (west -> east)   wind speed
      INTEGER,PARAMETER::PPWSPY = 049 ! V (south -> north) wind speed
! CFDCODES end

END MODULE cfdcodes_mod
