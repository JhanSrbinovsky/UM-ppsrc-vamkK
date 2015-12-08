! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module containing field header codes for Fieldcalc

MODULE FldCodes_mod

! Description:
! Stores the following for fieldcalc use:
!  1. Stash Codes
!  2. Met O8 Field codes
!  3. Met O8 Level codes
!  4. PP codes
!  5. PP vertical co-ordinate codes
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

IMPLICIT None

! ---- Fieldcalc Section Numbers ----
INTEGER, PARAMETER :: UMSectnNo =    20   !  => STASH codes 20xxx
INTEGER, PARAMETER :: ProcessNo =    10   !  Superstash Process

! ---- Input fields ----
! STASH codes:                             Variables:          Levels:
INTEGER, PARAMETER :: ST_Urho   =     2   ! Wind U-component    rho
INTEGER, PARAMETER :: ST_Vrho   =     3   ! Wind V-component    rho
INTEGER, PARAMETER :: ST_Htheta =    10   ! Specific humidity   theta
INTEGER, PARAMETER :: ST_Wtheta =   150   ! Wind W-component    theta
INTEGER, PARAMETER :: ST_Exnerrho=  255   ! Exner               rho
INTEGER, PARAMETER :: ST_Prho   =   407   ! Pressure            rho
INTEGER, PARAMETER :: ST_Ptheta =   408   ! Pressure            theta
INTEGER, PARAMETER :: ST_BlkCld =   266   ! Bulk Cloud Fraction theta
INTEGER, PARAMETER :: ST_Dust1  =   431   ! Dust MMR1           theta
INTEGER, PARAMETER :: ST_Dust2  =   432   ! Dust MMR2           theta
INTEGER, PARAMETER :: ST_Dust3  =   433   ! Dust MMR3           theta
INTEGER, PARAMETER :: ST_Dust4  =   434   ! Dust MMR4           theta
INTEGER, PARAMETER :: ST_Dust5  =   435   ! Dust MMR5           theta
INTEGER, PARAMETER :: ST_Dust6  =   436   ! Dust MMR6           theta
INTEGER, PARAMETER :: ST_Usurf  =  3219   ! Wind U-component    theta(1)
INTEGER, PARAMETER :: ST_Vsurf  =  3220   ! Wind V-component    theta(1)
INTEGER, PARAMETER :: ST_TScreen=  3236   ! Temperature         1.5m
INTEGER, PARAMETER :: ST_VisIncPpn=3281   ! Visibility (inc ppn)1.5m
INTEGER, PARAMETER :: ST_CldWat =  4205   ! Cloud Liquid Water  theta
INTEGER, PARAMETER :: ST_CldIce =  4206   ! Cloud Ice Content   theta
INTEGER, PARAMETER :: ST_ConCld =  5212   ! Conv Cloud Amount   theta
INTEGER, PARAMETER :: ST_GWSU   =  6201   ! Gravity Wave U      ?
INTEGER, PARAMETER :: ST_GWSV   =  6202   ! Gravity Wave V      ?
INTEGER, PARAMETER :: ST_Ttheta = 16004   ! Temperature         theta
INTEGER, PARAMETER :: ST_Ztheta = 16201   ! Geoptl Height       theta
INTEGER, PARAMETER :: ST_Zrho   = 16255   ! Geoptl Height       rho
INTEGER, PARAMETER :: ST_Dust   = 17257   ! Total dust concentration
                                          !  (microg/m3)        theta

INTEGER, PARAMETER :: ST_LSM    =    30   ! Land-sea mask       surface
INTEGER, PARAMETER :: ST_Orog   =    33   ! Orography           surface
INTEGER, PARAMETER :: ST_Pstar  =   409   ! Pressure            surface
INTEGER, PARAMETER :: ST_CClBP  =  5207   ! Conv Cloud Base Pressure
INTEGER, PARAMETER :: ST_CClTP  =  5208   ! Conv Cloud Top Pressure
INTEGER, PARAMETER :: ST_LCClBP =  5222   ! Lwst Conv Cloud Base Press
INTEGER, PARAMETER :: ST_LCClTP =  5223   ! Lwst Conv Cloud Top Press
INTEGER, PARAMETER :: ST_Ustd   = 15201   ! Wind U-component   std
INTEGER, PARAMETER :: ST_Vstd   = 15202   ! Wind V-component   std

! ---- Output diagnostics ----
! STASH codes:
INTEGER, PARAMETER :: ST_CPNRT  =  5205   ! Conv Precipitation rate
INTEGER, PARAMETER :: ST_LCClBI =  5224   ! Lwst Conv Cloud Base ICAO Ht
INTEGER, PARAMETER :: ST_LCClTI =  5225   ! Lwst Conv Cloud Top ICAO Ht
INTEGER, PARAMETER :: ST_Diverg = 20005   ! Divergence          std
INTEGER, PARAMETER :: ST_Vortic = 20006   ! Vorticity           std
INTEGER, PARAMETER :: ST_MWT    = 20007   ! Mountain Wave
                                          ! Turbulence Predictor
INTEGER, PARAMETER :: ST_PrSym  = 20014   ! Precip Symbol
INTEGER, PARAMETER :: ST_WXCode = 20015   ! Present Weather Code
INTEGER, PARAMETER :: ST_CAT    = 20016   ! CAT Predictor       std
INTEGER, PARAMETER :: ST_MxCAT  = 20017   ! Max CAT
INTEGER, PARAMETER :: ST_MxCATP = 20018   ! Max CAT Pressure
INTEGER, PARAMETER :: ST_MaxWU  = 20020   ! Max Wind U-component
INTEGER, PARAMETER :: ST_MaxWV  = 20021   ! Max Wind V-component
INTEGER, PARAMETER :: ST_MaxWP  = 20022   ! Max Wind Pressure
INTEGER, PARAMETER :: ST_MaxWI  = 20023   ! Max Wind ICAO Height
INTEGER, PARAMETER :: ST_TropP  = 20024   ! Tropopause Pressure
INTEGER, PARAMETER :: ST_TropT  = 20025   ! Tropopause Temperature
INTEGER, PARAMETER :: ST_TropZ  = 20026   ! Tropopause Height
INTEGER, PARAMETER :: ST_TropI  = 20027   ! Tropopause ICAO Height
INTEGER, PARAMETER :: ST_SnowPr = 20028   ! Snow Probability
INTEGER, PARAMETER :: ST_ContrL = 20029   ! Contrail Lower Lim
INTEGER, PARAMETER :: ST_ContrU = 20030   ! Contrail Upper Lim
INTEGER, PARAMETER :: ST_ThAdv  = 20031   ! Thermal Advection   std
INTEGER, PARAMETER :: ST_FreezZ = 20033   ! Freezing Level Geoptl Height
INTEGER, PARAMETER :: ST_FreezP = 20034   ! Freezing Level Pressure
INTEGER, PARAMETER :: ST_FreezI = 20035   ! Freezing Level ICAO Height
INTEGER, PARAMETER :: ST_Iso20Z = 20036   ! -20C Isotherm Geoptl Height
INTEGER, PARAMETER :: ST_Iso20P = 20037   ! -20C Isotherm Pressure
INTEGER, PARAMETER :: ST_Iso20I = 20038   ! -20C Isotherm ICAO Height
INTEGER, PARAMETER :: ST_CClBI  = 20039   ! Conv Cloud Base ICAO Height
INTEGER, PARAMETER :: ST_CClTI  = 20040   ! Conv Cloud Top ICAO Height
INTEGER, PARAMETER :: ST_MWBase = 20041   ! Max Wind Base - ICAO Hts
INTEGER, PARAMETER :: ST_MWTop  = 20042   ! Max Wind Top  - ICAO Hts
INTEGER, PARAMETER :: ST_MnIceP = 20043   ! Mean Icing Potential
INTEGER, PARAMETER :: ST_MxIceP = 20044   ! Max  Icing Potential
INTEGER, PARAMETER :: ST_MnICTP = 20045   ! Mean In-Cloud Turb Potential
INTEGER, PARAMETER :: ST_MxICTP = 20046   ! Max  In-Cloud Turb Potential
INTEGER, PARAMETER :: ST_MnCATPt= 20047   ! Mean CAT Potential
INTEGER, PARAMETER :: ST_MxCATPt= 20048   ! Max  CAT Potential
INTEGER, PARAMETER :: ST_CBHorE = 20049   ! CB Horizontal Extent
INTEGER, PARAMETER :: ST_P_CBB  = 20050   ! Pressure at CB Base
INTEGER, PARAMETER :: ST_P_CBT  = 20051   ! Pressure at CB Top
INTEGER, PARAMETER :: ST_P_ECBB = 20052   ! Pressure at Embedded CB Base
INTEGER, PARAMETER :: ST_P_ECBT = 20053   ! Pressure at Embedded CB Top
INTEGER, PARAMETER :: ST_I_CBB  = 20054   ! ICAO Ht at CB Base
INTEGER, PARAMETER :: ST_I_CBT  = 20055   ! ICAO Ht at CB Top
INTEGER, PARAMETER :: ST_I_ECBB = 20056   ! ICAO Ht at Embedded CB Base
INTEGER, PARAMETER :: ST_I_ECBT = 20057   ! ICAO Ht at Embedded CB Top
INTEGER, PARAMETER :: ST_DUSTC  = 20058   ! Surf dust conc
INTEGER, PARAMETER :: ST_DUSTY  = 20059   ! 2000-5000ft dust conc
INTEGER, PARAMETER :: ST_ZTD    = 20060   ! ZTD diagnostic
INTEGER, PARAMETER :: ST_Iso70Z = 20061   ! -70C Isotherm Geoptl Height
INTEGER, PARAMETER :: ST_Iso70P = 20062   ! -70C Isotherm Pressure
INTEGER, PARAMETER :: ST_Iso70I = 20063   ! -70C Isotherm ICAO Height
INTEGER, PARAMETER :: ST_TropP_GRIB2 = 20064   ! Tropopause Pressure
INTEGER, PARAMETER :: ST_TropT_GRIB2 = 20065   ! Tropopause Temperature
INTEGER, PARAMETER :: ST_TropZ_GRIB2 = 20066   ! Tropopause Height
INTEGER, PARAMETER :: ST_TropI_GRIB2 = 20067   ! Tropopause ICAO Height
INTEGER, PARAMETER :: ST_VisTot = 20070   ! 1.5m Vis incl. ppn and dust (m)
INTEGER, PARAMETER :: ST_VisDust= 20071   ! 1.5m Vis in dust only (m)

! MetO8 codes:
INTEGER, PARAMETER :: MO8_Z      =    2   ! Height
INTEGER, PARAMETER :: MO8_SnowPr =   27   ! Snow Probability
INTEGER, PARAMETER :: MO8_Iso20Z =   29   ! -20C Isotherm Geoptl Height
INTEGER, PARAMETER :: MO8_ThAdv  =   40   ! Thermal Advection
INTEGER, PARAMETER :: MO8_Iso20P =   41   ! -20C Isotherm Pressure
INTEGER, PARAMETER :: MO8_Iso20I =   42   ! -20C Isotherm ICAO Height
INTEGER, PARAMETER :: MO8_FreezI =   43   ! Freezing Level ICAO Height
INTEGER, PARAMETER :: MO8_ContrL =   44   ! Contrail Lower Lim
INTEGER, PARAMETER :: MO8_ContrU =   45   ! Contrail Upper Lim
INTEGER, PARAMETER :: MO8_MaxWU  =   46   ! Max Wind U-component
INTEGER, PARAMETER :: MO8_MaxWV  =   47   ! Max Wind V-component
INTEGER, PARAMETER :: MO8_MaxWI  =   48   ! Max Wind ICAO Height
INTEGER, PARAMETER :: MO8_TropI  =   49   ! Tropopause ICAO Height
INTEGER, PARAMETER :: MO8_FreezZ =   50   ! Freezing Level Geoptl Height
INTEGER, PARAMETER :: MO8_TropP  =   51   ! Tropopause Pressure
INTEGER, PARAMETER :: MO8_TropT  =   52   ! Tropopause Temperature
INTEGER, PARAMETER :: MO8_TropZ  =   53   ! Tropopause Height
INTEGER, PARAMETER :: MO8_MaxWP  =   54   ! Max Wind Pressure
INTEGER, PARAMETER :: MO8_FreezP =   59   ! Freezing Level Pressure
INTEGER, PARAMETER :: MO8_CPNRT  =   64   ! Convective Precipitation Rate
INTEGER, PARAMETER :: MO8_CClBI  =   87   ! Conv Cloud Base ICAO Height
INTEGER, PARAMETER :: MO8_CClTI  =   88   ! Conv Cloud Top ICAO Height
INTEGER, PARAMETER :: MO8_LCClBI =  204   ! Lwst Conv Cld Base ICAO Hght
INTEGER, PARAMETER :: MO8_LCClTI =  205   ! Lwst Conv Cld Top ICAO Hght
INTEGER, PARAMETER :: MO8_PrSym  =  302   ! Precip Symbol
INTEGER, PARAMETER :: MO8_WXCode =  303   ! Present Weather Code
INTEGER, PARAMETER :: MO8_Diverg =  304   ! Divergence
INTEGER, PARAMETER :: MO8_MWT    =  305   ! Mountain Wave
                                          ! Turbulence Predictor
INTEGER, PARAMETER :: MO8_CAT    =  307   ! CAT Predictor
INTEGER, PARAMETER :: MO8_MxCATP =  313   ! Max CAT Pressure
INTEGER, PARAMETER :: MO8_MxCAT  =  314   ! Max CAT
INTEGER, PARAMETER :: MO8_Vortic =  316   ! Vorticity
INTEGER, PARAMETER :: MO8_MxWBase = 334   ! Max Wind Base - ICAO Hts
INTEGER, PARAMETER :: MO8_MxWTop  = 335   ! Max Wind Top  - ICAO Hts
INTEGER, PARAMETER :: MO8_MnIceP  = 345   ! Mean Icing Potential
INTEGER, PARAMETER :: MO8_MxIceP  = 346   ! Max  Icing Potential
INTEGER, PARAMETER :: MO8_MnICTP  = 347   ! Mean In-Cloud Turb Potential
INTEGER, PARAMETER :: MO8_MxICTP  = 348   ! Max  In-Cloud Turb Potential
INTEGER, PARAMETER :: MO8_MnCATPt = 349   ! Mean CAT Potential
INTEGER, PARAMETER :: MO8_MxCATPt = 350   ! Max  CAT Potential
INTEGER, PARAMETER :: MO8_CBHorE  = 351   ! CB Horizontal Extent
INTEGER, PARAMETER :: MO8_P_CBB   = 352   ! Pressure at CB Base
INTEGER, PARAMETER :: MO8_P_CBT   = 353   ! Pressure at CB Top
INTEGER, PARAMETER :: MO8_P_ECBB  = 354   ! Pressure at Embedded CB Base
INTEGER, PARAMETER :: MO8_P_ECBT  = 355   ! Pressure at Embedded CB Top
INTEGER, PARAMETER :: MO8_I_CBB   = 356   ! ICAO Ht at CB Base
INTEGER, PARAMETER :: MO8_I_CBT   = 357   ! ICAO Ht at CB Top
INTEGER, PARAMETER :: MO8_I_ECBB  = 358   ! ICAO Ht at Embedded CB Base
INTEGER, PARAMETER :: MO8_I_ECBT  = 359   ! ICAO Ht at Embedded CB Top
INTEGER, PARAMETER :: MO8_DUSTC   = 381   ! Surf dust conc
INTEGER, PARAMETER :: MO8_DUSTY   = 382   ! 2000-5000ft dust conc
INTEGER, PARAMETER :: MO8_Iso70Z  = 473   ! -70C Isotherm Geoptl Height
INTEGER, PARAMETER :: MO8_Iso70P  = 474   ! -70C Isotherm Pressure
INTEGER, PARAMETER :: MO8_Iso70I  = 475   ! -70C Isotherm ICAO Height
INTEGER, PARAMETER :: MO8_ZTD     = 476   ! ZTD diagnostic
INTEGER, PARAMETER :: MO8_TropI_GRIB2 =  478   ! Tropopause ICAO Height
INTEGER, PARAMETER :: MO8_TropP_GRIB2 =  479   ! Tropopause Pressure
INTEGER, PARAMETER :: MO8_TropT_GRIB2 =  480   ! Tropopause Temperature
INTEGER, PARAMETER :: MO8_TropZ_GRIB2 =  481   ! Tropopause Height
INTEGER, PARAMETER :: MO8_VisTot  = 528   ! 1.5m Vis incl. ppn and dust (m)
INTEGER, PARAMETER :: MO8_VisDust = 529   ! 1.5m Vis in dust only (m)



! PP codes:
INTEGER, PARAMETER :: PP_Z       =     1  ! Height
INTEGER, PARAMETER :: PP_ICAOHt  =     4  ! ICAO Height
INTEGER, PARAMETER :: PP_P       =     8  ! Pressure
INTEGER, PARAMETER :: PP_T       =    16  ! Temperature
INTEGER, PARAMETER :: PP_Contr   =    35  ! Contrails
INTEGER, PARAMETER :: PP_U       =    56  ! Wind U-component
INTEGER, PARAMETER :: PP_V       =    57  ! Wind V-component
INTEGER, PARAMETER :: PP_Vortic  =    73  ! Vorticity
INTEGER, PARAMETER :: PP_Diverg  =    74  ! Divergence
INTEGER, PARAMETER :: PP_CPNRT   =    98  ! Convective Precipitation Rate
INTEGER, PARAMETER :: PP_ThAdv   =   189  ! Thermal Advection
INTEGER, PARAMETER :: PP_Prob    =   191  ! Snow Probability
INTEGER, PARAMETER :: PP_ZTD     =  1114  ! ZTD diagnostic
INTEGER, PARAMETER :: PP_PrSym   =  1564  ! Precip Symbol
INTEGER, PARAMETER :: PP_WXCode  =  1565  ! Present Weather Code
INTEGER, PARAMETER :: PP_MWT     =  1578  ! Mountain Wave Turbulence Predictor
INTEGER, PARAMETER :: PP_DUSTC   =  1643  ! Surf dust conc
INTEGER, PARAMETER :: PP_DUSTY   =  1644  ! 2000-5000ft dust conc
INTEGER, PARAMETER :: PP_VisTot  =  1649  ! 1.5m Vis incl. ppn and dust (m)
INTEGER, PARAMETER :: PP_VisDust =  1650  ! 1.5m Vis in dust only (m)
INTEGER, PARAMETER :: PP_CAT     =  1690  ! CAT Predictor
INTEGER, PARAMETER :: PP_MnIceP  =  1981  ! Mean Icing Potential
INTEGER, PARAMETER :: PP_MxIceP  =  1982  ! Max  Icing Potential
INTEGER, PARAMETER :: PP_MnICTP  =  1983  ! Mean In-Cloud Turb Potential
INTEGER, PARAMETER :: PP_MxICTP  =  1984  ! Max  In-Cloud Turb Potential
INTEGER, PARAMETER :: PP_MnCATPt =  1985  ! Mean CAT Potential
INTEGER, PARAMETER :: PP_MxCATPt =  1986  ! Max  CAT Potential
INTEGER, PARAMETER :: PP_CBHorE  =  1987  ! CB Horizontal Extent
INTEGER, PARAMETER :: PP_P_CBB   =  1988  ! Pressure at CB Base
INTEGER, PARAMETER :: PP_P_CBT   =  1989  ! Pressure at CB Top
INTEGER, PARAMETER :: PP_P_ECBB  =  1990  ! Pressure at Embedded CB Base
INTEGER, PARAMETER :: PP_P_ECBT  =  1991  ! Pressure at Embedded CB Top
INTEGER, PARAMETER :: PP_I_CBB   =  1992  ! ICAO Ht at CB Base
INTEGER, PARAMETER :: PP_I_CBT   =  1993  ! ICAO Ht at CB Top
INTEGER, PARAMETER :: PP_I_ECBB  =  1994  ! ICAO Ht at Embedded CB Base
INTEGER, PARAMETER :: PP_I_ECBT  =  1995  ! ICAO Ht at Embedded CB Top

! Vertical co-ordinate types:
INTEGER, PARAMETER :: VC_Height  = 1      ! Height
INTEGER, PARAMETER :: VC_Press   = 8      ! Pressure
INTEGER, PARAMETER :: VC_Turb    = 126    ! Turbulence level
INTEGER, PARAMETER :: VC_Surface = 129    ! Surface
INTEGER, PARAMETER :: VC_Trop    = 130    ! Tropopause
INTEGER, PARAMETER :: VC_MaxWind = 131    ! Max Wind Level
INTEGER, PARAMETER :: VC_Freez   = 132    ! Freezing Level
INTEGER, PARAMETER :: VC_Iso20   = 134    ! -20C Isotherm Level
INTEGER, PARAMETER :: VC_Upper   = 135    ! Upper Level
INTEGER, PARAMETER :: VC_Lower   = 136    ! Lower Level
INTEGER, PARAMETER :: VC_MnIceP  = 8      ! Mean Icing Potential
INTEGER, PARAMETER :: VC_MxIceP  = 8      ! Max  Icing Potential
INTEGER, PARAMETER :: VC_MnICTP  = 8      ! Mean In-Cloud Turb Potential
INTEGER, PARAMETER :: VC_MxICTP  = 8      ! Max  In-Cloud Turb Potential
INTEGER, PARAMETER :: VC_MnCATPt = 8      ! Mean CAT Potential
INTEGER, PARAMETER :: VC_MxCATPt = 8      ! Max  CAT Potential
INTEGER, PARAMETER :: VC_CBHorE  = 0      ! WAFC Cb Horizontal Extent
INTEGER, PARAMETER :: VC_P_CBB   = 138    ! WAFC Pressure at CB Base
INTEGER, PARAMETER :: VC_P_CBT   = 137    ! WAFC Pressure at CB Top
INTEGER, PARAMETER :: VC_P_ECBB  = 138    ! Pressure at Embedded CB Base
INTEGER, PARAMETER :: VC_P_ECBT  = 137    ! Pressure at Embedded CB Top
INTEGER, PARAMETER :: VC_I_CBB   = 136    ! ICAO Ht at CB Base
INTEGER, PARAMETER :: VC_I_CBT   = 135    ! ICAO Ht at CB Top
INTEGER, PARAMETER :: VC_I_ECBB  = 136    ! ICAO Ht at Embedded CB Base
INTEGER, PARAMETER :: VC_I_ECBT  = 135    ! ICAO Ht at Embedded CB Top
INTEGER, PARAMETER :: VC_CPNRT   = 129    ! Convective Precipitation Rate
INTEGER, PARAMETER :: VC_Iso70   = 2001   ! -70C Isotherm Level

! MetO8 levels:
INTEGER, PARAMETER :: LV_Special = 8888   ! Special level
INTEGER, PARAMETER :: LV_Surface = 9999   ! Surface

END MODULE FldCodes_mod
