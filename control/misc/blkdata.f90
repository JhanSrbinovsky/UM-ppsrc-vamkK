! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Block Data Subprogram : BLKDATA
!  
!   Purpose : Holds DATA for any variables that are in common blocks,
!             so that they are initialised in only one place.
!  
!  
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Misc
  BLOCK DATA BLKDATA
  USE Field_Types
  USE Atmos_Max_Sizes
  USE Control_Max_Sizes
  USE UM_ParParams
  USE entcoef_mod, ONLY: entcoef
 
  USE c_gwave_mod, ONLY: nsigma, amplitude_saturation, stress_saturation, &
                         beta_fix, frac_wl, lambdaz_min, lambdaz_max,     &
                         nsq_neutral, zav_converge, zav_iterate
  IMPLICIT NONE

!*L --------------------- Comdeck: CENVIR   ----------------------------
!
!   Purpose: COMDECK defining Character enviroment variables used
!            by portable IO to open and close files
!
!
! Type declarations
!
      CHARACTER(LEN=8) FT_ENVIRON(199)  ! Array holding enviroment variables
!                                  for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
!
!
!Common Blocks for character and integer arrays
!
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
!

!
!   Purpose: Data statements for character enviroment variables used
!            by portable IO to open/close files (links with CENVIR)
!
      DATA FT_ENVIRON/                                                  &
     &  'PPXREF  ','PPXREFU ','        ','STASHCTL','        ',         &
                                                                !  1- 5
     &  '        ','OUTPUT2 ','        ','        ','XHIST   ',         &
     &  'IHIST   ','THIST   ','ICECALVE','ERRFLAG ','CACHE1  ',         &
     &  'CACHE2  ','        ','ASWAP   ','OSWAP   ','AINITIAL',         &
     &  'ASTART  ','        ','APSUM1  ','APSTMP1 ','        ',         &
     &  '        ','AOMEAN  ','ATMANL  ','        ','OZONE   ',         &
     &  'SMCSNOWD','DSOILTMP','SOILTYPE','GENLAND ','SSTIN   ',         &
     &  'SICEIN  ','PERTURB ','CURNTIN ','        ','OINITIAL',         &
     &  'OSTART  ','        ','OPSUM1  ','OPSTMP1 ','        ',         &
     &  '        ','OCNANL  ','ATRACER ','OTRACER ','WFIN    ',         &
     &  'HFLUXIN ','PMEIN   ','ICEFIN  ','AIRTMP  ','SALINITY',         &
     &  'FLUXCORR','SWSPECTD','BAS_IND ','SLABHCON','PP0     ',         &
     &  'PP1     ','PP2     ','PP3     ','PP4     ','PP5     ',         &
     &  'PP6     ','PP7     ','PP8     ','PP9     ','OBS01   ',         &
                                                               !66-70
     &  'OBS02   ','OBS03   ','OBS04   ','OBS05   ','DUSTSOIL',         &
                                                               !71-75
     &  'BIOMASS ','RIVSTOR ','RIVCHAN ','RIVER2A ','LWSPECTD',         &
                                                               !76-80
     &  'SURGEOU1','SURGEOUT','PPSCREEN','PPSMC   ','WFOUT   ',         &
     &  'UARSOUT1','UARSOUT2','ICEFOUT ','MOSOUT  ','VERT_LEV',         &
     &  'SSTOUT  ','SICEOUT ','CURNOUT ','        ','DMSCONC ',         &
                                                               !91-95
     &  'OROG    ','TRANSP  ','OLABCIN ','OCNDEPTH',                    &
     &  'FOAMOUT1','FOAMOUT2','CXBKGERR','RFMOUT  ','FILE104 ',         &
     &  'FILE105 ','IDEALISE','TDF_dump','IAU_inc ','MURKFILE',         &
     &  'SULPEMIS','USRANCIL','USRMULTI','OUSRANCL','OUSRMULT',         &
     &  'SO2NATEM','CHEMOXID','AEROFCG ','CO2EMITS','TPPSOZON',         &
     &  'LANDFRAC','WLABCOU1','WLABCOU2','WLABCOU3','WLABCOU4',         &
                                                               !120-124
     &  'ALABCIN1','ALABCIN2','IOSCNTL ','OCFFEMIS','HORZGRID',         &
                                                               !125-129
     &  'SURFEMIS','AIRCREMS','STRATEMS','EXTRAEMS','RADONEMS',         &
                                                               !130-134
     &  'FRACINIT','VEGINIT ','DISTURB ','CACHED  ','SOOTEMIS',         &
                                                               !135-139
     &  'ALABCOU1','ALABCOU2','ALABCOU3','ALABCOU4','ALABCOU5',         &
                                                               !140-144
     &  'ALABCOU6','ALABCOU7','ALABCOU8','CARIOLO3','RPSEED  ',         & 
                                                               !145-149
     &  'PPVAR   ','PP10    ','ICFILE  ','VAR_GRID','ARCLBIOG',         &
                                                               !150-154
     &  'ARCLBIOM','ARCLBLCK','ARCLSSLT','ARCLSULP','ARCLDUST',         &
                                                               !155-159
     &  'ARCLOCFF','ARCLDLTA','TOPMEAN ','TOPSTDEV','PPMBC   ',         &
                                                               !160-164
     &  'UKCAPREC','UKCAACSW','UKCAACLW','UKCACRSW','UKCACRLW',         &
                                                               !165-169
     &  'UKCAFJXX','UKCAFJSC','UKCA2DO3','UKCA2CH4','UKCA2NOY',         &
                                                               !170-174
     &  'UKCA2PHO','UKCASTRD','UKCASTO3','UKCASTAR','UKCAFJAR',         &
                                                               !175-179
     &          20*'        '                                           &
     & /
! Text output file for STASH-related information is assigned to UNIT 200

      DATA LEN_FT_ENVIR/6,7,0,8,0, 0,7,0,0,5,                           &
                                                    !  1-10
     &                  5,5,8,7,6, 6,0,5,5,8,                           &
                                                    ! 11-20
     &                  6,0,6,7,0, 0,6,6,0,5,                           &
                                                    ! 21-30
     &                  8,8,8,7,5, 6,7,7,0,8,                           &
                                                    ! 31-40
     &                  6,0,6,7,0, 0,6,7,7,4,                           &
                                                    ! 41-50
     &                  7,5,6,6,8, 8,8,7,8,3,                           &
                                                    ! 51-60
     &                  3,3,3,3,3, 3,3,3,3,5,                           &
                                                    ! 61-70
     &                  5,5,5,5,8, 7,7,7,7,8,                           &
                                                    ! 71-80
     &                  8,8,8,5,5, 8,8,7,6,8,                           &
                                                    ! 81-90
     &                  6,7,7,0,7, 4,6,7,8,                             &
                                                    ! 91-99
     &                  8,8,8,6,7, 8,8,8,7,8,                           &
                                                    ! 100-109
     &                  8,8,8,8,8, 8,8,7,8,8,                           &
                                                    ! 110-119
     &                  8,8,8,8,8, 8,8,7,8,8,                           &
                                                    ! 120-129
     &                  8,8,8,8,8, 8,7,7,6,8,                           &
                                                    ! 130-139
     &                  8,8,8,8,8, 8,8,8,8,6,                           &
                                                    ! 140-149
     &                  5,4,6,8,8, 8,8,8,8,8,                           &
                                                    ! 150-159
     &                  8,8,7,8,5, 8,8,8,8,8,                           &
                                                    ! 160-169
     &                  8,8,8,8,8, 8,8,8,8,8,                           &
                                                    ! 170-179
     &                  20*0/                       ! 180-199


  END BLOCK DATA BLKDATA

