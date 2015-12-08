! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
MODULE aercmp3a_mod

IMPLICIT NONE
! AERCMP3A start
!     ------------------------------------------------------------------
!     MODULE TO SET INDICES OF AEROSOL COMPONENTS.
!
      ! maximum number of aerosol components
      ! N.B: this must be at least as large as the
      !      largest value in the list below
      INTEGER,PARAMETER:: NPD_AEROSOL_COMPONENT=32

      INTEGER,PARAMETER:: IP_WATER_SOLUBLE=1
      INTEGER,PARAMETER:: IP_DUST_LIKE=2
      INTEGER,PARAMETER:: IP_OCEANIC=3
      INTEGER,PARAMETER:: IP_SOOT=4
      INTEGER,PARAMETER:: IP_ASH=5
      INTEGER,PARAMETER:: IP_SULPHURIC=6
      INTEGER,PARAMETER:: IP_ACCUM_SULPHATE=10
      INTEGER,PARAMETER:: IP_AITKEN_SULPHATE=11
      INTEGER,PARAMETER:: IP_FRESH_SOOT=12
      INTEGER,PARAMETER:: IP_AGED_SOOT=13
      INTEGER,PARAMETER:: IP_SEASALT_FILM=15
      INTEGER,PARAMETER:: IP_SEASALT_JET=16
      INTEGER,PARAMETER:: IP_DUST_1=17
      INTEGER,PARAMETER:: IP_DUST_2=18
      INTEGER,PARAMETER:: IP_DUST_3=19
      INTEGER,PARAMETER:: IP_DUST_4=20
      INTEGER,PARAMETER:: IP_DUST_5=21
      INTEGER,PARAMETER:: IP_DUST_6=22
      INTEGER,PARAMETER:: IP_BIOMASS_1=23
      INTEGER,PARAMETER:: IP_BIOMASS_2=24
      INTEGER,PARAMETER:: IP_BIOGENIC=25
      INTEGER,PARAMETER:: IP_OCFF_FRESH=26
      INTEGER,PARAMETER:: IP_OCFF_AGED=27
      INTEGER,PARAMETER:: IP_DELTA=28
      INTEGER,PARAMETER:: IP_NITRATE=30
      INTEGER,PARAMETER:: ip_twobindust_1 = 31
      INTEGER,PARAMETER:: ip_twobindust_2 = 32

!     Aerosol data can be provided to radiaiton etc from different
!     sources, either prognostic, or different climatologies.
!     These can either effect radiative fluxes (RON), or are there
!     just for diagnostics (ROFF)
!     Simple land/sea climagology from Cusack(98):
      INTEGER,PARAMETER:: IP_AERSRC_CUSACK_RON=0
      INTEGER,PARAMETER:: IP_AERSRC_CUSACK_ROFF=10
!     CLASSIC prognostic aerosols:
      INTEGER,PARAMETER:: IP_AERSRC_CLASSIC_RON=1
      INTEGER,PARAMETER:: IP_AERSRC_CLASSIC_ROFF=11
!     Aerosol climatolgies, derived from CLASSIC:
      INTEGER,PARAMETER:: IP_AERSRC_ARCL_RON=2
      INTEGER,PARAMETER:: IP_AERSRC_ARCL_ROFF=12

! AERCMP3A end

END MODULE aercmp3a_mod
