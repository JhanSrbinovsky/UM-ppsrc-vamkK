! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Constants used by the aerosol climatology for NWP
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
MODULE arcl_mod

IMPLICIT NONE

! Available species:
!   1. Sulphate
!   2. Mineral dust
!   3. Sea-salt
!   4. Black-carbon (also named soot)
!   5. Biomass-burning
!   6. Fossil-fuel organic carbon
!   7. Delta aerosol
! Note: When adding a species, increase NPD_ARCL_SPECIES which 
!       must be equal to the largest IP_ARCL_X

  INTEGER, PARAMETER :: ip_arcl_sulp = 1
  INTEGER, PARAMETER :: ip_arcl_dust = 2
  INTEGER, PARAMETER :: ip_arcl_sslt = 3
  INTEGER, PARAMETER :: ip_arcl_blck = 4
  INTEGER, PARAMETER :: ip_arcl_biom = 5
  INTEGER, PARAMETER :: ip_arcl_ocff = 6
  INTEGER, PARAMETER :: ip_arcl_dlta = 7
  INTEGER, PARAMETER :: npd_arcl_species = 7

! List of components.
!   The number of components depends on the species:
!   1. Sulphate: 3 (accumulation, Aitken, dissolved)
!   2. Mineral dust: 6 size bins
!   3. Sea-salt: 2 (film and jet)
!   4. Black-carbon: 2 (fresh and aged)
!   5. Biomass-burning: 3 (fresh, aged, in-cloud)
!   6. Fossil-fuel organic carbon: 3 (fresh, aged, in-cloud)
!   7. Delta aerosol: 1

  INTEGER, PARAMETER :: IP_ARCL_SULP_AC = 1
  INTEGER, PARAMETER :: IP_ARCL_SULP_AK = 2
  INTEGER, PARAMETER :: IP_ARCL_SULP_DI = 3
  INTEGER, PARAMETER :: IP_ARCL_DUST_B1 = 4
  INTEGER, PARAMETER :: IP_ARCL_DUST_B2 = 5
  INTEGER, PARAMETER :: IP_ARCL_DUST_B3 = 6
  INTEGER, PARAMETER :: IP_ARCL_DUST_B4 = 7
  INTEGER, PARAMETER :: IP_ARCL_DUST_B5 = 8
  INTEGER, PARAMETER :: IP_ARCL_DUST_B6 = 9
  INTEGER, PARAMETER :: IP_ARCL_SSLT_FI = 10
  INTEGER, PARAMETER :: IP_ARCL_SSLT_JT = 11
  INTEGER, PARAMETER :: IP_ARCL_BLCK_FR = 12
  INTEGER, PARAMETER :: IP_ARCL_BLCK_AG = 13
  INTEGER, PARAMETER :: IP_ARCL_BIOM_FR = 14
  INTEGER, PARAMETER :: IP_ARCL_BIOM_AG = 15
  INTEGER, PARAMETER :: IP_ARCL_BIOM_IC = 16
  INTEGER, PARAMETER :: IP_ARCL_OCFF_FR = 17
  INTEGER, PARAMETER :: IP_ARCL_OCFF_AG = 18
  INTEGER, PARAMETER :: IP_ARCL_OCFF_IC = 19
  INTEGER, PARAMETER :: IP_ARCL_DLTA_DL = 20
  INTEGER, PARAMETER :: NPD_ARCL_COMPNTS = 20

! This must also match the number of components for each species:
  INTEGER, DIMENSION(npd_arcl_species), PARAMETER :: n_arcl_compnts_per_species &
   = (/                                                                         &
! IP_ARCL_SULP: Sulphate (accumulation, Aitken, dissolved)
  3,                                                                            &
! IP_ARCL_DUST: Six size bins
  6,                                                                            &
! IP_ARCL_SSLT: Sea-salt (film and jet)
  2,                                                                            &
! IP_ARCL_BLCK: Black-carbon (fresh and aged)
  2,                                                                            &
! IP_ARCL_BIOM: Biomass (fresh, aged, in-cloud)
  3,                                                                            &
! IP_ARCL_OCFF: Fossil-fuel org. carb. (fresh, aged, in-cld)
  3,                                                                            &
! IP_ARCL_DLTA: Delta aerosol
  1                                                                             &
  /)


END MODULE arcl_mod
