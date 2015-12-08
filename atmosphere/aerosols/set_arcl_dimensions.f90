! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine for setting internal dimensions and arrays associated
!  with the aerosol climatology for NWP.
!
! Purpose:
!   Translate model switches from CNTLATM to an internal array of
!   switches. Compute the number of aerosol species and components
!   that are requested.
!
! Method:
!   Straightforward.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols
!
! Description of code:
!   FORTRAN 90
!- ---------------------------------------------------------------------
MODULE set_arcl_dimensions_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE set_arcl_dimensions(                                                &
                                ! CNTLATM model switches
  l_use_arclbiom,                                                              &
  l_use_arclblck,                                                              &
  l_use_arclsslt,                                                              &
  l_use_arclsulp,                                                              &
  l_use_arcldust,                                                              &
  l_use_arclocff,                                                              &
  l_use_arcldlta,                                                              &
                                ! Number of aerosol species used
  n_arcl_species,                                                              &
                                ! Corresponding number of components
  n_arcl_compnts,                                                              &
                                ! Internal model switches
  l_use_arcl                                                                   &
  )

USE arcl_mod,        ONLY: npd_arcl_compnts, npd_arcl_species, ip_arcl_sulp,   &
                           ip_arcl_dust, ip_arcl_sslt, ip_arcl_blck,           &
                           ip_arcl_biom, ip_arcl_ocff, ip_arcl_dlta,           &
                           n_arcl_compnts_per_species

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!
! Arguments with intent(in)
!

!
! CNTLATM model switches
!
LOGICAL, INTENT(in) ::                                                         &
                                ! Biomass-burning aerosol climatology
  l_use_arclbiom,                                                              &
                                ! Black-carbon aerosol climatology
  l_use_arclblck,                                                              &
                                ! Sea-salt aerosol climatology
  l_use_arclsslt,                                                              &
                                ! Sulphate aerosol climatology
  l_use_arclsulp,                                                              &
                                ! Mineral dust aerosol climatology
  l_use_arcldust,                                                              &
                                ! Fossil-fuel organic carbon aerosol clim.
  l_use_arclocff,                                                              &
                                ! Delta aerosol climatology
  l_use_arcldlta

!
! Arguments with intent(out)
!

!
! Number of species used
!
INTEGER, INTENT(out) :: n_arcl_species

!
! Corresponding number of components
!
INTEGER, INTENT(out) :: n_arcl_compnts

!
! Internal model switches
!
LOGICAL, INTENT(out), DIMENSION(npd_arcl_species) :: l_use_arcl

!
! Local variables
!

INTEGER(kind=jpim), PARAMETER :: zhook_in  = 0
INTEGER(kind=jpim), PARAMETER :: zhook_out = 1
REAL(kind=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SET_ARCL_DIMENSIONS',zhook_in,zhook_handle)
n_arcl_species = 0
n_arcl_compnts = 0
l_use_arcl = .FALSE.

!
! For each CNTLATM model switch that is set, set the corresponding
! internal model switch, increase the number of requested species
! by 1, and increase the number of components by the corresponding
! amount.
!
IF (l_use_arclbiom) THEN

  l_use_arcl(ip_arcl_biom) = .TRUE.
  n_arcl_species = n_arcl_species + 1
  n_arcl_compnts = n_arcl_compnts +                                            &
    n_arcl_compnts_per_species(ip_arcl_biom)

END IF

IF (l_use_arclblck) THEN

  l_use_arcl(ip_arcl_blck) = .TRUE.
  n_arcl_species = n_arcl_species + 1
  n_arcl_compnts = n_arcl_compnts +                                            &
    n_arcl_compnts_per_species(ip_arcl_blck)

END IF

IF (l_use_arclsslt) THEN

  l_use_arcl(ip_arcl_sslt) = .TRUE.
  n_arcl_species = n_arcl_species + 1
  n_arcl_compnts = n_arcl_compnts +                                            &
    n_arcl_compnts_per_species(ip_arcl_sslt)

END IF

IF (l_use_arclsulp) THEN

  l_use_arcl(ip_arcl_sulp) = .TRUE.
  n_arcl_species = n_arcl_species + 1
  n_arcl_compnts = n_arcl_compnts +                                            &
    n_arcl_compnts_per_species(ip_arcl_sulp)

END IF

IF (l_use_arcldust) THEN

  l_use_arcl(ip_arcl_dust) = .TRUE.
  n_arcl_species = n_arcl_species + 1
  n_arcl_compnts = n_arcl_compnts +                                            &
    n_arcl_compnts_per_species(ip_arcl_dust)

END IF

IF (l_use_arclocff) THEN

  l_use_arcl(ip_arcl_ocff) = .TRUE.
  n_arcl_species = n_arcl_species + 1
  n_arcl_compnts = n_arcl_compnts +                                            &
    n_arcl_compnts_per_species(ip_arcl_ocff)

END IF

IF (l_use_arcldlta) THEN

  l_use_arcl(ip_arcl_dlta) = .TRUE.
  n_arcl_species = n_arcl_species + 1
  n_arcl_compnts = n_arcl_compnts +                                            &
    n_arcl_compnts_per_species(ip_arcl_dlta)

END IF

!
! Ensure that n_arcl_compnts is not zero, as it is used as
! array dimension in routines called later.
! n_arcl_compnts = 0 can only happen if n_arcl_species is also
! zero, in which case the value of n_arcl_compnts is unimportant.
!
n_arcl_compnts = MAX(1, n_arcl_compnts)

IF (lhook) CALL dr_hook('SET_ARCL_DIMENSIONS',zhook_out,zhook_handle)
RETURN

END SUBROUTINE set_arcl_dimensions
END MODULE set_arcl_dimensions_mod
