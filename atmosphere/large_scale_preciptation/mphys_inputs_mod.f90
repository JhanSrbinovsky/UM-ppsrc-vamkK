! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!+ data module for switches/options concerned with the precipitation scheme.
  ! Description:
  !   Module containing runtime options/data used by the precipitation scheme

  ! Method:
  !   Switches and associated data values used by the precipitation scheme
  !   are defined here and assigned default values. These may be overridden
  !   by namelist input.

  !   A description of what each switch or number refers to is provided
  !   with the namelist

  !   Any routine wishing to use these options may do so with the 'USE'
  !   statement.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Large Scale Precipitation

  ! Code Description:
  !   Language: FORTRAN 90


MODULE mphys_inputs_mod

  USE missing_data_mod, ONLY: rmdi, imdi

  IMPLICIT NONE

!===========================================================================
! INTEGER options set from RUN_PRECIP namelist
!===========================================================================

! Number of iterations of microphysics scheme
  INTEGER :: lsiter           = imdi ! loop inside of column
  INTEGER :: niter_bs         = imdi ! loop outside of column
!------------------------------------

!===========================================================================
! LOGICAL options set from RUN_PRECIP namelist
!===========================================================================

! Allow multiple iterations of microphysics scheme
  LOGICAL :: l_mcr_iter       = .FALSE.
!------------------------------------
! Share supersaturation
! between crystals and
! aggregates
  LOGICAL :: l_cry_agg_dep    = .FALSE.
!-----------------------------------

! Use iterative melting
  LOGICAL :: l_it_melting     = .FALSE.
!-----------------------------------

! Use improved warm rain microphysics scheme
  LOGICAL :: l_warm_new       = .FALSE.
!--------------------------------------

! Use 3b autoconversion rate and limit
  LOGICAL :: l_autoc_3b       = .FALSE.
!------------------------------------

! Hardwire
  LOGICAL :: l_autolim_3b     = .FALSE.
!------------------------------------

! Use generic
! ice psd
  LOGICAL :: l_psd            = .FALSE.
!------------------------------------

! Use global version
! (mid-lat selected if
! .false.)
  LOGICAL :: l_psd_global     = .FALSE.
!------------------------------------

! Activate the Hallett-Mossop
! process.
!
  LOGICAL :: l_hallett_mossop = .FALSE.
!------------------------------------

! Use murk aerosol to
! calculate the droplet
! number
  LOGICAL :: l_autoconv_murk  = .FALSE.
!------------------------------------

! Enable tapering of cloud droplets
! towards surface
  LOGICAL :: l_droplet_tpr    = .FALSE.
!------------------------------------

! Use the Clark et al (2008) aerosol
! scheme in MURK calculations
  LOGICAL :: l_clark_aero     = .FALSE.

! New variant of taper curve
! with variable aerosol at surface
  LOGICAL :: l_taper_new      = .FALSE.
!------------------------------------

! Use Abel & Shipway
! rain fall speeds
  LOGICAL :: l_rainfall_as    = .FALSE.
!------------------------------------

! Allow snow-rain collisions to produce
! graupel
  LOGICAL :: l_sr2graup       = .FALSE.
!------------------------------------

! Use Aerosol climatologies to generate drop number
  LOGICAL :: l_mcr_arcl       = .FALSE.
!-----------------------------------

! Include prognostic rain
  LOGICAL :: l_mcr_qrain      = .FALSE.
!-----------------------------------

! Include prognosic graupel
  LOGICAL :: l_mcr_qgraup     = .FALSE.
!-----------------------------------

! Prognostic rain lbcs active
  LOGICAL :: l_mcr_qrain_lbc  = .FALSE.
!-----------------------------------

! Prognostic graupel lbcs active
  LOGICAL :: l_mcr_qgraup_lbc = .FALSE.
!-----------------------------------

! Turns precipitation code on/off
  LOGICAL :: l_rain           = .FALSE.
!-----------------------------------

!===========================================================================
! REAL values set from RUN_PRECIP namelist
!===========================================================================

! Rain particle size distribution
! values
  REAL :: x1r                 = rmdi
  REAL :: x2r                 = rmdi
!------------------------------------

! Ice mass-diameter relationship
! values
  REAL :: ai                  = rmdi
  REAL :: bi                  = rmdi
  REAL :: aic                 = rmdi
  REAL :: bic                 = rmdi
!------------------------------------

! Best-Reynolds numbers for crystals
  REAL :: lsp_eic             = rmdi
  REAL :: lsp_fic             = rmdi
!------------------------------------


! Droplet taper height

  REAL :: z_peak_nd           = rmdi
!------------------------------------

! Droplet number at model level 1:

  REAL :: ndrop_surf          = rmdi
!------------------------------------
! Maximum droplet number assumed at model level 1:
! (currently not set by the UMUI but a candidate for
!  the be in future)

  REAL :: max_drop_surf = rmdi
!------------------------------------

! Cloud-rain correlation coefficient for inhomogeneity parametrization:
! i.e. when l_inhomog=.true.

  REAL :: c_r_correl = rmdi
!------------------------------------

! Aerosol climatology scaling factor, to account for spatial/temporal 
! inhomogeneity
  REAL :: arcl_inhom_sc  = rmdi
!------------------------------------

! Sub-grid inhomogeneity of qcl, nd and qrain
! for use when l_micro_kk = .true.
  REAL :: alpha_q = rmdi
  REAL :: alpha_n = rmdi
  REAL :: alpha_r = rmdi
!-----------------------------------

! Axial ratios for aggregates and crystals
  REAL :: ar  = rmdi
  REAL :: arc = rmdi
!-----------------------------------

! Maximum Temp for ice nuclei nucleation (deg C)
! Typically minus 10 C. 
  REAL :: tnuc = rmdi
!-----------------------------------

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

! Define the RUN_PRECIP namelist

  NAMELIST/run_precip/                                                  &
         l_cry_agg_dep, l_it_melting, l_autoc_3b, l_warm_new,           &
         l_autolim_3b,  l_psd, l_psd_global, l_hallett_mossop,          &
         l_autoconv_murk,                                               &
         x1r, x2r, c_r_correl,                                          &
         ai, bi, aic, bic, lsp_eic, lsp_fic, ar, arc, tnuc,             &
         lsiter, z_peak_nd, ndrop_surf,                                 &
         l_droplet_tpr, l_clark_aero,                                   &
         l_taper_new, max_drop_surf, l_rainfall_as,                     &
         l_mcr_iter, niter_bs, l_sr2graup, l_mcr_arcl, arcl_inhom_sc,   &
         l_mcr_qrain, l_mcr_qgraup, l_mcr_qrain_lbc, l_mcr_qgraup_lbc,  &
         l_rain

!===========================================================================
! LOGICAL options not set in namelist
!===========================================================================

! Second ice variable lbcs active
  LOGICAL :: l_mcr_qcf2_lbc   = .FALSE.
!-----------------------------------

! Include second ice variable
  LOGICAL :: l_mcr_qcf2       = .FALSE.
!-----------------------------------

END MODULE mphys_inputs_mod
