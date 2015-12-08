! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : RUNDATA

MODULE s_rundata

  USE scm_cntl_mod, ONLY: scm_nml

  USE s_maxdim, ONLY:                                                         &
    mx_rw_lng, mx_rw, mx_st_lv, mx_nlnd, mx_tr_lv, mx_tr_vars                 &
  , mx_wet_lv, mx_nlnd, mx_mod_lv, mx_o3_lv, mx_ntile

  USE scm_utils, ONLY:                                                        &
    rmdi, imdi, rw_lng, rw, nlnd, nml_nmod_lv, nmod_lv, nwet_lv, ntr_lv       &
  , ntr_var, o3_lv, ntile, interpolate, z_th, z_rh, nml_z_th, nml_z_rh        &
  , zhook_in, zhook_out, jpim, jprb, lhook, dr_hook

  USE s_interp_mod, ONLY: interp1d
  USE s_main_force, ONLY:                                                     &
    change_clim, dump_step, dump_days, resdump_days, runno_in, runno_out      &
  , ntrad1, ndayin, nminin, nsecin, timestep, min_trop_level, max_trop_level  &
  , co2start, co2end, co2rate, modelid, exname_in, exname_out                 &
  , scm_ntml           => ntml            &
  , scm_ntdsc          => ntdsc           &
  , scm_nbdsc          => nbdsc           &
  , scm_albsoil        => albsoil         &
  , scm_albobs_sw      => albobs_sw       &
  , scm_albobs_vis     => albobs_vis      &
  , scm_albobs_nir     => albobs_nir      &
  , scm_so2_em         => so2_em          &
  , scm_nh3_em         => nh3_em          &
  , scm_dms_em         => dms_em          &
  , scm_soot_em        => soot_em         &
  , scm_bmass_em       => bmass_em        &
  , scm_ocff_em        => ocff_em         &
  , scm_soot_hilem     => soot_hilem      &
  , scm_bmass_hilem    => bmass_hilem     &
  , scm_ocff_hilem     => ocff_hilem      &
  , scm_soot           => soot            &
  , scm_cort           => cort            &
  , scm_cord           => cord            &
  , scm_corvn          => corvn           &
  , scm_corw           => corw            &
  , scm_cclwp          => cclwp           &
  , scm_orog           => orog            &
  , scm_ddmfx          => ddmfx           &
  , scm_zh             => zh              &
  , scm_dOLR_rts       => dolr_rts        &
  , scm_fland_ctile    => fland_ctile     &
  , scm_land_alb       => land_alb        &
  , scm_sice_alb       => sice_alb        &
  , scm_tstar_land     => tstar_land      &
  , scm_tstar_sea      => tstar_sea       &
  , scm_tstar_sice     => tstar_sice      &
  , scm_co2_emits      => co2_emits       &
  , scm_co2flux        => co2flux         &
  , scm_sum_eng_fluxes => sum_eng_fluxes  &
  , scm_sum_moist_flux => sum_moist_flux  &
  , scm_tscrndcl_ssi   => tscrndcl_ssi    &
  , scm_tscrndcl_tile  => tscrndcl_tile   &
  , scm_tstbtrans      => tstbtrans       &
  , scm_ozone          => ozone           &
  , scm_co2            => co2             &
  , scm_so2            => so2             &
  , scm_dms            => dms             &
  , scm_nh3            => nh3             &
  , scm_so4_aitken     => so4_aitken      &
  , scm_so4_accu       => so4_accu        &
  , scm_so4_diss       => so4_diss        &
  , scm_aerosol        => aerosol         &
  , scm_aerosol_em     => aerosol_em      &
  , scm_free_tracers   => free_tracers    &
  , scm_soot_new       => soot_new        &
  , scm_soot_aged      => soot_aged       &
  , scm_soot_cld       => soot_cld        &
  , scm_bmass_new      => bmass_new       &
  , scm_bmass_aged     => bmass_aged      &
  , scm_bmass_cld      => bmass_cld       &
  , scm_ocff_new       => ocff_new        &
  , scm_ocff_aged      => ocff_aged       &
  , scm_ocff_cld       => ocff_cld        &
  , scm_nitr_acc       => nitr_acc        &
  , scm_nitr_diss      => nitr_diss

  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Allows SCM to read in RUNDATA namelist from forcing file, <scm_nml>.
!   RUNDATA contains data required for the run.
!
! Method:
!   Namelist RUNDATA is defined in this module and read in by contained
!   subroutine read_rundata.  Scalar variables are read directly to
!   s_main_force.  Array variables are declared as private single rank arrays
!   of fixed length (arrays cannot be varying length in namelists) and
!   reshaped/interpolated before being transferred to arrays of the
!   correct size/shape in s_main_force.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!=============================================================================
!
! Declarations:

  !--------------------------------------------------
  ! Fixed arrays to read in namelist before reshaping
  !--------------------------------------------------
  INTEGER, PRIVATE ::              &
    ntml  (mx_rw_lng*mx_rw) = imdi &! Top level of surface mixed layer
  , ntdsc (mx_rw_lng*mx_rw) = imdi &! Top level of decoupled SC layer
  , nbdsc (mx_rw_lng*mx_rw) = imdi  ! Bottom level of decoupled SC layer

  REAL, PRIVATE :: &
    albsoil    (mx_nlnd) = rmdi &
  , albobs_sw  (mx_nlnd) = rmdi &
  , albobs_vis (mx_nlnd) = rmdi &
  , albobs_nir (mx_nlnd) = rmdi

  ! Aerosols/Gases, emissions at ...
  REAL, PRIVATE ::                       &
    so2_em      (mx_rw_lng*mx_rw) = rmdi &! ... srf, so2
  , nh3_em      (mx_rw_lng*mx_rw) = rmdi &! ... srf, nh3
  , dms_em      (mx_rw_lng*mx_rw) = rmdi &! ... srf, dms
  , soot_em     (mx_rw_lng*mx_rw) = rmdi &! ... srf, soot
  , bmass_em    (mx_rw_lng*mx_rw) = rmdi &! ... srf, biomass
  , ocff_em     (mx_rw_lng*mx_rw) = rmdi &! ... srf, OC fossil-fuel
  , soot_hilem  (mx_rw_lng*mx_rw) = rmdi &! ... High lvl, soot
  , bmass_hilem (mx_rw_lng*mx_rw) = rmdi &! ... High lvl, biomass
  , ocff_hilem  (mx_rw_lng*mx_rw) = rmdi &! ... High lvl, OC fossil-fuel
  , soot        (mx_rw_lng*mx_rw) = rmdi  ! Snow soot content

  REAL, PRIVATE ::                                 &
    so2         (mx_rw_lng*mx_rw*mx_wet_lv) = rmdi &
  , so4_aitken  (mx_rw_lng*mx_rw*mx_wet_lv) = rmdi &
  , so4_accu    (mx_rw_lng*mx_rw*mx_wet_lv) = rmdi &
  , so4_diss    (mx_rw_lng*mx_rw*mx_wet_lv) = rmdi &
  , aerosol     (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , aerosol_em  (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , nh3         (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , dms         (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , soot_new    (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , soot_aged   (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , soot_cld    (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , bmass_new   (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , bmass_aged  (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , bmass_cld   (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , ocff_new    (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , ocff_aged   (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , ocff_cld    (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , nitr_acc    (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , nitr_diss   (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , co2         (mx_rw_lng*mx_rw*mx_mod_lv) = rmdi &
  , ozone       (mx_rw_lng*mx_rw*mx_o3_lv)  = rmdi &
  , free_tracers(mx_rw_lng*mx_rw*mx_tr_lv*mx_tr_vars) = rmdi


  ! Vertical correlation coeff.for ...
  REAL, PRIVATE ::                            &
    cort           (mx_rw_lng*mx_rw)  = rmdi  &! ... Temperature
  , cord           (mx_rw_lng*mx_rw)  = rmdi  &! ... Dew pt. depression
  , corvn          (mx_rw_lng*mx_rw)  = rmdi  &! ... Velocity vn
  , corw           (mx_rw_lng*mx_rw)  = rmdi   ! ... Vertical velocity

  REAL, PRIVATE ::                            &
    cclwp          (mx_rw_lng*mx_rw)  = rmdi  &! Condensed water path (kg/m2)
  , orog           (mx_rw_lng*mx_rw)  = rmdi  &! Orographic height (m)
  , ddmfx          (mx_rw_lng*mx_rw)  = rmdi  &! Downdraught mass flux
                                               ! at cloud-base
  , zh             (mx_rw_lng*mx_rw)  = rmdi  &! Top of boundary layer
                                               ! height above surface (m)
  , dOLR_rts       (mx_rw_lng*mx_rw)  = rmdi  &! TOA - surface upward LW
  , fland_ctile    (mx_rw_lng*mx_rw)  = rmdi  &! Land fraction on land points
  , land_alb       (mx_rw_lng*mx_rw)  = rmdi  &! Mean land albedo
  , sice_alb       (mx_rw_lng*mx_rw)  = rmdi  &! Sea-ice albedo
  , tstar_land     (mx_rw_lng*mx_rw)  = rmdi  &! Srf temperature, land (K)
  , tstar_sea      (mx_rw_lng*mx_rw)  = rmdi  &! Srf temperature, open sea (K)
  , tstar_sice     (mx_rw_lng*mx_rw)  = rmdi  &! Srf temperature, sea-ice (K)
  , co2_emits      (mx_rw_lng*mx_rw)  = rmdi  &
  , co2flux        (mx_rw_lng*mx_rw)  = rmdi  &
  , sum_eng_fluxes (mx_rw_lng*mx_rw)  = rmdi  &! Sum atmosphere fluxes
  , sum_moist_flux (mx_rw_lng*mx_rw)  = rmdi  &! Sum moist fluxes
  , tscrndcl_ssi   (mx_rw_lng*mx_rw)  = rmdi  &! Decoupled screen-level temp.
                                               ! over sea and sea-ice
  , tscrndcl_tile  (mx_nlnd*mx_ntile) = rmdi  &!
  , tstbtrans      (mx_rw_lng*mx_rw)  = rmdi   ! Time since the last
                                               ! transition to stability at
                                               ! the surface


  !---------------------------------------------------------------------------
  ! Allocatable arrays:
  ! Used when intepolating namelist forcing profiles from namelist resolution
  ! to SCM resolution
  !---------------------------------------------------------------------------
  REAL, ALLOCATABLE ::     &
    nml_so2        (:,:,:) &
  , nml_so4_aitken (:,:,:) &
  , nml_so4_accu   (:,:,:) &
  , nml_so4_diss   (:,:,:) &
  , nml_aerosol    (:,:,:) &
  , nml_aerosol_em (:,:,:) &
  , nml_nh3        (:,:,:) &
  , nml_dms        (:,:,:) &
  , nml_soot_new   (:,:,:) &
  , nml_soot_aged  (:,:,:) &
  , nml_soot_cld   (:,:,:) &
  , nml_bmass_new  (:,:,:) &
  , nml_bmass_aged (:,:,:) &
  , nml_bmass_cld  (:,:,:) &
  , nml_ocff_new   (:,:,:) &
  , nml_ocff_aged  (:,:,:) &
  , nml_ocff_cld   (:,:,:) &
  , nml_co2        (:,:,:) &
  , nml_nitr_acc   (:,:,:) &
  , nml_nitr_diss  (:,:,:) &
  , nml_ozone      (:,:,:) &
  , nml_free_tracers(:,:,:,:)


  !---------------------------------------------------------------------------
  ! Define namelist
  !---------------------------------------------------------------------------
  NAMELIST/rundata/                                                           &
    modelid, exname_in, exname_out, change_clim, dump_step, dump_days         &
  , resdump_days, runno_in, runno_out, ndayin, nminin, nsecin, ntrad1         &
  , min_trop_level, max_trop_level, timestep, co2start, co2end, co2rate, ntml &
  , nbdsc, ntdsc, albsoil, so2_em, nh3_em, dms_em, soot_em, bmass_em, ocff_em &
  , soot_hilem, bmass_hilem, ocff_hilem, soot, cort, cord, corvn, corw, cclwp &
  , orog, ddmfx, zh, dolr_rts, fland_ctile, land_alb, sice_alb, tstar_land    &
  , tstar_sea, tstar_sice, co2_emits, co2flux, sum_eng_fluxes, sum_moist_flux &
  , so2, so4_aitken, so4_accu, so4_diss, aerosol, aerosol_em, nh3, dms        &
  , soot_new, soot_cld, soot_aged, bmass_new, bmass_aged, bmass_cld, ocff_new &
  , ocff_aged, ocff_cld, nitr_acc, nitr_diss, co2, ozone, tscrndcl_ssi        &
  , tstbtrans, tscrndcl_tile, free_tracers, albobs_sw, albobs_vis, albobs_nir

  PRIVATE :: rundata

!=============================================================================
CONTAINS

  SUBROUTINE read_rundata

!    USE parkind1, ONLY: jpim, jprb
!    USE yomhook,  ONLY: lhook, dr_hook

    IMPLICIT NONE

    ! Dr Hook
    !=============================================================
!    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
!    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    ! Local variables
    INTEGER :: istatus
    INTEGER :: icode

    CHARACTER(LEN=11), PARAMETER :: routinename='read_rundata'

    IF (lhook) CALL dr_hook('READ_RUNDATA',zhook_in,zhook_handle)

    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(6,*) ' Error opening file on unit 10 from '//routinename
      WRITE(6,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) ' IOstat =', istatus
    END IF

    READ  (10, rundata)
    CLOSE (10)

    ! Reshape single level variables
    !==========================================================================
    IF (nlnd > 0) THEN
      IF (albsoil(1)        /= rmdi) scm_albsoil        &
           = RESHAPE( albsoil,       (/nlnd/))
      IF (albobs_sw(1)        /= rmdi) scm_albobs_sw    &
           = RESHAPE( albobs_sw,     (/nlnd/))
      IF (albobs_vis(1)        /= rmdi) scm_albobs_vis  &
           = RESHAPE( albobs_vis,    (/nlnd/))
      IF (albobs_nir(1)        /= rmdi) scm_albobs_nir  &
           = RESHAPE( albobs_nir,    (/nlnd/))
      IF (tscrndcl_tile(1)  /= rmdi) scm_tscrndcl_tile  &
           = RESHAPE( tscrndcl_tile, (/nlnd,ntile/))
    END IF

    IF (ntml(1)     /= rmdi) scm_ntml     = RESHAPE( ntml,     (/rw_lng,rw/))
    IF (ntdsc(1)    /= rmdi) scm_ntdsc    = RESHAPE( ntdsc,    (/rw_lng,rw/))
    IF (nbdsc(1)    /= rmdi) scm_nbdsc    = RESHAPE( nbdsc,    (/rw_lng,rw/))
    IF (so2_em(1)   /= rmdi) scm_so2_em   = RESHAPE( so2_em,   (/rw_lng,rw/))
    IF (nh3_em(1)   /= rmdi) scm_nh3_em   = RESHAPE( nh3_em,   (/rw_lng,rw/))
    IF (dms_em(1)   /= rmdi) scm_dms_em   = RESHAPE( dms_em,   (/rw_lng,rw/))
    IF (soot_em(1)  /= rmdi) scm_soot_em  = RESHAPE( soot_em,  (/rw_lng,rw/))
    IF (bmass_em(1) /= rmdi) scm_bmass_em = RESHAPE( bmass_em, (/rw_lng,rw/))
    IF (ocff_em(1)  /= rmdi) scm_ocff_em  = RESHAPE( ocff_em,  (/rw_lng,rw/))
    IF (soot(1)     /= rmdi) scm_soot     = RESHAPE( soot,     (/rw_lng,rw/))
    IF (cort(1)     /= rmdi) scm_cort     = RESHAPE( cort,     (/rw_lng,rw/))
    IF (cord(1)     /= rmdi) scm_cord     = RESHAPE( cord,     (/rw_lng,rw/))
    IF (corvn(1)    /= rmdi) scm_corvn    = RESHAPE( corvn,    (/rw_lng,rw/))
    IF (corw(1)     /= rmdi) scm_corw     = RESHAPE( corw,     (/rw_lng,rw/))
    IF (cclwp(1)    /= rmdi) scm_cclwp    = RESHAPE( cclwp,    (/rw_lng,rw/))
    IF (orog(1)     /= rmdi) scm_orog     = RESHAPE( orog,     (/rw_lng,rw/))
    IF (ddmfx(1)    /= rmdi) scm_ddmfx    = RESHAPE( ddmfx,    (/rw_lng,rw/))
    IF (zh(1)       /= rmdi) scm_zh       = RESHAPE( zh,       (/rw_lng,rw/))

    IF (soot_hilem(1)  /= rmdi) scm_soot_hilem                                &
                                      = RESHAPE( soot_hilem,  (/rw_lng,rw/))

    IF (bmass_hilem(1) /= rmdi) scm_bmass_hilem                               &
                                      = RESHAPE( bmass_hilem, (/rw_lng,rw/))

    IF (ocff_hilem(1)  /= rmdi) scm_ocff_hilem                                &
                                      = RESHAPE( ocff_hilem,    (/rw_lng,rw/))

    IF (dolr_rts(1)    /= rmdi) scm_dolr_rts                                  &
                                      = RESHAPE( dolr_rts,      (/rw_lng,rw/))

    IF (fland_ctile(1) /= rmdi) scm_fland_ctile                               &
                                      = RESHAPE( fland_ctile,   (/rw_lng,rw/))

    IF (land_alb(1)    /= rmdi) scm_land_alb                                  &
                                      = RESHAPE( land_alb,      (/rw_lng,rw/))

    IF (sice_alb(1)    /= rmdi) scm_sice_alb                                  &
                                      = RESHAPE( sice_alb,      (/rw_lng,rw/))

    IF (tstar_land(1)  /= rmdi) scm_tstar_land                                &
                                      = RESHAPE( tstar_land,    (/rw_lng,rw/))

    IF (tstar_sea(1)   /= rmdi) scm_tstar_sea                                 &
                                      = RESHAPE( tstar_sea,     (/rw_lng,rw/))

    IF (tstar_sice(1)  /= rmdi) scm_tstar_sice                                &
                                      = RESHAPE( tstar_sice,    (/rw_lng,rw/))

    IF (co2_emits(1)   /= rmdi) scm_co2_emits                                 &
                                      = RESHAPE( co2_emits,     (/rw_lng,rw/))

    IF (co2flux(1)     /= rmdi) scm_co2flux                                   &
                                      = RESHAPE( co2flux,       (/rw_lng,rw/))

    IF (sum_eng_fluxes(1) /= rmdi) scm_sum_eng_fluxes                         &
                                      = RESHAPE( sum_eng_fluxes,(/rw_lng,rw/))

    IF (sum_moist_flux(1) /= rmdi) scm_sum_moist_flux                         &
                                      = RESHAPE( sum_moist_flux,(/rw_lng,rw/))

    IF (tscrndcl_ssi(1)   /= rmdi) scm_tscrndcl_ssi                           &
                                      = RESHAPE( tscrndcl_ssi,  (/rw_lng,rw/))

    IF (tstbtrans(1)      /= rmdi) scm_tstbtrans                              &
                                      = RESHAPE( tstbtrans,     (/rw_lng,rw/))


    ! Reshape/Interpolate vertical profiles
    !==========================================================================
    IF (interpolate) THEN

      CALL alloc_rundata
      CALL interp_rundata
      CALL dealloc_rundata

    ELSE

      ! Reshape RUNDATA variables
      !==========================
      IF (ozone(1) /= rmdi)                                                   &
        scm_ozone  = RESHAPE( ozone, (/rw_lng,rw,o3_lv/) )

      IF (co2(1)   /= rmdi)                                                   &
        scm_co2    = RESHAPE( co2,   (/rw_lng,rw,nmod_lv/) )

      IF (so2(1)   /= rmdi)                                                   &
        scm_so2    = RESHAPE( so2,   (/rw_lng,rw,nwet_lv/) )

      IF (dms(1)   /= rmdi)                                                   &
        scm_dms    = RESHAPE( dms,   (/rw_lng,rw,nmod_lv/) )

      IF (nh3(1)   /= rmdi)                                                   &
        scm_nh3    = RESHAPE( nh3,   (/rw_lng,rw,nmod_lv/) )

      IF (so4_aitken(1) /= rmdi)                                              &
        scm_so4_aitken = RESHAPE( so4_aitken, (/rw_lng,rw,nwet_lv/) )

      IF (so4_accu(1)   /= rmdi)                                              &
        scm_so4_accu   = RESHAPE( so4_accu,   (/rw_lng,rw,nwet_lv/) )

      IF (so4_diss(1)   /= rmdi)                                              &
        scm_so4_diss   = RESHAPE( so4_diss,   (/rw_lng,rw,nwet_lv/) )

      IF (aerosol(1)    /= rmdi)                                              &
        scm_aerosol    = RESHAPE( aerosol,    (/rw_lng,rw,nmod_lv/) )

      IF (aerosol_em(1) /= rmdi)                                              &
        scm_aerosol_em = RESHAPE( aerosol_em, (/rw_lng,rw,nmod_lv/) )

      IF (free_tracers(1) /= rmdi)                                            &
        scm_free_tracers = RESHAPE(free_tracers,(/rw_lng,rw,ntr_lv,ntr_var/) )

      IF (soot_new(1)   /= rmdi)                                              &
        scm_soot_new   = RESHAPE( soot_new,   (/rw_lng,rw,nmod_lv/) )

      IF (soot_aged(1)  /= rmdi)                                              &
        scm_soot_aged  = RESHAPE( soot_aged,  (/rw_lng,rw,nmod_lv/) )

      IF (soot_cld(1)   /= rmdi)                                              &
        scm_soot_cld   = RESHAPE( soot_cld,   (/rw_lng,rw,nmod_lv/) )

      IF (bmass_new(1)  /= rmdi)                                              &
        scm_bmass_new  = RESHAPE( bmass_new,  (/rw_lng,rw,nmod_lv/) )

      IF (bmass_aged(1) /= rmdi)                                              &
        scm_bmass_aged = RESHAPE( bmass_aged, (/rw_lng,rw,nmod_lv/) )

      IF (bmass_cld(1)  /= rmdi)                                              &
        scm_bmass_cld  = RESHAPE( bmass_cld,  (/rw_lng,rw,nmod_lv/) )

      IF (ocff_new(1)   /= rmdi)                                              &
        scm_ocff_new   = RESHAPE( ocff_new,   (/rw_lng,rw,nmod_lv/) )

      IF (ocff_aged(1)  /= rmdi)                                              &
        scm_ocff_aged  = RESHAPE( ocff_aged,  (/rw_lng,rw,nmod_lv/) )

      IF (ocff_cld(1)   /= rmdi)                                              &
        scm_ocff_cld   = RESHAPE( ocff_cld,   (/rw_lng,rw,nmod_lv/) )

      IF (nitr_acc(1)   /= rmdi)                                              &
        scm_nitr_acc   = RESHAPE( nitr_acc,   (/rw_lng,rw,nmod_lv/) )

      IF (nitr_diss(1)  /= rmdi)                                              &
        scm_nitr_diss  = RESHAPE( nitr_diss,  (/rw_lng,rw,nmod_lv/) )

    END IF

    IF (lhook) CALL dr_hook('READ_RUNDATA',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_rundata

!-----------------------------------------------------------------------------

  SUBROUTINE interp_rundata

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    ! Local variables
    INTEGER :: i,j,k,n

    IF (lhook) CALL dr_hook('INTERP_RUNDATA',zhook_in,zhook_handle)

    !-------------------------------------------------------------------------
    ! Ozone
    !-------------------------------------------------------------------------
    IF (ozone(1) /= rmdi) THEN
      nml_ozone = RESHAPE( ozone, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,      nml_ozone(i,j,:)                                 &
            , z_th(1:o3_lv), scm_ozone(i,j,1:o3_lv) )
        END DO
      END DO

      ! Assume ozone is constant from lowest obs level to level 1
      DO k=1, o3_lv
        IF (z_th(k) <= nml_z_th(1)) scm_ozone(:,:,k) = nml_ozone(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! CO2
    !-------------------------------------------------------------------------
    IF (co2(1) /= rmdi) THEN
      nml_co2 = RESHAPE( co2, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_co2(i,j,:)                                 &
            , z_th(1:nmod_lv), scm_co2(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume co2 is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1)) scm_co2(:,:,k) = nml_co2(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! SO2
    !-------------------------------------------------------------------------
    IF (so2(1) /= rmdi) THEN
      nml_so2 = RESHAPE( so2, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_so2(i,j,:)                                 &
            , z_th(1:nwet_lv), scm_so2(i,j,1:nwet_lv) )
        END DO
      END DO

      ! Assume so2 is constant from lowest obs level to level 1
      DO k=1, nwet_lv
        IF (z_th(k) <= nml_z_th(1)) scm_so2(:,:,k) = nml_so2(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! DMS
    !-------------------------------------------------------------------------
    IF (dms(1) /= rmdi) THEN
      nml_dms = RESHAPE( dms, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_dms(i,j,:)                                 &
            , z_th(1:nmod_lv), scm_dms(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume dms is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1)) scm_dms(:,:,k) = nml_dms(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! NH3
    !-------------------------------------------------------------------------
    IF (nh3(1) /= rmdi) THEN
      nml_nh3 = RESHAPE( nh3, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_nh3(i,j,:)                                 &
            , z_th(1:nmod_lv), scm_nh3(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume nh3 is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1)) scm_nh3(:,:,k) = nml_nh3(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! SO4 Aitken
    !-------------------------------------------------------------------------
    IF (so4_aitken(1) /= rmdi) THEN
      nml_so4_aitken = RESHAPE( so4_aitken, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_so4_aitken(i,j,:)                          &
            , z_th(1:nwet_lv), scm_so4_aitken(i,j,1:nwet_lv) )
        END DO
      END DO

      ! Assume so4_aitken is constant from lowest obs level to level 1
      DO k=1, nwet_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_so4_aitken(:,:,k) = nml_so4_aitken(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! SO4 Accu
    !-------------------------------------------------------------------------
    IF (so4_accu(1) /= rmdi) THEN
      nml_so4_accu = RESHAPE( so4_accu, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_so4_accu(i,j,:)                            &
            , z_th(1:nwet_lv), scm_so4_accu(i,j,1:nwet_lv) )
        END DO
      END DO

      ! Assume so4_accu is constant from lowest obs level to level 1
      DO k=1, nwet_lv
        IF (z_th(k) <= nml_z_th(1)) scm_so4_accu(:,:,k) = nml_so4_accu(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! SO4 Diss
    !-------------------------------------------------------------------------
    IF (so4_diss(1) /= rmdi) THEN
      nml_so4_diss = RESHAPE( so4_diss, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_so4_diss(i,j,:)                            &
            , z_th(1:nwet_lv), scm_so4_diss(i,j,1:nwet_lv) )
        END DO
      END DO

      ! Assume so4_diss is constant from lowest obs level to level 1
      DO k=1, nwet_lv
        IF (z_th(k) <= nml_z_th(1)) scm_so4_diss(:,:,k) = nml_so4_diss(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Aerosol
    !-------------------------------------------------------------------------
    IF (aerosol(1) /= rmdi) THEN
      nml_aerosol = RESHAPE( aerosol, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_aerosol(i,j,:)                             &
            , z_th(1:nmod_lv), scm_aerosol(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume aerosol is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1)) scm_aerosol(:,:,k) = nml_aerosol(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Aerosol emissions
    !-------------------------------------------------------------------------
    IF (aerosol_em(1) /= rmdi) THEN
      nml_aerosol_em = RESHAPE( aerosol_em, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_aerosol_em(i,j,:)                          &
            , z_th(1:nmod_lv), scm_aerosol_em(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume aerosol_em is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_aerosol_em(:,:,k) = nml_aerosol_em(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Free tracers
    !-------------------------------------------------------------------------
    IF (free_tracers(1) /= rmdi) THEN
      nml_free_tracers =                                                      &
                 RESHAPE( free_tracers, (/rw_lng,rw,nml_nmod_lv,ntr_var/) )

      DO n=1, ntr_var
        DO j=1, rw
          DO i=1, rw_lng
            CALL interp1d                                                     &
              ( nml_z_th,       nml_free_tracers(i,j,:,n)                     &
              , z_th(1:ntr_lv), scm_free_tracers(i,j,1:ntr_lv,n) )
          END DO
        END DO

        ! Assume aerosol_em is constant from lowest obs level to level 1
        DO k=1, ntr_lv
          IF (z_th(k) <= nml_z_th(1))                                         &
            scm_free_tracers(:,:,k,n) = nml_free_tracers(:,:,1,n)
        END DO
      END DO   ! loop over n (ntr_var)
    END IF


    !-------------------------------------------------------------------------
    ! Soot - New
    !-------------------------------------------------------------------------
    IF (soot_new(1) /= rmdi) THEN
      nml_soot_new = RESHAPE( soot_new, (/rw_lng,rw,nml_nmod_lv/) )
      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_soot_new(i,j,:)                            &
            , z_th(1:nmod_lv), scm_soot_new(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume soot_new is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1)) scm_soot_new(:,:,k) = nml_soot_new(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Soot - Aged
    !-------------------------------------------------------------------------
    IF (soot_aged(1) /= rmdi) THEN
      nml_soot_aged = RESHAPE( soot_aged, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_soot_aged(i,j,:)                           &
            , z_th(1:nmod_lv), scm_soot_aged(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume soot_aged is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_soot_aged(:,:,k) = nml_soot_aged(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Soot - Cloud
    !-------------------------------------------------------------------------
    IF (soot_cld(1) /= rmdi) THEN
      nml_soot_cld = RESHAPE( soot_cld, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_soot_cld(i,j,:)                            &
            , z_th(1:nmod_lv), scm_soot_cld(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume soot_cld is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1)) scm_soot_cld(:,:,k) = nml_soot_cld(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Biomass - New
    !-------------------------------------------------------------------------
    IF (bmass_new(1) /= rmdi)  THEN
      nml_bmass_new = RESHAPE( bmass_new, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_bmass_new(i,j,:)                           &
            , z_th(1:nmod_lv), scm_bmass_new(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume bmass_new is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_bmass_new(:,:,k) = nml_bmass_new(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Biomass - Aged
    !-------------------------------------------------------------------------
    IF (bmass_aged(1) /= rmdi) THEN
      nml_bmass_aged = RESHAPE( bmass_aged, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_bmass_aged(i,j,:)                          &
            , z_th(1:nmod_lv), scm_bmass_aged(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume bmass_aged is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_bmass_aged(:,:,k) = nml_bmass_aged(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Biomass - Cloud
    !-------------------------------------------------------------------------
    IF (bmass_cld(1) /= rmdi) THEN
      nml_bmass_cld = RESHAPE( bmass_cld, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_bmass_cld(i,j,:)                           &
            , z_th(1:nmod_lv), scm_bmass_cld(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume bmass_cld is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_bmass_cld(:,:,k) = nml_bmass_cld(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Organic-Carbon Fossil-Fuels - New
    !-------------------------------------------------------------------------
    IF (ocff_new(1) /= rmdi) THEN
      nml_ocff_new = RESHAPE( ocff_new, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_ocff_new(i,j,:)                            &
            , z_th(1:nmod_lv), scm_ocff_new(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume ocff_new is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_ocff_new(:,:,k) = nml_ocff_new(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Organic-Carbon Fossil-Fuels - Aged
    !-------------------------------------------------------------------------
    IF (ocff_aged(1) /= rmdi) THEN
      nml_ocff_aged = RESHAPE( ocff_aged, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_ocff_aged(i,j,:)                           &
            , z_th(1:nmod_lv), scm_ocff_aged(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume ocff_aged is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_ocff_aged(:,:,k) = nml_ocff_aged(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Organic-Carbon Fossil-Fuels - Cloud
    !-------------------------------------------------------------------------
    IF (ocff_cld(1) /= rmdi) THEN
      nml_ocff_cld = RESHAPE( ocff_cld, (/rw_lng,rw,nml_nmod_lv/) )
      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_ocff_cld(i,j,:)                            &
            , z_th(1:nmod_lv), scm_ocff_cld(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume ocff_cld is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_ocff_cld(:,:,k) = nml_ocff_cld(:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! Nitrate - accumulated
    !-------------------------------------------------------------------------
    IF (nitr_acc(1) /= rmdi) THEN
      nml_nitr_acc = RESHAPE( nitr_acc, (/rw_lng,rw,nml_nmod_lv/) )
      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_nitr_acc(i,j,:)                            &
            , z_th(1:nmod_lv), scm_nitr_acc(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume nitr_acc is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_nitr_acc(:,:,k) = nml_nitr_acc(:,:,1)
      END DO
    END IF

    !-------------------------------------------------------------------------
    ! Nitrate - diss
    !-------------------------------------------------------------------------
    IF (nitr_diss(1) /= rmdi) THEN
      nml_nitr_diss = RESHAPE( nitr_diss, (/rw_lng,rw,nml_nmod_lv/) )
      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th,        nml_nitr_diss(i,j,:)                           &
            , z_th(1:nmod_lv), scm_nitr_diss(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Assume nitr_diss is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_nitr_diss(:,:,k) = nml_nitr_diss(:,:,1)
      END DO
    END IF


    IF (lhook) CALL dr_hook('INTERP_RUNDATA',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE interp_rundata

!-----------------------------------------------------------------------------

  SUBROUTINE alloc_rundata

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('ALLOC_RUNDATA',zhook_in,zhook_handle)

    ALLOCATE                                     &
      ( nml_so2        (rw_lng,rw,nml_nmod_lv)   &
      , nml_so4_aitken (rw_lng,rw,nml_nmod_lv)   &
      , nml_so4_accu   (rw_lng,rw,nml_nmod_lv)   &
      , nml_so4_diss   (rw_lng,rw,nml_nmod_lv)   &
      , nml_nh3        (rw_lng,rw,nml_nmod_lv)   &
      , nml_dms        (rw_lng,rw,nml_nmod_lv)   &
      , nml_aerosol    (rw_lng,rw,nml_nmod_lv)   &
      , nml_aerosol_em (rw_lng,rw,nml_nmod_lv)   &
      , nml_bmass_new  (rw_lng,rw,nml_nmod_lv)   &
      , nml_bmass_aged (rw_lng,rw,nml_nmod_lv)   &
      , nml_bmass_cld  (rw_lng,rw,nml_nmod_lv)   &
      , nml_soot_new   (rw_lng,rw,nml_nmod_lv)   &
      , nml_soot_aged  (rw_lng,rw,nml_nmod_lv)   &
      , nml_soot_cld   (rw_lng,rw,nml_nmod_lv)   &
      , nml_ocff_new   (rw_lng,rw,nml_nmod_lv)   &
      , nml_ocff_aged  (rw_lng,rw,nml_nmod_lv)   &
      , nml_ocff_cld   (rw_lng,rw,nml_nmod_lv)   &
      , nml_nitr_acc   (rw_lng,rw,nml_nmod_lv)   &
      , nml_nitr_diss  (rw_lng,rw,nml_nmod_lv)   &
      , nml_co2        (rw_lng,rw,nml_nmod_lv)   &
      , nml_ozone      (rw_lng,rw,nml_nmod_lv)   &
      , nml_free_tracers(rw_lng,rw,nml_nmod_lv,ntr_var) )

    IF (lhook) CALL dr_hook('ALLOC_RUNDATA',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE alloc_rundata

!-----------------------------------------------------------------------------

  SUBROUTINE dealloc_rundata

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('DEALLOC_RUNDATA',zhook_in,zhook_handle)

    DEALLOCATE           &
      ( nml_free_tracers &
      , nml_ozone        &
      , nml_co2          &
      , nml_nitr_diss    &
      , nml_nitr_acc     &
      , nml_ocff_cld     &
      , nml_ocff_aged    &
      , nml_ocff_new     &
      , nml_soot_cld     &
      , nml_soot_aged    &
      , nml_soot_new     &
      , nml_bmass_cld    &
      , nml_bmass_aged   &
      , nml_bmass_new    &
      , nml_aerosol_em   &
      , nml_aerosol      &
      , nml_dms          &
      , nml_nh3          &
      , nml_so4_diss     &
      , nml_so4_accu     &
      , nml_so4_aitken   &
      , nml_so2        )

    IF (lhook) CALL dr_hook('DEALLOC_RUNDATA',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE dealloc_rundata

!=============================================================================
END MODULE s_rundata
