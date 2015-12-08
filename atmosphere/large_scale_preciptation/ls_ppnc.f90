! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     SUBROUTINES LS_PPN and LS_PPNC------------------------------------
!    Purpose:
!            LS_PPN and LS_PPNC:
!             Calculate large-scale (dynamical) precipitation.
!             LS_PPNC is the gather/scatter routine which then
!             calls LSP_ICE.
!    Note: in all cases, level counters (incl subscripts) run from 1
!          (lowest model layer) to Q_LEVELS (topmost "wet" model
!          layer) - it is assumed that the bottom Q_LEVELS layers are
!          the "wet" layers.
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Documentation: UM Documentation Paper 26.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation
MODULE ls_ppnc_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE ls_ppnc( level, ix, n, n_iterations,                         &
 lsrain,lssnow,lssnow2,lsgraup, droplet_flux,                           &
 cf,cfl,cff,                                                            &
 qcf,qcl,t, qcf2,qrain,qgraup,                                          &
 n_drop_tpr, n_drop_out,                                                & 
 aerosol, land_fract,                                                   &
!--------------------------------------------------------
! Layer thicknesses and variables passed from layer above
!--------------------------------------------------------
 q, p_theta_levels, layer_thickness,                                    &
 deltaz, rhodz_dry, rhodz_moist,                                        &
 rhc_row_length, rhc_rows, bland, rhcrit,                               &
 vfall, vfall2, vfall_rain, vfall_graup,                                &
 frac_ice_above, cttemp, rainfrac, lsiter, niter_bs,                    &
 uk, vk, ukp1, vkp1, onelevel                                           & 
 )


  ! General atmosphere modules
  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE earth_constants_mod,   ONLY: g

  USE murk_inputs_mod, ONLY: l_murk
  
  ! Microphysics modules
  USE mphys_diags_mod,       ONLY: l_aggfr_diag,                        &
                                   l_psdep_diag, l_psaut_diag,          &
                                   l_psacw_diag, l_psacr_diag,          &
                                   l_psaci_diag, l_psmlt_diag,          &
                                   l_psmltevp_diag, l_praut_diag,       &
                                   l_pracw_diag, l_prevp_diag,          &
                                   l_pgaut_diag, l_pgacw_diag,          &
                                   l_pgacs_diag, l_pgmlt_diag,          &
                                   l_pifrw_diag, l_piprm_diag,          &
                                   l_pidep_diag, l_piacw_diag,          &
                                   l_piacr_diag, l_pimlt_diag,          &
                                   l_pimltevp_diag, l_pifall_diag,      &
                                   l_psfall_diag, l_prfall_diag,        &
                                   l_pgfall_diag, l_plset_diag,         &
                                   l_plevpset_diag, l_pifrr_diag,       &
                                   frac_agg, psdep, psaut,              &
                                   psacw, psacr, psaci, psmlt,          &
                                   psmltevp, praut, pracw, prevp,       &
                                   pgaut, pgacw, pgacs, pgmlt, pifrw,   &
                                   piprm, pidep, piacw, piacr, pimlt,   &
                                   pimltevp, pifall, psfall, prfall,    &
                                   pgfall, plset, plevpset, pifrr

  USE mphys_inputs_mod,      ONLY:  l_mcr_qcf2, l_mcr_qgraup,            &
                                   l_mcr_qrain

  USE cloud_inputs_mod,      ONLY:  l_pc2
                                  
  ! Grid bounds module

  USE atm_fields_bounds_mod, ONLY: qdims, pdims, tdims

  USE trignometric_mod,      ONLY: cos_theta_latitude,                  &
                                   fv_cos_theta_latitude

  USE level_heights_mod,     ONLY: r_theta_levels

  ! Dr Hook modules
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim

  ! Large scale precipitation modules
  USE lsp_ice_mod,           ONLY: lsp_ice
  USE lsp_scav_mod,          ONLY: lsp_scav

  IMPLICIT NONE

  INTEGER ::                                                            &
   level,                                                               &
                ! level number
   n,                                                                   &
                ! IN Number of points where pptn non-zero from above
!                    or where CF>CFMIN
   n_iterations,                                                        &
                     ! Number of iterations for iterative melting
   lsiter,                                                              &
                ! Number of total iterations of microphysics (inner loop)
   niter_bs,                                                            &
                ! Number of total iterations of microphysics (outer loop)
   rhc_row_length,rhc_rows,                                             &
   ix ( (qdims%i_end - qdims%i_start + 1) *                             &
        (qdims%j_end - qdims%j_start + 1) , 2 )
                                ! IN gather/scatter index
  REAL ::                                                               &
   cf (qdims%i_start:qdims%i_end, qdims%j_start:qdims%j_end),           &
                              ! IN Cloud fraction.
   cfl(qdims%i_start:qdims%i_end, qdims%j_start:qdims%j_end),           &
                              ! IN Cloud liquid fraction.
   cff(qdims%i_start:qdims%i_end, qdims%j_start:qdims%j_end),           &
                              ! IN Cloud ice fraction.
! CF, CFL and CFF are IN/OUT if the PC2 cloud scheme is in use.

    p_theta_levels( qdims%i_start:qdims%i_end,                          &
                    qdims%j_start:qdims%j_end ),                        &
    layer_thickness(qdims%i_start:qdims%i_end,                          &
                    qdims%j_start:qdims%j_end ),                        &
                                          ! IN thickness of layer (Pa)
    deltaz( qdims%i_start:qdims%i_end,                                  &
            qdims%j_start:qdims%j_end ),                                &
                                          ! IN thickness of layer (m)
    rhodz_dry(qdims%i_start:qdims%i_end,                                &
              qdims%j_start:qdims%j_end ),                              &
                                          ! Dry air density
                                          ! * layer thickness (kg m-2)
    rhodz_moist(qdims%i_start:qdims%i_end,                              &
                qdims%j_start:qdims%j_end ),                            &
                                          ! Moist air density
                                          ! * layer thickness (kg m-2)

    rhcrit(rhc_row_length,rhc_rows)
!                       IN Critical humidity for cloud formation.


  LOGICAL, INTENT(IN) ::                                                &
       bland(qdims%i_start:qdims%i_end,                                 &
             qdims%j_start:qdims%j_end )
                          !IN Land/sea mask
  REAL, INTENT(INOUT) ::                                                &
       q(qdims%i_start:qdims%i_end,                                     &
         qdims%j_start:qdims%j_end),                                    &
                            ! INOUT Specific humidity (kg water/kg air).
        qcf(qdims%i_start:qdims%i_end,                                  &
            qdims%j_start:qdims%j_end),                                 &
                            ! INOUT Cloud ice (kg per kg air).
        qcl(qdims%i_start:qdims%i_end,                                  &
            qdims%j_start:qdims%j_end),                                 &
                            ! INOUT Cloud liquid water (kg per kg air).
        qcf2(qdims%i_start:qdims%i_end,                                 &
             qdims%j_start:qdims%j_end),                                &
                            ! INOUT Cloud ice2 (kg per kg air).
        qrain(qdims%i_start:qdims%i_end,                                &
              qdims%j_start:qdims%j_end),                               &
                            ! INOUT Rain water (kg per kg air).
        qgraup(qdims%i_start:qdims%i_end,                               &
               qdims%j_start:qdims%j_end),                              &
                            ! INOUT Graupel water (kg per kg air).
        t(qdims%i_start:qdims%i_end,                                    &
          qdims%j_start:qdims%j_end),                                   &
                            ! INOUT Temperature (K).
        aerosol(qdims%i_start:qdims%i_end,                              &
                qdims%j_start:qdims%j_end),                             &
                            ! INOUT Aerosol (K).
        lsrain(qdims%i_start:qdims%i_end,                               &
               qdims%j_start:qdims%j_end),                              &
                            !INOUT Surface rainfall rate (kg m^-2 s^-1).
        lssnow(qdims%i_start:qdims%i_end,                               &
               qdims%j_start:qdims%j_end),                              &
                            !INOUT Surface snowfall rate (kg m^-2 s^-1).
        lssnow2(qdims%i_start:qdims%i_end,                              &
                qdims%j_start:qdims%j_end),                             &
                            !INOUT layer snowfall rate (kg m^-2 s^-1).
        lsgraup(qdims%i_start:qdims%i_end,                              &
                qdims%j_start:qdims%j_end),                             &
                            !INOUT layer graupelfall rate (kg m^-2 s^-1)
        droplet_flux(qdims%i_start:qdims%i_end,                         &
                     qdims%j_start:qdims%j_end),                        &
                            !INOUT water droplet flux / kg m^-2 s^-1
        cttemp(qdims%i_start:qdims%i_end,                               &
               qdims%j_start:qdims%j_end),                              &
                            ! INOUT Ice cloud top temperature (K)
        rainfrac(qdims%i_start:qdims%i_end,                             &
                 qdims%j_start:qdims%j_end),                            &
                            ! INOUT Rain fraction.
        frac_ice_above(qdims%i_start:qdims%i_end,                       &
                       qdims%j_start:qdims%j_end),                      &
                            ! INOUT Ice fraction from layer above
   vfall(qdims%i_start:qdims%i_end,                                     &
         qdims%j_start:qdims%j_end),                                    &
                                     ! INOUT fall velocity of ice (m per
   vfall2(qdims%i_start:qdims%i_end,                                    &
          qdims%j_start:qdims%j_end),                                   &
                                      ! INOUT fall vel. of rain (m/s)
   vfall_rain(qdims%i_start:qdims%i_end,                                &
              qdims%j_start:qdims%j_end),                               &
                                      ! INOUT fall vel. of rain (m/s)
   vfall_graup(qdims%i_start:qdims%i_end,                               &
               qdims%j_start:qdims%j_end)
                                      ! INOUT fall vel. of graupel (m/s)

  REAL,INTENT(IN) ::                                                    &
! U and V wind at levels k and k+1 on P grid, so defined using pdims. 
   uk  (pdims%i_start : pdims%i_end,                                    & 
        pdims%j_start : pdims%j_end),                                   & 
   ukp1(pdims%i_start : pdims%i_end,                                    & 
        pdims%j_start : pdims%j_end),                                   & 
   vk  (pdims%i_start : pdims%i_end,                                    & 
        pdims%j_start : pdims%j_end),                                   & 
   vkp1(pdims%i_start : pdims%i_end,                                    & 
        pdims%j_start : pdims%j_end)
! Level we are interested in for r_theta_level
  INTEGER, INTENT(IN) :: onelevel

  REAL ::                                                               &
       land_fract(qdims%i_start:qdims%i_end,                            &
                qdims%j_start:qdims%j_end)

  REAL ::                                                               &
                      !, Intent(IN)
        n_drop_tpr( qdims%i_start : qdims%i_end,                        &
                    qdims%j_start : qdims%j_end )
                               ! Tapered droplet number
  REAL ::                                                               &
                      ! INTENT(INOUT)
        n_drop_out( qdims%i_start : qdims%i_end,                        &
                    qdims%j_start : qdims%j_end )    
                      ! droplet number from autoconversion

!    Workspace usage ---------------------------------------------------

  REAL ::                                                               &
   cf_c(n),                                                             &
                        ! gathered Cloud fraction.
   q_c(n),                                                              &
                         ! gathered Specific humidity (kg water/kg air).
   qcf_c(n),                                                            &
                         ! gathered Cloud ice (kg per kg air).
   qcl_c(n),                                                            &
                         ! gathered Cloud liquid water (kg per kg air).
   qcf2_c(n),                                                           &
                          ! gathered cloud ice2 (kg per kg air).
   qrain_c(n),                                                          &
                           ! gathered rain (kg per kg air).
   qgraup_c(n),                                                         &
                            ! gathered graupel (kg per kg air).
   t_c(n),                                                              &
                         ! gathered Temperature (K).
   uk_c(n),                                                             &
                         ! gathered u wind on level k
   vk_c(n),                                                             &
                         ! gathered v wind on level k
   ukp1_c(n),                                                           &
                         ! gathered u wind on level k+1
   vkp1_c(n),                                                           &
                         ! gathered v wind on level k+1
   r_theta_levels_c(n),                                                 &
                         ! gathered distance from centre of Earth
   g_cos_theta_latitude_c(n),                                           &
                         ! gathered cos of latitude
   aero_c(n),                                                           &
                         ! gathered Aerosol.
   lsrain_c(n),                                                         &
                   !gathered Surface rainfall rate (kg per sq m per s).
   lssnow_c(n),                                                         &
                   !gathered Surface snowfall rate (kg per sq m per s).
   lssnow2_c(n),                                                        &
                    !gathered layer snowfall rate (kg per sq m per s).
   lsgraup_c(n),                                                        &
                    !gathered layer graupel fall rate (kg/sq m/s)
   droplet_flux_c(n),                                                   &
                         ! gathered water droplet flux / kg m-2 s-1
   cttemp_c(n),                                                         &
                              !gathered ice cloud top temperature.
   rainfrac_c(n),                                                       &
                              !gathered rain fraction.
   frac_ice_above_c(n),                                                 &
                              !gathered fraction of ice in layer above
   frac_agg_c(n),                                                       &
                         ! gathered aggregate fraction
   cfl_c(n),                                                            &
                         ! gathered Cloud liquid fraction.
   cff_c(n),                                                            &
                         ! gathered Cloud ice fraction.
   vfall_c(n),                                                          &
                         ! gathered fall velocity (m per s).
   vfall2_c(n),                                                         &
                          ! gathered fall velocity for qcf2 (m per s).
   vfall_rain_c(n),                                                     &
                          ! gathered fall velocity for qcf2 (m per s).
   vfall_graup_c(n),                                                    &
                          ! gathered fall velocity for qcf2 (m per s).
   rhc_c(n),                                                            &
                          ! gathered RH_crit value at points.
   n_drop_tpr_c(n),                                                     &
   n_drop_out_c(n),                                                     &
                          ! gathered droplet numbers
   land_fract_c(n)
                          ! gathered land fraction


  REAL ::                                                               &
   rhodz(n),                                                            &
                      ! WORK Used for air mass p.u.a. in successive
!                            layers.
   deltaz_c(n),                                                         &
                       ! Thickness of layer (m)
   rhodz_dry_c(n),                                                      &
                       ! Dry air density * layer thickness (kg m-2)
   rhodz_moist_c(n),                                                    &
                       ! Moist air density * layer thickness (kg m-2)

   p(n)           ! WORK Used for pressure at successive levels.
  LOGICAL :: bland_c(n)          ! gathered land/sea mask

! Microphysical process rate diagnostics (compressed arrays)
  REAL ::                                                               &
    psdep_c(n),                                                         &
                    ! Deposition of vapour to snow aggregates
    psaut_c(n),                                                         &
                    ! Autoconversion of aggregates from crystals
    psacw_c(n),                                                         &
                    ! Accretion of liq. water by snow aggregates
    psacr_c(n),                                                         &
                    ! Collection of rain by snow aggregates
    psaci_c(n),                                                         &
                    ! Collection of ice crystals by aggregates
    psmlt_c(n),                                                         &
                    ! Melting of snow aggregates
    psmltevp_c(n)  ! Evaporation of melting aggregates
  REAL ::                                                               &
    praut_c(n),                                                         &
                    ! Autoconversion of cloud drops to rain
    pracw_c(n),                                                         &
                    ! Accretion of liq. water by rain
    prevp_c(n)  ! Evaporation of rain
  REAL ::                                                               &
    pgaut_c(n),                                                         &
                    ! Autoconversion of graupel from aggregates
    pgacw_c(n),                                                         &
                    ! Accretion of liq. water by graupel
    pgacs_c(n),                                                         &
                    ! Collection of snow aggregates by graupel
    pgmlt_c(n)  ! Melting of graupel
  REAL ::                                                               &
    pifrw_c(n),                                                         &
                    ! Homogeneous freezing nucleation
    pifrr_c(n),                                                         &
                    ! Homogeneous freezing of rain
    piprm_c(n),                                                         &
                    ! Heterogeneous (primary) nucleation
    pidep_c(n),                                                         &
                    ! Deposition of vapour to ice crystals
    piacw_c(n),                                                         &
                    ! Accretion of liq. water by ice crystals
    piacr_c(n),                                                         &
                    ! Collection of rain by ice crystals
    pimlt_c(n),                                                         &
                    ! Melting of ice crystals
    pimltevp_c(n)  ! Evaporation of melting ice crystals
  REAL ::                                                               &
    pifall_c(n),                                                        &
                    ! Sedimentation of ice crystals
    psfall_c(n),                                                        &
                    ! Sedimentation of aggregates
    prfall_c(n),                                                        &
                    ! Sedimentation of rain
    pgfall_c(n) ! Sedimentation of graupel
  REAL ::                                                               &
    plset_c(n),                                                         &
                    ! Droplet settling of liquid water
    plevpset_c(n) ! Evaporated settled droplets

  INTEGER :: jj,i1,i2,ip

! Physical constants --------------------------------------------------
  REAL :: p1upong
  PARAMETER (                                                           &
   p1upong=1./g                                                         &
                           ! One upon g (sq seconds per m).
  )
!  Define local variables ----------------------------------------------
  INTEGER ::                                                            &
   multrhc,                                                             &
                     ! Zero if (rhc_row_length*rhc_rows) le 1, else 1
   multdrp
                     ! Zero for single profile of droplet number, else 1

  INTEGER :: i,ii,ij,                                                   &
                            ! Loop counters: I - horizontal field index;
            irhi,irhj,                                                  &
                            !  IRHI,IRHJ-indices for RHcrit.
            idpi,idpj
                            !  IDPI,IDPJ-indices for droplets

  INTEGER :: jblock   !Block size to be used in cache-blocked loop

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! needed for vatpoles for fv_cos_theta_latitude vs cos_theta_latitude
REAL, POINTER :: xx_cos_theta_latitude (:,:)

!-----------------------------------------------------------------------
!  Internal structure.
!  1. gather variables using index
!-----------------------------------------------------------------------
  IF (lhook) CALL dr_hook('LS_PPNC',zhook_in,zhook_handle)

IF ( l_vatpoles ) THEN
   xx_cos_theta_latitude => cos_theta_latitude
ELSE
   xx_cos_theta_latitude => fv_cos_theta_latitude
END IF ! vatpoles

!Set the block size.



  jblock = n


!Multiple RHC rows?

  IF ( (rhc_row_length * rhc_rows)  >   1) THEN
    multrhc = 1
  ELSE
    multrhc = 0
  END IF

  qcf2_c=0.0
  qrain_c=0.0
  qgraup_c=0.0
  aero_c=0.0

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(ii, ij,     &
!$OMP& irhi, irhj, i) 
  DO i = 1, n
    ii   =  ix(i,1)
    ij   =  ix(i,2)
    irhi = (multrhc * (ii - 1)) + 1
    irhj = (multrhc * (ij - 1)) + 1

    lsrain_c(i)          = lsrain(ii,ij)
    lssnow_c(i)          = lssnow(ii,ij)
    lssnow2_c(i)         = lssnow2(ii,ij)
    lsgraup_c(i)         = lsgraup(ii,ij)
    droplet_flux_c(i)    = droplet_flux(ii,ij)
    cttemp_c(i)          = cttemp(ii,ij)
    rainfrac_c(i)        = rainfrac(ii,ij)
    frac_ice_above_c(i)  = frac_ice_above(ii,ij)
    p(i)                 = p_theta_levels(ii,ij)
    rhodz(i)             = -p1upong*layer_thickness(ii,ij)
    deltaz_c(i)          = deltaz(ii,ij)
    uk_c(i)              = uk(ii,ij) 
    vk_c(i)              = vk(ii,ij) 
    ukp1_c(i)            = ukp1(ii,ij) 
    vkp1_c(i)            = vkp1(ii,ij)
    r_theta_levels_c(i)  = r_theta_levels(ii,ij,onelevel)
    g_cos_theta_latitude_c(i)  = xx_cos_theta_latitude(ii,ij)
    rhodz_dry_c(i)   = rhodz_dry(ii,ij)
    rhodz_moist_c(i) = rhodz_moist(ii,ij)
    bland_c(i)       = bland(ii,ij)
    cf_c(i)          = cf(ii,ij)
    cfl_c(i)         = cfl(ii,ij)
    cff_c(i)         = cff(ii,ij)
    qcf_c(i)         = qcf(ii,ij)
    qcl_c(i)         = qcl(ii,ij)
    q_c(i)           = q(ii,ij)
    t_c(i)           = t(ii,ij)
    n_drop_tpr_c(i)  = n_drop_tpr(ii, ij)

    IF (l_mcr_qcf2)   qcf2_c(i)   = qcf2(ii,ij)
    IF (l_mcr_qrain)  qrain_c(i)  = qrain(ii,ij)
    IF (l_mcr_qgraup) qgraup_c(i) = qgraup(ii,ij)
    IF (l_murk)       aero_c(i)   = aerosol(ii,ij)

    n_drop_out_c(i)  = n_drop_out(ii,ij)
    vfall_c(i)       = vfall(ii,ij)
    vfall2_c(i)      = vfall2(ii,ij)
    vfall_rain_c(i)  = vfall_rain(ii,ij)
    vfall_graup_c(i) = vfall_graup(ii,ij)
    rhc_c(i)         = rhcrit(irhi,irhj)

    land_fract_c(i) = land_fract(ii,ij)

        ! Process diagnostic arrays are initialised to zero in LSP_ICE
  END DO ! Loop over points
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
! ICE FORMATION/EVAPORATION/MELTING
! WATER CLOUD AND RAIN FORMATION/EVAPORATION
!-----------------------------------------------------------------------
! The call to LSP_ICE replaces the calls to LSP_EVAP, LSPFRMT
! and LSP_FORM.
! CFL_C contains cloud fraction for ice
! CFF_C contains cloud fraction for water

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(NONE)                      &
!$OMP& SHARED( n, p, rhodz,deltaz_c,rhodz_dry_c,rhodz_moist_c,          &
!$OMP& n_iterations, rhc_c, land_fract_c,                               &
!$OMP& qcf_c, qcl_c, q_c, qcf2_c, qrain_c, qgraup_c,                    &
!$OMP& n_drop_tpr_c, n_drop_out_c,                                      &
!$OMP& lsrain_c, vfall_rain_c, lssnow_c, vfall_c, lssnow2_c, vfall2_c,  &
!$OMP& lsgraup_c, vfall_graup_c, droplet_flux_c, frac_ice_above_c,      &
!$OMP& frac_agg_c,                                                      &
!$OMP& cttemp_c, rainfrac_c, t_c, cf_c, cfl_c, cff_c, bland_c,          &
!$OMP& psdep_c, psaut_c, psacw_c, psacr_c, psaci_c, psmlt_c,            &
!$OMP& psmltevp_c, praut_c, pracw_c, prevp_c, pgaut_c, pgacw_c, pgacs_c,&
!$OMP& pgmlt_c, pifrw_c, pifrr_c, piprm_c, pidep_c, piacw_c, piacr_c,   &
!$OMP& pimlt_c, pimltevp_c, pifall_c, psfall_c, prfall_c, pgfall_c,     &
!$OMP& plset_c, plevpset_c, lsiter, niter_bs, uk_c, vk_c, ukp1_c,       &
!$OMP& vkp1_c,  r_theta_levels_c, g_cos_theta_latitude_c, jblock )               &
!$OMP& PRIVATE(jj,i1,i2,ip)
  DO jj = 1, n, jblock
    !Cache-blocking applied to loop over points.
    !Set the starting index and the number of points to be passed through
    !to the called subroutine.
    i1 = jj
    i2 = MIN(jj+jblock-1,n)
    ip = i2-i1+1

    CALL lsp_ice( p(i1), rhodz(i1), deltaz_c(i1), rhodz_dry_c(i1),      &
      rhodz_moist_c(i1), n_iterations, ip, rhc_c(i1), land_fract_c(i1), &
      qcf_c(i1), qcl_c(i1), q_c(i1), qcf2_c(i1), qrain_c(i1),           &
      qgraup_c(i1),                                                     &
      n_drop_tpr_c(i1), n_drop_out_c(i1),                               &
      lsrain_c(i1),vfall_rain_c(i1), lssnow_c(i1),vfall_c(i1),          &
      lssnow2_c(i1),vfall2_c(i1), lsgraup_c(i1),                        &
      vfall_graup_c(i1), droplet_flux_c(i1),                            &
      frac_ice_above_c(i1), frac_agg_c(i1),                             &
      cttemp_c(i1), rainfrac_c(i1),                                     &
      t_c(i1),cf_c(i1),cfl_c(i1), cff_c(i1), bland_c(i1),               &
      psdep_c(i1),  psaut_c(i1),  psacw_c(i1), psacr_c(i1),             &
      psaci_c(i1),  psmlt_c(i1),  psmltevp_c(i1),                       &
      praut_c(i1),  pracw_c(i1),  prevp_c(i1),                          &
      pgaut_c(i1),  pgacw_c(i1),  pgacs_c(i1),  pgmlt_c(i1),            &
      pifrw_c(i1),  pifrr_c(i1),  piprm_c(i1),  pidep_c(i1),            &
      piacw_c(i1),  piacr_c(i1),  pimlt_c(i1),  pimltevp_c(i1),         &
      pifall_c(i1), psfall_c(i1), prfall_c(i1), pgfall_c(i1),           &
      plset_c(i1),  plevpset_c(i1), lsiter, niter_bs,                   &
      uk_c(i1), vk_c(i1), ukp1_c(i1), vkp1_c(i1),                       &
      r_theta_levels_c(i1), g_cos_theta_latitude_c(i1)                  &
      )

  END DO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  3.4 Lose aerosol by scavenging: call LSP_SCAV
!-----------------------------------------------------------------------

  IF (l_murk)  THEN
    CALL lsp_scav(n,lsrain_c,lssnow_c,droplet_flux_c,aero_c)
  END IF

!-----------------------------------------------------------------------
!  4  Scatter back arrays which will have been changed.

!-----------------------------------------------------------------------

! assumption that all pairs ii,ij generated for each i are distinct
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(ii, ij, i)
  DO i = 1, n

    ii         = ix(i,1)
    ij         = ix(i,2)
    t(ii,ij)   = t_c(i)
    q(ii,ij)   = q_c(i)
    qcf(ii,ij) = qcf_c(i)
    qcl(ii,ij) = qcl_c(i)

    IF (l_mcr_qcf2) qcf2(ii,ij)     = qcf2_c(i)
    IF (l_mcr_qrain) qrain(ii,ij)   = qrain_c(i)
    IF (l_mcr_qgraup) qgraup(ii,ij) = qgraup_c(i)
    IF (l_murk) aerosol(ii,ij)      = aero_c(i)

    n_drop_out(ii,ij)     = n_drop_out_c(i)
    lsrain(ii,ij)         = lsrain_c(i)
    lssnow(ii,ij)         = lssnow_c(i)
    lssnow2(ii,ij)        = lssnow2_c(i)
    lsgraup(ii,ij)        = lsgraup_c(i)
    droplet_flux(ii,ij)   = droplet_flux_c(i)
    cttemp(ii,ij)         = cttemp_c(i)
    rainfrac(ii,ij)       = rainfrac_c(i)
    frac_ice_above(ii,ij) = frac_ice_above_c(i)

    IF (l_pc2) THEN

      cff(ii,ij) = cff_c(i)
      cfl(ii,ij) = cfl_c(i)
      cf(ii,ij)  = cf_c(i)

    END IF  ! l_pc2

    vfall(ii,ij)       = vfall_c(i)
    vfall2(ii,ij)      = vfall2_c(i)
    vfall_rain(ii,ij)  = vfall_rain_c(i)
    vfall_graup(ii,ij) = vfall_graup_c(i)

        ! Only store process rates in array for diagnostic
        ! if a particular diagnostic is requested,
        ! otherwise overwriting will occur
        ! (space for the 3D array in MCR_CTL is only allocated
        ! if the diagnostic is active, to save memory)
    IF (l_aggfr_diag)    frac_agg(ii,ij, level) = frac_agg_c(i)
    IF (l_psdep_diag)    psdep(ii,ij, level)    = psdep_c(i)
    IF (l_psaut_diag)    psaut(ii,ij, level)    = psaut_c(i)
    IF (l_psacw_diag)    psacw(ii,ij, level)    = psacw_c(i)
    IF (l_psacr_diag)    psacr(ii,ij, level)    = psacr_c(i)
    IF (l_psaci_diag)    psaci(ii,ij, level)    = psaci_c(i)
    IF (l_psmlt_diag)    psmlt(ii,ij, level)    = psmlt_c(i)
    IF (l_psmltevp_diag) psmltevp(ii,ij, level) = psmltevp_c(i)
    IF (l_praut_diag)    praut(ii,ij, level)    = praut_c(i)
    IF (l_pracw_diag)    pracw(ii,ij, level)    = pracw_c(i)
    IF (l_prevp_diag)    prevp(ii,ij, level)    = prevp_c(i)
    IF (l_pgaut_diag)    pgaut(ii,ij, level)    = pgaut_c(i)
    IF (l_pgacw_diag)    pgacw(ii,ij, level)    = pgacw_c(i)
    IF (l_pgacs_diag)    pgacs(ii,ij, level)    = pgacs_c(i)
    IF (l_pgmlt_diag)    pgmlt(ii,ij, level)    = pgmlt_c(i)
    IF (l_pifrw_diag)    pifrw(ii,ij, level)    = pifrw_c(i)
    IF (l_piprm_diag)    piprm(ii,ij, level)    = piprm_c(i)
    IF (l_pidep_diag)    pidep(ii,ij, level)    = pidep_c(i)
    IF (l_piacw_diag)    piacw(ii,ij, level)    = piacw_c(i)
    IF (l_piacr_diag)    piacr(ii,ij, level)    = piacr_c(i)
    IF (l_pimlt_diag)    pimlt(ii,ij, level)    = pimlt_c(i)
    IF (l_pimltevp_diag) pimltevp(ii,ij, level) = pimltevp_c(i)
    IF (l_pifall_diag)   pifall(ii,ij, level)   = pifall_c(i)
    IF (l_psfall_diag)   psfall(ii,ij, level)   = psfall_c(i)
    IF (l_prfall_diag)   prfall(ii,ij, level)   = prfall_c(i)
    IF (l_pgfall_diag)   pgfall(ii,ij, level)   = pgfall_c(i)
    IF (l_plset_diag)    plset(ii,ij, level)    = plset_c(i)
    IF (l_plevpset_diag) plevpset(ii,ij, level) = plevpset_c(i)
    IF (l_pifrr_diag)    pifrr(ii, ij, level)   = pifrr_c(i)

  END DO ! Loop over points
!$OMP END PARALLEL DO


  IF (lhook) CALL dr_hook('LS_PPNC',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE ls_ppnc
END MODULE ls_ppnc_mod
