! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     SUBROUTINES LS_PPN and LS_PPNC------------------------------------
! Purpose:
!    LS_PPN and LS_PPNC:
!       Calculate large-scale (dynamical) precipitation.
!       LS_PPNC is the gather/scatter routine which then
!       calls LSP_ICE.
! Note: in all cases, level counters (incl subscripts) run from
!       1 (lowest model layer) to qdims%k_end
!       (topmost "wet" model layer)
!
! Programming standard: Unified Model Documentation Paper No 3
!
! Documentation: UM Documentation Paper 26.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation
MODULE ls_ppn_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE ls_ppn(                                                      &
  p_layer_boundaries, p_theta_levels, bland,                            &
!-----------------------
! primary fields and
! cloud fractions
!-----------------------
  cf, cfl, cff,                                                         &
  rhcrit,                                                               &
  lspice_dim1,lspice_dim2,lspice_dim3,                                  &
  rho_r2, q, qcf, qcl, t,                                               &
  qcf2, qrain, qgraup,                                                  &
!------------------------------------
! Wind field for lateral displacement 
! of falling ice by shear
!------------------------------------
  u_on_p, v_on_p,                                                       &
!-----------------------
! aerosol variables
!-----------------------
  sea_salt_film, sea_salt_jet,                                          &
  salt_dim1, salt_dim2, salt_dim3,                                      & 
  ukca_cdnc,                                                            &
  cdnc_dim1, cdnc_dim2, cdnc_dim3,                                      & 
  biogenic,                                                             &
  snow_depth, land_fract,                                               &
  so4_ait,                                                              &
  so4_acc,                                                              &
  so4_dis,                                                              &
  bmass_agd,                                                            &
  bmass_cld,                                                            &
  ocff_agd,                                                             &
  ocff_cld,                                                             &
  nitr_acc,                                                             &
  nitr_diss,                                                            &
  aerosol,                                                              &
  arcl,                                                                 &
!---------------------------
! Other variables for mphys
!---------------------------
  lsrain,lssnow,                                                        &
  lsrain3d, lssnow3d, rainfrac3d,                                       &
  n_drop_tpr, n_drop_3d,                                                &
  rhc_row_length, rhc_rows,                                             &
!------------------------------------------------------
! Variables for stochastic physics random parameters2
!------------------------------------------------------
  m_ci,                                                                 &
!-------------
! Error code
!-------------
  error          )

     ! General modules
  USE um_input_control_mod,  ONLY: l_mr_physics1

  USE arcl_mod,              ONLY: npd_arcl_compnts

     ! Microphysics modules
  USE mphys_inputs_mod,      ONLY: l_it_melting, l_mcr_iter, lsiter,    &
                                   niter_bs, l_warm_new, l_mcr_qrain,   &
                                   l_mcr_qcf2, l_mcr_qgraup
                                  
  USE mphys_bypass_mod,      ONLY: mphys_mod_top     

  USE mphys_constants_mod,   ONLY: m_ci_sav, l_calc_mp_const, iter_z,   &
                                   max_it_mlt

  USE mphys_diags_mod,       ONLY: l_point_diag, mphys_pts

     ! ENDGame modules
  USE level_heights_mod,     ONLY: r_theta_levels,  r_rho_levels,       &
                                   eta_theta_levels

  USE atm_fields_bounds_mod, ONLY: qdims, pdims_s, tdims,               &
                                   pdims

     ! Dr Hook modules
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim

     ! Large scale precipitation modules
  USE ls_ppnc_mod,           ONLY: ls_ppnc
  USE lsp_taper_ndrop_mod,   ONLY: lsp_taper_ndrop
  USE lspcon_mod,            ONLY: lspcon

  IMPLICIT NONE

  INTEGER, INTENT(IN) ::                                                &
    rhc_row_length, rhc_rows,                                           &
    lspice_dim1,lspice_dim2,lspice_dim3,                                &
! Dimensions for 3D arrays
    salt_dim1,                                                          &
                     ! Array dimensions for sea-salt arrays (equal
    salt_dim2,                                                          &
                     ! either to qdims%i_start:qdims%i_end,
                     ! qdims%j_start:qdims%j_end, and
                     ! 1:qdims%k_end, or
    salt_dim3,                                                          &
                     ! else 1,1,1, depending on L_SEASALT_CCN).
    cdnc_dim1, cdnc_dim2, cdnc_dim3
                     ! UKCA cloud drop number concentration dimensions

  REAL ::                                                               &
    cf( qdims%i_start:qdims%i_end,                                      &
        qdims%j_start:qdims%j_end,                                      &
        1:qdims%k_end ),                                                &
                                    ! IN Cloud fraction.
    p_theta_levels( qdims%i_start:qdims%i_end,                          &
                    qdims%j_start:qdims%j_end,                          &
                    1:qdims%k_end ),                                    &

! Note that the declaration beneath has been written to cope with the
! ENDGame modules, but as ENDGame hasn't yet got a qdims value for 0
! in the k dimension, this may need altering in the future

    p_layer_boundaries( tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end,                      &
                         0           :tdims%k_end ),                    &

    rhcrit( rhc_row_length, rhc_rows,                                   &
            1:qdims%k_end ),                                            &
                                                 ! IN Critical humidity
                                                 ! for cloud formation.
   cfl( qdims%i_start:qdims%i_end,                                      &
        qdims%j_start:qdims%j_end,                                      &
        1:qdims%k_end ),                                                &
                                        !IN Cloud liquid fraction.
   cff( qdims%i_start:qdims%i_end,                                      &
        qdims%j_start:qdims%j_end,                                      &
        1:qdims%k_end ),                                                &
                                        !IN Cloud ice fraction.
   rho_r2(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end,                                &
          pdims_s%k_start:pdims_s%k_end)
                                ! IN Air density * earth radius**2

  LOGICAL :: bland( tdims%i_start : tdims%i_end,                        &
                 tdims%j_start : tdims%j_end )                          
                                     ! IN Land/sea mask

  REAL, INTENT(INOUT) ::                                                &
   q( qdims%i_start:qdims%i_end,                                        &
      qdims%j_start:qdims%j_end,                                        &
      1:qdims%k_end ),                                                  &
                                         ! Specific humidity (kg water
   qcf( qdims%i_start:qdims%i_end,                                      &
        qdims%j_start:qdims%j_end,                                      &
        1:qdims%k_end ),                                                &
                                         ! Cloud ice (kg per kg air).
   qcl( qdims%i_start:qdims%i_end,                                      &
        qdims%j_start:qdims%j_end,                                      &
        1:qdims%k_end ),                                                &
                                         ! Cloud liquid water (kg per
   qcf2( qdims%i_start:qdims%i_end,                                     &
         qdims%j_start:qdims%j_end,                                     &
         1:qdims%k_end ),                                               &
                                         ! Ice (kg per kg air)
   qrain( qdims%i_start:qdims%i_end,                                    &
          qdims%j_start:qdims%j_end,                                    &
          1:qdims%k_end ),                                              &
                                         ! Rain (kg per kg air)
   qgraup( qdims%i_start:qdims%i_end,                                   &
           qdims%j_start:qdims%j_end,                                   &
           1:qdims%k_end ),                                             &
                                         ! Graupel (kg per kg air)
   t( qdims%i_start:qdims%i_end,                                        &
      qdims%j_start:qdims%j_end,                                        &
      1:qdims%k_end ),                                                  &
                                         ! Temperature (K).
   aerosol( qdims%i_start:qdims%i_end,                                  &
            qdims%j_start:qdims%j_end,                                  &
            1:qdims%k_end )
                                         ! 'Murk' tracer aerosol.

  REAL,INTENT(IN) ::                                                    &
    ! For calculating shear in falling ice cloud fraction calculation.
    u_on_p( pdims%i_start : pdims%i_end,                                & 
            pdims%j_start : pdims%j_end,                                & 
            pdims%k_start : pdims%k_end),                               & 
    v_on_p( pdims%i_start : pdims%i_end,                                & 
            pdims%j_start : pdims%j_end,                                & 
            pdims%k_start : pdims%k_end),                               &

    arcl(   qdims%i_start : qdims%i_end,                                &
            qdims%j_start : qdims%j_end,                                &
                        1 : qdims%k_end,                                &
            npd_arcl_compnts            )

  REAL,INTENT(INOUT) ::                                                 &
                                     !Sulphur Cycle tracers (mmr kg/kg)
     so4_ait( qdims%i_start:qdims%i_end,                                &
              qdims%j_start:qdims%j_end,                                &
              1:qdims%k_end ),                                          &
     so4_acc( qdims%i_start:qdims%i_end,                                &
              qdims%j_start:qdims%j_end,                                &
              1:qdims%k_end ),                                          &
     so4_dis( qdims%i_start:qdims%i_end,                                &
              qdims%j_start:qdims%j_end,                                &
              1:qdims%k_end ),                                          &

                                     !Biomass smoke tracers
     bmass_agd( qdims%i_start:qdims%i_end,                              &
                qdims%j_start:qdims%j_end,                              &
                1:qdims%k_end ),                                        &
     bmass_cld( qdims%i_start:qdims%i_end,                              &
                qdims%j_start:qdims%j_end,                              &
                1:qdims%k_end ),                                        &

                                     !Fossil-fuel organic carbon tracers
     ocff_agd( qdims%i_start:qdims%i_end,                               &
               qdims%j_start:qdims%j_end,                               &
               1:qdims%k_end ),                                         &
     ocff_cld( qdims%i_start:qdims%i_end,                               &
               qdims%j_start:qdims%j_end,                               &
               1:qdims%k_end ),                                         &

                                     !Ammonium nitrate tracers
     nitr_acc( qdims%i_start:qdims%i_end,                               &
                qdims%j_start:qdims%j_end,                              &
                1:qdims%k_end ),                                        &
      nitr_diss( qdims%i_start:qdims%i_end,                             &
                 qdims%j_start:qdims%j_end,                             &
                 1:qdims%k_end )

  REAL ::                                                               &
      snow_depth( qdims%i_start:qdims%i_end,                            &
                  qdims%j_start:qdims%j_end ),                          &
                                     ! IN Snow depth (m)
      land_fract( qdims%i_start:qdims%i_end,                            &
                  qdims%j_start:qdims%j_end )
                                     ! IN Land fraction

  REAL ::                                                               &
      sea_salt_film(salt_dim1, salt_dim2, salt_dim3),                   &
                                                         ! (m-3)
      sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)  ! (m-3)

  REAL ::                                                               &
      biogenic( qdims%i_start:qdims%i_end,                              &
                qdims%j_start:qdims%j_end,                              &
                1:qdims%k_end )
                                                  ! (m.m.r.)

  REAL :: ukca_cdnc(cdnc_dim1, cdnc_dim2, cdnc_dim3) 
!         CDNC from UKCA for 2nd indirect effect (m-3)

  REAL ::                                                               &
   lsrain( qdims%i_start:qdims%i_end,                                   &
           qdims%j_start:qdims%j_end ),                                 &
                               ! OUT Surface rainfall rate (kg / sq m /
   lssnow( qdims%i_start:qdims%i_end,                                   &
           qdims%j_start:qdims%j_end ),                                 &
                               ! OUT Surface snowfall rate (kg / sq m /
   lsgraup( qdims%i_start:qdims%i_end,                                  &
            qdims%j_start:qdims%j_end )
                               ! Graupel fall rate (kg/m2/s)
 
  REAL, INTENT(OUT) :: n_drop_tpr( qdims%i_start : qdims%i_end,         &
                                   qdims%j_start : qdims%j_end,         &
                                               1 : qdims%k_end ),       &
                       n_drop_3d(  qdims%i_start : qdims%i_end,         &
                                   qdims%j_start : qdims%j_end,         &
                                               1 : qdims%k_end )
                       ! Tapered droplet number and droplet number
                       ! from autoconversion scheme

  REAL ::                                                               &
      lsrain3d(lspice_dim1,lspice_dim2,lspice_dim3),                    &
                                                        ! OUT
!                           Rain rate out of each model layer
      lssnow3d(lspice_dim1,lspice_dim2,lspice_dim3),                    &
                                                        ! OUT
!                           Snow rate out of each model layer
      rainfrac3d(lspice_dim1,lspice_dim2,lspice_dim3) ! OUT
!                           rain fraction out of each model layer
! Variables for stochastic physics random parameters
  REAL,    INTENT(IN) :: m_ci   ! used to modify ice fall speed.

  INTEGER ::    error          ! OUT Return code - 0 if OK,
!                                                - 1 if bad arguments.

!    Workspace usage ---------------------------------------------------

  REAL ::  lsrain_mean(qdims%i_start:qdims%i_end,                       &
                       qdims%j_start:qdims%j_end )

  REAL ::  lssnow_mean(qdims%i_start:qdims%i_end,                       &
                       qdims%j_start:qdims%j_end )

  INTEGER ::                                                            &
   ix( (qdims%i_end - qdims%i_start + 1) *                              &
       (qdims%j_end - qdims%j_start + 1) , 2 ),                         &
                                 ! Index for compress/expand.
    n_iterations,                                                       &
                                 ! Number of iterations
    eta_ratio1               ! ratio of model level height to
                                 ! iterative melting height plus 1.

  REAL :: vfall( qdims%i_start:qdims%i_end,                             &
              qdims%j_start:qdims%j_end )
               ! snow fall velocity (m per s).
  REAL :: vfall2( qdims%i_start:qdims%i_end,                            &
               qdims%j_start:qdims%j_end )
               ! fall velocity for qcf2 (m/s)
  REAL :: lssnow2( qdims%i_start:qdims%i_end,                           &
                qdims%j_start:qdims%j_end )
               ! snowfall rate for qcf2
  REAL :: droplet_flux( qdims%i_start:qdims%i_end,                      &
                     qdims%j_start:qdims%j_end )
               ! water drop flux / kg m-2 s-1
  REAL :: vfall_rain( qdims%i_start:qdims%i_end,                        &
                   qdims%j_start:qdims%j_end )
               ! fall velocity for rain (m/s)
  REAL :: vfall_graup( qdims%i_start:qdims%i_end,                       &
                    qdims%j_start:qdims%j_end )
               ! fall vel. for graupel (m/s)
  REAL :: cttemp( qdims%i_start:qdims%i_end,                            &
               qdims%j_start:qdims%j_end )
  REAL :: rainfrac( qdims%i_start:qdims%i_end,                          &
                 qdims%j_start:qdims%j_end )
  REAL :: frac_ice_above( qdims%i_start:qdims%i_end,                    &
                       qdims%j_start:qdims%j_end )
               ! Cloud ice fraction passed
               ! in layer above
  REAL :: layer_thickness( qdims%i_start:qdims%i_end,                   &
                        qdims%j_start:qdims%j_end )
  REAL :: rho1( qdims%i_start:qdims%i_end,                              &
             qdims%j_start:qdims%j_end )
  REAL :: rho2( qdims%i_start:qdims%i_end,                              &
             qdims%j_start:qdims%j_end )
  REAL :: deltaz( qdims%i_start:qdims%i_end,                            &
               qdims%j_start:qdims%j_end,                               &
               1:qdims%k_end )
  REAL :: rhodz_dry( qdims%i_start:qdims%i_end,                         &
                  qdims%j_start:qdims%j_end,                            &
                  1:qdims%k_end )
  REAL :: rhodz_moist( qdims%i_start:qdims%i_end,                       &
                    qdims%j_start:qdims%j_end,                          &
                    1:qdims%k_end )
  REAL :: q_total( qdims%i_start:qdims%i_end,                           &
                qdims%j_start:qdims%j_end )

  REAL :: iter_eta                      ! eta value at which to
                                         ! start iterative melting

  REAL :: fqi_rain( qdims%i_start:qdims%i_end,                          &
                 qdims%j_start:qdims%j_end,                             &
                 1:qdims%k_end )
  REAL :: lamr( qdims%i_start:qdims%i_end,                              &
             qdims%j_start:qdims%j_end,                                 &
             1:qdims%k_end )


!  Physical constants -------------------------------------------------
  REAL :: cfmin
  PARAMETER (                                                           &
   cfmin=1.0e-3                                                         &
                           ! Used for LS_PPNC  compress.
  )
!  Define local variables ----------------------------------------------
  INTEGER :: i,k,j,it,                                                  &
                        ! Loop counters: I - horizontal field index;
!                                        K - vertical level index.
             kp1,                                                       & 
                        ! Index of level above: k=k+1, apart from  
                        ! when k=_dims%end when kp1=k.
             n          ! "nval" for WHEN routine.

  REAL :: work
                      ! work variable

  LOGICAL :: l_3ddiag ! Flag to determine if we want 3d diagnostics 
  
  REAL :: one_over_niter_bs  ! 1./niter_bs

!  Variables for Dr Hook:
 
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
  IF (lhook) CALL dr_hook('LS_PPN',zhook_in,zhook_handle)

  error=0

! Determine if we want 3d diagnosics
  l_3ddiag = (lspice_dim1  ==  (qdims%i_end - qdims%i_start + 1)  .AND. & 
     lspice_dim2  ==  (qdims%j_end - qdims%j_start + 1)  .AND.          & 
     lspice_dim3  ==  ( qdims%k_end ) )

!----------------------------------------------------------------------
! If multiple iterations is not selected, then ensure the number
! of all iterations is always exactly 1.
!----------------------------------------------------------------------

  IF ( .NOT. l_mcr_iter ) THEN

    lsiter   = 1
    niter_bs = 1

  END IF

!----------------------------------------------------------------------

! Define CX and CONSTP values

  IF (l_calc_mp_const .OR. m_ci /= m_ci_sav) THEN  

    ! If the microphysics constants are not set we can compute them
    ! and save to a module. However, if m_ci (random parameters) 
    ! changes some constants will change with it, so we
    ! need to recompute the constants in this case. Continuation
    ! runs (e.g. climate) will also need to recalculate the constants.


    CALL lspcon( m_ci )

  END IF

  

! ----------------------------------------------------------------------
! Calculate the (non-hydrostatic) layer thicknesses (deltaz) and air
! densities multiplied by deltaz (rhodz_moist and rhodz_dry).
! ----------------------------------------------------------------------
! We should note that this formulation, although better than the
! hydrostatic formulation, is still not entirely conservative. To ensure
! conservation we would need to rewrite the large-scale precipitation
! scheme to consider masses in terms of rho<q>, and
! not the current <rho>q formulation.

! We only need to calculate averages for the moist levels

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k,      &
!$OMP& q_total,rho1,rho2)
  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end

              ! Calculate densities at the boundaries of the layer
              ! by removing the r**2 term from rho_r2.
              ! Rho1 is the density at the lower boundary.
        rho1(i,j)= rho_r2(i,j,k)/( r_rho_levels(i,j,k) *                &
                               r_rho_levels(i,j,k) )

              ! Check whether there is a rho level above the current
              ! moist level.
        IF ( k  <   tdims%k_end ) THEN
                ! Rho2 is the density at the upper boundary.
          rho2(i,j)= rho_r2(i,j,k+1)/( r_rho_levels(i,j,k+1) *          &
                                 r_rho_levels(i,j,k+1) )

                ! Calculate the average value of rho across the layer
                ! multiplied by the layer thickness and the layer
                ! thickness.
          rhodz_moist(i,j,k) =                                          &
                        rho2(i,j) * ( r_theta_levels(i,j,k) -           &
                                           r_rho_levels(i,j,k) )        &
                     +  rho1(i,j) * ( r_rho_levels(i,j,k+1) -           &
                                      r_theta_levels(i,j,k) )
          deltaz(i,j,k) = r_rho_levels(i,j,k+1)                         &
                            - r_rho_levels(i,j,k)

          IF (k  ==  1) THEN  
            ! For the lowest layer we need to extend the lower
            ! boundary from the first rho level to the surface.
            ! The surface is the 0'th theta level.
            deltaz(i,j,1) = r_rho_levels(i,j,2)                         &
                              - r_theta_levels(i,j,0)
            rhodz_moist(i,j,1) = rhodz_moist(i,j,1)*deltaz(i,j,1)       &
                     / (r_rho_levels(i,j,2)-r_rho_levels(i,j,1))
          END IF  ! k  ==  1

        ELSE
          ! For a top layer higher than the highest rho level
          ! we can calculate a pseudo rho level. We will assume
          ! it has a similar density to the rho level below
          ! and that the intervening theta level is in the centre
          ! of the layer.
          deltaz(i,j,k) = 2.0*(r_theta_levels(i,j,k)                    &
                                -r_rho_levels(i,j,k))
          rhodz_moist(i,j,k) = rho1(i,j) * deltaz(i,j,k)

        END IF  ! k  <   tdims%k_end

        ! Calculate total moisture
        q_total(i,j) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)

        IF (l_mcr_qcf2) THEN
          q_total(i,j) = q_total(i,j) + qcf2(i,j,k)
        END IF  ! l_mcr_qcf2
          
        IF (l_mcr_qrain) THEN
          q_total(i,j) = q_total(i,j) + qrain(i,j,k)
        END IF  ! l_mcr_qrain

        IF (l_mcr_qgraup) THEN
          q_total(i,j) = q_total(i,j) + qgraup(i,j,k)
        END IF  ! l_mcr_qgraup


        ! Rho_r2 uses the moist density of air. If the mixing
        ! ratio framework is in place then we need to also know
        ! the dry density of air.
        IF (l_mr_physics1) THEN
          rhodz_dry(i,j,k) = rhodz_moist(i,j,k)                         &
               / (1.0 + q_total(i,j))
        ELSE
          rhodz_dry(i,j,k) = rhodz_moist(i,j,k)                         &
               * (1.0 - q_total(i,j)) 
        END IF  ! l_mr_physics1

      END DO  ! i
    END DO  ! j
  END DO  ! k
!$OMP END PARALLEL DO

!-------------------------------------------------------------
! Calculation of cloud droplet number. This is now calculated
! here for all models and not in lsp_autoc as was done in
! older versions of the Unified Model
!-------------------------------------------------------------

    CALL lsp_taper_ndrop(                                               &
                   ! (Full) Aerosol tracers
                            so4_ait, so4_acc, so4_dis, sea_salt_film,   &
                            biogenic, sea_salt_jet, bmass_agd,          &
                            bmass_cld, ocff_agd, ocff_cld,              &
                            nitr_acc, nitr_diss, arcl,                  &
                   ! Murk aerosol
                            aerosol,                                    &
                   ! CDNC from UKCA 
                            ukca_cdnc,                                  &
                            cdnc_dim1, cdnc_dim2, cdnc_dim3,            & 
                   ! Other parameters
                            rhodz_dry,  rhodz_moist,                    &
                            deltaz,                                     &
                            snow_depth, land_fract,                     &
                   ! Output parameters
                            n_drop_tpr                                  &
                                 )

!-----------------------------------------------------------------------
!  2. Loop round levels from top down (counting bottom level as level 1,
!     as is standard in the Unified model).
!-----------------------------------------------------------------------


!-----------------------------------------------
! Setup for iterative metlting
! Define constants outside of K loop
!-----------------------------------------------

      ! calculate level independent eta value
      ! iter_z is theta_height for level 13 in 38 level set.
      ! hence is now independent of number of levels

  iter_eta = iter_z / mphys_mod_top

!----------------------------------------------------------------------- 
! Internal structure.  
! 2a. Initialise outside of iterative loop. 
!----------------------------------------------------------------------- 
  lsrain_mean=0.0 
  lssnow_mean=0.0 
  lsrain3d=0.0 
  lssnow3d=0.0 
  rainfrac3d=0.0 

  one_over_niter_bs=1.0/niter_bs

  DO it = 1, niter_bs ! Substep outside of column 
!-----------------------------------------------------------------------
!  Internal structure - moved inside BS iter loop
!  2.b Initialise rain and snow to zero.
!   Initialise scavenged amounts of S Cycle tracers to 0 for full field
!-----------------------------------------------------------------------
  DO j = qdims%j_start, qdims%j_end
    DO i = qdims%i_start, qdims%i_end

      lsrain(i,j)=0.0
      lssnow(i,j)=0.0
      lssnow2(i,j)=0.0
      lsgraup(i,j)=0.0
      droplet_flux(i,j)=0.0
      cttemp(i,j)=0.0
      rainfrac(i,j)=0.0
      frac_ice_above(i,j)=0.0
      vfall(i,j)=0.0
      vfall2(i,j)=0.0
      vfall_rain(i,j)=0.0
      vfall_graup(i,j)=0.0

    END DO ! Loop over points,i
  END DO ! Loop over points,j

! Initialise n_drop_3d to zero before passing down code tree

  n_drop_3d(:,:,:) = 0.0

  DO k = qdims%k_end, 1, -1

    ! kp1 is index of model level above, unless we are the model  
    ! top in which case it is set to k. Used to calc vertical wind shear. 
    IF (k == qdims%k_end) THEN 
      kp1=k 
    ELSE 
      kp1=k+1 
    END IF

    IF (l_it_melting) THEN

           !-----------------------------------------------------------
           !      Calculate the number of fall and melting iterations
           !-----------------------------------------------------------
              !Version to remove dependency on bl_levels

      IF (eta_theta_levels(k) /= 0.0) THEN

                 ! We're not right at the surface, so can
                 ! calculate the number of iterations based
                 ! on eta value.

        eta_ratio1 =                                                    &
              INT( (iter_eta / eta_theta_levels(k) ) )+ 1

                    ! was eta_theta_levels_tmp(k+1)

        n_iterations = MIN(eta_ratio1, max_it_mlt)

      ELSE

                 ! eta_theta_levels = 0.0
                 ! need to avoid divide by 0.0 error
                 ! at this level, however, iterations should always
                 ! be the maximum number, so we can just set this.

        n_iterations = max_it_mlt

      END IF ! eta_theta_levels > 0

    ELSE     ! Do not use iterative melting

      n_iterations = 1

    END IF  ! l_it_melting

    n=0
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end

!-----------------------------------------------------------------------
!  2.5 Form INDEX IX to gather/scatter variables in LS_PPNC
!-----------------------------------------------------------------------

!  Set index where cloud fraction > CFMIN or where non-zero pptn
!  Note: whenimd is functionally equivalent to WHENILE (but autotasks).


        layer_thickness(i,j) =                                          &
         p_layer_boundaries(i,j,k) - p_layer_boundaries(i,j,k-1)

            ! Set up IF statement to determine whether to call the
            ! microphysics code for this grid box (i.e. if there is
            ! already condensate in the grid box or there is
            ! precipitation about to fall into the grid box)
        work = qcf(i,j,k)

            ! Include extra microphysics variables if in use
        IF (l_mcr_qcf2 )  work = work + qcf2(i,j,k) + lssnow2(i,j)
        IF (l_mcr_qrain)  work = work + qrain(i,j,k)
        IF (l_mcr_qgraup) work = work + qgraup(i,j,k)+lsgraup(i,j)
        work = work + droplet_flux(i,j)   ! droplet settling

        IF (cfl(i,j,k) > cfmin .OR.                                     &
           (lsrain(i,j)+lssnow(i,j)) > 0.0 .OR. work > 0.0) THEN
              ! include this grid box.
              ! Strictly speaking the CFL > CFMIN clause is too
              ! restrictive since ice nucleation does not require
              ! liquid water, but the code would be very messy.
          n = n + 1
          ix(n,1) = i
          ix(n,2) = j
              ! Note that mphys is done on this point if diagnostic is
              ! requested
          IF (l_point_diag) mphys_pts(i,j,k) = .TRUE.

        END IF

! set up rain fraction
        IF (l_mcr_qrain .AND. l_warm_new) THEN
          IF (rainfrac(i,j) == 0.0 .AND. qrain(i,j,k) > 0.0) THEN
! if cloud is present, assume this is a good proxy for the rain frac
            rainfrac(i,j) = MAXVAL(cf(i,j,k:qdims%k_end))
          END IF
          IF (rainfrac(i,j) == 0.0 .AND. qrain(i,j,k) > 0.0) THEN
! otherwise set to 0.5
            rainfrac(i,j) = 0.5
! set a lower limit for stability 
          ELSE IF (rainfrac(i,j) <= 0.01 .AND. qrain(i,j,k) > 0.0) THEN 
            rainfrac(i,j) = 0.01 
          END IF
        END IF

      END DO ! Loop over points,i 
    END DO ! Loop over points,j 

    IF (n > 0) THEN

      CALL ls_ppnc(k,ix,n,n_iterations,                                 &
                   lsrain,lssnow,lssnow2,lsgraup,droplet_flux,          &
                   cf(1,1,k),cfl(1,1,k),cff(1,1,k),                     &
                   qcf(1,1,k),qcl(1,1,k),t(1,1,k),                      &
                   qcf2(1,1,k),qrain(1,1,k),qgraup(1,1,k),              &
                   n_drop_tpr(1,1,k), n_drop_3d(1,1,k),                 &
                   aerosol(1,1,k), land_fract,                          &
                   q(1,1,k), p_theta_levels(1,1,k), layer_thickness,    &
                   deltaz(1,1,k), rhodz_dry(1,1,k), rhodz_moist(1,1,k), &
                   rhc_row_length, rhc_rows,                            &
                   bland, rhcrit(1,1,k),                                &
                   vfall, vfall2, vfall_rain, vfall_graup,              &
                   frac_ice_above,                                      & 
                   cttemp, rainfrac, lsiter, niter_bs,                  &
                   u_on_p(1,1,k),   v_on_p(1,1,k),                      &
                   u_on_p(1,1,kp1), v_on_p(1,1,kp1), k                  &
                    )
    END IF

! Copy rainfall and snowfall rates to 3D fields for diagnostic output

    IF ( l_3ddiag ) THEN

! Only copy rain and snow to 3D fields if arrays are dimensionalized.

      DO j = qdims%j_start, qdims%j_end

        DO i = qdims%i_start, qdims%i_end

          lsrain3d(i, j, k)   = lsrain3d(i, j, k) +                     &  
             (lsrain(i, j) + droplet_flux(i, j))*one_over_niter_bs

          lssnow3d(i, j, k)   = lssnow3d(i, j, k) +                     &
             (lssnow(i, j) + lssnow2(i, j) + lsgraup(i, j))             &
             *one_over_niter_bs

          rainfrac3d(i, j, k) = rainfrac3d(i, j, k) +                   &
             rainfrac(i, j)*one_over_niter_bs

        END DO

      END DO

    END IF

  END DO ! Loop over K

! If substepping outside of loop over K, then need to accumulate (mean)
! precip rates...
  DO j = qdims%j_start, qdims%j_end

    DO i = qdims%i_start, qdims%i_end

      lsrain_mean(i,j) = lsrain_mean(i,j)                             &
         + (lsrain(i, j) + droplet_flux(i, j))*one_over_niter_bs
      
! Add together ice crystals, snow aggregates and graupel
! for surface snow rate (kg/m2/s)
      IF (l_mcr_qcf2) THEN

        lssnow_mean(i, j) = lssnow_mean(i, j)                           &
                          + (lssnow(i, j)  + lssnow2(i, j)              &
                          + lsgraup(i, j)) * one_over_niter_bs
      ELSE

        lssnow_mean(i, j) = lssnow_mean(i, j)                           &
                          + ( lssnow(i, j) ) * one_over_niter_bs

      END IF

    END DO

  END DO


  END DO ! Outer substepping loop over it

! Recover meaned precip rates
    DO j = qdims%j_start, qdims%j_end

      DO i = qdims%i_start, qdims%i_end

        lssnow(i,j) = lssnow_mean(i,j)
        lsrain(i,j) = lsrain_mean(i,j) 

      END DO

    END DO

  IF (lhook) CALL dr_hook('LS_PPN',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE ls_ppn

END MODULE ls_ppn_mod
