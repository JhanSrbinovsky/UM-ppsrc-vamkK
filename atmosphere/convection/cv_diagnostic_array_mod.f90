! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Holds convective diagnostic arrays

MODULE cv_diagnostic_array_mod

IMPLICIT NONE
SAVE

! Description:
!   Module containing convective diagnostics passed between atmos_physic2
!   ni_cont_ctl and diagnostics_conv
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3  programming standards.
!
! Declarations:


REAL, ALLOCATABLE  ::              &
  uw_dp(:,:,:)        & ! x-stress for deep convection (kg m-1 s-2)
 ,vw_dp(:,:,:)        & ! y-stress for deep convection (kg m-1 s-2)
 ,uw_shall(:,:,:)     & ! x-stress for shallow convection (kg m-1 s-2)
 ,vw_shall(:,:,:)     & ! y-stress for shallow convection (kg m-1 s-2)
 ,uw_mid(:,:,:)       & ! x-stress for mid convection (kg m-1 s-2)
 ,vw_mid(:,:,:)       & ! y-stress for mid convection (kg m-1 s-2)
 ,up_flux(:,:,:)      & !  updraught mass flux
 ,up_flux_half(:,:,:) & ! up flux on half levels
 ,dwn_flux(:,:,:)     & ! downdraught mass flux
 ,entrain_up(:,:,:)   & ! fractional entrainment rate
 ,detrain_up(:,:,:)   & ! fractional detrainment rate
 ,entrain_dwn(:,:,:)  & ! fractional entrainment rate into downdraughts
 ,detrain_dwn(:,:,:)  & ! fractional detrainment rate from downdraughts
 ,deep_tops(:,:,:)      ! Frequency deep convection terminates on level

REAL, ALLOCATABLE  ::       &
  precip_deep(:,:)          & ! deep precipitation
 ,precip_shall(:,:)         & ! shallow precipitation
 ,precip_mid(:,:)           & ! mid precipitation
 ,precip_cong(:,:)          & ! congestus precipitation
 ,cape_out(:,:)             & ! CAPE - dilute value from convection scheme
 ,deep_ind(:,:)             & ! Deep indicator
 ,shallow_ind(:,:)          & ! shallow indicator
 ,congestus_ind(:,:)        & ! congestus indicator
 ,congestus_ind2(:,:)       & ! congestus indicator2
 ,mid_ind(:,:)              & ! mid indicator
 ,ntml_diag(:,:)            & ! NTML output diag
 ,ntpar_diag(:,:)           & ! NTPAR output diag
 ,freeze_diag(:,:)          & ! freezing level output diag
 ,kterm_diag(:,:)           & ! deep termination level diag
 ,wstar_up_diag(:,:)        & ! cumulus layer convective vel
 ,wstar_dn_diag(:,:)        & ! subcloud layer convective vel
 ,mb1_diag(:,:)             & ! cloud base mass flux 1
 ,mb2_diag(:,:)             & ! cloud base mass flux 2
 ,wqt_cb(:,:)               & ! w'qt cloud base flux
 ,wthetal_cb(:,:)           & ! w'thetal' cb flux
 ,wqt_inv(:,:)              & ! w'qt inversion base flux
 ,wthetal_inv(:,:)          & ! w'thetal' inversion flux
 ,sh_top(:,:)               & ! height of top of shallow conv (m)
 ,sh_base(:,:)              & ! height of base of shallow conv (m)
 ,cg_top(:,:)               & ! height of top of congestus conv (m)
 ,cg_base(:,:)              & ! height of base of congestus conv (m)
 ,cg_term(:,:)              & ! congestus conv termination level
 ,cape_ts_diag(:,:)         & ! cape timescale diagnostic deep (s) 
 ,ind_cape_reduced_diag(:,:)& ! indicate reduced cape timescale
 ,deep_cfl_limited(:,:)     & ! Indicator for deep points CFL limited
 ,mid_cfl_limited(:,:)        ! Indicator for mid points CFL limited


! Convective aerosol relative diagnostics

REAL, ALLOCATABLE ::     &
  conscav_dust(:,:,:)    & !col total scvngd for DIV1,DIV2,...,DIV6 dust
 ,conscav_so4ait(:,:)    & !column total scvngd so4ait
 ,conscav_so4acc(:,:)    & !column total scvngd so4acc
 ,conscav_so4dis(:,:)    & !column total scvngd so4dis
 ,conscav_agedsoot(:,:)  & !colm total scvngd agd soot
 ,conscav_agedbmass(:,:) & !colm total scvngd bmass
 ,conscav_agedocff(:,:)  & !colm total scvngd ocff
 ,conscav_nitracc(:,:)   & ! column total scavenged acc nitrate
 ,conscav_nitrdiss(:,:)    ! column total scavenged diss nitrate
                                    
! Diagnostics for S Cycle
REAL, ALLOCATABLE ::     &
  conwash_so2(:,:)       & !column total scvngd so2
 ,conwash_nh3(:,:)         !column total scvngd nh3

REAL, ALLOCATABLE ::                &
  qcl_incr_inhom_diag(:,:,:)        & ! qcl         increment for STASH
 ,qcf_incr_inhom_diag(:,:,:)        & ! qcf         increment for STASH
 ,bulk_cf_incr_inhom_diag(:,:,:)    & ! bcf   increment for STASH
 ,cf_liquid_incr_inhom_diag(:,:,:)  & ! cf_l  increment for STASH
 ,cf_frozen_incr_inhom_diag(:,:,:)    ! cf_f  increment for STASH

REAL, ALLOCATABLE ::               &
  T_incr_conv_only(:,:,:)          & ! temperature inc for conv without PC2 homo
 ,q_incr_conv_only(:,:,:)            ! humidity inc for conv without PC2 homo

REAL, ALLOCATABLE ::               &
  theta_diag(:,:,:)                & ! Holds Theta at end of convection
 ,q_diag(:,:,:)                      ! Holds q at the end of convection

REAL, ALLOCATABLE ::               &
  T_incr_diag_conv(:,:,:)          & ! temperature increment for conv
 ,q_incr_diag_conv(:,:,:)          & ! humidity increment for conv
 ,qcl_incr_diag_conv(:,:,:)        & ! qCL   increment for conv
 ,qcf_incr_diag_conv(:,:,:)        & ! qCF   increment for conv
 ,cf_liquid_incr_diag_conv(:,:,:)  & ! cf_l  increment for conv
 ,cf_frozen_incr_diag_conv(:,:,:)  & ! cf_f  increment for conv
 ,bulk_cf_incr_diag_conv(:,:,:)    & ! bcf   increment for conv
 ,u_incr_diag_conv(:,:,:)          & ! u wind  increment for conv
 ,v_incr_diag_conv(:,:,:)            ! v wind  increment for conv

REAL, ALLOCATABLE ::      &
  mf_deep(:,:,:)          & ! mass flux for deep convection
 ,mf_congest(:,:,:)       & ! mass flux for congestus convection
 ,mf_shall(:,:,:)         & ! mass flux for shallow convection
 ,mf_midlev(:,:,:)        & ! mass flux for mid-level convection
 ,dt_deep(:,:,:)          & ! dT for deep convection
 ,dt_congest(:,:,:)       & ! dT for congestus convection
 ,dt_shall(:,:,:)         & ! dT for shallow convection
 ,dt_midlev(:,:,:)        & ! dT for mid-level convection
 ,dq_deep(:,:,:)          & ! dq for deep convection
 ,dq_congest(:,:,:)       & ! dq for congestus convection
 ,dq_shall(:,:,:)         & ! dq for shallow convection
 ,dq_midlev (:,:,:)       & ! dq for mid-level convection
 ,du_deep(:,:,:)          & ! du for deep convection
 ,du_congest(:,:,:)       & ! du for congestus convection
 ,du_shall(:,:,:)         & ! du for shallow convection
 ,du_midlev(:,:,:)        & ! du for mid-level convection
 ,dv_deep(:,:,:)          & ! dv for deep convection
 ,dv_congest(:,:,:)       & ! dv for congestus convection
 ,dv_shall(:,:,:)         & ! dv for shallow convection
 ,dv_midlev(:,:,:)        & ! dv for mid-level convection
 ,wqt_flux_sh(:,:,:)      & ! w'qt flux  shallow
 ,wql_flux_sh(:,:,:)      & ! w'qt flux shallow
 ,wthetal_flux_sh(:,:,:)  & ! w'thetal' flux shallow
 ,wthetav_flux_sh(:,:,:)  & ! w'thetav' flux shallow
 ,dubydt_pout(:,:,:)      & ! du on p grid
 ,dvbydt_pout(:,:,:)      & ! dv on p grid
 ,conv_rain_3d(:,:,:)     & ! 3d conv rainfall rate
 ,conv_snow_3d(:,:,:)       ! 3d conv snowfall rate


END MODULE cv_diagnostic_array_mod
