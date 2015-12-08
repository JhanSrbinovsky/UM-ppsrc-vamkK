! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Initialise all the SW Diagnostics to False.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!    FORTRAN 77/90  with extensions listed in documentation.

!-----------------------------------------------------------------------

SUBROUTINE init_swdiag_logic(j_sw)

      USE sw_diag_mod
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

! arguments with intent in

      INTEGER, INTENT(IN) :: j_sw      ! call to SW radiation

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('INIT_SWDIAG_LOGIC',zhook_in,zhook_handle)


! Switches enabled as STASHflags: SW

! Fluxes and Heating Rates

        sw_diag(j_sw)%l_flux_up           = .FALSE.
        sw_diag(j_sw)%l_flux_down         = .FALSE.
        sw_diag(j_sw)%l_flux_up_clear     = .FALSE.
        sw_diag(j_sw)%l_flux_down_clear   = .FALSE.
        sw_diag(j_sw)%l_solar_out_toa     = .FALSE.
        sw_diag(j_sw)%l_solar_out_clear   = .FALSE.
        sw_diag(j_sw)%l_surface_down_flux = .FALSE.
        sw_diag(j_sw)%l_surf_down_clr     = .FALSE.
        sw_diag(j_sw)%l_surf_up_clr       = .FALSE.
        sw_diag(j_sw)%l_clear_hr          = .FALSE.
        sw_diag(j_sw)%l_net_flux_trop     = .FALSE.
        sw_diag(j_sw)%l_up_flux_trop      = .FALSE.
        sw_diag(j_sw)%l_flux_direct       = .FALSE.
        sw_diag(j_sw)%l_flux_diffuse      = .FALSE.

! UV-Fluxes

        sw_diag(j_sw)%l_uvflux_direct     = .FALSE.
        sw_diag(j_sw)%l_uvflux_up         = .FALSE.
        sw_diag(j_sw)%l_uvflux_net        = .FALSE.
        sw_diag(j_sw)%l_surf_uv           = .FALSE.
        sw_diag(j_sw)%l_surf_uv_clr       = .FALSE.

! Surface Albedo Diagnostics 

        sw_diag(j_sw)%l_direct_albedo     = .FALSE. 
        sw_diag(j_sw)%l_diffuse_albedo    = .FALSE. 
        sw_diag(j_sw)%l_vis_albedo_sc     = .FALSE. 
        sw_diag(j_sw)%l_nir_albedo_sc     = .FALSE. 

! Radiance

        sw_diag(j_sw)%l_toa_radiance      = .FALSE.

! Microphysical diagnostics

        sw_diag(j_sw)%re_strat_flag       = .FALSE.
        sw_diag(j_sw)%wgt_strat_flag      = .FALSE.
        sw_diag(j_sw)%lwp_strat_flag      = .FALSE.
        sw_diag(j_sw)%re_conv_flag        = .FALSE.
        sw_diag(j_sw)%wgt_conv_flag       = .FALSE.
        sw_diag(j_sw)%ntot_diag_flag      = .FALSE.
        sw_diag(j_sw)%strat_lwc_diag_flag = .FALSE.
        sw_diag(j_sw)%so4_ccn_diag_flag   = .FALSE.
        sw_diag(j_sw)%cond_samp_wgt_flag  = .FALSE.
        sw_diag(j_sw)%weighted_re_flag    = .FALSE.
        sw_diag(j_sw)%sum_weight_re_flag  = .FALSE.
        sw_diag(j_sw)%wgtd_warm_re_flag   = .FALSE.
        sw_diag(j_sw)%sum_wgt_warm_re_flag= .FALSE.
        sw_diag(j_sw)%seasalt_film_flag   = .FALSE.
        sw_diag(j_sw)%seasalt_jet_flag    = .FALSE.
        sw_diag(j_sw)%nc_diag_flag        = .FALSE.
        sw_diag(j_sw)%nc_weight_flag      = .FALSE.

! Diagnostics for MOSES II

        sw_diag(j_sw)%l_FlxSolBelow690nmSurf = .FALSE.
        sw_diag(j_sw)%l_FlxSeaBelow690nmSurf = .FALSE.

! Diagnostics for orography correction

        sw_diag(j_sw)%l_sol_bearing = .FALSE.
        sw_diag(j_sw)%l_orog_corr   = .FALSE.

! Extinction and absorptivity diagnostics

        sw_diag(j_sw)%l_cloud_extinction           = .FALSE.
        sw_diag(j_sw)%l_cloud_weight_extinction    = .FALSE.
        sw_diag(j_sw)%l_ls_cloud_extinction        = .FALSE.
        sw_diag(j_sw)%l_ls_cloud_weight_extinction = .FALSE.
        sw_diag(j_sw)%l_cnv_cloud_extinction       = .FALSE.
        sw_diag(j_sw)%l_cnv_cloud_weight_extinction= .FALSE.

      IF (lhook) CALL dr_hook('INIT_SWDIAG_LOGIC',zhook_out,zhook_handle)
      RETURN

END SUBROUTINE init_swdiag_logic
