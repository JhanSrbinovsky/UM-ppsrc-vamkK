! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Set the SW diagnostic flags to true if necessary

! Method:

! If a radiation diagnostic has been chosen in STASH then the
! flag of the corresponding radiation diagnostic in the structure
! SW_diag is set to true.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!-----------------------------------------------------------------------

SUBROUTINE set_swdiag_logic(sf,nitems,nsects, l_radiance, j_sw, i_off)

      USE lw_diag_mod
      USE sw_diag_mod
      USE rad_input_mod, ONLY: l_rad_perturb
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

! arguments with intent in

      INTEGER, INTENT(IN) :: nitems    ! item number
      INTEGER, INTENT(IN) :: nsects    ! section number
      INTEGER, INTENT(IN) :: j_sw      ! call to SW radiation
      INTEGER, INTENT(IN) :: i_off     ! offset for diagnostics

      LOGICAL, INTENT(IN) :: l_radiance    ! Flag for radiances

      LOGICAL, INTENT(IN) :: sf(0:nitems,0:nsects)
!        STASH Flags

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('SET_SWDIAG_LOGIC',zhook_in,zhook_handle)



! Switches enabled as STASHflags: SW

      IF (l_radiance.AND.(j_sw > 1)) THEN

        sw_diag(j_sw)%l_toa_radiance      = sf(297+i_off,1)

      ELSE

! Fluxes and Heating Rates

        sw_diag(j_sw)%l_flux_up           = sf(217+i_off,1)
        sw_diag(j_sw)%l_flux_down         = sf(218+i_off,1)
        sw_diag(j_sw)%l_solar_out_toa     = sf(208+i_off,1)
        sw_diag(j_sw)%l_surface_down_flux = sf(235+i_off,1)
        sw_diag(j_sw)%l_net_flux_trop     = sf(237+i_off,1)
        sw_diag(j_sw)%l_up_flux_trop      = sf(238+i_off,1)

! Diffuse & Direct Flux

        sw_diag(j_sw)%l_flux_direct       = sf(230+i_off,1)
        sw_diag(j_sw)%l_flux_diffuse      = sf(231+i_off,1)

! Diagnostics for photosynthetically active surface radiation

        sw_diag(j_sw)%l_FlxSolBelow690nmSurf = sf(259+i_off,1)
        sw_diag(j_sw)%l_FlxSeaBelow690nmSurf = sf(260+i_off,1)

! Diagnostics for orography correction

        sw_diag(j_sw)%l_orog_corr   = sf(295+i_off,1)

! Diagnostics for sea salt

        sw_diag(j_sw)%seasalt_film_flag   = sf(247+i_off,1)
        sw_diag(j_sw)%seasalt_jet_flag    = sf(248+i_off,1)

! Albedo scaling to obs diagnostics

        sw_diag(j_sw)%l_vis_albedo_sc        = sf(270+i_off,1) 
        sw_diag(j_sw)%l_nir_albedo_sc        = sf(271+i_off,1) 

      IF (.NOT.(l_rad_perturb.AND.(j_sw == 2))) THEN
!       For the incremental time-stepping scheme (l_rad_perturb) many
!       of the diagnostics are not calculated on the "cloud only"
!       radiation calls (j_sw==2).

! Clear-sky Fluxes and Heating Rates

        sw_diag(j_sw)%l_flux_up_clear     = sf(219+i_off,1)
        sw_diag(j_sw)%l_flux_down_clear   = sf(220+i_off,1)
        sw_diag(j_sw)%l_solar_out_clear   = sf(209+i_off,1)
        sw_diag(j_sw)%l_surf_down_clr     = sf(210+i_off,1)
        sw_diag(j_sw)%l_surf_up_clr       = sf(211+i_off,1)
        sw_diag(j_sw)%l_clear_hr          = sf(233+i_off,1)
        sw_diag(j_sw)%l_uvflux_direct     = sf(212+i_off,1)
        sw_diag(j_sw)%l_uvflux_up         = sf(213+i_off,1)
        sw_diag(j_sw)%l_uvflux_net        = sf(214+i_off,1)
        sw_diag(j_sw)%l_surf_uv           = sf(288+i_off,1)
        sw_diag(j_sw)%l_surf_uv_clr       = sf(289+i_off,1)

! Microphysical diagnostics

        sw_diag(j_sw)%re_strat_flag       = sf(221+i_off,1)
        sw_diag(j_sw)%wgt_strat_flag      = sf(223+i_off,1)
        sw_diag(j_sw)%lwp_strat_flag      = sf(224+i_off,1)
        sw_diag(j_sw)%re_conv_flag        = sf(225+i_off,1)
        sw_diag(j_sw)%wgt_conv_flag       = sf(226+i_off,1)
        sw_diag(j_sw)%ntot_diag_flag      = sf(241+i_off,1)
        sw_diag(j_sw)%strat_lwc_diag_flag = sf(242+i_off,1)
        sw_diag(j_sw)%so4_ccn_diag_flag   = sf(243+i_off,1)
        sw_diag(j_sw)%cond_samp_wgt_flag  = sf(244+i_off,1)
        sw_diag(j_sw)%weighted_re_flag    = sf(245+i_off,1)
        sw_diag(j_sw)%sum_weight_re_flag  = sf(246+i_off,1)
        sw_diag(j_sw)%wgtd_warm_re_flag   = sf(254+i_off,1)
        sw_diag(j_sw)%sum_wgt_warm_re_flag= sf(255+i_off,1)
        sw_diag(j_sw)%nc_diag_flag        = sf(280+i_off,1)
        sw_diag(j_sw)%nc_weight_flag      = sf(281+i_off,1)

! Diagnostics for orography correction

        sw_diag(j_sw)%l_sol_bearing = sf(292+i_off,1)

! Extinction and absorptivity diagnostics

        sw_diag(j_sw)%l_cloud_extinction           = sf(262+i_off,1)
        sw_diag(j_sw)%l_cloud_weight_extinction    = sf(263+i_off,1)
        sw_diag(j_sw)%l_ls_cloud_extinction        = (sf(264+i_off,1)  &
           .OR.(lw_diag(j_sw)%l_isccp_weights))
        sw_diag(j_sw)%l_ls_cloud_weight_extinction = (sf(265+i_off,1)  &
           .OR.(lw_diag(j_sw)%l_isccp_weights))
        sw_diag(j_sw)%l_cnv_cloud_extinction       = (sf(266+i_off,1)  &
           .OR.(lw_diag(j_sw)%l_isccp_weights))
        sw_diag(j_sw)%l_cnv_cloud_weight_extinction= (sf(267+i_off,1)  &
           .OR.(lw_diag(j_sw)%l_isccp_weights))

! Surface Albedo Diagnostics 
        sw_diag(j_sw)%l_direct_albedo        = sf(268+i_off,1) 
        sw_diag(j_sw)%l_diffuse_albedo       = sf(269+i_off,1) 

      END IF ! .not.(l_rad_perturb.and.(j_sw == 2))

      END IF

IF (lhook) CALL dr_hook('SET_SWDIAG_LOGIC',zhook_out,zhook_handle)
RETURN

END SUBROUTINE set_swdiag_logic
