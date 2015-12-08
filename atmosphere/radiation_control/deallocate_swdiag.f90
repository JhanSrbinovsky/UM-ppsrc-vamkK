! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Deallocate the SW diagnostics

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!    FORTRAN 77/90  with extensions listed in documentation.

!-----------------------------------------------------------------------

SUBROUTINE deallocate_swdiag(j_sw)

  USE sw_diag_mod
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: j_sw              ! call to SW radiation

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


! Deallocate SW Diagnostics

  IF (lhook) CALL dr_hook('DEALLOCATE_SWDIAG',zhook_in,zhook_handle)

  IF (ASSOCIATED(sw_diag(j_sw)%flux_up)) &
    DEALLOCATE(sw_diag(j_sw)%flux_up)
  IF (ASSOCIATED(sw_diag(j_sw)%flux_down)) &
    DEALLOCATE(sw_diag(j_sw)%flux_down)
  IF (ASSOCIATED(sw_diag(j_sw)%flux_up_clear)) &
    DEALLOCATE(sw_diag(j_sw)%flux_up_clear)
  IF (ASSOCIATED(sw_diag(j_sw)%flux_down_clear)) &
    DEALLOCATE(sw_diag(j_sw)%flux_down_clear)
  IF (ASSOCIATED(sw_diag(j_sw)%solar_out_toa)) &
    DEALLOCATE(sw_diag(j_sw)%solar_out_toa)
  IF (ASSOCIATED(sw_diag(j_sw)%solar_out_clear)) &
    DEALLOCATE(sw_diag(j_sw)%solar_out_clear)
  IF (ASSOCIATED(sw_diag(j_sw)%surface_down_flux)) &
    DEALLOCATE(sw_diag(j_sw)%surface_down_flux)
  IF (ASSOCIATED(sw_diag(j_sw)%surf_down_clr)) &
    DEALLOCATE(sw_diag(j_sw)%surf_down_clr)
  IF (ASSOCIATED(sw_diag(j_sw)%surf_up_clr)) &
    DEALLOCATE(sw_diag(j_sw)%surf_up_clr)
  IF (ASSOCIATED(sw_diag(j_sw)%net_flux_trop)) &
    DEALLOCATE(sw_diag(j_sw)%net_flux_trop)
  IF (ASSOCIATED(sw_diag(j_sw)%up_flux_trop)) &
    DEALLOCATE(sw_diag(j_sw)%up_flux_trop)
  IF (ASSOCIATED(sw_diag(j_sw)%clear_hr)) &
    DEALLOCATE(sw_diag(j_sw)%clear_hr)
! Radiances
  IF (ASSOCIATED(sw_diag(j_sw)%toa_radiance)) &
    DEALLOCATE(sw_diag(j_sw)%toa_radiance)
! UV-Fluxes
  IF (ASSOCIATED(sw_diag(j_sw)%uvflux_direct)) &
    DEALLOCATE(sw_diag(j_sw)%uvflux_direct)
  IF (ASSOCIATED(sw_diag(j_sw)%uvflux_up)) &
    DEALLOCATE(sw_diag(j_sw)%uvflux_up)
  IF (ASSOCIATED(sw_diag(j_sw)%uvflux_net)) &
    DEALLOCATE(sw_diag(j_sw)%uvflux_net)
  IF (ASSOCIATED(sw_diag(j_sw)%surf_uv)) &
    DEALLOCATE(sw_diag(j_sw)%surf_uv)
  IF (ASSOCIATED(sw_diag(j_sw)%surf_uv_clr)) &
    DEALLOCATE(sw_diag(j_sw)%surf_uv_clr)
! Surface Albedos 
  IF (ASSOCIATED(sw_diag(j_sw)%direct_albedo)) & 
    DEALLOCATE(sw_diag(j_sw)%direct_albedo) 
  IF (ASSOCIATED(sw_diag(j_sw)%diffuse_albedo)) & 
    DEALLOCATE(sw_diag(j_sw)%diffuse_albedo) 
  IF (ASSOCIATED(sw_diag(j_sw)%vis_albedo_sc)) & 
    DEALLOCATE(sw_diag(j_sw)%vis_albedo_sc) 
  IF (ASSOCIATED(sw_diag(j_sw)%nir_albedo_sc)) & 
    DEALLOCATE(sw_diag(j_sw)%nir_albedo_sc) 
! Direct and diffuse downward SW Fluxes
  IF (ASSOCIATED(sw_diag(j_sw)%flux_direct)) &
    DEALLOCATE(sw_diag(j_sw)%flux_direct)
  IF (ASSOCIATED(sw_diag(j_sw)%flux_diffuse)) &
    DEALLOCATE(sw_diag(j_sw)%flux_diffuse)
! Microphysical diagnostics
  IF (ASSOCIATED(sw_diag(j_sw)%re_strat)) &
    DEALLOCATE(sw_diag(j_sw)%re_strat)
  IF (ASSOCIATED(sw_diag(j_sw)%wgt_strat)) &
    DEALLOCATE(sw_diag(j_sw)%wgt_strat)
  IF (ASSOCIATED(sw_diag(j_sw)%lwp_strat)) &
    DEALLOCATE(sw_diag(j_sw)%lwp_strat)
  IF (ASSOCIATED(sw_diag(j_sw)%re_conv)) &
    DEALLOCATE(sw_diag(j_sw)%re_conv)
  IF (ASSOCIATED(sw_diag(j_sw)%wgt_conv)) &
    DEALLOCATE(sw_diag(j_sw)%wgt_conv)
  IF (ASSOCIATED(sw_diag(j_sw)%ntot_diag)) &
    DEALLOCATE(sw_diag(j_sw)%ntot_diag)
  IF (ASSOCIATED(sw_diag(j_sw)%strat_lwc_diag)) &
    DEALLOCATE(sw_diag(j_sw)%strat_lwc_diag)
  IF (ASSOCIATED(sw_diag(j_sw)%so4_ccn_diag)) &
    DEALLOCATE(sw_diag(j_sw)%so4_ccn_diag)
  IF (ASSOCIATED(sw_diag(j_sw)%cond_samp_wgt)) &
    DEALLOCATE(sw_diag(j_sw)%cond_samp_wgt)
  IF (ASSOCIATED(sw_diag(j_sw)%weighted_re)) &
    DEALLOCATE(sw_diag(j_sw)%weighted_re)
  IF (ASSOCIATED(sw_diag(j_sw)%sum_weight_re)) &
    DEALLOCATE(sw_diag(j_sw)%sum_weight_re)
  IF (ASSOCIATED(sw_diag(j_sw)%weighted_warm_re)) &
    DEALLOCATE(sw_diag(j_sw)%weighted_warm_re)
  IF (ASSOCIATED(sw_diag(j_sw)%sum_weight_warm_re)) &
    DEALLOCATE(sw_diag(j_sw)%sum_weight_warm_re)
  IF (ASSOCIATED(sw_diag(j_sw)%nc_diag)) &
    DEALLOCATE(sw_diag(j_sw)%nc_diag)
  IF (ASSOCIATED(sw_diag(j_sw)%nc_weight)) &
    DEALLOCATE(sw_diag(j_sw)%nc_weight)
! Diagnostics for MOSES II
  IF (ASSOCIATED(sw_diag(j_sw)%FlxSolBelow690nmSurf)) &
    DEALLOCATE(sw_diag(j_sw)%FlxSolBelow690nmSurf)
  IF (ASSOCIATED(sw_diag(j_sw)%FlxSeaBelow690nmSurf)) &
    DEALLOCATE(sw_diag(j_sw)%FlxSeaBelow690nmSurf)
! Orography correction diagnostics:
  IF (ASSOCIATED(sw_diag(j_sw)%orog_corr)) &
    DEALLOCATE(sw_diag(j_sw)%orog_corr)
  IF (ASSOCIATED(sw_diag(j_sw)%sol_bearing)) &
    DEALLOCATE(sw_diag(j_sw)%sol_bearing)

  IF (lhook) CALL dr_hook('DEALLOCATE_SWDIAG',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE deallocate_swdiag
