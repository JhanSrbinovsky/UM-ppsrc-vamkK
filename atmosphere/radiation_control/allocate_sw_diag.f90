! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Allocate the required memory space for the radiation
!          diagnostics.

! Method: If a radiation diagnostic is required the correct amount
!         of memory space is allocated otherise a minimal amount
!         of space is allocated.
!         All radiation diagnostics are initialised to zero

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!    FORTRAN 77/90  with extensions listed in documentation.

!-----------------------------------------------------------------------

SUBROUTINE allocate_sw_diag(row_length, rows, model_levels,     &
                            cloud_levels, ntiles, j_sw)

! Modules
   USE spec_sw_lw
   USE sw_diag_mod

   USE yomhook, ONLY: lhook, dr_hook
   USE parkind1, ONLY: jprb, jpim
   USE missing_data_mod, ONLY: rmdi
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: row_length       ! Length of rows
   INTEGER, INTENT(IN) :: rows             ! Number of rows
   INTEGER, INTENT(IN) :: model_levels     ! Number of model levels
   INTEGER, INTENT(IN) :: cloud_levels     ! Number of cloud levels
   INTEGER, INTENT(IN) :: ntiles           ! Number of land surface tiles 
   INTEGER, INTENT(IN) :: j_sw             ! call to SW radiation

   INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
   INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
   REAL(KIND=jprb)               :: zhook_handle

   IF (lhook) CALL dr_hook('ALLOCATE_SW_DIAG',zhook_in,zhook_handle)


   IF (sw_diag(j_sw)%l_flux_up) THEN
     ALLOCATE(sw_diag(j_sw)%flux_up(row_length,rows,model_levels+1))
   ELSE
     ALLOCATE(sw_diag(j_sw)%flux_up(1,1,1))
   END IF
   sw_diag(j_sw)%flux_up = 0.0

   IF (sw_diag(j_sw)%l_flux_down) THEN
     ALLOCATE(sw_diag(j_sw)%flux_down(row_length,rows,model_levels+1))
   ELSE
     ALLOCATE(sw_diag(j_sw)%flux_down(1,1,1))
   END IF
   sw_diag(j_sw)%flux_down = 0.0

   IF (sw_diag(j_sw)%l_flux_up_clear) THEN
     ALLOCATE(sw_diag(j_sw)%flux_up_clear(row_length,rows, &
                                          model_levels+1))
   ELSE
     ALLOCATE(sw_diag(j_sw)%flux_up_clear(1,1,1))
   END IF
   sw_diag(j_sw)%flux_up_clear = 0.0

   IF (sw_diag(j_sw)%l_flux_down_clear) THEN
     ALLOCATE(sw_diag(j_sw)%flux_down_clear(row_length,rows, &
                                            model_levels+1))
   ELSE
     ALLOCATE(sw_diag(j_sw)%flux_down_clear(1,1,1))
   END IF
   sw_diag(j_sw)%flux_down_clear = 0.0

   IF (sw_diag(j_sw)%l_solar_out_toa) THEN
      ALLOCATE(sw_diag(j_sw)%solar_out_toa(row_length,rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%solar_out_toa(1,1))
   END IF
   sw_diag(j_sw)%solar_out_toa = 0.0

   IF (sw_diag(j_sw)%l_solar_out_clear) THEN
      ALLOCATE(sw_diag(j_sw)%solar_out_clear(row_length,rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%solar_out_clear(1,1))
   END IF
   sw_diag(j_sw)%solar_out_clear = 0.0

   IF (sw_diag(j_sw)%l_surface_down_flux) THEN
      ALLOCATE(sw_diag(j_sw)%surface_down_flux(row_length,rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%surface_down_flux(1,1))
   END IF
   sw_diag(j_sw)%surface_down_flux = 0.0

   IF (sw_diag(j_sw)%l_surf_down_clr) THEN
      ALLOCATE(sw_diag(j_sw)%surf_down_clr(row_length,rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%surf_down_clr(1,1))
   END IF
   sw_diag(j_sw)%surf_down_clr = 0.0

   IF (sw_diag(j_sw)%l_surf_up_clr) THEN
      ALLOCATE(sw_diag(j_sw)%surf_up_clr(row_length,rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%surf_up_clr(1,1))
   END IF
   sw_diag(j_sw)%surf_up_clr = 0.0

   IF (sw_diag(j_sw)%l_net_flux_trop) THEN
      ALLOCATE(sw_diag(j_sw)%net_flux_trop(row_length,rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%net_flux_trop(1,1))
   END IF
   sw_diag(j_sw)%net_flux_trop = 0.0

   IF (sw_diag(j_sw)%l_up_flux_trop) THEN
      ALLOCATE(sw_diag(j_sw)%up_flux_trop(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%up_flux_trop(1,1))
   END IF
   sw_diag(j_sw)%up_flux_trop = 0.0

   IF (sw_diag(j_sw)%l_clear_hr) THEN
      ALLOCATE(sw_diag(j_sw)%clear_hr(row_length, rows, model_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%clear_hr(1,1,1))
   END IF
   sw_diag(j_sw)%clear_hr = 0.0

!  Radiance
   IF (sw_diag(j_sw)%l_toa_radiance) THEN
      ALLOCATE(sw_diag(j_sw)%toa_radiance(row_length, rows,       &
                                      sw_spectrum(j_sw)%n_band))
   ELSE
      ALLOCATE(sw_diag(j_sw)%toa_radiance(1,1,1))
   END IF
   sw_diag(j_sw)%toa_radiance = 0.0

! Direct and Diffuse Downward Flux
   IF (sw_diag(j_sw)%l_flux_direct) THEN
      ALLOCATE(sw_diag(j_sw)%flux_direct(row_length, rows,        &
                                         model_levels+1))
   ELSE
      ALLOCATE(sw_diag(j_sw)%flux_direct(1,1,1))
   END IF
   sw_diag(j_sw)%flux_direct = 0.0

   IF (sw_diag(j_sw)%l_flux_diffuse) THEN
      ALLOCATE(sw_diag(j_sw)%flux_diffuse(row_length, rows,       &
                                          model_levels+1))
   ELSE
      ALLOCATE(sw_diag(j_sw)%flux_diffuse(1,1,1))
   END IF
   sw_diag(j_sw)%flux_diffuse = 0.0

! UV-Fluxes
   IF (sw_diag(j_sw)%l_uvflux_direct) THEN
      ALLOCATE(sw_diag(j_sw)%uvflux_direct(row_length, rows,      &
                                          model_levels+1))
   ELSE
      ALLOCATE(sw_diag(j_sw)%uvflux_direct(1,1,1))
   END IF
   sw_diag(j_sw)%uvflux_direct = 0.0

   IF (sw_diag(j_sw)%l_uvflux_up) THEN
      ALLOCATE(sw_diag(j_sw)%uvflux_up(row_length, rows,          &
                                          model_levels+1))
   ELSE
      ALLOCATE(sw_diag(j_sw)%uvflux_up(1,1,1))
   END IF
   sw_diag(j_sw)%uvflux_up = 0.0

   IF (sw_diag(j_sw)%l_uvflux_net) THEN
      ALLOCATE(sw_diag(j_sw)%uvflux_net(row_length, rows,         &
                                          model_levels+1))
   ELSE
      ALLOCATE(sw_diag(j_sw)%uvflux_net(1,1,1))
   END IF
   sw_diag(j_sw)%uvflux_net = 0.0

   IF (sw_diag(j_sw)%l_surf_uv) THEN
      ALLOCATE(sw_diag(j_sw)%surf_uv(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%surf_uv(1,1))
   END IF
   sw_diag(j_sw)%surf_uv = 0.0

   IF (sw_diag(j_sw)%l_surf_uv_clr) THEN
      ALLOCATE(sw_diag(j_sw)%surf_uv_clr(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%surf_uv_clr(1,1))
   END IF
   sw_diag(j_sw)%surf_uv_clr = 0.0

! Surface Albedo Diagnostics 
   IF (sw_diag(j_sw)%l_direct_albedo) THEN 
      ALLOCATE(sw_diag(j_sw)%direct_albedo(row_length, rows,      & 
                             sw_spectrum(j_sw)%n_band )) 
   ELSE 
      ALLOCATE(sw_diag(j_sw)%direct_albedo(1,1,1)) 
   END IF 
   sw_diag(j_sw)%direct_albedo = rmdi 
 
   IF (sw_diag(j_sw)%l_diffuse_albedo) THEN 
      ALLOCATE(sw_diag(j_sw)%diffuse_albedo(row_length, rows,     & 
                             sw_spectrum(j_sw)%n_band )) 
   ELSE 
      ALLOCATE(sw_diag(j_sw)%diffuse_albedo(1,1,1)) 
   END IF 
   sw_diag(j_sw)%diffuse_albedo = rmdi 
 
   IF (sw_diag(j_sw)%l_vis_albedo_sc) THEN 
      ALLOCATE(sw_diag(j_sw)%vis_albedo_sc(row_length, rows,      & 
                             ntiles )) 
   ELSE 
      ALLOCATE(sw_diag(j_sw)%vis_albedo_sc(1,1,1)) 
   END IF 
   sw_diag(j_sw)%vis_albedo_sc = rmdi 
 
   IF (sw_diag(j_sw)%l_nir_albedo_sc) THEN 
      ALLOCATE(sw_diag(j_sw)%nir_albedo_sc(row_length, rows,      & 
                             ntiles )) 
   ELSE 
     ALLOCATE(sw_diag(j_sw)%nir_albedo_sc(1,1,1)) 
   END IF 
   sw_diag(j_sw)%nir_albedo_sc = rmdi 


! Microphysical diagnostics
   IF (sw_diag(j_sw)%re_strat_flag) THEN
     ALLOCATE(sw_diag(j_sw)%re_strat(row_length, rows, cloud_levels))
   ELSE
     ALLOCATE(sw_diag(j_sw)%re_strat(1,1,1))
   END IF
   sw_diag(j_sw)%re_strat = 0.0

   IF (sw_diag(j_sw)%wgt_strat_flag) THEN
     ALLOCATE(sw_diag(j_sw)%wgt_strat(row_length, rows, cloud_levels))
   ELSE
     ALLOCATE(sw_diag(j_sw)%wgt_strat(1,1,1))
   END IF
   sw_diag(j_sw)%wgt_strat = 0.0

   IF (sw_diag(j_sw)%lwp_strat_flag) THEN
     ALLOCATE(sw_diag(j_sw)%lwp_strat(row_length, rows, cloud_levels))
   ELSE
     ALLOCATE(sw_diag(j_sw)%lwp_strat(1,1,1))
   END IF
   sw_diag(j_sw)%lwp_strat = 0.0

   IF (sw_diag(j_sw)%re_conv_flag) THEN
      ALLOCATE(sw_diag(j_sw)%re_conv(row_length, rows, cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%re_conv(1,1,1))
   END IF
   sw_diag(j_sw)%re_conv = 0.0

   IF (sw_diag(j_sw)%wgt_conv_flag) THEN
      ALLOCATE(sw_diag(j_sw)%wgt_conv(row_length, rows,cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%wgt_conv(1,1,1))
   END IF
   sw_diag(j_sw)%wgt_conv = 0.0

   IF (sw_diag(j_sw)%ntot_diag_flag) THEN
      ALLOCATE(sw_diag(j_sw)%ntot_diag(row_length, rows, cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%ntot_diag(1,1,1))
   END IF
   sw_diag(j_sw)%ntot_diag = 0.0

   IF (sw_diag(j_sw)%strat_lwc_diag_flag) THEN
      ALLOCATE(sw_diag(j_sw)%strat_lwc_diag(row_length, rows,      &
                                            cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%strat_lwc_diag(1,1,1))
   END IF
   sw_diag(j_sw)%strat_lwc_diag = 0.0

   IF (sw_diag(j_sw)%so4_ccn_diag_flag) THEN
      ALLOCATE(sw_diag(j_sw)%so4_ccn_diag(row_length, rows,        &
                                          cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%so4_ccn_diag(1,1,1))
   END IF
   sw_diag(j_sw)%so4_ccn_diag = 0.0

   IF (sw_diag(j_sw)%cond_samp_wgt_flag) THEN
      ALLOCATE(sw_diag(j_sw)%cond_samp_wgt(row_length, rows,       &
                                           cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%cond_samp_wgt(1,1,1))
   END IF
   sw_diag(j_sw)%cond_samp_wgt = 0.0

   IF (sw_diag(j_sw)%weighted_re_flag) THEN
      ALLOCATE(sw_diag(j_sw)%weighted_re(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%weighted_re(1,1))
   END IF
   sw_diag(j_sw)%weighted_re = 0.0

   IF (sw_diag(j_sw)%sum_weight_re_flag) THEN
      ALLOCATE(sw_diag(j_sw)%sum_weight_re(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%sum_weight_re(1,1))
   END IF
   sw_diag(j_sw)%sum_weight_re = 0.0

   IF (sw_diag(j_sw)%wgtd_warm_re_flag) THEN
      ALLOCATE(sw_diag(j_sw)%weighted_warm_re(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%weighted_warm_re(1,1))
   END IF
   sw_diag(j_sw)%weighted_warm_re = 0.0

   IF (sw_diag(j_sw)%sum_wgt_warm_re_flag) THEN
      ALLOCATE(sw_diag(j_sw)%sum_weight_warm_re(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%sum_weight_warm_re(1,1))
   END IF
   sw_diag(j_sw)%sum_weight_warm_re = 0.0

   IF (sw_diag(j_sw)%nc_diag_flag) THEN
     ALLOCATE(sw_diag(j_sw)%nc_diag(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%nc_diag(1,1))
   END IF
   sw_diag(j_sw)%nc_diag = 0.0

   IF (sw_diag(j_sw)%nc_weight_flag) THEN
      ALLOCATE(sw_diag(j_sw)%nc_weight(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%nc_weight(1,1))
   END IF
   sw_diag(j_sw)%nc_weight = 0.0

! Diagnostics for MOSES II:
   IF (sw_diag(j_sw)%l_FlxSolBelow690nmSurf) THEN
      ALLOCATE(sw_diag(j_sw)%FlxSolBelow690nmSurf(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%FlxSolBelow690nmSurf(1,1))
   END IF
   sw_diag(j_sw)%FlxSolBelow690nmSurf = 0.0

   IF (sw_diag(j_sw)%l_FlxSeaBelow690nmSurf) THEN
      ALLOCATE(sw_diag(j_sw)%FlxSeaBelow690nmSurf(row_length, rows))
   ELSE
      ALLOCATE(sw_diag(j_sw)%FlxSeaBelow690nmSurf(1,1))
   END IF
   sw_diag(j_sw)%FlxSeaBelow690nmSurf = 0.0

! Orography correction diagnostics:
   IF (sw_diag(j_sw)%l_orog_corr) THEN
     ALLOCATE(sw_diag(j_sw)%orog_corr(row_length, rows))
   ELSE
     ALLOCATE(sw_diag(j_sw)%orog_corr(1,1))
   END IF
   sw_diag(j_sw)%orog_corr = 1.0

   IF (sw_diag(j_sw)%l_sol_bearing) THEN
     ALLOCATE(sw_diag(j_sw)%sol_bearing(row_length, rows))
   ELSE
     ALLOCATE(sw_diag(j_sw)%sol_bearing(1,1))
   END IF
   sw_diag(j_sw)%sol_bearing = 0.0

! Extinction diagnostics:
   IF (sw_diag(j_sw)%l_cloud_extinction) THEN
      ALLOCATE(sw_diag(j_sw)%cloud_extinction(row_length,rows,      &
                                              cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%cloud_extinction(1,1,1))
   END IF
   sw_diag(j_sw)%cloud_extinction = 0.0

   IF (sw_diag(j_sw)%l_cloud_weight_extinction) THEN
      ALLOCATE(sw_diag(j_sw)%cloud_weight_extinction(row_length,    &
                                       rows,cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%cloud_weight_extinction(1,1,1))
   END IF
   sw_diag(j_sw)%cloud_weight_extinction = 0.0

   IF (sw_diag(j_sw)%l_ls_cloud_extinction) THEN
      ALLOCATE(sw_diag(j_sw)%ls_cloud_extinction(row_length,        &
                                      rows, cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%ls_cloud_extinction(1,1,1))
   END IF
   sw_diag(j_sw)%ls_cloud_extinction = 0.0

   IF (sw_diag(j_sw)%l_ls_cloud_weight_extinction) THEN
      ALLOCATE(sw_diag(j_sw)%ls_cloud_weight_extinction(            &
                                row_length,rows,cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%ls_cloud_weight_extinction(1,1,1))
   END IF
   sw_diag(j_sw)%ls_cloud_weight_extinction = 0.0

   IF (sw_diag(j_sw)%l_cnv_cloud_extinction) THEN
      ALLOCATE(sw_diag(j_sw)%cnv_cloud_extinction(row_length,       &
                                           rows,cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%cnv_cloud_extinction(1,1,1))
   END IF
   sw_diag(j_sw)%cnv_cloud_extinction = 0.0

   IF (sw_diag(j_sw)%l_cnv_cloud_weight_extinction) THEN
      ALLOCATE(sw_diag(j_sw)%cnv_cloud_weight_extinction(           &
                               row_length,rows, cloud_levels))
   ELSE
      ALLOCATE(sw_diag(j_sw)%cnv_cloud_weight_extinction(1,1,1))
   END IF
   sw_diag(j_sw)%cnv_cloud_weight_extinction = 0.0

   IF (lhook) CALL dr_hook('ALLOCATE_SW_DIAG',zhook_out,zhook_handle)

END SUBROUTINE allocate_sw_diag
