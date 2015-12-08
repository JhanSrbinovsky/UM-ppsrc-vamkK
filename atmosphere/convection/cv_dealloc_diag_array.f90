! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! DEALLOCATE arrays for convective diagnostics
!
SUBROUTINE cv_dealloc_diag_array( )


USE cv_diagnostic_array_mod , ONLY:                                     &
        precip_deep, precip_shall ,precip_mid ,precip_cong ,cape_out    &
       ,deep_ind ,shallow_ind ,congestus_ind ,congestus_ind2 ,mid_ind   &
       ,ntml_diag ,ntpar_diag ,freeze_diag ,kterm_diag                  &
       ,wstar_up_diag ,wstar_dn_diag ,mb1_diag ,mb2_diag ,wqt_cb        &
       ,wthetal_cb ,wqt_inv ,wthetal_inv ,sh_top ,sh_base ,cg_top       &
       ,cg_base ,cg_term ,cape_ts_diag ,ind_cape_reduced_diag           &
       ,conscav_dust ,conscav_so4ait ,conscav_so4acc ,conscav_so4dis    &
       ,conscav_agedsoot ,conscav_agedbmass ,conscav_agedocff           &
       ,conscav_nitracc ,conscav_nitrdiss ,conwash_so2 ,conwash_nh3     &
       ,uw_dp, vw_dp ,uw_shall ,vw_shall ,uw_mid ,vw_mid ,up_flux       &
       ,dwn_flux ,entrain_up ,detrain_up ,entrain_dwn ,detrain_dwn      &
       ,up_flux_half ,T_incr_diag_conv ,q_incr_diag_conv                &
       ,qcl_incr_diag_conv ,qcf_incr_diag_conv                          & 
       ,cf_liquid_incr_diag_conv ,cf_frozen_incr_diag_conv              &
       ,bulk_cf_incr_diag_conv,u_incr_diag_conv ,v_incr_diag_conv       &
       ,theta_diag, q_diag                                              &
       ,mf_deep ,mf_congest ,mf_shall ,mf_midlev                        &
       ,dt_deep ,dt_congest ,dt_shall ,dt_midlev                        &
       ,dq_deep ,dq_congest ,dq_shall ,dq_midlev                        &
       ,du_deep ,du_congest ,du_shall ,du_midlev                        &
       ,dv_deep ,dv_congest ,dv_shall ,dv_midlev                        &
       ,wqt_flux_sh ,wql_flux_sh ,wthetal_flux_sh ,wthetav_flux_sh      &
       ,dubydt_pout ,dvbydt_pout ,conv_rain_3d ,conv_snow_3d            &
       ,t_incr_conv_only,q_incr_conv_only                               &
       ,qcl_incr_inhom_diag ,qcf_incr_inhom_diag                        &
       ,bulk_cf_incr_inhom_diag ,cf_liquid_incr_inhom_diag              &
       ,cf_frozen_incr_inhom_diag, deep_cfl_limited, mid_cfl_limited    &
       ,deep_tops

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   DEALLOCATE arrays required by convection for a full model timestep on the
!   first convection substep of a model timestep.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.

! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------




! Required by Dr hook 
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('CV_DEALLOC_DIAG_ARRAY',zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Clear up allocatable arrays

  IF(ALLOCATED (qcl_incr_inhom_diag)) THEN
    DEALLOCATE (qcl_incr_inhom_diag)
  END IF

  IF(ALLOCATED (qcf_incr_inhom_diag)) THEN
     DEALLOCATE (qcf_incr_inhom_diag)
  END IF

  IF(ALLOCATED (bulk_cf_incr_inhom_diag)) THEN
     DEALLOCATE (bulk_cf_incr_inhom_diag)
  END IF

  IF(ALLOCATED (cf_liquid_incr_inhom_diag)) THEN
    DEALLOCATE (cf_liquid_incr_inhom_diag)
  END IF

  IF(ALLOCATED (cf_frozen_incr_inhom_diag)) THEN
    DEALLOCATE (cf_frozen_incr_inhom_diag)
  END IF

DEALLOCATE(deep_tops)
DEALLOCATE(vw_mid)
DEALLOCATE(uw_mid)
DEALLOCATE(vw_shall)
DEALLOCATE(uw_shall)
DEALLOCATE(vw_dp)
DEALLOCATE(uw_dp)
 
DEALLOCATE(detrain_dwn)
DEALLOCATE(entrain_dwn)
DEALLOCATE(detrain_up)
DEALLOCATE(entrain_up)

DEALLOCATE(dwn_flux)
DEALLOCATE(up_flux)
DEALLOCATE(up_flux_half)

DEALLOCATE(dv_midlev)
DEALLOCATE(dv_shall)
DEALLOCATE(dv_congest)
DEALLOCATE(dv_deep)
DEALLOCATE(du_midlev)
DEALLOCATE(du_shall)
DEALLOCATE(du_congest)
DEALLOCATE(du_deep)
DEALLOCATE(dq_midlev)
DEALLOCATE(dq_shall)
DEALLOCATE(dq_congest)
DEALLOCATE(dq_deep)
DEALLOCATE(dt_midlev)
DEALLOCATE(dt_shall)
DEALLOCATE(dt_congest)
DEALLOCATE(dt_deep)
DEALLOCATE(mf_midlev)
DEALLOCATE(mf_shall)
DEALLOCATE(mf_congest)
DEALLOCATE(mf_deep)

DEALLOCATE(dvbydt_pout)
DEALLOCATE(dubydt_pout)

DEALLOCATE(wthetav_flux_sh)
DEALLOCATE(wthetal_flux_sh)
DEALLOCATE(wql_flux_sh)
DEALLOCATE(wqt_flux_sh)

DEALLOCATE(conv_snow_3d)
DEALLOCATE(conv_rain_3d)

DEALLOCATE(q_incr_conv_only)
DEALLOCATE(T_incr_conv_only)

DEALLOCATE(v_incr_diag_conv)
DEALLOCATE(u_incr_diag_conv)
DEALLOCATE(bulk_cf_incr_diag_conv)
DEALLOCATE(cf_frozen_incr_diag_conv)
DEALLOCATE(cf_liquid_incr_diag_conv)
DEALLOCATE(qcf_incr_diag_conv)
DEALLOCATE(qcl_incr_diag_conv)
DEALLOCATE(q_incr_diag_conv)
DEALLOCATE(T_incr_diag_conv)

DEALLOCATE(q_diag)
DEALLOCATE(theta_diag)

! 2d arrays
! Convective aerosol relative diagnostics


 DEALLOCATE(conscav_dust)

 DEALLOCATE(conscav_so4ait )
 DEALLOCATE(conscav_so4acc )
 DEALLOCATE(conscav_so4dis )
 DEALLOCATE(conscav_agedsoot )
 DEALLOCATE(conscav_agedbmass )
 DEALLOCATE(conscav_agedocff )
 DEALLOCATE(conscav_nitracc )
 DEALLOCATE(conscav_nitrdiss )
                                    
 DEALLOCATE(conwash_so2 )
 DEALLOCATE(conwash_nh3 )

 DEALLOCATE(precip_deep)
 DEALLOCATE(precip_shall)
 DEALLOCATE(precip_mid) 
 DEALLOCATE(precip_cong)
 DEALLOCATE(cape_out)
 DEALLOCATE(deep_ind)
 DEALLOCATE(shallow_ind)
 DEALLOCATE(congestus_ind)
 DEALLOCATE(congestus_ind2)
 DEALLOCATE(mid_ind)
 DEALLOCATE(ntml_diag) 
 DEALLOCATE(ntpar_diag)
 DEALLOCATE(freeze_diag)
 DEALLOCATE(kterm_diag)
 DEALLOCATE(wstar_up_diag)
 DEALLOCATE(wstar_dn_diag)
 DEALLOCATE(mb1_diag)
 DEALLOCATE(mb2_diag)
 DEALLOCATE(wqt_cb )
 DEALLOCATE(wthetal_cb ) 
 DEALLOCATE(wqt_inv )
 DEALLOCATE(wthetal_inv )
 DEALLOCATE(sh_top )
 DEALLOCATE(sh_base )
 DEALLOCATE(cg_top )
 DEALLOCATE(cg_base )
 DEALLOCATE(cg_term )
 DEALLOCATE(cape_ts_diag )
 DEALLOCATE(ind_cape_reduced_diag )
 DEALLOCATE(deep_cfl_limited)
 DEALLOCATE(mid_cfl_limited)

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook('CV_DEALLOC_DIAG_ARRAY',zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE cv_dealloc_diag_array
