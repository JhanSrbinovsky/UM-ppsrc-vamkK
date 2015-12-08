! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sisl_resetcon_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_sisl_resetcon(                                          &
         row_length, rows, n_rows, model_levels, halo_i, halo_j,      &
         offx, offy,  l_datastart, l_slice,l_test_tracer, l_cartesian,&
         l_shallow, l_const_grav,eccentricity, height_domain,         &
         tprofile_number,trefer_number,dtheta_dz1,z_top_of_model,     &
         t_surf_ref_in, p_surf_ref_in, thetav,rho,t_surface_in,       &
         p_surface_in,p_star,exner, l_inc_solver,                     &
         m_v, m_r, m_gr, m_cl, m_cf, m_cf2,                           &
         f1_comp, f2_comp, f3_comp, ih )

USE atm_fields_bounds_mod
USE atmos_constants_mod
USE conversions_mod
USE earth_constants_mod
USE domain_params
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE proc_info_mod,         ONLY : me, global_row_length,          &
                                  global_rows, model_domain
USE timestep_mod,          ONLY : timestep,timestep_number
USE level_heights_mod,     ONLY : eta_theta_levels, &
                                  eta_rho_levels,   &
                                  xi3_at_theta=>r_theta_levels,        &
                                  xi3_at_rho=>r_rho_levels,&
                                  xi3_at_u=>r_at_u,&
                                  xi3_at_v=>r_at_v,&
                                  xi3_at_u_w=>r_at_u_w,&
                                  xi3_at_v_w=>r_at_v_w
USE eg_helmholtz_mod
USE eg_set_helmholtz_mod
USE eg_dry_static_adj_ref_pro_mod
USE integrity_mod
USE horiz_grid_mod
USE ref_pro_mod
USE metric_terms_mod
USE Field_Types
USE PrintStatus_mod
USE eg_alpha_mod
USE eg_swap_bounds_mod

USE tprofile_mod, ONLY: tp_dthetadz, tp_isothermal, tp_bruntv,        &
                        tp_bv_isoth, tp_dyn_core, tp_dyn_core_lam,    &
                        tp_namelist, tp_dump
IMPLICIT NONE
!
! Description: Resets the reference profile to last timestep state
!  
!
! Method: Reduced implementation of eg_SISL_setcon
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

REAL, INTENT(IN) :: ih

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,        &
                       l_datastart(3), halo_i, halo_j, offx,          &
                       offy, tprofile_number,           &
                       trefer_number

! Loop index bounds for arrays defined on p, u, v points respectively


LOGICAL, INTENT(IN) :: l_slice, l_test_tracer, l_inc_solver,          &
                       l_cartesian, l_shallow, l_const_grav

REAL, INTENT(IN) :: dtheta_dz1(3)

REAL, INTENT(IN) ::  z_top_of_model

! Surface temperature and pressure for refernce profile
REAL, INTENT(IN) :: t_surf_ref_in, p_surf_ref_in

REAL  t_surf_ref2d(1-offx:row_length+offx,                            &
                   1-offy:rows+offy)                                  &
,     p_surf_ref2d(1-offx:row_length+offx,                            &
                   1-offy:rows+offy)

! Surface temperature and pressure for initial data
REAL, INTENT(IN) :: t_surface_in, p_surface_in

REAL  t_surface2d(1-offx:row_length+offx,                             &
                  1-offy:rows+offy)                                   &
,     p_surface2d(1-offx:row_length+offx,                             &
                  1-offy:rows+offy)


REAL, INTENT(IN)  ::                                                  &
  f1_comp (row_length,rows),                                          &
  f2_comp (row_length,rows),                                          &
  f3_comp (row_length,rows)

REAL, INTENT(IN) :: eccentricity
REAL, INTENT(IN) :: height_domain


REAL, INTENT(INOUT) ::                                                &
  thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL, INTENT(IN)    ::                                                &
  m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
  m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),        &
  m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),       &
  m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),       &
  m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),       &
  m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

!Output Arrays from this routine
REAL, INTENT(INOUT) ::                                                &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1)


! Local variables

INTEGER :: i,j,k, option_bubble

REAL :: dtheta_dz1_ref(3)

! Tolerance for initial and reference surface pressures to be considered
! to differ:
REAL, PARAMETER :: p_tol=0.1

! Tolerance for initial and reference surface temperatures to be
! considered to differ:
REAL, PARAMETER :: t_tol=0.01


REAL, INTENT(IN) :: p_star(1-offx:row_length+offx, 1-offy:rows+offy)


REAL kp2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER filter

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SISL_RESETCON',zhook_in,zhook_handle)

reference_profile_changed = .TRUE.

kp2    = (1.0-kappa)/kappa

!-------------------------------------------------------------------------
! Initialise theta, rho, exner fields
!-------------------------------------------------------------------------

! If we are running with a different background state we need
! to recompute the reference profiles.
! Might as well always recompute!

If (PrintStatus >= PrStatus_Normal) THEN 
   WRITE(6,*)  'EG_SISL_Resetcon: calculate reference profile'
END IF

dtheta_dz1_ref      = 0.0 ! Reset the theta gradient to zero

DO k = 1, model_levels
   thetav_ref_pro(:,:,k) = thetav(:,:,k)
   exner_ref_pro(:,:,k)  = exner (:,:,k)
   rho_ref_pro(:,:,k)    = rho(:,:,k)
END DO
k = model_levels+1
exner_ref_pro(:,:,0)  = p_star(:,:)
exner_ref_pro(:,:,k)  = exner(:,:,k)
thetav_ref_pro(:,:,0) = thetav(:,:,0)





 
    CALL eg_dry_static_adj_ref_pro(thetav_ref_pro)

  DO k = pdims%k_start, pdims%k_end 
    DO j = pdims%j_start, pdims%j_end 
      DO i = pdims%i_start, pdims%i_end 

!   recompute density to satisfy equation of state
          rho_ref_pro (i,j,k) = p_zero/(R*                            &
                          (intw_w2rho(k,1)*thetav_ref_pro(i,j,k)      &
                         + intw_w2rho(k,2)*thetav_ref_pro(i,j,k-1) )) &
                         *(exner_ref_pro(i,j,k))**((1.0-kappa)/kappa)

      END DO
    END DO
  END DO


! DEPENDS ON: swap_bounds
    CALL swap_bounds(exner_ref_pro,row_length,rows,                   &
                            model_levels+2,                           &
                            offx,offy,fld_type_p,.FALSE.)
    CALL eg_swap_bounds(thetav_ref_pro,tdims_s,fld_type_p,.FALSE.)

    CALL eg_swap_bounds(rho_ref_pro,pdims_s,fld_type_p,.FALSE.)

!----------------------------------------------------------------------
! Compute Helmholtz coefficients
!----------------------------------------------------------------------

CALL eg_set_helmholtz (                                               &
       row_length, rows, n_rows, model_levels, halo_i, halo_j,        &
       offx, offy,  l_slice, l_test_tracer,                           &
       l_cartesian, l_shallow, l_inc_solver,                          &
       f1_comp, f2_comp, f3_comp, ih )



  IF (integrity_test)                                                 &
    CALL update_hash_m(exner_ref_pro,   SIZE(exner_ref_pro),  'piref',&
                       thetav_ref_pro,  SIZE(thetav_ref_pro), 'tvref',&
                       rho_ref_pro,     SIZE(rho_ref_pro),    'r_ref',&
                       thetav_ref_eta,  SIZE(thetav_ref_eta), 'tvree',&
                       rho_ref_eta,     SIZE(rho_ref_eta),    'rreta')

IF (lhook) CALL dr_hook('EG_SISL_RESETCON',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_sisl_resetcon
END MODULE eg_sisl_resetcon_mod
