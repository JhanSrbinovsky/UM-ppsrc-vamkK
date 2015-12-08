! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_full_wind_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_sl_full_wind(                                           &
                row_length, rows, n_rows, model_levels, halo_i,       &
                halo_j, offx, offy,  g_i_pe,depart_scheme,            &
                depart_order, high_order_scheme, monotone_scheme,     &
                depart_high_order_scheme, depart_monotone_scheme,     &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_cartesian,                                          &
                l_shallow, l_rk_dps,                                  &
                l_high, l_mono,                                       &
                l_depart_high, l_depart_mono, lam_max_cfl,            &
                etadot, u_np1, v_np1, w_np1, etadot_np1,              &
                u, v, w, r_u, r_v, r_w, r_u_d, r_v_d, r_w_d,          &
                error_code )


USE parkind1,          ONLY : jpim, jprb       !DrHook
USE yomhook,           ONLY : lhook, dr_hook   !DrHook
USE proc_info_mod,     ONLY : mype=>me,                               &
                              nproc=>n_proc,                          &
                              nproc_x=>n_procx,                       &
                              nproc_y=>n_procy,                       &
                              global_row_length, global_rows,         &
                              at_extremity,gc_proc_row_group,         &
                              gc_proc_col_group,model_domain

USE timestep_mod,      ONLY : timestep
USE level_heights_mod, ONLY : eta_theta_levels, eta_rho_levels,       &
                              xi3_at_theta=>r_theta_levels,           &
                              xi3_at_rho=>r_rho_levels,               &
                              xi3_at_u=>r_at_u, xi3_at_v=>r_at_v,     &
                              xi3_at_u_w=>r_at_u_w,                   &
                              xi3_at_v_w=>r_at_v_w
USE eg_v_at_poles_mod
USE eg_sl_wind_u_mod
USE eg_sl_wind_v_mod
USE eg_sl_wind_w_mod
USE integrity_mod
USE atm_fields_bounds_mod

USE UM_ParParams
USE departure_pts_mod
USE eg_swap_bounds_mod

IMPLICIT NONE
!
! Description:
!   Find departure points on u,v,w-grid 
!   and apply rotation on timelevel n
!   quantities R_u, R_v, R_w.
!  
!
! Method: ENDGame formulation version 1.01,
!         section 7.2.
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Model dimensions

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,        &
                       first_constant_r_rho_level

! MPP options

INTEGER, INTENT(IN) ::                                                &
  halo_i,                                                             &
                     ! Size of halo in i.
  halo_j,                                                             &
                     ! Size of halo in j.
  offx,                                                               &
                     ! Size of small halo in i
  offy,                                                               &
                     ! Size of small halo in j.
  g_i_pe(1-halo_i:global_row_length+halo_i),                          &
                     ! processor on my processor-row
                     ! holding a given value in i direction
  lam_max_cfl(2)     ! Max CFL for a LAM allowed near the
                     ! boundaries

! Loop index bounds for arrays defined on p, u, v points respectively

! Integer parameters for advection

INTEGER, INTENT(IN) ::                                                &
  high_order_scheme,                                                  &
                     ! a code saying which high order scheme to
                     ! use. 1 = tensor tri-cubic lagrange order
                     ! (j,i,k) no other options available at
                     ! present.
  monotone_scheme,                                                    &
                     ! a code saying which monotone scheme to use.
                     ! 1 = tri-linear order (j,i,k)
                     ! no other options available at present.
  depart_scheme,                                                      &
                     ! code saying which departure point scheme to
                     ! use.
  depart_order,                                                       &
                     ! for the chosen departure point scheme how
                     ! many iterations/terms to use.
  depart_high_order_scheme,                                           &
                     ! code choosing high order
                     ! interpolation scheme used in
                     ! Departure_point routine
  depart_monotone_scheme,                                             &
                     ! code choosing monotone
                     ! interpolation scheme used in
                     ! Departure_point routine
  interp_vertical_search_tol,                                         &
                      ! used in interpolation code.
  check_bottom_levels ! used in interpolation code, and is
                      ! the number of levels to check to see
                      ! if the departure point lies inside the
                      ! orography.

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                &
  l_high,                                                             &
                   ! True, if high order interpolation required.
  l_mono,                                                             &
                   ! True, if interpolation required to be monotone.
  l_depart_high,                                                      &
                   ! True if high order interpolation scheme to be
                   ! used in Departure scheme
  l_depart_mono
                   ! True if monotone interpolation scheme to be
                   ! used in Departure scheme

LOGICAL, INTENT(IN) :: l_cartesian, l_shallow, l_rk_dps

REAL, INTENT(INOUT) ::                                                &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,                 &
         model_levels),                                               &
  v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,               &
        model_levels),                                                &
  w(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,                 &
         0:model_levels),                                             &
  u_np1(-offx:row_length-1+offx,1-offy:rows+offy, model_levels),      &
  v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy, model_levels),    &
  w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,                 &
             0:model_levels)

! Timelevel n arrival point quantities

REAL, INTENT(INOUT) ::                                                &
  r_u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),         &
  r_v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),       &
  r_w(row_length,rows,0:model_levels)


INTEGER :: error_code   ! Non-zero on exit if error detected.

! Timelevel n departure point quantities

REAL, INTENT(OUT) ::                                                  &
  r_u_d(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),       &
  r_v_d(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),     &
  r_w_d(row_length,rows,0:model_levels)

! Halo-ed copies of R_u, R_v, R_w for interpolation

REAL ::                                                               &
  rwork_u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,           &
           model_levels),                                             &
  rwork_v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,         &
           model_levels),                                             &
  rwork_w(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,            &
            0:model_levels)

! Local variables

INTEGER :: i,j,k

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SL_FULL_WIND',zhook_in,zhook_handle)

IF (integrity_test)                                                   &
  CALL  check_hash_m(R_u,             SIZE(R_u),             'R_u__', &
                     R_v,             SIZE(R_v),             'R_v__', &
                     R_w,             SIZE(R_w),             'R_w__')

IF( model_domain == mt_global ) THEN
   IF( at_extremity(psouth) ) THEN

      CALL eg_v_at_poles(r_u,r_v, 1.0, udims%j_start, vdims%j_start,&
                         udims_s,vdims_s)

   END IF
   IF( at_extremity(pnorth) ) THEN

      CALL eg_v_at_poles(r_u,r_v, -1.0, udims%j_end, vdims%j_end,&
                         udims_s,vdims_s)

   END IF
END IF

! set workspace

k = 0
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      rwork_w(i,j,k) = r_w(i,j,k)
    END DO
  END DO

!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE)                &
!$OMP& SHARED(model_levels,udims,vdims,pdims,r_u,r_v,r_w,     &
!$OMP& rwork_u,rwork_v,rwork_w) SCHEDULE(STATIC)
DO k=1, model_levels

  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end
      rwork_u(i,j,k) = r_u(i,j,k)
    END DO
  END DO

  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
      rwork_v(i,j,k) = r_v(i,j,k)
    END DO
  END DO

  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      rwork_w(i,j,k) = r_w(i,j,k)
    END DO
  END DO

END DO
!$OMP END PARALLEL DO

CALL eg_swap_bounds( rwork_u,udims_l,fld_type_u,.TRUE.)
CALL eg_swap_bounds( rwork_v,vdims_l,fld_type_v,.TRUE.)
CALL eg_swap_bounds( rwork_w,wdims_l,fld_type_p,.FALSE.)

CALL eg_sl_wind_w(                                                    &
                row_length, rows, n_rows, model_levels, g_i_pe,       &
                model_domain, depart_scheme, l_rk_dps,                &
                depart_order, high_order_scheme, monotone_scheme,     &
                depart_high_order_scheme, depart_monotone_scheme,     &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_high, l_mono,                                       &
                l_depart_high, l_depart_mono,                         &
                l_cartesian,l_shallow, lam_max_cfl,                   &
                etadot, u_np1, v_np1, w_np1, etadot_np1,              &
                u, v, w, rwork_u, rwork_v, rwork_w, r_w_d,            &
                error_code )

CALL eg_sl_wind_u(                                                    &
                row_length, rows, n_rows, model_levels,  g_i_pe,      &
                model_domain, depart_scheme, l_rk_dps,                &
                depart_order, high_order_scheme, monotone_scheme,     &
                depart_high_order_scheme,depart_monotone_scheme,      &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_high, l_mono,                                       &
                l_depart_high, l_depart_mono,                         &
                l_cartesian,l_shallow, lam_max_cfl,                   &
                etadot, u_np1, v_np1, w_np1, etadot_np1,              &
                u, v, w, rwork_u, rwork_v, rwork_w, r_u_d,            &
                error_code )

CALL eg_sl_wind_v(                                                    &
                row_length, rows, n_rows, model_levels,  g_i_pe,      &
                model_domain, depart_scheme, l_rk_dps,                &
                depart_order, high_order_scheme, monotone_scheme,     &
                depart_high_order_scheme,depart_monotone_scheme,      &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_high, l_mono,                                       &
                l_depart_high, l_depart_mono,                         &
                l_cartesian,l_shallow, lam_max_cfl,                   &
                etadot, u_np1, v_np1, w_np1, etadot_np1,              &
                u, v, w, rwork_u, rwork_v, rwork_w, r_v_d,error_code )



CALL eg_swap_bounds( r_u_d,udims_s,fld_type_u,.TRUE.)
CALL eg_swap_bounds( r_v_d,vdims_s,fld_type_v,.TRUE.)


IF (INTEGRITY_TEST)                                                   &
  CALL update_hash_m(R_u_d,           SIZE(R_u_d),           'R_u_d', &
                     R_v_d,           SIZE(R_v_d),           'R_v_d', &
                     R_w_d,           SIZE(R_w_d),           'R_w_d', &
                     R_u,             SIZE(R_u),             'R_u__', &
                     R_v,             SIZE(R_v),             'R_v__', &
                     R_w,             SIZE(R_w),             'R_w__', &
                     depart_xi1_u,    SIZE(depart_xi1_u),    'dxi1u', &
                     depart_xi2_u,    SIZE(depart_xi2_u),    'dxi2u', &
                     depart_xi3_u,    SIZE(depart_xi3_u),    'dxi3u', &
                     depart_xi1_v,    SIZE(depart_xi1_v),    'dxi1v', &
                     depart_xi2_v,    SIZE(depart_xi2_v),    'dxi2v', &
                     depart_xi3_v,    SIZE(depart_xi3_v),    'dxi3v', &
                     depart_xi1_w,    SIZE(depart_xi1_w),    'dxi1w', &
                     depart_xi2_w,    SIZE(depart_xi2_w),    'dxi2w', &
                     depart_xi3_w,    SIZE(depart_xi3_w),    'dxi3w')

IF (lhook) CALL dr_hook('EG_SL_FULL_WIND',zhook_out,zhook_handle)

END SUBROUTINE eg_sl_full_wind
END MODULE eg_sl_full_wind_mod
