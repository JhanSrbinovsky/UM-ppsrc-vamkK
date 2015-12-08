! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_wind_v_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_sl_wind_v(                                              &
                row_length, rows, n_rows, model_levels,g_i_pe,        &
                model_domain,depart_scheme, l_rk_dps,                 &
                depart_order, high_order_scheme, monotone_scheme,     &
                depart_high_order_scheme, depart_monotone_scheme,     &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_high, l_mono, l_depart_high, l_depart_mono,         &
                l_cartesian, l_shallow, lam_max_cfl,                  &
                etadot, u_np1, v_np1, w_np1, etadot_np1,              &
                u, v, w, rwork_u, rwork_v, rwork_w,                   &
                r_v_d, error_code )

USE timestep_mod,      ONLY : timestep
USE um_parvars,        ONLY : offx, offy, halo_i, halo_j
USE proc_info_mod,     ONLY : datastart=>l_datastart, mype=>me,       &
                              nproc=>n_proc, nproc_x=>n_procx,        &
                              nproc_y=>n_procy,                       &
                              global_row_length, global_rows,         &
                              at_extremity,gc_proc_col_group,         &
                              gc_proc_row_group

USE level_heights_mod, ONLY : eta_theta_levels, eta_rho_levels,       &
                              xi3_at_theta=>r_theta_levels,           &
                              xi3_at_rho=>r_rho_levels,               &
                              xi3_at_u=>r_at_u,                       &
                              xi3_at_v=>r_at_v

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE eg_interpolation_eta_mod
USE departure_point_eta_mod
USE eg_dep_pnt_cart_eta_mod
USE horiz_grid_mod
USE atm_fields_bounds_mod

USE ereport_mod, ONLY : ereport
USE Field_Types

USE departure_pts_mod
USE eg_parameters_mod, ONLY : interp_dpt_pt
USE conversions_mod,   ONLY : pi

IMPLICIT NONE
!
! Description: Find v-grid departure point and apply rotation on R_v.
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

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,    &
                       model_domain, first_constant_r_rho_level

! MPP options

INTEGER, INTENT(IN) ::                                             &
  g_i_pe(1-halo_i:global_row_length+halo_i),                       &
                     ! processor on my processor-row
                     ! holding a given value in i direction
  lam_max_cfl(2)     ! Max CFL for a LAM allowed near the
                     ! boundaries

! Integer parameters for advection

INTEGER, INTENT(IN) ::                                            &
  high_order_scheme,                                              &
                     ! a code saying which high order scheme to
                     ! use. 1 = tensor tri-cubic lagrange order
                     ! (j,i,k) no other options available at
                     ! present.
  monotone_scheme,                                                &
                     ! a code saying which monotone scheme to use.
                     ! 1 = tri-linear order (j,i,k)
                     ! no other options available at present.
  depart_scheme,                                                  &
                     ! code saying which departure point scheme to
                     ! use.
  depart_order,                                                   &
                     ! for the chosen departure point scheme how
                     ! many iterations/terms to use.
  depart_high_order_scheme,                                       &
                     ! code choosing high order
                     ! interpolation scheme used in Depart routine
  depart_monotone_scheme,                                         &
                     ! code choosing monotone
                     ! interpolation scheme used in Depart routine
  interp_vertical_search_tol,                                     &
                      ! used in interpolation code.
  check_bottom_levels ! used in interpolation code, and is
                      ! the number of levels to check to see
                      ! if the departure point lies inside the
                      ! orography.

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                            &
  l_high,                                                         &
                   ! True, if high order interpolation required.
  l_mono,                                                         &
                   ! True, if interpolation required to be monotone.
  l_depart_high,                                                  &
                   ! True if high order interpolation scheme to be
                   ! used in Depart scheme
  l_depart_mono
                   ! True if monotone interpolation scheme to be
                   ! used in Depart scheme


LOGICAL, INTENT(IN) :: l_cartesian, l_shallow, l_rk_dps

REAL, INTENT(IN) ::                                                &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),  &
  u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,              &
         model_levels),                                            &
  v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,            &
        model_levels),                                             &
  w(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,              &
         0:model_levels),                                          &
  u_np1(-offx:row_length-1+offx,1-offy:rows+offy, model_levels),   &
  v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy, model_levels), &
  w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
  etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,              &
             0:model_levels)

! Halo-ed copies of R_u, R_v, R_w for interpolation

REAL, INTENT(IN) ::                                            &
  rwork_u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,    &
           model_levels),                                      &
  rwork_v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,  &
           model_levels),                                      &
  rwork_w(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,     &
            0:model_levels)

INTEGER :: error_code   ! Non-zero on exit if error detected.

! Timelevel n (rotated) v-departure point quantity

REAL, INTENT(OUT) ::                          &
  r_v_d(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels)

! Local variables

! row_length, number of rows and indexing offset for
! the dep point type computed

INTEGER :: dep_row_len, dep_rows, off_i, off_j, off_k, offz,     &
           number_of_inputs

INTEGER :: i, j, k

REAL    :: rm_21, rm_22, rm_23, temp
! Temporary local
REAL :: temp_csalpha


! Timelevel n R quantities interpolated on a v-point

REAL ::                                        &
  ru_vd(row_length,0:n_rows-1,model_levels),   &
  rv_vd(row_length,0:n_rows-1,model_levels),   &
  rw_vd(row_length,0:n_rows-1,model_levels)


REAL :: alpha

! End of header

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SL_WIND_V',zhook_in,zhook_handle)

dep_row_len = vdims%i_end - vdims%i_start + 1
dep_rows    = vdims%j_end - vdims%j_start + 1

off_i = 0
off_j = 1
off_k = 0
offz  = 0

number_of_inputs = 1
!      alpha = alpha_v
alpha = 0.5

IF ( .NOT. interp_dpt_pt ) THEN
  IF( model_domain /= mt_global ) THEN
      CALL eg_dep_pnt_cart_eta(                                         &
                 row_length, rows, n_rows, model_levels, halo_i,        &
                 halo_j, offx, offy, mype, nproc, nproc_x, nproc_y,     &
                 global_row_length, global_rows, datastart,             &
                 at_extremity, g_i_pe, gc_proc_row_group,               &
                 gc_proc_col_group, model_domain,                       &
                 fld_type_v, dep_row_len, dep_rows, off_i, off_j, off_k,&
                 offz,depart_scheme, depart_order, l_rk_dps,            &
                 depart_high_order_scheme,                              &
                 depart_monotone_scheme, first_constant_r_rho_level,    &
                 interp_vertical_search_tol, check_bottom_levels,       &
                 l_depart_high,                                         &
                 l_depart_mono, lam_max_cfl, alpha, timestep,           &
                 eta_theta_levels, eta_rho_levels,                      &
                 xi3_at_theta, xi3_at_rho, xi3_at_u, xi3_at_v,          &
                 etadot, u_np1, v_np1, w_np1,                           &
                 etadot_np1, u, v, w, depart_xi1_v, depart_xi2_v,       &
                 depart_xi3_v )
  
  ELSE
  
      CALL departure_point_eta(                                         &
                 row_length, rows, n_rows, model_levels, g_i_pe,        &
                 fld_type_v, dep_row_len, dep_rows,                     &
                 off_i, off_j, off_k, offz,l_rk_dps, depart_order,      &
                 depart_high_order_scheme, depart_monotone_scheme,      &
                 first_constant_r_rho_level, l_depart_high,             &
                 l_depart_mono, l_shallow, lam_max_cfl,                 &
                 alpha, eta_theta_levels, eta_rho_levels,               &
                 xi3_at_theta, xi3_at_rho, xi3_at_u, xi3_at_v,          &
                 etadot, u_np1, v_np1,etadot_np1, u, v, w, depart_xi1_v,&
                 depart_xi2_v, depart_xi3_v )
  END IF
END IF

IF( .NOT. l_cartesian ) THEN
  CALL eg_interpolation_eta(                                          &
                     eta_rho_levels,                                  &
                     fld_type_u, number_of_inputs,                    &
                     row_length, rows, model_levels,                  &
                     rows,                                            &
                     row_length, n_rows, model_levels,                &
                     high_order_scheme, monotone_scheme,              &
                     model_domain, l_high, l_mono,                    &
                     depart_xi3_v, depart_xi1_v,                      &
                     depart_xi2_v, mype, nproc, nproc_x, nproc_y,     &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     0, 0, offx, offy,                                &
                     error_code,                                      &
                     rwork_u, ru_vd)
ENDIF

  CALL eg_interpolation_eta(                                          &
                     eta_rho_levels,fld_type_v, number_of_inputs,     &
                     row_length, n_rows, model_levels,                &
                     rows,row_length, n_rows, model_levels,           &
                     high_order_scheme, monotone_scheme,              &
                     model_domain, l_high, l_mono,                    &
                     depart_xi3_v, depart_xi1_v,                      &
                     depart_xi2_v, mype, nproc, nproc_x, nproc_y,     &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     0, 0, offx, offy,error_code,                     &
                     rwork_v, rv_vd)

IF( .NOT. (l_cartesian .OR. l_shallow) ) THEN
  CALL eg_interpolation_eta(                                          &
                     eta_theta_levels,                                &
                     fld_type_w, number_of_inputs,                    &
                     row_length, rows, model_levels+1,                &
                     rows,                                            &
                     row_length, n_rows, model_levels,                &
                     high_order_scheme, monotone_scheme,              &
                     model_domain, l_high, l_mono,                    &
                     depart_xi3_v, depart_xi1_v,                      &
                     depart_xi2_v, mype, nproc, nproc_x, nproc_y,     &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     0, 0, offx, offy,error_code,                     &
                     rwork_w, rw_vd)
ENDIF

IF( .NOT. l_cartesian ) THEN
  IF (l_shallow) THEN
!$OMP PARALLEL DO PRIVATE(i,j,k,temp,temp_csalpha,rm_21,rm_22)        &
!$OMP& DEFAULT(NONE) SHARED( model_levels,vdims,xi1_p,depart_xi1_v,   &
!$OMP& xi2_v,depart_xi2_v,r_v_d,ru_vd,rv_vd) SCHEDULE(STATIC)
   DO k=1, model_levels
      DO j=vdims%j_start, vdims%j_end
         DO i=vdims%i_start, vdims%i_end

! Coded for spericals only at the moment

            temp  = COS(xi1_p(i) - depart_xi1_v(i,j,k))
            temp_csalpha = SIN(xi2_v(j))*SIN(depart_xi2_v(i,j,k))     &
                           + COS(xi2_v(j))*COS(depart_xi2_v(i,j,k)) * &
                           temp

            rm_21 = -SIN(xi1_p(i) - depart_xi1_v(i,j,k)) *            &
                    (SIN(xi2_v(j)) + SIN(depart_xi2_v(i,j,k))) /      &
                    (1.0 + temp_csalpha)
            rm_22 = (COS(xi2_v(j)) * COS(depart_xi2_v(i,j,k)) +       &
                      temp * (1.0 + SIN(xi2_v(j)) *                   &
                      SIN(depart_xi2_v(i,j,k)))) /                    &
                      (1.0 + temp_csalpha)

            r_v_d(i,j,k) = rm_21*ru_vd(i,j,k) + rm_22*rv_vd(i,j,k)
         END DO
      END DO
   END DO
!$OMP END PARALLEL DO

 ELSE ! If not L_shallow

!$OMP PARALLEL DO PRIVATE(i,j,k,temp,rm_21,rm_22,rm_23) DEFAULT(NONE)   &
!$OMP& SHARED(model_levels,vdims,xi2_v,xi1_p,depart_xi1_v,depart_xi2_v, &
!$OMP& r_v_d,ru_vd,rv_vd,rw_vd) SCHEDULE(STATIC)
   DO k=1, model_levels
      DO j=vdims%j_start, vdims%j_end
         DO i=vdims%i_start, vdims%i_end

! Coded for spericals only at the moment

            temp  = SIN(xi2_v(j))                                     &
                   *COS(xi1_p(i)-depart_xi1_v(i,j,k))

            rm_21 =-SIN(xi2_v(j))                                     &
                   *SIN(xi1_p(i)-depart_xi1_v(i,j,k))
            rm_22 = COS(xi2_v(j))*COS(depart_xi2_v(i,j,k))            &
                   +SIN(depart_xi2_v(i,j,k))*temp
            rm_23 = COS(xi2_v(j))*SIN(depart_xi2_v(i,j,k))            &
                   -COS(depart_xi2_v(i,j,k))*temp

            r_v_d(i,j,k) = rm_21*ru_vd(i,j,k) +                       &
                           rm_22*rv_vd(i,j,k) +                       &
                           rm_23*rw_vd(i,j,k)
         END DO
      END DO
   END DO 
!$OMP END PARALLEL DO

   END IF  !End of if shallow

ELSE ! If cartesian then

!$OMP PARALLEL DO PRIVATE(i,j,k) DEFAULT(NONE) SHARED(model_levels,   &
!$OMP& vdims,r_v_d,rv_vd) SCHEDULE(STATIC)
   DO k=1, model_levels
      DO j=vdims%j_start, vdims%j_end
         DO i=vdims%i_start, vdims%i_end

            r_v_d(i,j,k) = rv_vd(i,j,k)

         END DO
      END DO
   END DO
!$OMP END PARALLEL DO
END IF

IF (lhook) CALL dr_hook('EG_SL_WIND_V',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_sl_wind_v
END MODULE eg_sl_wind_v_mod
