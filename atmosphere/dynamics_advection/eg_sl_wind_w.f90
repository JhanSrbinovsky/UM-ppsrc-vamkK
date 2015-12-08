! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sl_wind_w_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_sl_wind_w(                                              &
                row_length, rows, n_rows, model_levels, g_i_pe,       &
                model_domain, depart_scheme, l_rk_dps,                &
                depart_order, high_order_scheme, monotone_scheme,     &
                depart_high_order_scheme, depart_monotone_scheme,     &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_high, l_mono,                                       &
                l_depart_high, l_depart_mono,                         &
                l_cartesian, l_shallow, lam_max_cfl,                  &
                etadot, u_np1, v_np1, w_np1, etadot_np1,              &
                u, v, w, rwork_u, rwork_v, rwork_w, r_w_d,error_code )



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
                              xi3_at_v=>r_at_v,                       &
                              xi3_at_u_w=>r_at_u_w,                   &
                              xi3_at_v_w=>r_at_v_w

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


IMPLICIT NONE
!
! Description: Find w-grid departure point and apply rotation on R_w.
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
                       model_domain, first_constant_r_rho_level

! MPP options
INTEGER, INTENT(IN) ::                                                &
  g_i_pe(1-halo_i:global_row_length+halo_i),                          &
                     ! processor on my processor-row
                     ! holding a given value in i direction
  lam_max_cfl(2)     ! Max CFL for a LAM allowed near the
                     ! boundaries


! Integer parameters for advection

INTEGER, INTENT(IN) ::                                               &
  high_order_scheme,                                                 &
                     ! a code saying which high order scheme to
                     ! use. 1 = tensor tri-cubic lagrange order
                     ! (j,i,k) no other options available at
                     ! present.
  monotone_scheme,                                                   &
                     ! a code saying which monotone scheme to use.
                     ! 1 = tri-linear order (j,i,k)
                     ! no other options available at present.
  depart_scheme,                                                     &
                     ! code saying which departure point scheme to
                     ! use.
  depart_order,                                                      &
                     ! for the chosen departure point scheme how
                     ! many iterations/terms to use.
  depart_high_order_scheme,                                          &
                     ! code choosing high order
                     ! interpolation scheme used in Depart routine
  depart_monotone_scheme,                                            &
                     ! code choosing monotone
                     ! interpolation scheme used in Depart routine
  interp_vertical_search_tol,                                        &
                      ! used in interpolation code.
  check_bottom_levels ! used in interpolation code, and is
                      ! the number of levels to check to see
                      ! if the departure point lies inside the
                      ! orography.

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                               &
  l_high,                                                            &
                   ! True, if high order interpolation required.
  l_mono,                                                            &
                   ! True, if interpolation required to be monotone.
  l_depart_high,                                                     &
                   ! True if high order interpolation scheme to
                   ! be used in Departure scheme
  l_depart_mono
                   ! True if monotone interpolation scheme to
                   ! be used in Departure scheme


LOGICAL, INTENT(IN) :: l_cartesian
LOGICAL, INTENT(IN) :: l_shallow, l_rk_dps
REAL, INTENT(INOUT)  ::                                               &
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


! Halo-ed copies of R_u, R_v, R_w for interpolation

REAL, INTENT(IN) ::                                                   &
  rwork_u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,           &
           model_levels),                                             &
  rwork_v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,         &
           model_levels),                                             &
  rwork_w(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,            &
            0:model_levels)

INTEGER :: error_code   ! Non-zero on exit if error detected.

! Timelevel n (rotated) w-departure point quantity

REAL, INTENT(OUT) :: r_w_d(row_length,rows,0:model_levels)

! Local variables

! row_length, number of rows and indexing offset for
! the dep point type computed

INTEGER :: dep_row_len, dep_rows, off_i, off_j, off_k, offz,          &
           number_of_inputs

INTEGER :: i, j, k

REAL    :: rm_31, rm_32, rm_33, temp

! Timelevel n R quantities interpolated on a w-point

REAL ::                                                               &
  ru_wd(row_length,rows,0:model_levels),                              &
  rv_wd(row_length,rows,0:model_levels),                              &
  rw_wd(row_length,rows,0:model_levels),                              &
  work_u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,            &
         0:model_levels+1),                                           &
  work_v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,          &
         0:model_levels+1),                                           &
  xi3_rho_ext(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,        &
              0:model_levels+1),                                      &
              ! vertical coordinate at rho points
  xi3_u_ext(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,         &
            0:model_levels+1),                                        &
              ! vertical coordinate at u points
  xi3_v_ext(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,       &
            0:model_levels+1),                                        &
  eta_rho_ext(0:model_levels+1)


REAL :: alpha


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SL_WIND_W',zhook_in,zhook_handle)


dep_row_len = pdims%i_end - pdims%i_start + 1
dep_rows    = pdims%j_end - pdims%j_start + 1
off_i = 0
off_j = 0
off_k = 1
offz  = 1

!      alpha = alpha_w
alpha = 0.5

! Copy u,v components into work_u, work_v and extend them
! below/above bottom/top boundary.

! Same for height fields
! Loop over extended halos to avoid swap_bounds.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SCHEDULE(STATIC)        &
!$OMP& SHARED(model_levels,udims,halo_j,halo_i,work_u,u,               &
!$OMP& xi3_u_ext,xi3_at_u,work_v,v,vdims,eta_rho_ext,                  &
!$OMP& eta_rho_levels,pdims,xi3_rho_ext,xi3_at_rho,xi3_v_ext,xi3_at_v)
DO k = 1, model_levels
  DO j = udims%j_start-halo_j, udims%j_end+halo_j
    DO i = udims%i_start-halo_i, udims%i_end+halo_i
      work_u(i,j,k) = u(i,j,k)
    END DO
  END DO

  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      xi3_u_ext(i,j,k) = xi3_at_u(i,j,k)
    END DO
  END DO

  DO j = vdims%j_start-halo_j, vdims%j_end+halo_j
    DO i = vdims%i_start-halo_i, vdims%i_end+halo_i
      work_v(i,j,k) = v(i,j,k)
    END DO
  END DO

  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      xi3_v_ext(i,j,k) = xi3_at_v(i,j,k)
    END DO
  END DO


  eta_rho_ext(k) = eta_rho_levels(k)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      xi3_rho_ext(i,j,k) = xi3_at_rho(i,j,k)
    END DO
  END DO

END DO
!$OMP END PARALLEL DO

DO j = udims%j_start-halo_j, udims%j_end+halo_j
  DO i = udims%i_start-halo_i, udims%i_end+halo_i
    work_u(i,j,0) = u(i,j,1)
    work_u(i,j,model_levels+1) = u(i,j,model_levels)
  END DO
END DO

DO j = vdims%j_start-halo_j, vdims%j_end+halo_j
  DO i = vdims%i_start-halo_i, vdims%i_end+halo_i
    work_v(i,j,0) = v(i,j,1)
    work_v(i,j,model_levels+1) = v(i,j,model_levels)
  END DO
END DO

DO j = udims%j_start, udims%j_end
  DO i = udims%i_start, udims%i_end
    xi3_u_ext(i,j,0) = xi3_at_u_w(i,j,0)
    xi3_u_ext(i,j,model_levels+1) =                                   &
              intw_p2u(i,1)*xi3_at_theta(i,j,model_levels) +          &
              intw_p2u(i,2)*xi3_at_theta(i+1,j,model_levels)
  END DO
END DO

DO j = vdims%j_start, vdims%j_end
  DO i = vdims%i_start, vdims%i_end
    xi3_v_ext(i,j,0) = xi3_at_v_w(i,j,0)
    xi3_v_ext(i,j,model_levels+1) =                                   &
                 intw_p2v(j,1)*xi3_at_theta(i,j,model_levels) +       &
                 intw_p2v(j,2)*xi3_at_theta(i,j+1,model_levels)
  END DO
END DO


eta_rho_ext(0)              = eta_theta_levels(0)
eta_rho_ext(model_levels+1) = eta_theta_levels(model_levels)

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    xi3_rho_ext(i,j,0) = xi3_at_theta(i,j,0)
    xi3_rho_ext(i,j,model_levels+1) = xi3_at_theta(i,j,model_levels)
  END DO
END DO

number_of_inputs = 1

IF( model_domain /= mt_global ) THEN
    CALL eg_dep_pnt_cart_eta(                                         &
               row_length, rows, n_rows, model_levels, halo_i,        &
               halo_j, offx, offy, mype, nproc, nproc_x, nproc_y,     &
               global_row_length, global_rows, datastart,             &
               at_extremity, g_i_pe, gc_proc_row_group,               &
               gc_proc_col_group, model_domain,                       &
               fld_type_w, dep_row_len, dep_rows, off_i ,off_j, off_k,&
               offz,depart_scheme, depart_order, l_rk_dps,            &
               depart_high_order_scheme,                              &
               depart_monotone_scheme, first_constant_r_rho_level,    &
               interp_vertical_search_tol, check_bottom_levels,       &
               l_depart_high,                                         &
               l_depart_mono, lam_max_cfl, alpha, timestep,           &
               eta_theta_levels, eta_rho_ext,                         &
               xi3_at_theta, xi3_rho_ext, xi3_u_ext, xi3_v_ext,       &
               etadot, u_np1, v_np1, w_np1,                           &
               etadot_np1, work_u, work_v, w,                         &
               depart_xi1_w, depart_xi2_w, depart_xi3_w )
ELSE

    CALL departure_point_eta(                                         &
               row_length, rows, n_rows, model_levels, g_i_pe,        &
               fld_type_w, dep_row_len, dep_rows,                     &
               off_i ,off_j, off_k,offz, l_rk_dps,depart_order,       &
               depart_high_order_scheme, depart_monotone_scheme,      &
               first_constant_r_rho_level, l_depart_high,             &
               l_depart_mono, l_shallow, lam_max_cfl, alpha,          &
               eta_theta_levels, eta_rho_ext, xi3_at_theta,           &
               xi3_rho_ext, xi3_u_ext, xi3_v_ext,                     &
               etadot, u_np1, v_np1, etadot_np1,                      &
               work_u, work_v, w, depart_xi1_w, depart_xi2_w,         &
               depart_xi3_w )
END IF

! Copy Rwork_u, v components into work_u, work_v and extend them
! below/above bottom/top boundary.

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(model_levels,   &
!$OMP& udims,halo_i,halo_j,work_u,rwork_u,vdims,work_v,rwork_v)       &
!$OMP& SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = udims%j_start-halo_j, udims%j_end+halo_j
    DO i = udims%i_start-halo_i, udims%i_end+halo_i
      work_u(i,j,k) = rwork_u(i,j,k)
    END DO
  END DO
!END DO

!DO k=1, model_levels
  DO j = vdims%j_start-halo_j, vdims%j_end+halo_j
    DO i = vdims%i_start-halo_i, vdims%i_end+halo_i
      work_v(i,j,k) = rwork_v(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

DO j = udims%j_start-halo_j, udims%j_end+halo_j
  DO i = udims%i_start-halo_i, udims%i_end+halo_i
    work_u(i,j,0) = rwork_u(i,j,1)
    work_u(i,j,model_levels+1) = rwork_u(i,j,model_levels)
  END DO
END DO

DO j = vdims%j_start-halo_j, vdims%j_end+halo_j
  DO i = vdims%i_start-halo_i, vdims%i_end+halo_i
    work_v(i,j,0) = rwork_v(i,j,1)
    work_v(i,j,model_levels+1) = rwork_v(i,j,model_levels)
  END DO
END DO

IF( .NOT. (l_shallow .OR. l_cartesian) ) THEN
CALL eg_interpolation_eta(                                            &
                     eta_rho_ext,                                     &
                     fld_type_u, number_of_inputs,                    &
                     row_length, rows, model_levels+2,                &
                     rows,                                            &
                     row_length, rows, model_levels+1,                &
                     high_order_scheme, monotone_scheme,              &
                     model_domain, l_high, l_mono,                    &
                     depart_xi3_w, depart_xi1_w, depart_xi2_w,        &
                     mype, nproc, nproc_x, nproc_y,                   &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     0, 0, offx, offy,error_code,                     &
                     work_u, ru_wd)

CALL eg_interpolation_eta(                                            &
                     eta_rho_ext,fld_type_v, number_of_inputs,        &
                     row_length, n_rows, model_levels+2,              &
                     rows,                                            &
                     row_length, rows, model_levels+1,                &
                     high_order_scheme, monotone_scheme,              &
                     model_domain, l_high, l_mono,                    &
                     depart_xi3_w, depart_xi1_w,                      &
                     depart_xi2_w, mype, nproc, nproc_x, nproc_y,     &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     0, 0, offx, offy,error_code,                     &
                     work_v, rv_wd)
ENDIF

CALL eg_interpolation_eta(                                            &
                     eta_theta_levels,                                &
                     fld_type_w,                                      &
                     number_of_inputs,                                &
                     row_length, rows, model_levels+1,                &
                     rows,                                            &
                     row_length, rows, model_levels+1,                &
                     high_order_scheme, monotone_scheme,              &
                     model_domain, l_high, l_mono,                    &
                     depart_xi3_w, depart_xi1_w,                      &
                     depart_xi2_w, mype, nproc, nproc_x, nproc_y,     &
                     halo_i, halo_j,                                  &
                     global_row_length, datastart, at_extremity,      &
                     g_i_pe, gc_proc_row_group, gc_proc_col_group,    &
                     0, 0, offx, offy,error_code,                     &
                     rwork_w, rw_wd)


! If spherical geometry, Rotate Ru_wd, Rv_wd, Rw_wd
IF( .NOT. l_cartesian ) THEN

  IF (l_shallow) THEN

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(model_levels,   &
!$OMP& pdims,r_w_d,rw_wd) SCHEDULE(STATIC)
    DO k = 0, model_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          r_w_d(i,j,k) = rw_wd(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  ELSE   !If deep then:

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,temp,rm_31,rm_32,rm_33)   &
!$OMP& SHARED(model_levels,pdims,xi2_p,xi1_p,depart_xi1_w,depart_xi2_w, &
!$OMP& r_w_d,ru_wd,rv_wd,rw_wd) SCHEDULE(STATIC)
    DO k=0, model_levels
      DO j=pdims%j_start, pdims%j_end
        DO i=pdims%i_start, pdims%i_end

! Coded only for spherical at the moment

          temp  = COS(xi2_p(j))                                       &
                 *COS(xi1_p(i)-depart_xi1_w(i,j,k))

          rm_31 = COS(xi2_p(j))                                       &
                 *SIN(xi1_p(i)-depart_xi1_w(i,j,k))
          rm_32 = SIN(xi2_p(j))                                       &
                 *COS(depart_xi2_w(i,j,k))                            &
                 -SIN(depart_xi2_w(i,j,k))*temp
          rm_33 = SIN(xi2_p(j))                                       &
                 *SIN(depart_xi2_w(i,j,k))                            &
                 +COS(depart_xi2_w(i,j,k))*temp

          r_w_d(i,j,k) = rm_31*ru_wd(i,j,k) +                         &
                         rm_32*rv_wd(i,j,k) +                         &
                         rm_33*rw_wd(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END IF !End of If shallow

ELSE !Else, if cartesian then:

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(model_levels,   &
!$OMP& pdims,r_w_d,rw_wd) SCHEDULE(STATIC)
  DO k=0, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        r_w_d(i,j,k) = rw_wd(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

IF (lhook) CALL dr_hook('EG_SL_WIND_W',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_sl_wind_w
END MODULE eg_sl_wind_w_mod
