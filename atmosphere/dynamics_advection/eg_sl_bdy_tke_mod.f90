! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SL advection of boundary layer tke fields
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Dynamics Advection

MODULE eg_sl_bdy_tke_mod
IMPLICIT NONE 
 
! Description:  
!  
! Code description:
! Language: Fortran 95.  
! This code is written to UMDP3 standards.  
CONTAINS
SUBROUTINE eg_sl_bdy_tke(
                row_length, rows, n_rows, model_levels,               &
                halo_i, halo_j, offx, offy,datastart, g_i_pe,         &
                eg_version_num, outerno,                              &
                depart_scheme, depart_order, high_order_scheme,       &
                monotone_scheme, ritchie_high_order_scheme,           &
                ritchie_monotone_scheme,                              &
                first_constant_r_rho_level,                           &
                interp_vertical_search_tol, check_bottom_levels,      &
                l_code_test, l_cartesian, l_sl_halo_reprod,           &
                l_free_slip, l_high, l_mono,                          &
                l_ritchie_high, l_ritchie_mono,                       &
!               flag to switch physics fields:
                l_pc2,                                                &
!
                alpha_theta,                                          &
                deta_xi3_theta, delta_xi1, delta_xi2,                 &
                base_xi1, base_xi2, intw_u2p, intw_v2p,               &
                intw_p2u, intw_p2v, intw_rho2w, intw_w2rho,           &
! reciprocals needed by interpolation
! Reciprocals needed by interpolation subroutine
                          recip_lambda_p_m, recip_lambda_p_0,          &
                          recip_lambda_p_p, recip_lambda_p_p2,         &
                          recip_lambda_u_m, recip_lambda_u_0,          &
                          recip_lambda_u_p, recip_lambda_u_p2,         &
                          recip_phi_p_m, recip_phi_p_0,                &
                          recip_phi_p_p, recip_phi_p_p2,               &
                          recip_phi_v_m, recip_phi_v_0,                &
                          recip_phi_v_p, recip_phi_v_p2,               &
                cos_theta_latitude, sec_theta_latitude,               &
                sin_theta_latitude, cos_v_latitude,                   &
                sec_v_latitude, sin_v_latitude,                       &
                tan_theta_latitude, tan_v_latitude,                   &
                cos_theta_longitude, sin_theta_longitude,             &
                h1_p, h1_xi1_u, h1_xi2_v, h1_p_eta, h2_p,             &
                h2_xi1_u, h2_xi2_v, h2_p_eta, h3_p, h3_xi1_u,         &
                h3_xi2_v, h3_p_eta, u_np1, v_np1, w_np1,              &
                e_trb, tsq, qsq, cov,                                 &
                depart_xi1_w, depart_xi2_w, depart_xi3_w, etadot,     &
                error_code )


USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

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
                              xi3_at_u=>r_at_u, xi3_at_v=>r_at_v

USE mym_option_mod, ONLY:                                             &
                            bdy_tke, l_adv_turb_field, tke_levels

USE atm_fields_bounds_mod
USE integrity_mod
USE horiz_grid_mod
USE ref_pro_mod

USE Field_Types

IMPLICIT NONE

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


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
  offy
                     ! Size of small halo in j.



INTEGER, INTENT(IN) ::                                                &
  datastart(3),                                                       &
                     ! First gridpoints held by this processor.
  g_i_pe(1-halo_i:global_row_length+halo_i)
                     ! processor on my processor-row
                     ! holding a given value in i direction


! Loop index bounds for arrays defined on p, u, v points respectively

INTEGER, INTENT(IN) :: pdims%i_start, pdims%i_end, pdims%j_start, pdims%j_end,        &
                       udims%i_start, udims%i_end, udims%j_start, udims%j_end,        &
                       vdims%i_start, vdims%i_end, vdims%j_start, vdims%j_end

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
  ritchie_high_order_scheme,                                          &
                     ! code choosing high order
                     ! interpolation scheme used in Ritchie routine
  ritchie_monotone_scheme,                                            &
                     ! code choosing monotone
                     ! interpolation scheme used in Ritchie routine
  interp_vertical_search_tol,                                         &
                      ! used in interpolation code.
  check_bottom_levels ! used in interpolation code, and is
                      ! the number of levels to check to see
                      ! if the departure point lies inside the
                      ! orography.

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                &
  l_free_slip,                                                        &
                   ! True, free-slip lower boundary condition
  l_high,                                                             &
                   ! True, if high order interpolation required.
  l_mono,                                                             &
                   ! True, if interpolation required to be monotone.
  l_ritchie_high,                                                     &
                   ! True if high order interpolation scheme to be
                   ! used in Ritchie scheme
  l_ritchie_mono,                                                     &
                   ! True if monotone interpolation scheme to be
                   ! used in Ritchie scheme
  l_sl_halo_reprod,                                                   &
                   ! if true then sl code bit reproducible with
                   ! any sensible halo size
  l_pc2

LOGICAL, INTENT(IN) :: l_code_test, l_cartesian

INTEGER, INTENT(IN) :: outerno,  eg_version_num
                                           ! Outer loop iteration number
                                           ! ENDGAME formulation version

REAL, INTENT(IN) :: alpha_theta  ! time-weight for theta eqn
                                 ! and model timestep



REAL, INTENT(IN) ::                                                   &
  deta_xi3_theta(1-offx:row_length+offx,1-offy:rows+offy,             &
                 0:model_levels)

! Horizontal grid coordinates & interpolation weights
REAL, INTENT(IN) ::                                                   &
  intw_u2p(1:row_length,2),                                           &
  intw_v2p(1:rows,2),                                                 &
  intw_p2u(0:row_length,2),                                           &
  intw_p2v(0:n_rows,2),                                               &
  intw_rho2w(model_levels,2),                                         &
  intw_w2rho(model_levels,2)

! Grid base and spacing
REAL :: delta_xi1, delta_xi2, base_xi1, base_xi2

! Trigonometric functions
REAL, INTENT(IN) ::                                                   &
  cos_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy),        &
  sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy),        &
  sin_theta_latitude(row_length, rows),                               &
  tan_theta_latitude(row_length, rows),                               &
  cos_v_latitude(1-offx:row_length+offx, -offy:n_rows-1+offy),        &
  sec_v_latitude(1-offx:row_length+offx, -offy:n_rows-1+offy),        &
  sin_v_latitude(row_length, 0:n_rows-1),                             &
  tan_v_latitude(row_length, 0:n_rows-1),                             &
  cos_theta_longitude(row_length, rows),                              &
  sin_theta_longitude(row_length, rows)


! Metric terms

REAL, INTENT(IN) ::                                                   &
  h1_p(1-offx:row_length+offx,1-offy:rows+offy,model_levels),         &
  h1_xi1_u(-offx:row_length-1+offx,1-offy:rows+offy,                  &
           model_levels),                                             &
  h1_xi2_v(1-offx:row_length+offx,-offy:n_rows-1+offy,                &
           model_levels),                                             &
  h1_p_eta(1-offx:row_length+offx,1-offy:rows+offy,                   &
           0:model_levels),                                           &
  h2_p(1-offx:row_length+offx,1-offy:rows+offy,model_levels),         &
  h2_xi1_u(-offx:row_length-1+offx,1-offy:rows+offy,                  &
           model_levels),                                             &
  h2_xi2_v(1-offx:row_length+offx,-offy:n_rows-1+offy,                &
           model_levels),                                             &
  h2_p_eta(1-offx:row_length+offx,1-offy:rows+offy,                   &
           0:model_levels),                                           &
  h3_p(1-offx:row_length+offx,1-offy:rows+offy,model_levels),         &
  h3_xi1_u(-offx:row_length-1+offx,1-offy:rows+offy,                  &
           model_levels),                                             &
  h3_xi2_v(1-offx:row_length+offx,-offy:n_rows-1+offy,                &
           model_levels),                                             &
  h3_p_eta(1-offx:row_length+offx,1-offy:rows+offy,                   &
           0:model_levels)



REAL, INTENT (INOUT) :: e_trb(1-halo_i:row_length+halo_i,                 &
                                1-halo_j:rows+halo_j,model_levels)
REAL, INTENT (INOUT) :: tsq(1-halo_i:row_length+halo_i,                   &
                                1-halo_j:rows+halo_j,model_levels)
REAL, INTENT (INOUT) :: qsq(1-halo_i:row_length+halo_i,                   &
                                1-halo_j:rows+halo_j,model_levels)
REAL, INTENT (INOUT) :: cov(1-halo_i:row_length+halo_i,                   &
                                1-halo_j:rows+halo_j,model_levels)

REAL :: super_array(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,      &
                0:model_levels, moisture_array_size )

REAL :: super_array_out(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,  &
                0:model_levels, moisture_array_size )

INTEGER count

      IF (lhook) CALL dr_hook('eg_sl_bdy_tke',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 4.0  Calculate turbulent variables at departure point.
!              Trajectory limit near ground already done in section 2.
!              Trajectory limit at top required if
!              model_levels > tke_levels - 1.
! ----------------------------------------------------------------------

      IF (Error_Code  ==  0 .AND.                                       &
                            bdy_tke > 0 .AND. l_adv_turb_field) THEN
! ----------------------------------------------------------------------
! Section 4.1  Put values into super arrays
! ----------------------------------------------------------------------

!       e_trb(:,:,K) denotes valus on the (K-1)th theta level.
        DO k = 1, tke_levels - 1
          DO j = 1, rows
            DO i = 1, row_length
              super_array(i,j,k,1) = e_trb(i,j,k + 1)
            END DO
          END DO
        END DO
        k = 0
          DO j = 1, rows
            DO i = 1, row_length
              super_array(i,j,k,1) = e_trb(i,j,k + 1)
            END DO
          END DO

        count = 1

!       (:,:,K) denotes valus on the (K-1)th theta level.
        IF (bdy_tke == 3) THEN
          DO k = 1, tke_levels - 1
            DO j = 1, rows
              DO i = 1, row_length
                super_array(i,j,k,2) = tsq(i,j,k + 1)
                super_array(i,j,k,3) = qsq(i,j,k + 1)
                super_array(i,j,k,4) = cov(i,j,k + 1)
              END DO
            END DO
          END DO
          count = 4
        END IF

        CALL Swap_Bounds(                                             &
             super_array, row_length, rows,                           &
             count*(model_levels+1),                                  &
             halo_i, halo_j, fld_type_p, .false. )


        DO k = 1, tke_levels - 1
          DO j = 1, rows
            DO i = 1, row_length
              e_trb(i,j,k + 1) = super_array_out(i,j,k,1)
            END DO
          END DO
        END DO
        IF (bdy_tke == 3) THEN
          DO k = 1, tke_levels - 1
            DO j = 1, rows
              DO i = 1, row_length
                tsq(i,j,k + 1) = super_array_out(i,j,k,2)
                qsq(i,j,k + 1) = super_array_out(i,j,k,3)
                cov(i,j,k + 1) = super_array_out(i,j,k,4)
              END DO
            END DO
          END DO
        END IF
      END IF ! error_code=0 .AND. ! bdy_tke > 0 .AND. l_adv_turb_field

! End of routine.
      IF (lhook) CALL dr_hook('eg_sl_bdy_tke',zhook_out,zhook_handle)

END SUBROUTINE eg_sl_bdy_tke
END MODULE eg_sl_thermo_mod
