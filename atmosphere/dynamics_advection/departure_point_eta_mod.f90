! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Find departure point 
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics

 MODULE departure_point_eta_mod
 IMPLICIT NONE
 
 CONTAINS
 SUBROUTINE departure_point_eta(                                        &
                  row_length, rows, n_rows, model_levels,g_i_pe,        &
                  pnt_type, dep_row_len, dep_rows,                      &
                  off_i, off_j, off_k,offz, l_rk_dps,depart_order,      &
                  depart_high_order_scheme,                             &
                  depart_monotone_scheme, first_constant_r_rho_level,   &
                  l_depart_high,l_depart_mono, l_shallow,               &
                  lam_max_cfl, alpha, eta_theta_levels, eta_rho_levels, &
                  xi3_at_theta, xi3_at_rho, xi3_at_u, xi3_at_v,         &
                  etadot, u_np1, v_np1, etadot_np1,                     &
                  u, v, w, depart_xi1,depart_xi2,depart_xi3) 

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE conversions_mod
USE earth_constants_mod
USE eg_interpolation_eta_mod
USE eg_adjust_vert_bound_mod
USE horiz_grid_mod

USE domain_params
USE UM_ParParams
USE UM_parvars, ONLY:  offx,offy,halo_i,halo_j, datastart,              &
                       gc_proc_row_group,gc_proc_col_group,at_extremity,&
                       nproc, nproc_x, nproc_y
USE proc_info_mod,  ONLY:  model_domain,global_row_length,global_rows,mype=>me
USE Field_Types
USE atm_fields_bounds_mod
USE dynamics_input_mod,    ONLY : l_regular

USE eg_check_sl_domain_mod
USE eg_parameters_mod, ONLY : l_slice,interp_dpt_pt
USE atm_step_local,    ONLY : cycleno
USE timestep_mod,      ONLY : timestep
USE departure_pts_mod, ONLY:  depart_xi1_w_disp,  depart_xi2_w_disp
IMPLICIT NONE
!
! Description: 
!   Find u,v,w,rho departure point using a local cartesian transform method.
!
! Method: 
!   ENDGame formulation version 3.02 (John Thuburns modification)
!
! Code Description:
!   Language: FORTRAN 90.
!   This code is written to UMDP3 v8.0 programming standards.
!
! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
 
! Model dimensions

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,          &
                       pnt_type, first_constant_r_rho_level


INTEGER, INTENT(IN) :: dep_row_len, dep_rows, off_i, off_j,             &
                       off_k, offz

! MPP options
INTEGER, INTENT(IN) ::                                                  &
        g_i_pe(1-halo_i:global_row_length+halo_i),                      &
                           ! processor on my processor-row 
                           ! holding a given value in i direction
        LAM_max_cfl(2)     ! Max CFL for a LAM allowed near the 
                           ! boundaries 
 
! Integer parameters for advection

INTEGER, INTENT(IN) ::                                                  &
        depart_order,                                                   &
                           ! for the chosen departure point scheme how
                           ! many iterations/terms to use.
        depart_high_order_scheme,                                       &
                           ! code choosing high order
                           ! interpolation scheme used in Depart routine
        depart_monotone_scheme
                           ! code choosing monotone
                           ! interpolation scheme used in Depart routine

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                  &
        l_depart_high,                                                  &
                         ! True if high order interpolation scheme to 
                         ! be used in Depart scheme                         
        l_depart_mono,                                                  &
                         ! True if monotone interpolation scheme to     
                         ! be used in Depart scheme
        l_shallow

LOGICAL, INTENT(IN) :: l_rk_dps ! Flag to enable RK2 scheme

REAL, INTENT(IN) :: alpha

REAL, INTENT(IN) ::                                                     &
        eta_theta_levels(0:model_levels),                               &
        eta_rho_levels(1-offz:model_levels+offz)

REAL, INTENT(IN) ::                                                     &
        xi3_at_theta(1-halo_i:row_length+halo_i,                        &
                     1-halo_j:rows+halo_j,0:model_levels),              &
                    ! vertical coordinate at theta points
        xi3_at_rho(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,     &
                   1-offz:model_levels+offz),                           &
                    ! vertical coordinate at rho points. Vert offsetting 
                    ! to allow use of extended height when needed 
        xi3_at_u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,      &
                   1-offz:model_levels+offz),                           &
                    ! vertical coordinate at u points
        xi3_at_v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,    &
                   1-offz:model_levels+offz),                           &
                    ! vertical coordinate at v points
        etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
        u(-halo_i:row_length-1+halo_i,1-halo_j:rows+halo_j,             &
                   1-offz:model_levels+offz),                           &
                    ! Use vertical offsetting to allow use of extended   
                    ! fields for w-point interpolation
        v(1-halo_i:row_length+halo_i,-halo_j:n_rows-1+halo_j,           &
                   1-offz:model_levels+offz),                           &
                    ! Use vertical offsetting to allow use of extended   
                    ! fields for w-point interpolation
        w(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,              &
                   0:model_levels),                                     &
        u_np1(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),   &
        v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels), &
        etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,             &
                   0:model_levels)

! Departure point coordinates

REAL, DIMENSION(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,       &
                1-off_k:model_levels), INTENT(OUT) ::                   &
                      depart_xi1, depart_xi2, depart_xi3

! Local variables

INTEGER :: uv_mlevels_in, uv_mlevels_out, w_mlevels_in,                 &
           w_mlevels_out, uv_first_constant_rho,                        &
           number_of_inputs, i, j, k, kk, extra_i, extra_j,             &
           error_code

REAL :: beta, dxi1, dxi2,                                               &
              max_xi1, min_xi1, max_xi2, min_xi2

REAL, DIMENSION(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,       &
                1-off_k:model_levels) ::                                &
               int_wind_u, int_wind_v, int_wind_w

! Components of rotation matrix

REAL :: RM_11, RM_12, Q_33, gam

! Arrival and departure velocities

REAL :: u_a, v_a, w_a, u_d, v_d, w_d
REAL :: snxi2_d, csxi2_d, csxi1_d, snxi1_d
REAL :: x_d, y_d, z_d, Is, rad

! Temporary local

REAL    :: temp_Csalpha
INTEGER :: jv_begin, jv_stop, dep_its


REAL :: ddepxi2dy,pm_pi 

! End of header

IF (lhook) CALL dr_hook('DEPARTURE_POINT_ETA',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 1.0 Initialisation 
!-----------------------------------------------------------------------

dep_its = depart_order
IF( L_RK_DPS ) dep_its = 1

jv_begin = vdims%j_start
jv_stop  = vdims%j_end
IF( .NOT. L_SLICE ) THEN
  IF( model_domain == mt_Global ) THEN
    IF( at_extremity(PSouth) ) jv_begin = jv_begin + 1
    IF( at_extremity(PNorth) ) jv_stop  = jv_stop  - 1
  END IF
END IF

pm_pi = pi
IF( at_extremity(PSouth) ) pm_pi = -pm_pi

Is = 1.0
IF( l_shallow ) Is = 0.0
      
IF ( L_regular ) THEN
   dxi1 = delta_xi1
   dxi2 = delta_xi2
ELSE  ! variable resolution
   dxi1 = xi1_u(udims%i_end)-xi1_u(udims%i_end-1)
   dxi2 = xi2_v(vdims%i_end)-xi2_v(vdims%i_end-1)
END IF ! L_regular

! Set Lateral boundary limits.

   min_xi1 = (2-halo_i)*dxi1
   max_xi1 = 2.*pi + (halo_i-2)*dxi1

   min_xi2 = -pi*.5
   max_xi2 = pi*.5

number_of_inputs = 1 
beta = 1.0 - alpha

!-----------------------------------------------------------------------
! 2.0 Compute a departure point for a u, v, w, rho point.
!-----------------------------------------------------------------------


IF( L_RK_DPS .AND. cycleno < 3 ) THEN
   SELECT CASE(pnt_type)
      CASE(fld_type_u)

!$OMP PARALLEL DO PRIVATE(i,j,k,u_a,v_a,w_a,rad,x_d,y_d,z_d) DEFAULT(NONE)     &
!$OMP&  SHARED(model_levels,udims,u_np1,intw_p2u,intw_v2p,intw_w2rho,v_np1,    &
!$OMP&  etadot_np1,IS,xi3_at_u,timestep,xi1_u,csxi2_p,snxi2_p,                 &
!$OMP&  eta_rho_levels,depart_xi1,depart_xi2, depart_xi3) SCHEDULE(STATIC)
         DO k = 1, model_levels
            DO j = udims%j_start, udims%j_end
               DO i = udims%i_start, udims%i_end

! Arrival points

                  u_a    = u_np1(i,j,k)
                  v_a    = intw_p2u(i,1)*(intw_v2p(j,1)*v_np1(i,j-1,k)         &
                                         +intw_v2p(j,2)*v_np1(i,j,k) )         &
                          +intw_p2u(i,2)*(intw_v2p(j,1)*v_np1(i+1,j-1,k)       &
                                         +intw_v2p(j,2)*v_np1(i+1,j,k) )
                  w_a    = intw_p2u(i,1)*(intw_w2rho(k,1)*etadot_np1(i,j,k)    &
                                         +intw_w2rho(k,2)*etadot_np1(i,j,k-1)) &
                          +intw_p2u(i,2)*(intw_w2rho(k,1)*etadot_np1(i+1,j,k)  &
                                         +intw_w2rho(k,2)*etadot_np1(i+1,j,k-1))

                  rad    = 1.0/(Earth_radius                                   &
                                   + Is*(xi3_at_u(i,j,k) - Earth_radius))
                  u_a    = u_a*rad
                  v_a    = v_a*rad

! Update departure points

                  x_d = -timestep*u_a
                  y_d = -timestep*v_a
                  z_d = SQRT(1.0 - x_d**2 - y_d**2)

                  depart_xi1(i,j,k) = xi1_u(i) +                               &
                                      ATAN2(x_d,z_d*csxi2_p(j)-y_d*snxi2_p(j))
                  depart_xi2(i,j,k) = ASIN(y_d*csxi2_p(j)+z_d*snxi2_p(j))
                  depart_xi3(i,j,k) = eta_rho_levels(k) - timestep*w_a
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO

      CASE(fld_type_v)

!$OMP PARALLEL DO PRIVATE(i,j,k,u_a,v_a,w_a,rad,x_d,y_d,z_d) DEFAULT(NONE)     &
!$OMP&  SHARED(model_levels,jv_begin,jv_stop,vdims,intw_p2v,intw_u2p,u_np1,    &
!$OMP&  intw_w2rho,etadot_np1,Is,xi3_at_v,timestep,v_np1,                      &
!$OMP&  xi1_p,csxi2_v,snxi2_v,depart_xi1,depart_xi2,depart_xi3,                &
!$OMP&  eta_rho_levels) SCHEDULE(STATIC)
         DO k = 1, model_levels
            DO j = jv_begin, jv_stop
               DO i = vdims%i_start, vdims%i_end

! Arrival points

                  u_a    = intw_p2v(j,1)*(intw_u2p(i,1)*u_np1(i-1,j,k)         &
                                         +intw_u2p(i,2)*u_np1(i,j,k))          &
                          +intw_p2v(j,2)*(intw_u2p(i,1)*u_np1(i-1,j+1,k)       &
                                         +intw_u2p(i,2)*u_np1(i,j+1,k))
                  v_a    = v_np1(i,j,k)
                  w_a    = intw_p2v(j,1)*(intw_w2rho(k,1)*etadot_np1(i,j,k)    &
                                         +intw_w2rho(k,2)*etadot_np1(i,j,k-1)) &
                          +intw_p2v(j,2)*(intw_w2rho(k,1)*etadot_np1(i,j+1,k)  &
                                         +intw_w2rho(k,2)*etadot_np1(i,j+1,k-1))

                  rad    = 1.0/(Earth_radius                                   &
                                   + Is*(xi3_at_v(i,j,k) - Earth_radius))

                  u_a = u_a*rad
                  v_a = v_a*rad

! Update departure points

                  x_d = -timestep*u_a
                  y_d = -timestep*v_a
                  z_d = SQRT(1.0 - x_d**2 - y_d**2)

                  depart_xi1(i,j,k) = xi1_p(i) +                               &
                                      ATAN2(x_d,z_d*csxi2_v(j)-y_d*snxi2_v(j))
                  depart_xi2(i,j,k) = ASIN(y_d*csxi2_v(j)+z_d*snxi2_v(j))
                  depart_xi3(i,j,k) = eta_rho_levels(k) - timestep*w_a
               END DO
            END DO
         END DO                                       
!$OMP END PARALLEL DO

      CASE(fld_type_w)

!$OMP PARALLEL DO PRIVATE(i,j,k,u_a,v_a,w_a,rad,x_d,y_d,z_d,ddepxi2dy)         &
!$OMP&  DEFAULT(NONE) SHARED(model_levels,pdims,intw_u2p,u_np1,intw_v2p,       &
!$OMP&  v_np1,Is,xi3_at_theta,timestep,csxi2_p,snxi2_p,intw_rho2w,etadot_np1,  &
!$OMP&  depart_xi3,xi1_p, depart_xi2_w_disp,                                   &
!$OMP&  depart_xi1_w_disp,depart_xi1,depart_xi2,pm_pi,xi2_p,eta_theta_levels)  &
!$OMP& SCHEDULE(STATIC)
         DO k = 0, model_levels

           IF( k == 0 ) THEN
             DO j = pdims%j_start, pdims%j_end
                DO i = pdims%i_start, pdims%i_end

! Arrival points

                  u_a    = intw_u2p(i,1)*u_np1(i-1,j,k+1) +                       &
                           intw_u2p(i,2)*u_np1(i,j,k+1)
                  v_a    = intw_v2p(j,1)*v_np1(i,j-1,k+1) +                       &
                           intw_v2p(j,2)*v_np1(i,j,k+1)

                  rad    = 1.0/(Earth_radius                                      &
                                   + Is*(xi3_at_theta(i,j,k) - Earth_radius))

                  u_a    = u_a*rad
                  v_a    = v_a*rad

! Update departure points

                  x_d = -timestep*u_a
                  y_d = -timestep*v_a
                  z_d = SQRT(1.0 - x_d**2 - y_d**2)

                  depart_xi1_w_disp(i,j,k) = ATAN2(x_d,z_d*csxi2_p(j)-       &
                                                       y_d*snxi2_p(j))
                  depart_xi1(i,j,k) = xi1_p(i) + depart_xi1_w_disp(i,j,k)
                                  
                  depart_xi2(i,j,k) = ASIN(y_d*csxi2_p(j)+z_d*snxi2_p(j))

                  ddepxi2dy = (csxi2_p(j)-y_d*snxi2_p(j)/z_d)                &
                              /COS(depart_xi2(i,j,k))

                  IF(ddepxi2dy > 0.0) THEN
                    depart_xi2_w_disp(i,j,k) = depart_xi2(i,j,k) - xi2_p(j)
                  ELSE
                    depart_xi2_w_disp(i,j,k) = pm_pi- depart_xi2(i,j,k)      &
                                                             - xi2_p(j)
                  END IF

                  depart_xi3(i,j,k) = eta_theta_levels(k)
                END DO
             END DO
           ELSE IF( k == model_levels ) THEN
             DO j = pdims%j_start, pdims%j_end
               DO i = pdims%i_start, pdims%i_end

! Arrival points

                  u_a    = intw_u2p(i,1)*u_np1(i-1,j,k) +                    &
                           intw_u2p(i,2)*u_np1(i,j,k)
                  v_a    = intw_v2p(j,1)*v_np1(i,j-1,k) +                    &
                           intw_v2p(j,2)*v_np1(i,j,k)

                  rad    = 1.0/(Earth_radius                                 &
                                   + Is*(xi3_at_theta(i,j,k) - Earth_radius))

                  u_a = u_a*rad
                  v_a = v_a*rad

! Update departure points

                  x_d = -timestep*u_a
                  y_d = -timestep*v_a
                  z_d = SQRT(1.0 - x_d**2 - y_d**2)

                  depart_xi1_w_disp(i,j,k) = ATAN2(x_d,z_d*csxi2_p(j)-       &
                                                      y_d*snxi2_p(j))
                  depart_xi1(i,j,k) = xi1_p(i) + depart_xi1_w_disp(i,j,k)

                  depart_xi2(i,j,k) = ASIN(y_d*csxi2_p(j)+z_d*snxi2_p(j))

                  ddepxi2dy = (csxi2_p(j)-y_d*snxi2_p(j)/z_d)                 &
                              /COS(depart_xi2(i,j,k))

                  IF(ddepxi2dy > 0.0) THEN
                    depart_xi2_w_disp(i,j,k) = depart_xi2(i,j,k) - xi2_p(j)
                  ELSE
                    depart_xi2_w_disp(i,j,k) = pm_pi- depart_xi2(i,j,k)       &
                                                    - xi2_p(j)
                  END IF

                  depart_xi3(i,j,k) = eta_theta_levels(k)
               END DO
             END DO

           ELSE
            DO j = pdims%j_start, pdims%j_end
               DO i = pdims%i_start, pdims%i_end

! Arrival points

                  u_a = intw_rho2w(k,1)*(intw_u2p(i,1)*u_np1(i-1,j,k+1)        &
                                        +intw_u2p(i,2)*u_np1(i,j,k+1))         &
                       +intw_rho2w(k,2)*(intw_u2p(i,1)*u_np1(i-1,j,k)          &
                                        +intw_u2p(i,2)*u_np1(i,j,k))
                  v_a = intw_rho2w(k,1)*(intw_v2p(j,1)*v_np1(i,j-1,k+1)        &
                                        +intw_v2p(j,2)*v_np1(i,j,k+1))         &
                       +intw_rho2w(k,2)*(intw_v2p(j,1)*v_np1(i,j-1,k)          &
                                        +intw_v2p(j,2)*v_np1(i,j,k))
                  w_a = etadot_np1(i,j,k)

                  rad    = 1.0/(Earth_radius                                   &
                                   + Is*(xi3_at_theta(i,j,k) - Earth_radius))

                  u_a = u_a*rad
                  v_a = v_a*rad

! Update departure points

                  x_d = -timestep*u_a
                  y_d = -timestep*v_a
                  z_d = SQRT(1.0 - x_d**2 - y_d**2)

                  depart_xi1_w_disp(i,j,k) = ATAN2( x_d,z_d*csxi2_p(j)         &
                                                   -y_d*snxi2_p(j))
                  depart_xi1(i,j,k) = xi1_p(i) + depart_xi1_w_disp(i,j,k)

                  depart_xi2(i,j,k) = ASIN(y_d*csxi2_p(j)+z_d*snxi2_p(j))

                  ddepxi2dy = (csxi2_p(j)-y_d*snxi2_p(j)/z_d)                  &
                              /COS(depart_xi2(i,j,k))

                  IF(ddepxi2dy > 0.0) THEN
                    depart_xi2_w_disp(i,j,k) = depart_xi2(i,j,k) - xi2_p(j)
                  ELSE
                    depart_xi2_w_disp(i,j,k) = pm_pi- depart_xi2(i,j,k)        &
                                                    - xi2_p(j)
                  END IF

                  depart_xi3(i,j,k) = eta_theta_levels(k) - timestep*w_a
               END DO
            END DO
           ENDIF
       END DO                                       
!$OMP END PARALLEL DO

      CASE(fld_type_p)

!$OMP PARALLEL DO PRIVATE(i,j,k,u_a,v_a,w_a,rad,x_d,y_d,z_d) DEFAULT(NONE)    &
!$OMP& SHARED(model_levels,pdims,intw_u2p,u_np1,intw_v2p,v_np1,intw_w2rho,    &
!$OMP& etadot_np1,Is,xi3_at_rho,timestep,xi1_p,csxi2_p,                       &
!$OMP& snxi2_p,eta_rho_levels,depart_xi1,depart_xi2,depart_xi3)               &
!$OMP& SCHEDULE(STATIC)
         DO k = 1, model_levels
            DO j = pdims%j_start, pdims%j_end
               DO i = pdims%i_start, pdims%i_end

! Arrival points

                  u_a    = intw_u2p(i,1)*u_np1(i-1,j,k) +                      &
                           intw_u2p(i,2)*u_np1(i,j,k)
                  v_a    = intw_v2p(j,1)*v_np1(i,j-1,k) +                      &
                           intw_v2p(j,2)*v_np1(i,j,k)
                  w_a    = intw_w2rho(k,1)*etadot_np1(i,j,k) +                 &
                           intw_w2rho(k,2)*etadot_np1(i,j,k-1)

                  rad    = 1.0/(Earth_radius                                   &
                                   + Is*(xi3_at_rho(i,j,k) - Earth_radius))

                  u_a    = u_a*rad
                  v_a    = v_a*rad

! Update departure points

                  x_d = -timestep*u_a
                  y_d = -timestep*v_a
                  z_d = sqrt(1.0 - x_d**2 - y_d**2)

                  depart_xi1(i,j,k) = xi1_p(i) +                               &
                                      ATAN2(x_d,z_d*csxi2_p(j)-y_d*snxi2_p(j))
                  depart_xi2(i,j,k) = ASIN(y_d*csxi2_p(j)+z_d*snxi2_p(j))
                  depart_xi3(i,j,k) = eta_rho_levels(k) - timestep*w_a
               END DO
            END DO
         END DO                                       
!$OMP END PARALLEL DO

! Check the values are in range

   END SELECT

   CALL eg_check_sl_domain( model_domain, depart_xi2, depart_xi1,        &
                      dep_row_len, dep_rows, model_levels+off_k,         &
                      max_xi1, min_xi1, max_xi2, min_xi2)

   IF(.NOT.interp_dpt_pt .OR. cycleno > 1 ) THEN
      CALL eg_adjust_vert_bound(                                         &
            row_length, rows, model_levels, g_i_pe,                      &
            pnt_type, dep_row_len,                                       &
            dep_rows, off_i, off_j, off_k, offz,                         &
            etadot, etadot_np1, depart_xi1, depart_xi2,                  &
            depart_xi3)
   END IF

END IF

IF( L_RK_DPS .AND. cycleno == 1 ) THEN
   IF (lhook) CALL dr_hook('DEPARTURE_POINT_ETA',zhook_out,zhook_handle)
   RETURN
END IF


SELECT CASE (pnt_type)
   CASE(fld_type_u)

! Settings for levels used in Interpolation subroutine
      uv_mlevels_in  = model_levels                                  
      uv_mlevels_out = model_levels                                  
      w_mlevels_in   = model_levels + 1                                 
      w_mlevels_out  = model_levels   
      uv_first_constant_rho = first_constant_r_rho_level

   CASE(fld_type_v)

! Settings for levels used in Interpolation subroutine
      uv_mlevels_in  = model_levels                                  
      uv_mlevels_out = model_levels                                  
      w_mlevels_in   = model_levels  + 1                                
      w_mlevels_out  = model_levels                     
      uv_first_constant_rho = first_constant_r_rho_level             

   CASE(fld_type_w)

! Settings for levels used in Interpolation subroutine
      uv_mlevels_in  = model_levels  + 2                                
      uv_mlevels_out = model_levels  + 1                              
      w_mlevels_in   = model_levels  + 1
      w_mlevels_out  = model_levels  + 1               
      uv_first_constant_rho = first_constant_r_rho_level + 2

   CASE(fld_type_p)

! Settings for levels used in Interpolation subroutine
      uv_mlevels_in  = model_levels                                  
      uv_mlevels_out = model_levels                                  
      w_mlevels_in   = model_levels + 1                                
      w_mlevels_out  = model_levels                     
      uv_first_constant_rho = first_constant_r_rho_level

END SELECT


!-----------------------------------------------------------------------
! 
! 2.2: Iterate algorithm
! 
!-----------------------------------------------------------------------

DO kk=1, dep_its

   CALL eg_interpolation_eta(                                            &
                           eta_rho_levels, fld_type_u, number_of_inputs, &
                           row_length, rows, uv_mlevels_in, rows,        &
                           dep_row_len, dep_rows, uv_mlevels_out,        &
                           depart_high_order_scheme,                     &
                           depart_monotone_scheme,model_domain,          &
                           L_depart_high, L_depart_mono,                 &
                           depart_xi3, depart_xi1, depart_xi2,           &
                           mype, nproc, nproc_x, nproc_y,                &
                           halo_i, halo_j,                               &
                           global_row_length, datastart, at_extremity,   &
                           g_i_pe, gc_proc_row_group, gc_proc_col_group, &
                           0, 0,offx,offy,error_code,                    &
                           u,int_wind_u)

   CALL eg_interpolation_eta(                                            &
                           eta_rho_levels,fld_type_v, number_of_inputs,  &
                           row_length, n_rows, uv_mlevels_in,            &
                           rows,                                         &
                           dep_row_len, dep_rows, uv_mlevels_out,        &
                           depart_high_order_scheme,                     &
                           depart_monotone_scheme,model_domain,          &
                           L_depart_high, L_depart_mono,                 &
                           depart_xi3, depart_xi1, depart_xi2,           &
                           mype, nproc, nproc_x, nproc_y,                &
                           halo_i, halo_j,                               &
                           global_row_length, datastart, at_extremity,   &
                           g_i_pe, gc_proc_row_group, gc_proc_col_group, &
                           0,0,offx,offy,error_code,                     &
                           v,int_wind_v)

   CALL eg_interpolation_eta(                                            &
                           eta_theta_levels,fld_type_w,                  &
                           number_of_inputs,                             &
                           row_length, rows, w_mlevels_in,               &
                           rows,                                         &
                           dep_row_len, dep_rows, w_mlevels_out,         &
                           depart_high_order_scheme,                     &
                           depart_monotone_scheme,model_domain,          &
                           L_depart_high, L_depart_mono,                 &
                           depart_xi3, depart_xi1, depart_xi2,           &
                           mype, nproc, nproc_x, nproc_y,                &
                           halo_i, halo_j,                               &
                           global_row_length, datastart, at_extremity,   &
                           g_i_pe, gc_proc_row_group, gc_proc_col_group, &
                           0, 0,offx,offy,error_code,                    &
                           w,int_wind_w)

   SELECT CASE(pnt_type)
      CASE(fld_type_u)

! Compute u-departure point using current iteration estimates

!$OMP PARALLEL DO PRIVATE(i,j,k,u_a,v_a,w_a,u_d,v_d,w_d,rad,x_d,               &
!$OMP& y_d,z_d,rm_11,rm_12,                                                    &
!$OMP& gam,q_33,temp_Csalpha,snxi1_d,csxi1_d,snxi2_d,csxi2_d) DEFAULT(NONE)    &
!$OMP& SHARED(model_levels,udims,depart_xi1,xi1_u,depart_xi2,snxi2_p,csxi2_p,  &
!$OMP& timestep,u_np1,intw_p2u,intw_v2p,v_np1,intw_w2rho,etadot_np1,           &
!$OMP& Is,xi3_at_u,int_wind_u,int_wind_v,int_wind_w,depart_xi3,eta_rho_levels, &
!$OMP& alpha,beta) SCHEDULE(STATIC)
         DO k = 1, model_levels
            DO j = udims%j_start, udims%j_end
               DO i = udims%i_start, udims%i_end

                  snxi1_d = SIN(depart_xi1(i,j,k)-xi1_u(i))
                  csxi1_d = COS(depart_xi1(i,j,k)-xi1_u(i))
                  snxi2_d = SIN(depart_xi2(i,j,k))
                  csxi2_d = COS(depart_xi2(i,j,k))

                  q_33         = snxi2_d*snxi2_p(j)+csxi2_d*csxi2_p(j)*csxi1_d
                  gam          = (1.0 + q_33)
                  temp_Csalpha = 1.0/gam
                  gam          = 0.5*gam*timestep

! Components of rotation matrix

                  rm_11 = (csxi2_d*csxi2_p(j)                                  &
                               +(1.0+snxi2_p(j)*snxi2_d)*csxi1_d)              &
                         *temp_Csalpha
                  rm_12 =-Snxi1_d*(Snxi2_p(j) + Snxi2_d)*temp_Csalpha

! Arrival points

                  u_a    = u_np1(i,j,k)
                  v_a    = intw_p2u(i,1)*(intw_v2p(j,1)*v_np1(i,j-1,k)         &
                                         +intw_v2p(j,2)*v_np1(i,j,k) )         &
                          +intw_p2u(i,2)*(intw_v2p(j,1)*v_np1(i+1,j-1,k)       &
                                         +intw_v2p(j,2)*v_np1(i+1,j,k) )
                  w_a    = intw_p2u(i,1)*(intw_w2rho(k,1)*etadot_np1(i,j,k)    &
                                         +intw_w2rho(k,2)*etadot_np1(i,j,k-1)) &
                          +intw_p2u(i,2)*(intw_w2rho(k,1)*etadot_np1(i+1,j,k)  &
                                         +intw_w2rho(k,2)*etadot_np1(i+1,j,k-1))

                  rad    = 1.0/(Earth_radius                                   &
                                   + Is*(xi3_at_u(i,j,k) - Earth_radius))
                  u_a    = u_a *rad
                  v_a    = v_a *rad

! Rotated velocities at departure points

                  u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
                  v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)
                  w_d = int_wind_w(i,j,k)

! Update departure points

                  x_d = -gam*(alpha*u_a + beta*u_d)
                  y_d = -gam*(alpha*v_a + beta*v_d)
                  z_d = SQRT(1.0 - x_d**2 - y_d**2)

                  depart_xi1(i,j,k) = xi1_u(i) +                               &
                                      ATAN2(x_d,z_d*csxi2_p(j)-y_d*snxi2_p(j))
                  depart_xi2(i,j,k) = ASIN(y_d*csxi2_p(j)+z_d*snxi2_p(j))

                  depart_xi3(i,j,k) = eta_rho_levels(k)                        &
                                     -timestep*(alpha*w_a + beta*w_d)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO

      CASE(fld_type_v)

!$OMP PARALLEL DO PRIVATE(i,j,k,u_a,v_a,w_a,u_d,v_d,w_d,rad,x_d,              &
!$OMP& y_d,z_d,rm_11,rm_12,                                                   & 
!$OMP& gam,q_33,temp_Csalpha,snxi1_d,csxi1_d,snxi2_d,csxi2_d) DEFAULT(NONE)   &
!$OMP& SHARED(model_levels,vdims,depart_xi1,depart_xi2,                       &
!$OMP& timestep,u_np1,intw_w2rho,etadot_np1, v_np1,                           &
!$OMP& Is,int_wind_u,int_wind_v,int_wind_w, xi3_at_v,alpha,beta,              &
!$OMP& intw_u2p,eta_rho_levels,depart_xi3,jv_begin,jv_stop, csxi2_v,          &
!$OMP& intw_p2v,xi1_p,snxi2_v) SCHEDULE(STATIC)
         DO k = 1, model_levels
            DO j = jv_begin, jv_stop
               DO i = vdims%i_start, vdims%i_end

                  snxi1_d = SIN(depart_xi1(i,j,k)-xi1_p(i))
                  csxi1_d = COS(depart_xi1(i,j,k)-xi1_p(i))
                  snxi2_d = SIN(depart_xi2(i,j,k))
                  csxi2_d = COS(depart_xi2(i,j,k))

                  q_33         = snxi2_d*snxi2_v(j)+csxi2_d*csxi2_v(j)*csxi1_d
                  gam          = (1.0 + q_33)
                  temp_Csalpha = 1.0/gam
                  gam          = 0.5*gam*timestep

! Components of rotation matrix
                  rm_11 = (csxi2_d*csxi2_v(j)                                  &
                            +(1.0+snxi2_v(j)*snxi2_d)*csxi1_d)*temp_Csalpha
                  rm_12 =-snxi1_d*(snxi2_v(j) + snxi2_d)*temp_Csalpha

! Arrival points

                  u_a    = intw_p2v(j,1)*(intw_u2p(i,1)*u_np1(i-1,j,k)         &
                                         +intw_u2p(i,2)*u_np1(i,j,k))          &
                          +intw_p2v(j,2)*(intw_u2p(i,1)*u_np1(i-1,j+1,k)       &
                                         +intw_u2p(i,2)*u_np1(i,j+1,k))
                  v_a    = v_np1(i,j,k)
                  w_a    = intw_p2v(j,1)*(intw_w2rho(k,1)*etadot_np1(i,j,k)    &
                                         +intw_w2rho(k,2)*etadot_np1(i,j,k-1)) &
                          +intw_p2v(j,2)*(intw_w2rho(k,1)*etadot_np1(i,j+1,k)  &
                                         +intw_w2rho(k,2)*etadot_np1(i,j+1,k-1))

                  rad    = 1.0/(Earth_radius                                   &
                                   + Is*(xi3_at_v(i,j,k) - Earth_radius))

                  u_a = u_a *rad
                  v_a = v_a *rad

! Rotated velocities at departure points

                  u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
                  v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k) 
                  w_d = int_wind_w(i,j,k)

! Update departure points

                  x_d = -gam*(alpha*u_a + beta*u_d)
                  y_d = -gam*(alpha*v_a + beta*v_d)
                  z_d = SQRT(1.0 - x_d**2 - y_d**2)

                  depart_xi1(i,j,k) = xi1_p(i) +                               &
                                      ATAN2(x_d,z_d*csxi2_v(j)-y_d*snxi2_v(j))
                  depart_xi2(i,j,k) = ASIN(y_d*csxi2_v(j)+z_d*snxi2_v(j))
                  depart_xi3(i,j,k) = eta_rho_levels(k)                        &
                                     - timestep*(alpha*w_a + beta*w_d)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO

      CASE(fld_type_w)

! Compute w-departure point using current iteration estimates 
! Use no vertical shear assumption for horizontal momentum

         k = 0
         DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end

               snxi1_d = SIN(depart_xi1(i,j,k)-xi1_p(i))
               csxi1_d = COS(depart_xi1(i,j,k)-xi1_p(i))
               snxi2_d = SIN(depart_xi2(i,j,k))
               csxi2_d = COS(depart_xi2(i,j,k))

               q_33         = snxi2_d*snxi2_p(j)+csxi2_d*csxi2_p(j)*csxi1_d
               gam          = (1.0 + q_33)
               temp_Csalpha = 1.0/gam
               gam          = 0.5*gam*timestep

! Components of rotation matrix
               rm_11 = (csxi2_d*csxi2_p(j) + (1+snxi2_p(j)*snxi2_d)*csxi1_d)   &
                       *temp_Csalpha
               rm_12 =-snxi1_d*(snxi2_p(j) + snxi2_d)*temp_Csalpha

! Arrival points

               u_a    = intw_u2p(i,1)*u_np1(i-1,j,k+1) +                       &
                        intw_u2p(i,2)*u_np1(i,j,k+1)
               v_a    = intw_v2p(j,1)*v_np1(i,j-1,k+1) +                       &
                        intw_v2p(j,2)*v_np1(i,j,k+1)

               rad    = 1.0/(Earth_radius                                      &
                                   + Is*(xi3_at_theta(i,j,k) - Earth_radius))

               u_a    = u_a *rad
               v_a    = v_a *rad

! Rotated velocities at departure points

               u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
               v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)

! Update departure points

               x_d = -gam*(alpha*u_a + beta*u_d)
               y_d = -gam*(alpha*v_a + beta*v_d)
               z_d = SQRT(1.0 - x_d**2 - y_d**2)

               depart_xi1_w_disp(i,j,k) = ATAN2(x_d,z_d*csxi2_p(j)           &
                                           -y_d*snxi2_p(j))
               depart_xi1(i,j,k) = xi1_p(i) + depart_xi1_w_disp(i,j,k)

               depart_xi2(i,j,k) = ASIN(y_d*csxi2_p(j)+z_d*snxi2_p(j))

               ddepxi2dy = (csxi2_p(j)-y_d*snxi2_p(j)/z_d)                   &
                              /COS(depart_xi2(i,j,k))

               IF(ddepxi2dy > 0.0) THEN
                 depart_xi2_w_disp(i,j,k) = depart_xi2(i,j,k) - xi2_p(j)
               ELSE
                 depart_xi2_w_disp(i,j,k) = pm_pi- depart_xi2(i,j,k) - xi2_p(j)
               END IF

               depart_xi3(i,j,k) = eta_theta_levels(k)

            END DO
         END DO

!$OMP PARALLEL DO PRIVATE(i,j,k,u_a,v_a,w_a,u_d,v_d,w_d,rad,x_d,y_d,z_d,     &
!$OMP& rm_11,rm_12, gam,q_33,temp_Csalpha,snxi1_d,csxi1_d,                   &
!$OMP& snxi2_d,csxi2_d,ddepxi2dy) DEFAULT(NONE)                              &
!$OMP& SHARED(model_levels,pdims,depart_xi1,depart_xi2,snxi2_p,csxi2_p,      &
!$OMP& timestep,u_np1,intw_v2p,v_np1,etadot_np1,                             &
!$OMP& Is,int_wind_u,int_wind_v,int_wind_w, depart_xi2_w_disp,               &
!$OMP& depart_xi1_w_disp, beta,alpha,xi2_p,eta_theta_levels, depart_xi3,     &
!$OMP& pm_pi, intw_rho2w,xi1_p,xi3_at_theta ,intw_u2p) SCHEDULE(STATIC)
         DO k = 1, model_levels-1
            DO j = pdims%j_start, pdims%j_end
               DO i = pdims%i_start, pdims%i_end

                  snxi1_d = SIN(depart_xi1(i,j,k)-xi1_p(i))
                  csxi1_d = COS(depart_xi1(i,j,k)-xi1_p(i))
                  snxi2_d = SIN(depart_xi2(i,j,k))
                  csxi2_d = COS(depart_xi2(i,j,k))

                  q_33         = snxi2_d*snxi2_p(j)+csxi2_d*csxi2_p(j)*csxi1_d
                  gam          = (1.0 + q_33)
                  temp_Csalpha = 1.0/gam
                  gam          = 0.5*gam*timestep

! Components of rotation matrix
                  rm_11 = (csxi2_d*csxi2_p(j) +                                &
                               (1+snxi2_p(j)*snxi2_d)*csxi1_d)*temp_Csalpha
                  rm_12 =-snxi1_d*(snxi2_p(j) + snxi2_d)*temp_Csalpha

! Arrival points

                  u_a = intw_rho2w(k,1)*(intw_u2p(i,1)*u_np1(i-1,j,k+1)        &
                                        +intw_u2p(i,2)*u_np1(i,j,k+1))         &
                       +intw_rho2w(k,2)*(intw_u2p(i,1)*u_np1(i-1,j,k)          &
                                        +intw_u2p(i,2)*u_np1(i,j,k))
                  v_a = intw_rho2w(k,1)*(intw_v2p(j,1)*v_np1(i,j-1,k+1)        &
                                        +intw_v2p(j,2)*v_np1(i,j,k+1))         &
                       +intw_rho2w(k,2)*(intw_v2p(j,1)*v_np1(i,j-1,k)          &
                                        +intw_v2p(j,2)*v_np1(i,j,k))
                  w_a = etadot_np1(i,j,k)

                  rad    = 1.0/(Earth_radius                                   &
                                   + Is*(xi3_at_theta(i,j,k) - Earth_radius))

                  u_a = u_a *rad
                  v_a = v_a *rad

! Rotated velocities at departure points

                  u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
                  v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)
                  w_d = int_wind_w(i,j,k)

! Update departure points

                  x_d = -gam*(alpha*u_a + beta*u_d)
                  y_d = -gam*(alpha*v_a + beta*v_d)
                  z_d = SQRT(1.0 - x_d**2 - y_d**2)

                  depart_xi1_w_disp(i,j,k) = ATAN2(x_d,z_d*csxi2_p(j)          &
                                               -y_d*snxi2_p(j))
                  depart_xi1(i,j,k) = xi1_p(i) + depart_xi1_w_disp(i,j,k)

                  depart_xi2(i,j,k) = ASIN(y_d*csxi2_p(j)+z_d*snxi2_p(j))

                  ddepxi2dy = (csxi2_p(j)-y_d*snxi2_p(j)/z_d)                  &
                              /COS(depart_xi2(i,j,k))

                  IF(ddepxi2dy > 0.0) THEN
                    depart_xi2_w_disp(i,j,k) = depart_xi2(i,j,k) - xi2_p(j)
                  ELSE
                    depart_xi2_w_disp(i,j,k) = pm_pi- depart_xi2(i,j,k)        &
                                                       - xi2_p(j)
                  END IF

                  depart_xi3(i,j,k) = eta_theta_levels(k)                      &
                                      -timestep*(alpha*w_a + beta*w_d)

               END DO
            END DO
         END DO
!$OMP END PARALLEL DO

! no vertical momentum shear above top rho level
         k = model_levels
         DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end

               snxi1_d = SIN(depart_xi1(i,j,k)-xi1_p(i))
               csxi1_d = COS(depart_xi1(i,j,k)-xi1_p(i))
               snxi2_d = SIN(depart_xi2(i,j,k))
               csxi2_d = COS(depart_xi2(i,j,k))

               q_33         = snxi2_d*snxi2_p(j)+csxi2_d*csxi2_p(j)*csxi1_d
               gam          = (1.0 + q_33)
               temp_Csalpha = 1.0/gam
               gam          = 0.5*gam*timestep

! Components of rotation matrix
               rm_11 = (csxi2_d*csxi2_p(j)                                     &
                               + (1+snxi2_p(j)*snxi2_d)*csxi1_d)*temp_Csalpha
               rm_12 =-Snxi1_d*(snxi2_p(j) + snxi2_d)*temp_Csalpha

! Arrival points

               u_a    = intw_u2p(i,1)*u_np1(i-1,j,k) +                         &
                        intw_u2p(i,2)*u_np1(i,j,k)
               v_a    = intw_v2p(j,1)*v_np1(i,j-1,k) +                         &
                        intw_v2p(j,2)*v_np1(i,j,k)

               rad    = 1.0/(Earth_radius                                      &
                                   + Is*(xi3_at_theta(i,j,k) - Earth_radius))

               u_a = u_a *rad
               v_a = v_a *rad

! Rotated velocities at departure points

               u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
               v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)

! Update departure points

               x_d = -gam*(alpha*u_a + beta*u_d)
               y_d = -gam*(alpha*v_a + beta*v_d)
               z_d = SQRT(1.0 - x_d**2 - y_d**2)

               depart_xi1_w_disp(i,j,k) = ATAN2(x_d,z_d*csxi2_p(j)-          &
                                                y_d*snxi2_p(j))
               depart_xi1(i,j,k) = xi1_p(i) + depart_xi1_w_disp(i,j,k)
               depart_xi2(i,j,k) = ASIN(y_d*csxi2_p(j)+z_d*snxi2_p(j))

               ddepxi2dy = (csxi2_p(j)-y_d*snxi2_p(j)/z_d)                   &
                              /COS(depart_xi2(i,j,k))

               IF(ddepxi2dy > 0.0) THEN
                 depart_xi2_w_disp(i,j,k) = depart_xi2(i,j,k) - xi2_p(j)
               ELSE
                 depart_xi2_w_disp(i,j,k) = pm_pi- depart_xi2(i,j,k) - xi2_p(j)
               END IF

               depart_xi3(i,j,k) = eta_theta_levels(k)

            END DO
         END DO
                                                
      CASE(fld_type_p)

!$OMP PARALLEL DO PRIVATE(i,j,k,u_a,v_a,w_a,u_d,v_d,w_d,rad,                  &
!$OMP& x_d,y_d,z_d,rm_11,rm_12,                                               &
!$OMP& gam,q_33,temp_Csalpha,snxi1_d,csxi1_d,snxi2_d,csxi2_d) DEFAULT(NONE)   &
!$OMP& SHARED(model_levels,pdims,depart_xi1,depart_xi2,snxi2_p,csxi2_p,       &
!$OMP&   timestep,u_np1,intw_v2p,v_np1,intw_w2rho,etadot_np1,                 &
!$OMP&  Is,int_wind_u,int_wind_v,int_wind_w,alpha,beta,eta_rho_levels,        &
!$OMP&  depart_xi3,xi3_at_rho,xi1_p,intw_u2p) SCHEDULE(STATIC)
         DO k = 1, model_levels
            DO j = pdims%j_start, pdims%j_end
               DO i = pdims%i_start, pdims%i_end

                  snxi1_d = SIN(depart_xi1(i,j,k)-xi1_p(i))
                  csxi1_d = COS(depart_xi1(i,j,k)-xi1_p(i))
                  snxi2_d = SIN(depart_xi2(i,j,k))
                  csxi2_d = COS(depart_xi2(i,j,k))

                  q_33         = snxi2_d*snxi2_p(j)+csxi2_d*csxi2_p(j)*csxi1_d
                  gam          = (1.0 + q_33)
                  temp_Csalpha = 1.0/gam
                  gam          = 0.5*gam*timestep

! Components of rotation matrix
                  rm_11 = (csxi2_d*csxi2_p(j) + (1.0+snxi2_p(j)*snxi2_d)       &
                          *csxi1_d)*temp_Csalpha
                  rm_12 =-snxi1_d*(snxi2_p(j) + snxi2_d)*temp_Csalpha

! Arrival points

                  u_a    = intw_u2p(i,1)*u_np1(i-1,j,k) +                      &
                           intw_u2p(i,2)*u_np1(i,j,k)
                  v_a    = intw_v2p(j,1)*v_np1(i,j-1,k) +                      &
                           intw_v2p(j,2)*v_np1(i,j,k)
                  w_a    = intw_w2rho(k,1)*etadot_np1(i,j,k) +                 &
                           intw_w2rho(k,2)*etadot_np1(i,j,k-1)

                 rad    = 1.0/(Earth_radius                                    &
                                   + Is*(xi3_at_rho(i,j,k) - Earth_radius))

                  u_a    = u_a *rad
                  v_a    = v_a *rad

! Rotated velocities at departure points

                  u_d = rm_11*int_wind_u(i,j,k) + rm_12*int_wind_v(i,j,k)
                  v_d =-rm_12*int_wind_u(i,j,k) + rm_11*int_wind_v(i,j,k)
                  w_d = int_wind_w(i,j,k)

! Update departure points

                  x_d = -gam*(alpha*u_a + beta*u_d)
                  y_d = -gam*(alpha*v_a + beta*v_d)
                  z_d = sqrt(1.0 - x_d**2 - y_d**2)

                  depart_xi1(i,j,k) = xi1_p(i) +                               &
                                      ATAN2(x_d,z_d*csxi2_p(j)-y_d*snxi2_p(j))
                  depart_xi2(i,j,k) = ASIN(y_d*csxi2_p(j)+z_d*snxi2_p(j))
                  depart_xi3(i,j,k) = eta_rho_levels(k)                        &
                                     -timestep*(alpha*w_a + beta*w_d)

               END DO
            END DO
         END DO
!$OMP END PARALLEL DO                                      

   END SELECT


   CALL eg_check_sl_domain( model_domain, depart_xi2, depart_xi1,           &
                         dep_row_len, dep_rows, model_levels+off_k,         &
                         max_xi1, min_xi1, max_xi2, min_xi2)

   IF(.NOT.interp_dpt_pt .OR. kk < dep_its) THEN
      CALL eg_adjust_vert_bound(                                            &
               row_length, rows, model_levels, g_i_pe,                      &
               pnt_type, dep_row_len,                                       &
               dep_rows, off_i, off_j, off_k, offz,                         &
               etadot, etadot_np1, depart_xi1, depart_xi2,                  &
               depart_xi3)
   END IF

END DO ! dep_its iterations

IF (lhook) CALL dr_hook('DEPARTURE_POINT_ETA',zhook_out,zhook_handle)

RETURN
END SUBROUTINE DEPARTURE_POINT_ETA
END MODULE
