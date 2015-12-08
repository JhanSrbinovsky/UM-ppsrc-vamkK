! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE eg_diff_ctl
      SUBROUTINE eg_diff_ctl(                                           &
                             L_Backwards,                               &
                             timestep, pos_timestep, neg_timestep,      &
                             theta, moist, u, v, w, rho,                &
                             exner_rho_levels, exner_theta_levels,      &
                             r_theta_levels, r_rho_levels,              &
                             r_at_u, r_at_v,                            &
                             eta_theta_levels, eta_rho_levels,          &
                             offx, offy, halo_i, halo_j,                &
                             row_length, rows, n_rows,                  &
                             model_levels, flat_level,                  &
                             L_tardiff_q, w_conv_limit,                 &
                             tardiffq_factor, tardiffq_test,            &
                             tardiffq_start, tardiffq_end,              &
                             theta_star, moist_star, S_u, S_v, S_w )

! Purpose:
!          Subroutine to interface to diffusion code
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

      USE horiz_grid_mod, ONLY: xi1_p, xi1_u, xi2_p, xi2_v,             &
                                csxi2_p, csxi2_v
      USE turb_diff_mod, ONLY: L_subfilter_horiz, l_subfilter_vert,     &
          turb_startlev_horiz,turb_endlev_horiz,                        &
          turb_startlev_vert,turb_endlev_vert
      USE turb_diff_ctl_mod, ONLY: visc_m, visc_h, max_diff
      USE proc_info_mod, ONLY: model_domain

      USE swapable_field_mod, ONLY : swapable_field_pointer_type

      USE UM_ParVars, ONLY: at_extremity, PSouth,                       &
                            fld_type_p, fld_type_u, fld_type_v
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      LOGICAL                                                           &
        L_vertical_diffusion                                            &
      , L_tardiff_q                                                     &
      , L_Backwards

      INTEGER                                                           &
        row_length                                                      &
                         ! number of point on a row.
      , rows                                                            &
                         ! number of rows.
      , n_rows                                                          &
                         ! number of v rows.
      , model_levels                                                    &
                         ! number of model levels.
      , flat_level                                                      &
                         ! first flat level
      , offx                                                            &
                     ! Size of small halo in i
      , offy                                                            &
                     ! Size of small halo in j.
      , halo_i                                                          &
                       ! Size of halo in i direction.
      , halo_j         ! Size of halo in j direction.

      INTEGER                                                           &
        tardiffq_test                                                   &
      , tardiffq_start                                                  &
      , tardiffq_end

      REAL                                                              &
       timestep                                                         &
                            ! atmosphere model timestep
      ,pos_timestep                                                     &
                    ! = +timestep.
      ,neg_timestep                                                     &
                    ! = -timestep.
      ,w_conv_limit

      REAL                                                              &
           ! primary model variables
        u(0:row_length+1, 0:rows+1, model_levels)               &
      , v(0:row_length+1, 0:n_rows+1, model_levels)             &
      , w(0:row_length+1, 0:rows+1, 0:model_levels)             &
      , rho(0:row_length+1, 0:rows+1, model_levels)             &
      , exner_rho_levels(0:row_length+1,                        &
                         0:rows+1, model_levels)                &
      , exner_theta_levels(0:row_length+1,                      &
                           0:rows+1, model_levels)              &
      , moist (1-offx:row_length+offx,                          &
               1-offy:rows+offy, 0:model_levels)                &
      , theta (0:row_length+1, 0:rows+1, 0:model_levels)

      REAL                                                              &
           ! vertical co-ordinate arrays.
        r_theta_levels (1-halo_i:row_length+halo_i,                     &
                        1-halo_j:rows+halo_j, 0:model_levels)           &
      , r_rho_levels (1-halo_i:row_length+halo_i,                       &
                      1-halo_j:rows+halo_j, model_levels)               &
      , r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
                model_levels)                                           &
      , r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
                model_levels)                                           &
      , eta_theta_levels (0:model_levels)                               &
      , eta_rho_levels (model_levels)


      REAL :: tardiffq_factor  ! targeted diffusion coefficient


      REAL, ALLOCATABLE ::                                              &
       coeff_u(:,:,:)                                               &
                             ! visc_m interpolated onto
      ,coeff_v (:,:,:)                                              &
                             ! appropriate points
      ,coeff_th(:,:,:)                                              &
                             ! for use in the turb_diff_*
      ,work_th(:,:,:)                                              &
                             ! for use in the turb_diff_*
      ,coeff_centre(:,:,:)                                          &
                             ! diffusion routines.                      
      ,w_coeff_u(:,:,:)                                             &
                             ! horizontal diffusion coefficients for    
      ,w_coeff_v(:,:,:)        ! diffusion of w

! local
      INTEGER :: i, j, k, l
      INTEGER :: info

      REAL                                                              &
       vert_diffusion_coeff_wind                                        &
      ,vert_diffusion_coeff_theta                                       &
      ,vert_diffusion_coeff_q                                           &
      , vdiffuv_factor                                                  &
      , ramp_lat_radians

!   variables with intent INOUT
!   For ENDGame, these are Fast physics fields
      REAL                                                              &
        S_u(0:row_length+1, 0:rows+1, model_levels)                     &
      , S_v(0:row_length+1, 0:n_rows+1, model_levels)                   &
      , S_w(row_length, rows, 0:model_levels)                           &
      , moist_star(0:row_length+1, 0:rows+1, 0:model_levels)            &
      , theta_star(0:row_length+1, 0:rows+1, 0:model_levels)

! Local variables
      INTEGER :: levels
      INTEGER :: i_start, i_stop    ! i start/end points
      INTEGER :: j_start, j_stop    ! j start/end points

      INTEGER :: i_field

      REAL                                                              &
! Interpolation weights
        w1(0:row_length+1)                                              &
     ,  w2(0:rows+1)                                                    &
     ,  weight ! vertical weight

      REAL ::  delta_z(0:row_length+1, 0:rows+1, model_levels)
      REAL ::  recip_r_squared_delz(row_length, rows, model_levels)
      REAL ::  r_theta_uv(0:row_length+1, 0:rows+1, flat_level)

      TYPE(swapable_field_pointer_type) :: fields_to_swap(2)
                                           ! multivariate swapbounds
   
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 1.0  Horizontal Diffusion
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('EG_DIFF_CTL',zhook_in,zhook_handle)

      i_start = 1
      i_stop = row_length
      j_start = 1
      j_stop = rows
      IF (at_extremity(PSouth)) j_start = 2
      
      IF ( L_subfilter_horiz) THEN

        ALLOCATE( coeff_u(0:row_length+1, 0:rows+1, model_levels))
        ALLOCATE( coeff_v(0:row_length+1, 0:rows+1, model_levels))
        ALLOCATE( coeff_th(0:row_length+1, 0:rows+1, model_levels))
        ALLOCATE( work_th(0:row_length+2, 0:rows+2, model_levels))
        ALLOCATE( coeff_centre(0:row_length+1, 0:rows+1, model_levels))
        ALLOCATE( w_coeff_u(0:row_length+1, 0:rows+1, model_levels) )
        ALLOCATE( w_coeff_v(0:row_length+1, 0:rows+1, model_levels) ) 
        coeff_u(:,:,:) = 0.0
        coeff_v(:,:,:) = 0.0
        coeff_th(:,:,:) = 0.0    
        work_th(:,:,:) = 0.0    
        coeff_centre(:,:,:) = 0.0
        w_coeff_u(:,:,:) = 0.0
        w_coeff_v(:,:,:) = 0.0

    DO k = 1, model_levels - 2

      IF (k >= turb_startlev_horiz .AND.                        &
        k <= turb_endlev_horiz) THEN

        DO j = 1, rows
          DO i = 1, row_length
            visc_h(i,j,k) = MIN(visc_h(i,j,k), max_diff(i,j))
            visc_m(i,j,k) = MIN(visc_m(i,j,k), max_diff(i,j))
          END DO
        END DO

      ELSE

        DO j = 1, rows
          DO i = 1, row_length
            visc_h(i,j,k) = 0.0
            visc_m(i,j,k) = 0.0
          END DO
        END DO

      END IF

    END DO

  ! Set visc_m and visc_h at the top two levels

    DO k = model_levels - 1, model_levels

      IF (turb_startlev_horiz <= k .AND. &
                  turb_endlev_horiz   >= k) THEN

        DO j = 1, rows
          DO i = 1, row_length
            visc_h(i,j,k) = visc_h(i,j,model_levels-2)
            visc_m(i,j,k) = visc_m(i,j,model_levels-2)
          END DO
        END DO

      ELSE

        DO j = 1, rows
          DO i = 1, row_length
            visc_h(i,j,k) = 0.0
            visc_m(i,j,k) = 0.0
          END DO
        END DO

      END IF

    END DO

    i_field = 0

    i_field = i_field + 1
    fields_to_swap(i_field) % field      => visc_m(:,:,:)
    fields_to_swap(i_field) % field_type =  fld_type_p
    fields_to_swap(i_field) % levels     =  model_levels
    fields_to_swap(i_field) % rows       =  rows
    fields_to_swap(i_field) % vector     =  .FALSE.

    i_field = i_field + 1
    fields_to_swap(i_field) % field      => visc_h(:,:,:)
    fields_to_swap(i_field) % field_type =  fld_type_p
    fields_to_swap(i_field) % levels     =  model_levels
    fields_to_swap(i_field) % rows       =  rows
    fields_to_swap(i_field) % vector     =  .FALSE.

  ! DEPENDS ON: swap_bounds_mv
    CALL swap_bounds_mv( fields_to_swap, i_field,                   &
                                 row_length, halo_i, halo_j)

    ELSE  IF (l_subfilter_vert) THEN

    DO k = 1, model_levels - 2

      IF (k < turb_startlev_vert .OR. k > turb_endlev_vert) THEN
        DO j = 1, rows
          DO i = 1, row_length
            visc_h(i,j,k) = 0.0
            visc_m(i,j,k) = 0.0
          END DO
        END DO
      END IF

    END DO

  ! Set visc_m and visc_h at the top two levels

    DO k = model_levels - 1, model_levels
      IF (turb_startlev_vert <= k .AND. turb_endlev_vert >= k) THEN

        DO j = 1, rows
          DO i = 1, row_length
            visc_h(i,j,k) = visc_h(i,j,model_levels-2)
            visc_m(i,j,k) = visc_m(i,j,model_levels-2)
          END DO
        END DO
      ELSE
        DO j = 1, rows
          DO i = 1, row_length
            visc_h(i,j,k) = 0.0
            visc_m(i,j,k) = 0.0
          END DO
        END DO
      END IF
    END DO

    i_field = 0

    i_field = i_field + 1
    fields_to_swap(i_field) % field      => visc_m(:,:,:)
    fields_to_swap(i_field) % field_type =  fld_type_p
    fields_to_swap(i_field) % levels     =  model_levels
    fields_to_swap(i_field) % rows       =  rows
    fields_to_swap(i_field) % vector     =  .FALSE.

    i_field = i_field + 1
    fields_to_swap(i_field) % field      => visc_h(:,:,:)
    fields_to_swap(i_field) % field_type =  fld_type_p
    fields_to_swap(i_field) % levels     =  model_levels
    fields_to_swap(i_field) % rows       =  rows
    fields_to_swap(i_field) % vector     =  .FALSE.

! DEPENDS ON: swap_bounds_mv
    CALL swap_bounds_mv( fields_to_swap, i_field,                   &
                                 row_length, halo_i, halo_j)

      END IF !  L_subfilter_horiz

       IF ( L_subfilter_horiz) THEN

!        Run with a positive timestep IF integrating backwards.
         IF (L_Backwards) timestep = pos_timestep
!
!Interpolate visc_m onto appropriate points for the diffusion routines.
! At this point the diffusion coefficients are on w (theta) points.
!
! Horizontal weights
!
         DO i = 0, row_length+1
           w1(i) = (xi1_p(i+1) - xi1_u(i)) / (xi1_p(i+1) - xi1_p(i))
         END DO

         DO j = 0, rows+1
           w2(j) = (xi2_p(j+1) - xi2_v(j)) / (xi2_p(j+1) - xi2_p(j))
         END DO

         DO j = 0, rows + 2
           DO i = 0, row_length + 2
             work_th(i,j,1) = visc_m(i,j,1)
           END DO
         END DO

         DO k = 2, model_levels
           DO j = 0, rows + 2
             DO i = 0, row_length + 2
               weight = ( r_rho_levels(i,j,k) -                         &
                           r_theta_levels(i,j,k-1) ) /                  &
                       (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
               work_th(i,j,k) = ( weight * visc_m(i,j,k) +              &
                                  (1. - weight) * visc_m(i,j,k-1) )
             END DO
           END DO
         END DO ! k = 2, model_levels

         DO k = 1, model_levels

           DO j = 1, rows          
             DO i = 0, row_length + 1
               coeff_u(i,j,k) = w1(i) * visc_h(i,j,k) +                 &
                                 (1. - w1(i)) * visc_h(i+1,j,k)
               w_coeff_u(i,j,k) = w1(i) * visc_m(i,j,k) +               &
                                 (1. - w1(i)) * visc_m(i+1,j,k)     
             END DO
           END DO

           DO j = 0, rows + 1
             DO i = 1, row_length
               coeff_v(i,j,k) = w2(j) * visc_h(i,j,k) +                 &
                                (1. - w2(j)) * visc_h(i,j+1,k)
               w_coeff_v(i,j,k) = w2(j) * visc_m(i,j,k) +               &
                                    (1. - w2(j)) * visc_m(i,j+1,k)  
             END DO
           END DO

           DO j = 0, rows + 1
             DO i = 0, row_length + 1
               coeff_centre(i,j,k) = w2(j) *                            &
                                         ( w1(i) * work_th(i,j,k) +     &
                                    (1. - w1(i)) * work_th(i+1,j,k) )   &
                                              + (1. - w2(j)) *          &
                                        ( w1(i) * work_th(i,j+1,k) +    &
                                   (1. - w1(i)) * work_th(i+1,j+1,k) )
               coeff_th(i,j,k) = work_th(i,j,k)
             END DO
           END DO

         END DO !  k = 1, model_levels

! calculate dr/dz about theta levels
      DO k = 1, flat_level - 1
        DO j = 0, rows + 1
          DO i = 0, row_length + 1
            delta_z(i,j,k) = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
          END DO
        END DO
      END DO !   k = 1, flat_level - 1
      DO k = flat_level, model_levels
        DO j = 0, rows + 1
          DO i = 0, row_length + 1
            delta_z(i,j,k) = 1.0
          END DO
        END DO
      END DO !  k = flat_level, model_levels
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at flat levels and then no need to take special action in later loops
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            recip_r_squared_delz(i,j,k) = 1.0 / ( r_theta_levels(i,j,k) &
                                              * r_theta_levels(i,j,k) * &
                                                       delta_z(i,j,k) ) 
          END DO
        END DO
      END DO !   k = 1, model_levels

!  CALL diffusion routines with (3D) coefficients from the subgrid
!  turbulence scheme.

! DEPENDS ON: eg_turb_diff
           CALL eg_turb_diff(                                           &
                             theta(0,0,1), fld_type_p,                  &
                             row_length, rows, rows, n_rows,            &
                             model_levels, model_levels,                &
                             offx, offy, halo_i, halo_j,                &
                             offx, offy, offx, offy,                    &
                             xi1_p, xi1_u, xi2_p, xi2_v,                &
                             csxi2_p, csxi2_v,                          &
                             coeff_u, coeff_v,                          &
                             delta_z, recip_r_squared_delz ,            &
                             1, rows, .FALSE.,                          &
                             theta_star(1-offx, 1-offy,1) )

!DEPENDS ON: eg_turb_diff
           CALL eg_turb_diff(                                           &
                             moist(1-offx,1-offy,1), fld_type_p,        &
                             row_length, rows, rows, n_rows,            &
                             model_levels, model_levels,                &
                             offx, offy, halo_i, halo_j,                &
                             offx, offy, offx, offy,                    &
                             xi1_p, xi1_u, xi2_p, xi2_v,                &
                             csxi2_p, csxi2_v,                          &
                             coeff_u, coeff_v,                          &
                             delta_z, recip_r_squared_delz ,            &
                             1, rows, .FALSE.,                          &
                             moist_star(1-offx, 1-offy,1) )

! DEPENDS ON: eg_turb_diff
           CALL eg_turb_diff(                                           &
                             w(0,0,1), fld_type_p,                      &
                             row_length, rows, rows, n_rows,            &
                             model_levels, model_levels,                &
                             offx, offy, halo_i, halo_j,                &
                             0, 0, offx, offy,                          &
                             xi1_p, xi1_u, xi2_p, xi2_v,                &
                             csxi2_p, csxi2_v,                          &
                             w_coeff_u, w_coeff_v,                      &
                             delta_z, recip_r_squared_delz ,            &
                             1, rows, .FALSE., S_w(1,1,1) )

! calculate dr/dz about theta levels
      DO k = 1, flat_level
        DO j = 0, rows + 1
          DO i = 0, row_length + 1
            r_theta_uv(i,j,k) = w1(i) * r_theta_levels(i,j,k) +         &
                               (1. - w1(i))* r_theta_levels(i+1,j,k)
          END DO
        END DO
      END DO !  k = 1, flat_level - 1
      DO k = 1, flat_level - 1
        DO j = 0, rows + 1
          DO i = 0, row_length + 1
            delta_z(i,j,k) = r_theta_uv(i,j,k+1) - r_theta_uv(i,j,k)
          END DO
        END DO
      END DO !   k = 1, flat_level - 1
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at flat levels and then no need to take special action in later loops
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            recip_r_squared_delz(i,j,k) = 1.0 / ( r_at_u(i,j,k) *       &
                                                  r_at_u(i,j,k) *       &
                                                 delta_z(i,j,k) ) 
          END DO
        END DO
      END DO !   k = 1, model_levels

! DEPENDS ON: eg_turb_diff
           CALL eg_turb_diff(                                           &
                             u, fld_type_u,                             &
                             row_length, rows, rows, n_rows,            &
                             model_levels, model_levels,                &
                             offx, offy, halo_i, halo_j,                &
                             offx, offy, offx, offy,                    &
                             xi1_u, xi1_p, xi2_p, xi2_v,                &
                             csxi2_p, csxi2_v,                          &
                             coeff_th, coeff_centre,                    &
                             delta_z, recip_r_squared_delz,             &
                             1, rows, .TRUE., S_u )

! calculate dr/dz about theta levels
      DO k = 1, flat_level
        DO j = 0, rows + 1
          DO i = 0, row_length + 1
            r_theta_uv(i,j,k) = w2(j) * r_theta_levels(i,j,k) +         &
                               (1. - w2(j)) * r_theta_levels(i,j+1,k)
          END DO
        END DO
      END DO !  k = 1, flat_level - 1
      DO k = 1, flat_level - 1
        DO j = 0, rows + 1
          DO i = 0, row_length + 1
            delta_z(i,j,k) = r_theta_uv(i,j,k+1) - r_theta_uv(i,j,k)
          END DO
        END DO
      END DO !   k = 1, flat_level - 1
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at flat levels and then no need to take special action in later loops
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            recip_r_squared_delz(i,j,k) = 1.0 / ( r_at_v(i,j,k) *       &
                                                  r_at_v(i,j,k) *       &
                                                 delta_z(i,j,k) ) 
          END DO
        END DO
      END DO !   k = 1, model_levels

! DEPENDS ON: eg_turb_diff
           CALL eg_turb_diff(                                           &
                             v, fld_type_v,                             &
                             row_length, rows, n_rows, rows,            &
                             model_levels, model_levels,                &
                             offx, offy, halo_i, halo_j,                &
                             offx, offy, offx, offy,                    &
                             xi1_p, xi1_u, xi2_v, xi2_p,                &
                             csxi2_v, csxi2_p,                          &
                             coeff_centre, coeff_th,                    &
                             delta_z, recip_r_squared_delz,             &
                             j_start, j_stop, .TRUE., S_v )

!        Go back to negative timestep IF integrating backwards.
         IF (L_Backwards) timestep = neg_timestep

       END IF  !  L_subfilter_horiz

! ----------------------------------------------------------------------
! Section 2.0  Vertical Diffusion
! ----------------------------------------------------------------------


       IF (L_subfilter_horiz) THEN
         DEALLOCATE (coeff_u)
         DEALLOCATE (coeff_v)
         DEALLOCATE (coeff_th)
         DEALLOCATE (work_th)
         DEALLOCATE (coeff_centre)
         DEALLOCATE (w_coeff_u)
         DEALLOCATE (w_coeff_v) 
       END IF

      IF (lhook) CALL dr_hook('EG_DIFF_CTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE eg_diff_ctl
