! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine Calc_PV

      Subroutine Calc_PV                                                &
     &                  (u, v, theta, rho,                              &
     &                   r_theta_levels, r_rho_levels,                  &
     &                   r_at_u, r_at_v,                                &
     &                   sec_v_latitude, tan_v_latitude,                &
     &                   sec_theta_latitude, f3_at_v,                   &
     &                   delta_lambda, delta_phi,                       &
     &                   row_length, rows, n_rows, model_levels,        &
     &                   off_x, off_y, halo_i, halo_j,                  &
     &                   at_extremity,                                  &
     &                   pv)

! Description:
!          Calculates potential vorticity at PV points; ie. midway
!          horizontally between v points on rho levels.
!
! Method:
!          Discretisation of equation 7 in UMDP 13.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE atm_fields_bounds_mod, ONLY:                                  &
          udims, vdims, udims_s, vdims_s, tdims_s, pdims_s,             &
          udims_l, vdims_l

      USE Field_Types
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Domain_params
      USE um_input_control_mod,  ONLY: model_domain
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
        at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Parameters
      Integer                                                           &
         PNorth,                                                        &
                      ! North processor address in the neighbor array
         PEast,                                                         &
                      ! East processor address in the neighbor array
         PSouth,                                                        &
                      ! South processor address in the neighbor array
         PWest,                                                         &
                      ! West processor address in the neighbor array
         NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
         PNorth   = 1,                                                  &
         PEast    = 2,                                                  &
         PSouth   = 3,                                                  &
         PWest    = 4,                                                  &
         NoDomain = -1)

      Integer                                                           &
        row_length                                                      &
                        ! number of points on a row
      , rows                                                            &
                        ! number of rows of data
      , n_rows                                                          &
                        ! number of rows of data on a v row
      , model_levels                                                    &
                        ! number of levels of data
      , off_x                                                           &
      , off_y                                                           &
      , halo_i                                                          &
      , halo_j

      Real                                                              &
        delta_lambda                                                    &
      , delta_phi

      Real                                                              &
         theta(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)                           &
        ,rho(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end,                             &
             pdims_s%k_start:pdims_s%k_end)                             &  
        ,u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end, &
           udims_s%k_start:udims_s%k_end)                               &
        ,v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end, &
           vdims_s%k_start:vdims_s%k_end)                               &
      , r_theta_levels(1-halo_i:row_length+halo_i,                      &
                         1-halo_j:rows+halo_j,0:model_levels)           &
      , r_rho_levels(1-halo_i:row_length+halo_i,                        &
                       1-halo_j:rows+halo_j, model_levels)              &
      ,r_at_u (udims_l%i_start:udims_l%i_end,udims_l%j_start:udims_l%j_end, &
           udims_l%k_start:udims_l%k_end)                               &
      ,r_at_v (vdims_l%i_start:vdims_l%i_end,vdims_l%j_start:vdims_l%j_end, &
           vdims_l%k_start:vdims_l%k_end)

      Real                                                              &
        sec_v_latitude (vdims_s%i_start:vdims_s%i_end,                  &
                        vdims_s%j_start:vdims_s%j_end)                  &
      , tan_v_latitude(vdims%i_start:vdims%i_end,                       &
                       vdims%j_start:vdims%j_end)                       &
      , sec_theta_latitude (tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end)

      Real                                                              &
           ! Coriolis term
        f3_at_v (vdims_s%i_start:vdims_s%i_end,                         &
                 vdims_s%j_start:vdims_s%j_end)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
        pv (udims%i_start:udims%i_end,                                  &
            vdims%j_start:vdims%j_end, model_levels)

! Local variables

      Integer                                                           &
       i, j, k

      Real                                                              &
        recip_delta_lambda                                              &
      , recip_delta_phi                                                 &
      , weight1                                                         &
      , weight2

      Real                                                              &
        r_at_pv(udims%i_start:udims%i_end,                              &
                vdims%j_start:vdims%j_end, model_levels)                &
      , r_at_pv_on_theta_levels(udims%i_start:udims%i_end,              &
                                vdims%j_start:vdims%j_end,model_levels) &
      , dtheta_dr(tdims_s%i_start:tdims_s%i_end,                        &
                  tdims_s%j_start:tdims_s%j_end)                        &
      , du_dr(udims_s%i_start:udims_s%i_end,                            &
              udims_s%j_start:udims_s%j_end)                            &
      , dv_dr(vdims_s%i_start:vdims_s%i_end,                            &
              vdims_s%j_start:vdims_s%j_end)                            &
      , dtheta_dx(udims_s%i_start:udims_s%i_end,                        &
                  udims_s%j_start:udims_s%j_end)                        &
      , dtheta_dx_ave(udims%i_start:udims%i_end,                        &
                      vdims%j_start:vdims%j_end,model_levels)           &
      , du_dr_ave(udims%i_start:udims%i_end,                            &
                      vdims%j_start:vdims%j_end)                        &
      , dv_dr_ave(udims%i_start:udims%i_end,                            &
                      vdims%j_start:vdims%j_end)                        &
      , dtheta_dy(vdims_s%i_start:vdims_s%i_end,                        &
                  vdims_s%j_start:vdims_s%j_end)                        &
      , dtheta_dy_ave(udims%i_start:udims%i_end,                        &
                      vdims%j_start:vdims%j_end,model_levels)           &
      , dtheta_dr_ave(udims%i_start:udims%i_end,                        &
                      vdims%j_start:vdims%j_end)                        &
      , x_term(udims%i_start:udims%i_end,                               &
               vdims%j_start:vdims%j_end)                               &
      , y_term(udims%i_start:udims%i_end,                               &
               vdims%j_start:vdims%j_end)                               &
      , z_term(udims%i_start:udims%i_end,                               &
               vdims%j_start:vdims%j_end)                               &
      , r_at_u_on_theta_levels(udims_s%i_start:udims_s%i_end,           &
                               udims_s%j_start:udims_s%j_end)           &
      , density(pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end)                          &
      , density_ave(udims%i_start:udims%i_end,                          &
                    vdims%j_start:vdims%j_end)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 0. Initialisation
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CALC_PV',zhook_in,zhook_handle)
      recip_delta_lambda = 1./ delta_lambda
      recip_delta_phi = 1./ delta_phi

! ----------------------------------------------------------------------
! Section 1. Calculate horizontal part of terms
! ----------------------------------------------------------------------

! calculate dtheta/dx and dtheta/dy
      DO k = 1, model_levels
        DO j = udims_s%j_start,udims_s%j_end
          DO i = udims%i_start,udims%i_end
            dtheta_dx(i,j) = (theta(i+1,j,k) - theta(i,j,k)) *          &
                                recip_delta_lambda *                    &
                              sec_theta_latitude(i,j) * 2.0 /           &
                               (r_theta_levels(i+1,j,k) +               &
                                r_theta_levels(i,j,k) )           
          END DO
        END DO
        
        IF (.NOT. l_vatpoles) THEN

        IF (model_domain  /=  mt_bi_cyclic_LAM) THEN
          IF (at_extremity(PSouth)) THEN
            DO i = udims%i_start,udims%i_end
              dtheta_dx(i,udims%j_start) = 0.
            END DO
          END IF
          IF (at_extremity(PNorth)) THEN
            DO i = udims%i_start,udims%i_end
              dtheta_dx(i,udims%j_end) = 0.
            END DO
          END IF
        END IF
                
        END IF  ! vatpoles
  
                 
        DO j = vdims%j_start,vdims%j_end
          DO i = vdims_s%i_start,vdims_s%i_end 
            dtheta_dy(i,j) = (theta(i,j+1,k) - theta(i,j,k)) *          &
                                recip_delta_phi * 2. /                  &
                               (r_theta_levels(i,j+1,k) +               &
                                r_theta_levels(i,j,k) )
          END DO
        END DO

        IF ( l_vatpoles ) THEN

        IF (model_domain  /=  mt_bi_cyclic_LAM) THEN
          IF (at_extremity(PSouth)) THEN
            DO i = vdims%i_start,vdims%i_end 
              dtheta_dy(i,vdims%j_start) = 0.
            END DO
          END IF
          IF (at_extremity(PNorth)) THEN
            DO i = vdims%i_start,vdims%i_end 
              dtheta_dy(i,vdims%j_end) = 0.
            END DO
          END IF
        END IF
                
        END IF  ! vatpoles



        DO j = vdims%j_start,vdims%j_end
          DO i = udims%i_start,udims%i_end
            dtheta_dx_ave(i,j,k) = 0.5 * (dtheta_dx(i,j+1) +            &
                                          dtheta_dx(i,j) )
          END DO
        END DO

        DO j = vdims%j_start,vdims%j_end
          DO i = udims%i_start,udims%i_end
            dtheta_dy_ave(i,j,k) = 0.5 * (dtheta_dy(i+1,j) +            &
                                          dtheta_dy(i,j) )
          END DO
        END DO

        DO j = udims_s%j_start,udims_s%j_end 
          DO i =udims%i_start,udims%i_end
            r_at_u_on_theta_levels(i,j) =                               &
                            (r_theta_levels(i+1,j,k) +                  &
                             r_theta_levels(i,j,k) ) * 0.5
          END DO
        END DO
                               
                       

        DO j = vdims%j_start,vdims%j_end
          DO i =udims%i_start,udims%i_end
            r_at_pv(i,j,k) = 0.5 * (r_at_u(i,j+1,k) + r_at_u(i,j,k))
            r_at_pv_on_theta_levels(i,j,k) = 0.5 *                      &
                                   ( r_at_u_on_theta_levels(i,j+1) +    &
                                     r_at_u_on_theta_levels(i,j) )
          END DO
        END DO

! end loop over levels
      END DO

! loop over model levels for rest of code.
      DO k = 1, model_levels

! ----------------------------------------------------------------------
! Section 1.1 Calculate x_term and y_term.
! ----------------------------------------------------------------------

        IF (k  ==  1) THEN
          DO j = vdims%j_start,vdims%j_end
            DO i = udims%i_start,udims%i_end
              x_term(i,j) = - dtheta_dx_ave(i,j,k)
              y_term(i,j) = dtheta_dy_ave(i,j,k)
            END DO
          END DO

        ELSE

          DO j = vdims%j_start,vdims%j_end
            DO i = udims%i_start,udims%i_end
              weight1 = (r_at_pv_on_theta_levels(i,j,k) -               &
                         r_at_pv(i,j,k) ) /                             &
                        (r_at_pv_on_theta_levels(i,j,k) -               &
                         r_at_pv_on_theta_levels(i,j,k-1) )
              weight2 = 1.0 - weight1

              x_term(i,j) = - ( weight2 * dtheta_dx_ave(i,j,k-1) +      &
                                weight1 * dtheta_dx_ave(i,j,k) )
              y_term(i,j) = ( weight2 * dtheta_dy_ave(i,j,k-1) +        &
                              weight1 * dtheta_dy_ave(i,j,k) )
            END DO
          END DO

        END IF

! ----------------------------------------------------------------------
! Section 1.2 Calculate z_term.
! ----------------------------------------------------------------------
        DO j = vdims%j_start,vdims%j_end  
          DO i = udims%i_start,udims%i_end
          
            z_term(i,j) = f3_at_v(i,j) + (                              &
                        (v(i+1,j,k)-v(i,j,k)) * recip_delta_lambda *    &
                         sec_v_latitude((i-udims%i_start+1),j) +        &
                         0.5 * (u(i,j+1,k) + u(i,j,k)) *                &
                         tan_v_latitude((i-udims%i_start+1),j) -        &
                         (u(i,j+1,k) - u(i,j,k)) * recip_delta_phi      &
                         ) / r_at_pv(i,j,k)
         
          END DO
        END DO


! ----------------------------------------------------------------------
! Section 2. Multiply horizontal terms by vertical terms and form
!            full pv.
! ----------------------------------------------------------------------

        IF (k  ==  1) THEN
! Since theta gradient is un-defined use gradient from levels above.
          DO j = tdims_s%j_start,tdims_s%j_end
            DO i = tdims_s%i_start,tdims_s%i_end
              dtheta_dr(i,j) = (theta(i,j,k+1) - theta(i,j,k)) /        &
                               (r_theta_levels(i,j,k+1) -               &
                                r_theta_levels(i,j,k) )
            END DO
          END DO

! Now average to required location.
          DO j = vdims%j_start,vdims%j_end
            DO i = udims%i_start,udims%i_end
              dtheta_dr_ave(i,j) = 0.25 * (dtheta_dr(i,j) +             &
                                           dtheta_dr(i+1,j) +           &
                                           dtheta_dr(i,j+1) +           &
                                           dtheta_dr(i+1,j+1) )
            END DO
          END DO

        ELSE IF (k  /=  2) THEN
! no code required at level 2 since the value is the same as level 1.

! calculate theta gradient.
          DO j = tdims_s%j_start,tdims_s%j_end
            DO i = tdims_s%i_start,tdims_s%i_end
              dtheta_dr(i,j) = (theta(i,j,k) - theta(i,j,k-1)) /        &
                               (r_theta_levels(i,j,k) -                 &
                                r_theta_levels(i,j,k-1) )
            END DO
          END DO

! Now average to required location.
          DO j = vdims%j_start,vdims%j_end
            DO i = udims%i_start,udims%i_end
              dtheta_dr_ave(i,j) = 0.25 * (dtheta_dr(i,j) +             &
                                           dtheta_dr(i+1,j) +           &
                                           dtheta_dr(i,j+1) +           &
                                           dtheta_dr(i+1,j+1) )
            END DO
          END DO

        END IF

        IF (k  ==  1) THEN
          DO j =  vdims_s%j_start,vdims_s%j_end
            DO i = vdims_s%i_start,vdims_s%i_end
              dv_dr(i,j) = (v(i,j,k+1) - v(i,j,k)) /                    &
                           (r_at_v(i,j,k+1) -                           &
                            r_at_v(i,j,k) )
            END DO
          END DO
          DO j = udims_s%j_start,udims_s%j_end 
            DO i = udims_s%i_start,udims_s%i_end
              du_dr(i,j) = (u(i,j,k+1) - u(i,j,k)) /                    &
                           (r_at_u(i,j,k+1) -                           &
                            r_at_u(i,j,k) )
            END DO
          END DO

        ELSE IF (k  ==  model_levels) THEN

          DO j = vdims_s%j_start,vdims_s%j_end
            DO i = vdims_s%i_start,vdims_s%i_end 
              dv_dr(i,j) = (v(i,j,k) - v(i,j,k-1)) /                    &
                           (r_at_v(i,j,k) -                             &
                            r_at_v(i,j,k-1) )
            END DO
          END DO
          DO j = udims_s%j_start,udims_s%j_end 
            DO i = udims_s%i_start,udims_s%i_end
              du_dr(i,j) = (u(i,j,k) - u(i,j,k-1)) /                    &
                           (r_at_u(i,j,k) -                             &
                            r_at_u(i,j,k-1) )
            END DO
          END DO

        ELSE

          DO j = vdims_s%j_start,vdims_s%j_end
            DO i = vdims_s%i_start,vdims_s%i_end 
              dv_dr(i,j) = (v(i,j,k+1) - v(i,j,k-1)) /                  &
                           (r_at_v(i,j,k+1) -                           &
                            r_at_v(i,j,k-1) )
            END DO
          END DO
          DO j = udims_s%j_start,udims_s%j_end 
            DO i = udims_s%i_start,udims_s%i_end
              du_dr(i,j) = (u(i,j,k+1) - u(i,j,k-1)) /                  &
                           (r_at_u(i,j,k+1) -                           &
                            r_at_u(i,j,k-1) )
            END DO
          END DO

        END IF

! now average quantities
        DO j = vdims%j_start,vdims%j_end
          DO i = udims%i_start,udims%i_end
            du_dr_ave(i,j) = 0.5 * (du_dr(i,j+1) + du_dr(i,j))
          END DO
        END DO

        DO j = vdims%j_start,vdims%j_end
          DO i = udims%i_start,udims%i_end
            dv_dr_ave(i,j) = 0.5 * (dv_dr(i+1,j) + dv_dr(i,j))
          END DO
        END DO

! convert rho to true density by removing factor of r squared.
        DO j = pdims_s%j_start,pdims_s%j_end
          DO i = pdims_s%i_start,pdims_s%i_end
            density(i,j) = rho(i,j,k) / (r_rho_levels(i,j,k) *          &
                                         r_rho_levels(i,j,k) )
          END DO
        END DO

        DO j = vdims%j_start,vdims%j_end 
          DO i =udims%i_start,udims%i_end
            density_ave(i,j) = 0.25 * (density(i,j+1) + density(i+1,j) +&
                                       density(i+1,j+1) + density(i,j) )
          END DO
        END DO

! Calculate full PV.
        DO j = vdims%j_start,vdims%j_end
          DO i = udims%i_start,udims%i_end
            pv(i,j,k) = (x_term(i,j) * dv_dr_ave(i,j) +                 &
                         y_term(i,j) * du_dr_ave(i,j) +                 &
                         z_term(i,j) * dtheta_dr_ave(i,j) )             &
                        / density_ave(i,j)
          END DO
        END DO

        IF ( l_vatpoles ) THEN
! until we decide how to appropriately resolve PV at the poles, set to zero.
        IF (model_domain  ==  mt_global) THEN
          IF (at_extremity(PSouth)) THEN
            DO i = udims%i_start,udims%i_end 
              pv(i,vdims%j_start,k) = 0.0
            END DO
          END IF
          IF (at_extremity(PNorth)) THEN
            DO i = udims%i_start,udims%i_end 
              pv(i,vdims%j_end,k) = 0.0
            END DO
          END IF
        END IF
                
        END IF  ! vatpoles


! end loop over model levels
      END DO

! end of routine

      IF (lhook) CALL dr_hook('CALC_PV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Calc_PV

