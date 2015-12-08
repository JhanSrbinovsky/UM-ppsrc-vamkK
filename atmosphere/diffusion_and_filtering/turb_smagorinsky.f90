! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine turb_Smagorinsky
!
      SUBROUTINE turb_Smagorinsky(                                      &
                                  u, v, w, z0, timestep,                &
                                  row_length, rows, n_rows,             &
                                  model_levels, bl_levels,              &
                                  r_theta_levels, r_rho_levels,         &
                                  cos_theta_latitude,                   &
                                  cos_v_latitude,                       &
                                  delta_lambda, delta_phi )

      USE turb_diff_mod, ONLY: L_subfilter_blend, diff_factor, mix_factor
      USE turb_diff_ctl_mod, ONLY: visc_m, visc_h, rneutml_sq, max_diff,&
                                   shear, delta_smag
      USE atmos_constants_mod, ONLY: vkman
      USE earth_constants_mod, ONLY: earth_radius
      USE proc_info_mod, ONLY: global_row_length, global_rows,          &
                               model_domain, me, n_proc
      USE global_2d_sums_mod, ONLY: global_2d_sums

! Description: Calculates coefficients for use in subgrid turbulence
!              scheme
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      IMPLICIT NONE

! halo information

! Variables with Intent (In)

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, bl_levels        ! number of bl levels.

      Real                                                              &
     &  u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)     &
     &, z0(row_length,rows)                                             &
                             ! roughness length
     &, timestep                                                        &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                    1-halo_j:rows+halo_j, 0:model_levels)         &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                    1-halo_j:rows+halo_j, model_levels)           &
     &, delta_lambda                                                    &
     &, delta_phi                                                       &
     &, cos_theta_latitude (1-Offx:row_length+Offx,1-Offy:rows+Offy)    &
     &, cos_v_latitude (1-Offx:row_length+Offx,1-Offy:n_rows+Offy)

!Local parameters

      Integer                                                           &
        i, j, k                                                         &
                                  ! Loop indices
      , j_start, j_stop

      Real                                                              &
     &  delta_x                                                         &
        ! Sea level longitude spacing in m
     &, delta_y                                                         &
        ! Sea level latitude spacing in m
     &, delta_z(row_length, rows, model_levels)                         &
        ! Backward difference dz at rho points
     &, delta_zn(0:row_length+1,0:rows+1, model_levels)                 &
        ! Backward difference dz at theta points
        ! Note: k=1 is special - uses surface
        ! Note: delta_zn(k) is correct for theta(k-1)
     &, delta_zn_u                                                      &
        ! Average dz at theta points to u points (theta levels)
     &, delta_zn_v                                                      &
        ! Average dz at theta points to v points (theta levels)
     &, rdz(row_length, rows, model_levels)                             &
        ! 1/delta_z
     &, rdzn_u(0:row_length, rows, model_levels)                        &
        ! 1/delta_zn_u
     &, rdzn_v(row_length, 0:rows, model_levels)                        &
        ! 1/delta_zn_v
     &, smallp                                                          &
                          ! A small number
     &, rneutml           ! neutral mixing length

      Real                                                              &
     &  dx_rho                                                          &
!               ! dx on rho points
     &, dx_theta_u                                                      &
!               ! dx above u points on theta levels
     &, dx_rho_centre                                                   &
!               ! dx in centre of grid box on rho levels
     &, dy_rho                                                          &
!               ! dy on rho points
     &, dy_theta_v                                                      &
!               ! dy above v points on theta levels
     &, dy_rho_centre                                                   &
!               ! dy in centre of grid box on rho levels
     &, cx_rho                                                          &
!               ! reciprocal of dx_rho
     &, cy_rho                                                          &
!               ! reciprocal of dy_rho
     &, cx_theta_u(0:row_length, rows, model_levels)                    &
!               ! reciprocal of dx_theta_u
     &, cx_rho_centre(0:row_length, 0:rows, model_levels)               &
!               ! reciprocal of dx_rho_centre
     &, cy_theta_v(row_length, 0:rows, model_levels)                    &
!               ! reciprocal of dy_theta_v
     &, cy_rho_centre(0:row_length, 0:rows, model_levels)               &
!               ! reciprocal of dy_rho_centre
     &, cx2_rho(row_length, rows, model_levels)                         &
!               ! square of cx_rho
     &, cy2_rho(row_length, rows, model_levels)                         &
!               ! square of cy_rho
     &, rmlmax(row_length, rows)                                        &
!               ! basic mixing length (lambda0)
     &, z(row_length, rows, model_levels)

      Parameter (smallp=1.e-14)

      Real                                                              &
     &  weight_pl                                                       &
     &, weight_min                                                      &
     &, weight_2pl                                                      &
     &, weight_2min                                                     &
     &, ssq11                                                           &
                      ! ii component of shear stress
     &, ssq22                                                           &
                      ! jj component of shear stress
     &, ssq33                                                           &
                      ! kk component of shear stress
     &, ssq13                                                           &
                      ! ik component of shear stress
     &, ssq23                                                           &
                      ! jk component of shear stress
     &, ssq12                                                           &
                      ! ij component of shear stress
     &, sum(row_length, rows, model_levels-1)   ! sum of Sij (ssqij)

       REAL                                                             &
        sum_s(2*model_levels)                                         &
      , sum_n(2*model_levels)                                         &
      , temp_s(row_length, 2*model_levels)                            &
      , temp_n(row_length, 2*model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
     
!----------------------------------------------------------------------
! Initialise variables
!----------------------------------------------------------------------
     IF (lhook) CALL dr_hook('TURB_SMAGORINSKY',zhook_in,zhook_handle)

     j_start = 1
     j_stop = rows
     IF (model_domain  ==  mt_global ) THEN
       IF (at_extremity(PSouth)) j_start = 2
       IF (at_extremity(PNorth)) j_stop = rows - 1
     END IF

     Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            z(i,j,k) = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
          End Do
        End Do
     End Do

!----------------------------------------------------------------------
! Calculate grid spacings
!----------------------------------------------------------------------

      Do k = 1, model_levels
      
        Do j = j_start, j_stop
          Do i = 0, row_length
            dx_theta_u =                                                &
     &         0.5*(r_theta_levels(i,j,k) + r_theta_levels(i+1,j,k))    &
     &       * delta_lambda * cos_theta_latitude(i,j)
            cx_theta_u(i,j,k) =  1. / dx_theta_u
          End Do
        End Do
        
        Do j = j_start - 1, j_stop

          Do i = 1, row_length
            dy_theta_v =                                                &
     &         0.5*(r_theta_levels(i,j,k) + r_theta_levels(i,j+1,k))    &
     &       * delta_phi
            cy_theta_v(i,j,k) =  1. / dy_theta_v           
          End Do

          Do i = 0, row_length
            dx_rho_centre = 0.25*(                                      &
     &         r_rho_levels(i,j+1,k) + r_rho_levels(i+1,j+1,k)          &
     &        +r_rho_levels(i,j,k) + r_rho_levels(i+1,j,k) )            &
     &                * delta_lambda * cos_v_latitude(i,j)
            cx_rho_centre(i,j,k) =  1. / dx_rho_centre            
            dy_rho_centre = 0.25*(                                      &
     &         r_rho_levels(i,j+1,k) + r_rho_levels(i+1,j+1,k)          &
     &        +r_rho_levels(i,j,k) + r_rho_levels(i+1,j,k) )            &
     &        * delta_phi
            cy_rho_centre(i,j,k) = 1. / dy_rho_centre 
          End Do
        End Do

         Do j = j_start, j_stop
          Do i = 1, row_length
            dx_rho = r_rho_levels(i,j,k) * delta_lambda * cos_theta_latitude(i,j)
            cx_rho =  1. / dx_rho
            cx2_rho(i,j,k) = cx_rho * cx_rho
          End Do
        End Do

        Do j = 1, rows
          Do i = 1, row_length
            dy_rho = r_rho_levels(i,j,k) * delta_phi
            cy_rho =  1. / dy_rho            
            cy2_rho(i,j,k) = cy_rho * cy_rho            
          End Do
        End Do

      End Do   !  k = 1, model_levels

      IF (model_domain  ==  mt_global ) THEN
          If(at_extremity(PSouth))then
            Do k = 1, model_levels
              Do i = 1, row_length
                temp_s(i,k) = cx_theta_u(i,2,k)
                temp_s(i,k+model_levels) = cx2_rho(i,2,k)
              End Do
            End Do  !  k = 1, model_levels

            CALL global_2d_sums(temp_s, row_length, 1, 0, 0,            &
                              2*model_levels, sum_s, gc_proc_row_group)

            Do k = 1, model_levels
              Do i = 1, row_length
                cx_theta_u(i,1,k) = sum_s(k) / global_row_length
                cx2_rho(i,1,k) = sum_s(k+model_levels) /               &
                                               global_row_length
              End Do
            End Do  ! k = 1, model_levels
          End If
          If(at_extremity(PNorth))then
            Do k = 1, model_levels
              Do i = 1, row_length
                temp_n(i,k) = cx_theta_u(i,rows-1,k)
                temp_n(i,k+model_levels) = cx2_rho(i,rows-1,k)
              End Do
            End Do  !  k = 1, model_levels

            CALL global_2d_sums(temp_n, row_length, 1, 0, 0,            &
                              2*model_levels, sum_n, gc_proc_row_group)

            Do k = 1, model_levels
              Do i = 1, row_length
                cx_theta_u(i,rows,k) = sum_n(k) / global_row_length
                cx2_rho(i,rows,k) = sum_n(k+model_levels) /               &
                                               global_row_length
              End Do
            End Do  ! k = 1, model_levels
          End If
      END IF

! Vertical grid spacings used in calculation of rate of strain term

! Note: delta_zn(k) is at theta(k-1)
      Do j = 0, rows + 1
        Do i = 0, row_length + 1
          delta_zn(i,j,1) = r_rho_levels(i,j,1) - r_theta_levels(i,j,0)
            END DO
          END DO  

      DO k = 2, model_levels
        DO j = 0, rows + 1
          DO i = 0, row_length + 1
            delta_zn(i,j,k) = r_rho_levels(i,j,k)- r_rho_levels(i,j,k-1)
          END DO
        END DO  
      END DO ! k = 2, model_levels

      Do k = 1, model_levels

        Do j = 1, rows
          Do i = 1, row_length
            delta_z(i,j,k) = r_theta_levels(i,j,k) -                    &
                             r_theta_levels(i,j,k-1)
            rdz(i,j,k) = 1. / delta_z(i,j,k)
          End Do
        End Do

! Note: delta_zn_u/v(k) is at theta(k-1)

        Do j = 1, rows
          Do i = 0, row_length
            delta_zn_u = 0.5 * (delta_zn(i,j,k) + delta_zn(i+1,j,k))
            rdzn_u(i,j,k) = 1. / delta_zn_u
          End Do
        End Do
        
        Do j = 0, rows
          Do i = 1, row_length
            delta_zn_v = 0.5 * (delta_zn(i,j,k) + delta_zn(i,j+1,k))
            rdzn_v(i,j,k) = 1. / delta_zn_v
          End Do
        End Do
                
      End Do  !  k = 1, model_levels

!--------------------------------------------------------
! As in the LEM code
! _Now calculate half-squared strain rate SSQ on w-points
! _CX=1./DX, CY=1./DY, RDZ(K)=1./DZ(K),
!       RDZN(K) =1./DZN(K)
! _SSQ= 0.5*^^DU_I/DX_J+DU_J/DX_I^^**2
! _SSQIJ= (DU_I/DX_J+DU_J/DX_I)**2
! _Hence SSQ= SUM(I,J) {0.5*(SSQIJ)}
! _Note that for a simple shear S, SSQ=S**2
!   (as in the notation of Mason and Callen 1986)
!--------------------------------------------------------

      Do k = 1, model_levels - 1

        Do j = j_start, j_stop
          Do i = 1, row_length                
! Vertical weights
          weight_pl = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k))/    &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
          weight_min = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/ &
     &                 (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
          weight_2pl = delta_z(i,j,k)/(delta_z(i,j,k)+delta_z(i,j,k+1))
          weight_2min=delta_z(i,j,k+1)/(delta_z(i,j,k)+delta_z(i,j,k+1))

! ssq11 = 2 * backward difference (du/dx)^2 averaged to w point

          ssq11 =  2.0*(                                                &
     &      weight_pl*cx2_rho(i,j,k+1)*(u(i,j,k+1)-u(i-1,j,k+1))**2 +   &
     &      weight_min*cx2_rho(i,j,k)*(u(i,j,k)-u(i-1,j,k))**2   )

! ssq22 = 2 * backward difference (dv/dy)^2 averaged to w point  

          ssq22 = 2.0*(                                                 &
     &      weight_pl*cy2_rho(i,j,k+1)*(v(i,j,k+1)-v(i,j-1,k+1))**2 +   &
     &      weight_min*cy2_rho(i,j,k)*(v(i,j,k)-v(i,j-1,k))**2  )

! ssq33 =  2 * backward difference (dw/dz)^2 averaged to w point 
   
          ssq33 =  2.0*(                                                &
     &      weight_2pl*((w(i,j,k)-w(i,j,k-1))*rdz(i,j,k))**2 +          &
     &      weight_2min*((w(i,j,k+1)-w(i,j,k))*rdz(i,j,k+1))**2  )

! ssq13 =  (du/dz+dw/dx)^2 on u/theta averaged u:(i,i-1) to w point 
! Note: rdzn_u/v(k) is at theta(k-1)
! ** Needs rdzn_u(i-1:i,j,k+1),cx_theta_u(i-1:i,j,k))
! ssq31 = ssq13

          ssq13= (                                                      &
     &       ((u(i,j,k+1)-u(i,j,k))*rdzn_u(i,j,k+1) +                   &
     &        (w(i+1,j,k)-w(i,j,k))*cx_theta_u(i,j,k))**2 +             &
     &       ((u(i-1,j,k+1)-u(i-1,j,k))*rdzn_u(i-1,j,k+1) +             &
     &        (w(i,j,k)-w(i-1,j,k))*cx_theta_u(i-1,j,k))**2             &
     &            )*0.5      ! _averaging ssq13 over 2 points

! ssq23 =  (dw/dy+dv/dz)^2 on v/theta averaged v:(j-1,j) to w point 
! Note: rdzn_u/v(k) is at theta(k-1)
!**Needs rdzn_v(i,j-1:j,k+1),cy_theta_v(i,j-1,k)
! ssq32 = ssq23

          ssq23 = (                                                     &
     &       ((w(i,j,k)-w(i,j-1,k))*cy_theta_v(i,j-1,k) +               &
     &        (v(i,j-1,k+1)-v(i,j-1,k))*rdzn_v(i,j-1,k+1))**2 +         &
     &       ((w(i,j+1,k)-w(i,j,k))*cy_theta_v(i,j,k) +                 &
     &        (v(i,j,k+1)-v(i,j,k))*rdzn_v(i,j,k+1))**2                 &
     &             )*0.5        ! _averaging ssq23 over 2 points
!
! define points e(i,j) to be at u(i,j+1/2) and v(i+1/2,j)
! so a backward difference du/dy is at e(i-1,j)
! and a backward difference dv/dx is at e(i-1,j) 
!
! ssq12 = (du/dy+dv/dx)^2 on pi/rho averaged 
!         rho:(k:k+1)  u:(i-1:i) v:(j,j+1)  
!**Needs cy_rho_centre(i-1:i,j-1:j,k:k+1),
!        cx_rho_centre(i-1:i,j-1:j,k:k+1))
!   
! ssq21 = ssq12
          ssq12 = 0.25*weight_pl* (                                     &
!Level k+1
! (du/dy+dv/dx)^2 on rho level at  u(i-1,j+1/2), v(i-1/2,j)=e(i-1,j)
     &     ((u(i-1,j+1,k+1)-u(i-1,j,k+1))*cy_rho_centre(i-1,j,k+1)+     &
     &      (v(i,j,k+1)-v(i-1,j,k+1))*cx_rho_centre(i-1,j,k+1))**2      &
! (du/dy+dv/dx)^2 on rho level at  u(i,j+1/2), v(i+1/2,j)=e(i,j)        
     &    +((u(i,j+1,k+1)-u(i,j,k+1))*cy_rho_centre(i,j,k+1) +          &
     &      (v(i+1,j,k+1)-v(i,j,k+1))*cx_rho_centre(i,j,k+1))**2        &
! (du/dy+dv/dx)^2 on rho level at u(i-1,j-1/2), v(i-1/2,j-1)=e(i-1,j-1)
     &    +((u(i-1,j,k+1)-u(i-1,j-1,k+1))*cy_rho_centre(i-1,j-1,k+1) +  &
     &      (v(i,j-1,k+1)-v(i-1,j-1,k+1))*cx_rho_centre(i-1,j-1,k+1))**2&
! (du/dy+dv/dx)^2 on rho level at  u(i,j-1/2), v(i+1/2,j-1)=e(i,j-1)
     &    +((u(i,j,k+1)-u(i,j-1,k+1))*cy_rho_centre(i,j-1,k+1) +        &
     &      (v(i+1,j-1,k+1)-v(i,j-1,k+1))*cx_rho_centre(i,j-1,k+1))**2  &
     &        ) +  0.25*weight_min * (                                  &
!Level k
! (du/dy+dv/dx)^2 on rho level at  u(i-1,j+1/2), v(i-1/2,j)=e(i-1,j)
     &     ((u(i-1,j+1,k)-u(i-1,j,k))*cy_rho_centre(i-1,j,k)+           &
     &      (v(i,j,k)-v(i-1,j,k))*cx_rho_centre(i-1,j,k) )**2           &
! (du/dy+dv/dx)^2 on rho level at  u(i,j+1/2), v(i+1/2,j)=e(i,j)         
     &    +((u(i,j+1,k)-u(i,j,k))*cy_rho_centre(i,j,k) +                &
     &      (v(i+1,j,k)-v(i,j,k))*cx_rho_centre(i,j,k))**2              &
! (du/dy+dv/dx)^2 on rho level at  u(i-1,j-1/2), v(i-1/2,j-1)=e(i,j-1)
     &    +((u(i-1,j,k)-u(i-1,j-1,k))*cy_rho_centre(i-1,j-1,k) +        &
     &      (v(i,j-1,k)-v(i-1,j-1,k))*cx_rho_centre(i-1,j-1,k))**2      &
! (du/dy+dv/dx)^2 on rho level at  u(i,j-1/2), v(i+1/2,j-1)=e(i,j-1)
     &    +((u(i,j,k)-u(i,j-1,k))*cy_rho_centre(i,j-1,k) +              &
     &      (v(i+1,j-1,k)-v(i,j-1,k))*cx_rho_centre(i,j-1,k))**2        &
     &        )      ! _averaging s12 over 8 points

! Note: ssq superficially looks as if a factor of 0.5 is missing. 
! It is not.
          sum(i,j,k) = ssq11 + ssq22 + ssq33 +                          &
                               ssq13 + ssq23 + ssq12 + smallp

          End Do   ! on i
        End Do    ! on j
 
      End Do ! k = 1, model_levels - 1

      IF (model_domain  ==  mt_global ) THEN

          If(at_extremity(PSouth))then
            Do k = 1, model_levels - 1
              Do i = 1, row_length
                temp_s(i,k) = sum(i,2,k)
              End Do
            End Do  !  k = 1, model_levels - 1

            CALL global_2d_sums(temp_s, row_length, 1, 0, 0,            &
                                model_levels - 1, sum_s,                &
                                gc_proc_row_group )

            Do k = 1, model_levels - 1
              Do i = 1, row_length
                sum(i,1,k) = sum_s(k) / global_row_length
              End Do
            End Do  ! k = 1, model_levels - 1
          End If

          If(at_extremity(PNorth))then
            Do k = 1, model_levels - 1
              Do i = 1, row_length
                temp_n(i,k) = sum(i,rows-1,k)
              End Do
            End Do  !  k = 1, model_levels - 1

            CALL global_2d_sums(temp_n, row_length, 1, 0, 0,            &
                                model_levels - 1, sum_n,                &
                                gc_proc_row_group )

            Do k = 1, model_levels - 1
              Do i = 1, row_length
                sum(i,rows,k) = sum_n(k) / global_row_length
              End Do
            End Do  ! k = 1, model_levels - 1
          End If

      END IF  !  model_domain  ==  mt_global

      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
           shear(i,j,k) = SQRT(sum(i,j,k)) ! already been *0.5 in ssq
          End Do
        End Do
      End Do ! k = 1, model_levels - 1

! Horizontal grid spacings used in calculation of mixing length and
! maximum diffusivity.
! delta_lambda and delta_phi are in radians
!
      delta_y = Earth_radius * delta_phi 
      Do j = 1, rows
        Do i = 1, row_length
          delta_x = Earth_radius * delta_lambda * cos_v_latitude(i,j)
          IF (l_subfilter_blend) THEN
            delta_smag(i,j) = MIN(delta_x, delta_y)
          ELSE
            delta_smag(i,j) = MAX(delta_x, delta_y)
          END IF
          rmlmax(i,j) = mix_factor * delta_smag(i,j)
!
! maximum diffusion coefficient allowed in this run
!
          max_diff(i,j) = diff_factor /                                 &
                            ( (1.0/(delta_x*delta_x) +                  &
                               1.0/(delta_y*delta_y) ) * 4.0 * timestep)
        End Do
      End Do
                  
      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
            
! APL: ought to combine these two lines but will presumably change answers!
            rneutml=SQRT(1./                                            &
     &        (1./(VKMAN*(z(i,j,k)+z0(i,j)))**2+1./RMLMAX(i,j)**2) )
     
            rneutml_sq(i,j,k)=rneutml*rneutml
            visc_m(i,j,k) = shear(i,j,k)
            visc_h(i,j,k) = shear(i,j,k)
!
! visc is now S and still needs to be multiplied by
! stability functions and lambda^2* (done in BL code)
!
          End Do
        End Do
      End Do ! k = 1, model_levels - 1

      IF (lhook) CALL dr_hook('TURB_SMAGORINSKY',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE TURB_Smagorinsky

