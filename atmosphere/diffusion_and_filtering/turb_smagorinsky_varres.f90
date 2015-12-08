! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine turb_smagorinsky_varres
!
SUBROUTINE turb_smagorinsky_varres(                                     &
                                   u, v, w, z0, timestep,               &
                                   row_length, rows, n_rows,            &
                                   model_levels, bl_levels,             &
                                   r_theta_levels, r_rho_levels,        &
                                   cos_theta_latitude,                  &
                                   cos_v_latitude,                      &
                                   lambda_p, phi_p, lambda_u, phi_v,    &
                                   dlambda_p, dphi_p, dlambda_u, dphi_v,&
                                   recip_dlamp, recip_dphip,            &
                                   recip_dlamu, recip_dphiv )

      USE turb_diff_mod, ONLY: L_subfilter_blend, diff_factor, mix_factor
      USE turb_diff_ctl_mod, ONLY: visc_m, visc_h, rneutml_sq, max_diff,&
                                   shear, delta_smag
      USE atmos_constants_mod, ONLY: vkman
      USE earth_constants_mod, ONLY: earth_radius

! Description: Calculates coefficients for use in subgrid turbulence
!              scheme when using variable resolution option.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
USE vectlib_mod, ONLY :  oneover_v, sqrt_v

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars

IMPLICIT NONE

! halo information


! Variables with Intent (In)

INTEGER, INTENT(IN) ::                                            &
  row_length                                                      &
                         ! number of points on a row.
, rows                                                            &
                         ! number of rows.
, n_rows                                                          &
                         ! number of v rows.
, model_levels                                                    &
                         ! number of model levels.
, bl_levels              ! number of bl levels.

REAL, INTENT(IN) ::                                               &
  u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)     &
, z0(row_length,rows)                                             &
                             ! roughness length
, timestep                                                        &
, r_theta_levels (1-halo_i:row_length+halo_i,                     &
                    1-halo_j:rows+halo_j, 0:model_levels)         &
, r_rho_levels (1-halo_i:row_length+halo_i,                       &
                    1-halo_j:rows+halo_j, model_levels)           &
, cos_theta_latitude (1-Offx:row_length+Offx,1-Offy:rows+Offy)    &
, cos_v_latitude (1-Offx:row_length+Offx,1-Offy:n_rows+Offy)

!VarRes horizontal co-ordinate information
REAL, INTENT(IN) ::                                               &
  lambda_p(1-halo_i:row_length+halo_i)                            &
, phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
, lambda_u(1-halo_i:row_length+halo_i)                            &
, phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
, dlambda_p(1-halo_i:row_length+halo_i)                           &
, dlambda_u(1-halo_i:row_length+halo_i)                           &
, dphi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)        &
, dphi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)      &
, recip_dlamp(1-halo_i:row_length+halo_i)                         &
, recip_dphip(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
, recip_dlamu(1-halo_i:row_length+halo_i)                         &
, recip_dphiv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

!Local parameters

INTEGER :: i,j,k ! Loop counters
INTEGER row_by_rl ! row times row_length for vectlib routines

REAL ::                                                           &
  delta_x                                                         &
        ! Sea level longitude spacing in m
, delta_y                                                         &
        ! Sea level latitude spacing in m
, delta_z(row_length, rows, model_levels)                         &
        ! Backward difference dz at rho points
, delta_zn(0:row_length+1,0:rows+1, model_levels)                 &
        ! Backward difference dz at theta points
        ! Note: k=1 is special - uses surface
        ! Note: delta_zn(k) is correct for theta(k-1)
, delta_zn_u                                                      &
        ! Average dz at theta points to u points (theta levels)
, delta_zn_v                                                      &
        ! Average dz at theta points to v points (theta levels)
, rdz(row_length, rows, model_levels)                             &
        ! 1/delta_z
, rdzn_u(0:row_length, rows, model_levels)                        &
        ! 1/delta_zn_u
, rdzn_v(row_length, 0:rows, model_levels)                        &
        ! 1/delta_zn_v
, r_ave 

REAL ::                                                           &
  dx_rho                                                          &
!               ! dx on rho points
, dx_theta_u                                                      &
!               ! dx above u points on theta levels
, dx_rho_centre                                                   &
!               ! dx in centre of grid box on rho levels
, dy_rho                                                          &
!               ! dy on rho points
, dy_theta_v                                                      &
!               ! dy above v points on theta levels
, dy_rho_centre                                                   &
!               ! dy in centre of grid box on rho levels
, cx_rho                                                          &
!               ! reciprocal of dx_rho
, cy_rho                                                          &
!               ! reciprocal of dy_rho
, cx_theta_u(0:row_length, rows, model_levels)                    &
!               ! reciprocal of dx_theta_u
, cx_rho_centre(0:row_length, 0:rows, model_levels)               &
!               ! reciprocal of dx_rho_centre
, cy_theta_v(row_length, 0:rows, model_levels)                    &
!               ! reciprocal of dy_theta_v
, cy_rho_centre(0:row_length, 0:rows, model_levels)               &
!               ! reciprocal of dy_rho_centre
, cx2_rho(row_length, rows, model_levels)                         &
!               ! square of cx_rho
, cy2_rho(row_length, rows, model_levels)                         &                    
!               ! square of cy_rho
, rmlmax(row_length, rows)                                        &  
, temp(row_length, rows)                                          &
, rmlmax_sq_inv(row_length, rows)                                 &
!               ! basic mixing length (lambda0)
, z(row_length, rows, model_levels)

REAL, PARAMETER :: smallp = 1.e-14  ! A small number

REAL ::                                                           &
  weight_pl                                                       &
, weight_min                                                      &
, weight_2pl                                                      &
, weight_2min                                                     &
, w1(0:row_length)                                                &
, w2(0:row_length)                                                &
, w3(0:row_length, 0:rows)                                        &
, w4(0:row_length, 0:rows)                                        &
, w5(row_length)                                                  &
, w6(row_length)                                                  &
, w7(row_length, rows)                                            &
, w8(row_length, rows)                                            &
, ssq11                                                           &                                      
                      ! ii component of shear stress
, ssq22                                                           &
                      ! jj component of shear stress
, ssq33                                                           &
                      ! kk component of shear stress
, ssq13                                                           &
                      ! ik component of shear stress
, ssq23                                                           &
                      ! jk component of shear stress
, ssq12                                                           &
                      ! ij component of shear stress
, sum(row_length, rows, model_levels-1)   
                      ! sum of Sij (ssqij)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('TURB_SMAGORINSKY_VARRES',zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Initialise variables
!----------------------------------------------------------------------
      sum(:,:,:) = 0.0

!$OMP PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED(model_levels,halo_j,rows,row_length,halo_i,visc_m,      &
!$OMP& visc_h,z,r_theta_levels,w,w1,w2,w3,w4,w5,w6,w7,w8,dlambda_p,   &
!$OMP& cos_theta_latitude,r_rho_levels,cx_theta_u,dphi_p,             & 
!$OMP& cx_rho_centre,cy_rho_centre,dlambda_u,cx2_rho,recip_dphip,     &
!$OMP& cy2_rho, delta_z,rdz,delta_zn,rdzn_u,rdzn_v,u,v,shear,sum,     &
!$OMP& diff_factor,mix_factor,max_diff,rneutml_sq,rmlmax_sq_inv,      &
!$OMP& rmlmax,z0,cy_theta_v,phi_v,recip_dlamu,lambda_p,lambda_u,      &
!$OMP& recip_dphiv,timestep,dphi_v,cos_v_latitude,recip_dlamp,phi_p,  &
!$OMP& delta_smag,l_subfilter_blend)                                  &
!$OMP& PRIVATE(i,j,k,dx_theta_u,r_ave,dx_rho_centre,dy_theta_v,cx_rho,&
!$OMP& cy_rho,delta_zn_u,delta_zn_v,weight_pl,weight_min,weight_2pl,  &
!$OMP& weight_2min,ssq11,ssq22,ssq33,ssq13,ssq23,ssq12, row_by_rl,    &
!$OMP& delta_y,delta_x,dy_rho_centre,dy_rho,dx_rho,temp)

row_by_rl = rows*row_length

!$OMP DO SCHEDULE(STATIC)
DO k= 1, model_levels
  DO j = 1, rows
    DO i = 1, row_length
      z(i,j,k) = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
    END DO
  END DO   
END DO
!$OMP END DO NOWAIT

! Horizontal weights
! w1(i) = weight for point p(i) averaged to u(i)
! w2(i) = weight for point p(i+1) averaged to u(i)
! Needed i=(0:row_length)

!$OMP DO SCHEDULE(STATIC)
DO i=0, row_length
  w1(i) = (lambda_p(i+1) - lambda_u(i)) * recip_dlamp(i)
  w2(i) = (lambda_u(i) - lambda_p(i)) * recip_dlamp(i)
END DO
!$OMP END DO NOWAIT

! Vertical weights
! w3(i,j) = weight for point p(j+1) averaged to v(j)
! w4(i,j) = weight for point p(j) averaged to v(j)
! Needed j=(0:rows)

!$OMP DO SCHEDULE(STATIC)
DO j = 0, rows
  DO i = 0, row_length
    w3(i,j) = (phi_v(i,j) - phi_p(i,j)) * recip_dphip(i,j)  
    w4(i,j) = (phi_p(i,j+1)- phi_v(i,j)) * recip_dphip(i,j)
  END DO
END DO
!$OMP END DO NOWAIT

! Horizontal weights
! w5(i) = weight for point u(i) averaged to p(i)
! w6(i) = weight for point u(i-1) averaged to p(i)
! Needed i=(1:row_length)

!$OMP DO SCHEDULE(STATIC)
DO i=1, row_length
  w5(i) = ( lambda_p(i) - lambda_u(i-1) ) * recip_dlamu(i-1)
  w6(i) = ( lambda_u(i) - lambda_p(i) ) * recip_dlamu(i-1)
END DO
!$OMP END DO NOWAIT

! Vertical weights
! w7(i,j) = weight for point v(j) averaged to p(j)
! w8(i,j) = weight for point v(j-1) averaged to p(j)
! Needed j=(1:rows)

!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    w7(i,j) = ( phi_p(i,j) - phi_v(i,j-1) ) * recip_dphiv(i,j-1)
    w8(i,j) = ( phi_v(i,j) - phi_p(i,j) ) * recip_dphiv(i,j-1)
  END DO
END DO
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels

  DO j = 1, rows        
    DO i = 0, row_length    

      dx_theta_u =                                                &
         (w1(i)*r_theta_levels(i,j,k) +                           &
          w2(i)*r_theta_levels(i+1,j,k))                          &
       * dlambda_p(i)*cos_theta_latitude(i,j)

      cx_theta_u(i,j,k) =  1./dx_theta_u
      
    END DO    
  END DO
  
  DO j = 0, rows  
    DO i = 1, row_length

      dy_theta_v =                                                &
         (w3(i,j)*r_theta_levels(i,j,k) +                         &
          w4(i,j)*r_theta_levels(i,j+1,k))                        &
       * dphi_p(i,j)
     
      cy_theta_v(i,j,k) =  1./dy_theta_v

    END DO    
  END DO
  
  DO j = 0, rows
    DO i = 0, row_length
    
      r_ave=                                                      &
        w4(i,j)*(w2(i)*r_rho_levels(i+1,j,k) +                    & 
                 w1(i)*r_rho_levels(i,j,k)  ) +                   &
        w3(i,j)*(w2(i)*r_rho_levels(i+1,j+1,k) +                  &
                 w1(i)*r_rho_levels(i,j+1,k) )

      dx_rho_centre =  r_ave * dlambda_p(i)*cos_v_latitude(i,j)

      cx_rho_centre(i,j,k) =  1./dx_rho_centre

      dy_rho_centre = r_ave * dphi_p(i,j)

      cy_rho_centre(i,j,k) = 1./dy_rho_centre


    END DO
  END DO

  DO j = 1, rows
    DO i = 1, row_length

! Needs checking
      dx_rho = r_rho_levels(i,j,k)* dlambda_u(i-1)                &
                      *cos_theta_latitude(i,j)

! Needs checking
      dy_rho = r_rho_levels(i,j,k) * dphi_v(i,j-1)

      cx_rho = 1./dx_rho
      cy_rho =  1./dy_rho

      cx2_rho(i,j,k) = cx_rho*cx_rho
      cy2_rho(i,j,k) = cy_rho*cy_rho
      
    END DO
  END DO

END DO   !k
!$OMP END DO NOWAIT

! Vertical grid spacings used in calculation of rate of strain term

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels

  DO j = 1, rows
    DO i = 1, row_length
      delta_z(i,j,k) = r_theta_levels(i,j,k)                      &
                -r_theta_levels(i,j,k-1)
      rdz(i,j,k) = 1./delta_z(i,j,k)
    END DO
  END DO

! Note: delta_zn(k) is at theta(k-1)
  IF (k  ==  1) THEN
    DO j = 0, rows+1
      DO i = 0, row_length+1
        delta_zn(i,j,k) = r_rho_levels(i,j,1)                     &                
                         -r_theta_levels(i,j,0)
      END DO
    END DO  
  ELSE
! Note: delta_zn(k) is at theta(k-1)
    DO j = 0, rows+1
      DO i = 0, row_length+1
        delta_zn(i,j,k) = r_rho_levels(i,j,k)                     &
                        - r_rho_levels(i,j,k-1)
      END DO
    END DO  
  END IF

END DO  ! k
!$OMP END DO 
! Note: delta_zn_u/v(k) is at theta(k-1)

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels

  DO j = 1, rows
    DO i = 0, row_length
      delta_zn_u =                                                &
             w1(i)*delta_zn(i,j,k)+w2(i)*delta_zn(i+1,j,k)
      rdzn_u(i,j,k) = 1./delta_zn_u
    END DO
  END DO
  
  DO j = 0, rows
    DO i = 1, row_length
      delta_zn_v =                                                &
             w4(i,j)*delta_zn(i,j,k)+w3(i,j)*delta_zn(i,j+1,k)
      rdzn_v(i,j,k) = 1./delta_zn_v
    END DO
  END DO
  
END DO
!$OMP END DO

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

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels -1

  DO j = 1, rows

    DO i = 1, row_length
          
! Vertical weights

      weight_pl = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k))/    &
                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
      weight_min = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/ &
                 (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
      weight_2pl = delta_z(i,j,k)/(delta_z(i,j,k)+delta_z(i,j,k+1))
      weight_2min=delta_z(i,j,k+1)/(delta_z(i,j,k)+delta_z(i,j,k+1))

! ssq11 = 2 * backward difference (du/dx)^2 averaged to w point

      ssq11 =  2.0*(                                                &
        weight_pl*cx2_rho(i,j,k+1)*(u(i,j,k+1)-u(i-1,j,k+1))**2 +   &
        weight_min*cx2_rho(i,j,k)*(u(i,j,k)-u(i-1,j,k))**2   )

! ssq22 = 2 * backward difference (dv/dy)^2 averaged to w point  

      ssq22 = 2.0*(                                                 &
        weight_pl*cy2_rho(i,j,k+1)*(v(i,j,k+1)-v(i,j-1,k+1))**2 +   &
        weight_min*cy2_rho(i,j,k)*(v(i,j,k)-v(i,j-1,k))**2  )

! ssq33 =  2 * backward difference (dw/dz)^2 averaged to w point 
   
      ssq33 =  2.0*(                                                &
        weight_2pl*((w(i,j,k)-w(i,j,k-1))*rdz(i,j,k))**2 +          &
        weight_2min*((w(i,j,k+1)-w(i,j,k))*rdz(i,j,k+1))**2  )

! ssq13 =  (du/dz+dw/dx)^2 on u/theta averaged u:(i,i-1) to w point 
! Note: rdzn_u/v(k) is at theta(k-1)
! ** Needs rdzn_u(i-1:i,j,k+1),cx_theta_u(i-1:i,j,k))
! ssq31 = ssq13

      ssq13=                                                        &
        w5(i)*( (u(i,j,k+1)-u(i,j,k))*rdzn_u(i,j,k+1)               &
             +(w(i+1,j,k)-w(i,j,k))*cx_theta_u(i,j,k) )**2          &
       +w6(i)*( (u(i-1,j,k+1)-u(i-1,j,k))*rdzn_u(i-1,j,k+1)         &
             +(w(i,j,k)-w(i-1,j,k))*cx_theta_u(i-1,j,k) )**2              
                        ! weighted averaging ssq13 over 2 points

! ssq23 =  (dw/dy+dv/dz)^2 on v/theta averaged v:(j-1,j) to w point 
! Note: rdzn_u/v(k) is at theta(k-1)
!**Needs rdzn_v(i,j-1:j,k+1),cy_theta_v(i,j-1,k)
! ssq32 = ssq23

      ssq23 =                                                       &
        w8(i,j)*( (w(i,j,k)-w(i,j-1,k))*cy_theta_v(i,j-1,k)         &
             +(v(i,j-1,k+1)-v(i,j-1,k))*rdzn_v(i,j-1,k+1))**2       &
       +w7(i,j)*( (w(i,j+1,k)-w(i,j,k))*cy_theta_v(i,j,k)           &
             +(v(i,j,k+1)-v(i,j,k))*rdzn_v(i,j,k+1))**2
                       !  weighted averaging ssq23 over 2 points
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
      ssq12 = weight_pl* (                                              &
!Level k+1
        w7(i,j)*w6(i)*                                                  &
! (du/dy+dv/dx)^2 on rho level at  u(i-1,j+1/2), v(i-1/2,j)=e(i-1,j)         
        ((u(i-1,j+1,k+1)-u(i-1,j,k+1))*cy_rho_centre(i-1,j,k+1)         &
        +(v(i,j,k+1)-v(i-1,j,k+1))*cx_rho_centre(i-1,j,k+1))**2         &
       +w7(i,j)*w5(i)*                                                  &
! (du/dy+dv/dx)^2 on rho level at  u(i,j+1/2), v(i+1/2,j)=e(i,j)         
        ((u(i,j+1,k+1)-u(i,j,k+1))*cy_rho_centre(i,j,k+1)               &
        +(v(i+1,j,k+1)-v(i,j,k+1))*cx_rho_centre(i,j,k+1))**2           &
       +w8(i,j)*w6(i)*                                                  &
! (du/dy+dv/dx)^2 on rho level at u(i-1,j-1/2), v(i-1/2,j-1)=e(i-1,j-1)          
        ((u(i-1,j,k+1)-u(i-1,j-1,k+1))*cy_rho_centre(i-1,j-1,k+1)       &
        +(v(i,j-1,k+1)-v(i-1,j-1,k+1))*cx_rho_centre(i-1,j-1,k+1))**2   &
       +w8(i,j)*w5(i)*                                                  &
! (du/dy+dv/dx)^2 on rho level at  u(i,j-1/2), v(i+1/2,j-1)=e(i,j-1)         
        ((u(i,j,k+1)-u(i,j-1,k+1))*cy_rho_centre(i,j-1,k+1)             &
        +(v(i+1,j-1,k+1)-v(i,j-1,k+1))*cx_rho_centre(i,j-1,k+1))**2     &
        ) +  weight_min * (                                             &
!Level k 
        w7(i,j)*w6(i)*                                                  &
! (du/dy+dv/dx)^2 on rho level at  u(i-1,j+1/2), v(i-1/2,j)=e(i-1,j)         
        ((u(i-1,j+1,k)-u(i-1,j,k))*cy_rho_centre(i-1,j,k)               &
        +(v(i,j,k)-v(i-1,j,k))*cx_rho_centre(i-1,j,k) )**2              &
       +w7(i,j)*w5(i)*                                                  &
! (du/dy+dv/dx)^2 on rho level at  u(i,j+1/2), v(i+1/2,j)=e(i,j)         
        ((u(i,j+1,k)-u(i,j,k))*cy_rho_centre(i,j,k)                     &
        +(v(i+1,j,k)-v(i,j,k))*cx_rho_centre(i,j,k))**2                 &
       +w8(i,j)*w6(i)*                                                  &
! (du/dy+dv/dx)^2 on rho level at  u(i-1,j-1/2), v(i-1/2,j-1)=e(i,j-1)         
        ((u(i-1,j,k)-u(i-1,j-1,k))*cy_rho_centre(i-1,j-1,k)             &
       +(v(i,j-1,k)-v(i-1,j-1,k))*cx_rho_centre(i-1,j-1,k))**2          &
       +w8(i,j)*w5(i)*                                                  &
! (du/dy+dv/dx)^2 on rho level at  u(i,j-1/2), v(i+1/2,j-1)=e(i,j-1)         
        ((u(i,j,k)-u(i,j-1,k))*cy_rho_centre(i,j-1,k)                   &
        +(v(i+1,j-1,k)-v(i,j-1,k))*cx_rho_centre(i,j-1,k))**2           &
       )      ! weighted averaging s12 over 8 points

! Note: the 2nd term looks like it has a factor of 0.5 missing 
! it is not.
      sum(i,j,k) = sum(i,j,k) + ssq11+ssq22+ssq33+ssq13+ssq23+ssq12+smallp

    END DO   ! on i
  END DO    ! on j
  
  CALL sqrt_v(row_by_rl, sum(1,1,k), shear(1,1,k))

END DO ! k = 1, model_levels -1
!$OMP END DO NOWAIT


! Horizontal grid spacings used in calculation of mixing length and
! maximum diffusivity.
! delta_lambda and delta_phi are in radians
!

!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    
    delta_y = Earth_radius*dphi_v(i,j-1)  
    delta_x = Earth_radius*dlambda_u(i-1)*cos_v_latitude(i,j)
    IF (l_subfilter_blend) THEN
      delta_smag(i,j) = MIN(delta_x, delta_y)
    ELSE
      delta_smag(i,j) = MAX(delta_x, delta_y)
    END IF
    rmlmax(i,j) = mix_factor * delta_smag(i,j)
!
! maximum diffusion coefficient allowed in this run
!
    max_diff(i,j) = diff_factor /                                       &
                        ( (1.0/(delta_x*delta_x) +                      &
                           1.0/(delta_y*delta_y)) * 4.0 * timestep)
    
  END DO  ! on i
END DO  ! on j
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    rmlmax_sq_inv(i,j) = 1./(rmlmax(i,j)**2)
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels-1 
  
  DO j = 1, rows
    DO i = 1, row_length
      temp(i,j)= (vkman*(z(i,j,k)+z0(i,j)))**2
    END DO
  END DO

  CALL oneover_v(row_by_rl, temp, temp)

  DO j = 1, rows
    DO i = 1, row_length
      temp(i,j)= temp(i,j) + rmlmax_sq_inv(i,j)
    END DO
  END DO

  CALL oneover_v(row_by_rl, temp, temp)

  DO j = 1, rows
    DO i = 1, row_length
      rneutml_sq(i,j,k) = temp(i,j)
      visc_m(i,j,k) = shear(i,j,k)
      visc_h(i,j,k) = shear(i,j,k)
    END DO
  END DO
!
! visc is now S and still needs to be multiplied by
! stability functions and lambda^2 (done in BL code) 
!

END DO  !  k = 1, model_levels -1
!$OMP END DO NOWAIT

!$OMP END PARALLEL

IF (lhook) CALL dr_hook('TURB_SMAGORINSKY_VARRES',zhook_out,zhook_handle)
RETURN
END SUBROUTINE turb_smagorinsky_varres

