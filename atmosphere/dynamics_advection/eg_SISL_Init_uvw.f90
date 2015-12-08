! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_sisl_init_uvw_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_sisl_init_uvw(                                          &
         row_length, rows, n_rows, model_levels, ih, g_theta,         &
         u, v, w, thetav, rho, m_v, m_cl, m_cf, m_r,                  &
         m_gr, m_cf2, exner, exner_star, r_u, r_v, r_w,               &
         l_call_from_solver, l_call_from_f1sp,psi_w_surf,psi_w_lid)

USE level_heights_mod, ONLY : eta_theta_levels, eta_rho_levels,        &
                              xi3_at_theta=>r_theta_levels
USE um_parvars,        ONLY : offx, offy, halo_i, halo_j, Pnorth, Psouth
USE eg_alpha_mod,      ONLY : alpha_u, alpha_v, alpha_w
USE timestep_mod,      ONLY : timestep
USE proc_info_mod,     ONLY : model_domain, at_extremity
USE metric_terms_mod,  ONLY : h1_xi1_u, h1_xi2_v, h1_p_eta, h2_xi1_u, &
                              h2_xi2_v, h2_p_eta, h3_xi1_u, h3_xi2_v, & 
                              h3_p_eta, deta_xi3, deta_xi3_theta,     &
                              deta_xi3_u,  deta_xi3_v
USE coriolis_mod,      ONLY : f1_star, f2_star, f3_star,              &
                              f1_comp, f2_comp, f3_comp
USE parkind1,          ONLY : jpim, jprb       !DrHook
USE yomhook,           ONLY : lhook, dr_hook   !DrHook
USE dynamics_input_mod, ONLY : l_simple_coriolis
USE atmos_constants_mod
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE Field_Types
USE eg_swap_bounds_mod

IMPLICIT NONE
!
! Description: computes time level n arrival point quantities:
!              Ru, Rv, Rw, Rtheta, Rrho, Rm_v, Rm_cl, Rm_cf
!  
!
! Method: ENDGame formulation version 1.01,
!         section 11 (solution procedure), paragraph 1.

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

INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels

REAL                                                                    &
  psi_w_surf(row_length,rows),                                          &
  psi_w_lid (row_length,rows)


! Loop index bounds for arrays defined on p, u, v points respectively

! SI time weights     & hydrostatic switch

REAL, INTENT(IN) :: ih

! Gravity arrays

REAL, INTENT(IN) ::                                                   &
  g_theta (1-offx:row_length+offx, 1-offy:rows+offy,                  &
           0:model_levels)

REAL, INTENT(IN) ::                                                   &
  u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),           &
  v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),         &
  w(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),          &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
  exner_star(1-offx:row_length+offx,1-offy:rows+offy),                &
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels ),         &
  m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),       &
  m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),      &
  m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),      &
  m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),       &
  m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels ),      &
  m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels )


LOGICAL, INTENT(IN) :: l_call_from_solver, l_call_from_f1sp

! Timelevel n arrival point quantities

REAL, INTENT(INOUT) ::                                                &
  r_u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),         &
  r_v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),       &
  r_w(row_length,rows,0:model_levels)

REAL, INTENT(IN)    ::                                                &
  thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

! Local variables

INTEGER :: i,j,k, kp1

REAL :: beta_u_dt, beta_v_dt, beta_w_dt
REAL :: p_grad_coeff, thetav_ave

! Allocate temps to aid loop vectorization and readability

REAL ::                                                               &
  rho_wet(1-offx:row_length+offx,1-offy:rows+offy,model_levels),      &
  theta_wet(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),  &
  mix_fact(1-offx:row_length+offx,1-offy:rows+offy),                  &
  ustar(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),       &
  vstar(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),     &
  wstar(row_length,rows,0:model_levels),                              &
  work1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  work2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
  dxi1_u(udims%i_start:udims%i_end),                                          &
  dxi1_p(pdims%i_start:pdims%i_end),                                          &
  dxi2_v(vdims%j_start:vdims%j_end),                                          &
  dxi2_p(pdims%j_start:pdims%j_end),                                          &
  deta_w(0:model_levels),                                             &
  deta_rho(model_levels)

REAL temp



! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_SISL_INIT_UVW',zhook_in,zhook_handle)
 

! Calculate total rho and virtual potential temperature at all
! grid points including haloes - saves the fill!

k = 0
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

      mix_fact(i,j) = 1.0 + m_v(i,j,k) + m_cf(i,j,k) + m_cl(i,j,k)    &
                          + m_r(i,j,k) + m_gr(i,j,k) + m_cf2(i,j,k)
      theta_wet(i,j,k) = thetav(i,j,k)/mix_fact(i,j)

  END DO
END DO


DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

         temp = 1.0 + m_v(i,j,k) + m_cf(i,j,k) + m_cl(i,j,k)          &
                    + m_r(i,j,k) + m_gr(i,j,k) + m_cf2(i,j,k)
         theta_wet(i,j,k) = thetav(i,j,k)/temp
         rho_wet(i,j,k)   = rho(i,j,k)*( intw_w2rho(k,1)*temp         &
                                        +intw_w2rho(k,2)*mix_fact(i,j))
         mix_fact(i,j) = temp

    END DO
  END DO
END DO

CALL eg_swap_bounds( rho_wet, pdims_s, fld_type_p, .FALSE. )
CALL eg_swap_bounds( theta_wet, tdims_s, fld_type_p, .FALSE. )

IF(l_call_from_solver.OR.l_call_from_f1sp) THEN
  beta_u_dt = alpha_u*timestep
  beta_v_dt = alpha_v*timestep
  beta_w_dt = alpha_w*timestep
ELSE
  beta_u_dt = (1.0-alpha_u)*timestep
  beta_v_dt = (1.0-alpha_v)*timestep
  beta_w_dt = (1.0-alpha_w)*timestep
END IF

DO i=pdims%i_start, pdims%i_end
  dxi1_p(i) = xi1_u(i)-xi1_u(i-1)
END DO

DO i=udims%i_start, udims%i_end
  dxi1_u(i) = xi1_p(i+1)-xi1_p(i)
END DO

DO j=vdims%j_start, vdims%j_end
  dxi2_v(j) = xi2_p(j+1)-xi2_p(j)
END DO

DO j=pdims%j_start, pdims%j_end
  dxi2_p(j) = xi2_v(j)-xi2_v(j-1)
END DO

k = 0
  deta_w(k) = eta_rho_levels(k+1) - eta_theta_levels(k)

DO k=1, model_levels-1
  deta_w(k) = eta_rho_levels(k+1) - eta_rho_levels(k)
  deta_rho(k)   = eta_theta_levels(k) - eta_theta_levels(k-1)
END DO

k = model_levels
deta_rho(k) = eta_theta_levels(k) - eta_theta_levels(k-1)

IF( .NOT. l_simple_coriolis ) THEN
!----------------------------------------------------------------------
! Compute u-Coriolis acceleration term (Omega x U)u eqns:
!   (7.10) ENDGame formulation vn1.01.
!----------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,vdims,vstar,dxi1_p,deta_rho,h1_xi2_v,      &
!$OMP& h3_xi2_v,deta_xi3_v,intw_p2v,rho_wet,v)
  DO k=1, model_levels
    DO j=vdims%j_start, vdims%j_end
      DO i=vdims%i_start, vdims%i_end
        vstar(i,j,k) = dxi1_p(i)*deta_rho(k)*h1_xi2_v(i,j,k)*           &
                       h3_xi2_v(i,j,k)*deta_xi3_v(i,j,k)*               &
                       ( intw_p2v(j,1)*rho_wet(i,j,k) +                 &
                       intw_p2v(j,2)*rho_wet(i,j+1,k) )*v(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

!----------------------------------------------------------------------
! Compute v-Coriolis acceleration term (Omega x U)v eqns:
!   (7.11) ENDGame formulation vn1.01.
!----------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,udims,ustar,dxi2_p,deta_rho,h2_xi1_u,      &
!$OMP& h3_xi1_u,deta_xi3_u,intw_p2u,rho_wet,u)
  DO k=1, model_levels
    DO j=udims%j_start, udims%j_end
      DO i=udims%i_start, udims%i_end
        ustar(i,j,k) = dxi2_p(j)*deta_rho(k)*h2_xi1_u(i,j,k)*           &
                       h3_xi1_u(i,j,k)*deta_xi3_u(i,j,k)*               &
                       ( intw_p2u(i,1)*rho_wet(i,j,k) +                 &
                         intw_p2u(i,2)*rho_wet(i+1,j,k) )*u(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO


  CALL eg_swap_bounds( ustar, udims_s, fld_type_u, .TRUE. )
  CALL eg_swap_bounds( vstar, vdims_s, fld_type_v, .TRUE.)

! Bottom BC: rho(k=0) = rho(k=1/2)

  k = 0
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      wstar(i,j,k) = dxi1_p(i)*dxi2_p(j)*h1_p_eta(i,j,k)*             &
                     h2_p_eta(i,j,k)*rho_wet(i,j,k+1)*w(i,j,k)
    END DO
  END DO

! Top BC:
  k = model_levels
  wstar(:,:,k) = 0.0

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,pdims,dxi1_p,dxi2_p,h1_p_eta,h2_p_eta,     &
!$OMP& intw_rho2w,rho_wet,w,wstar)
  DO k=1, model_levels-1
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        wstar(i,j,k) = dxi1_p(i)*dxi2_p(j)*h1_p_eta(i,j,k)*             &
                 h2_p_eta(i,j,k)*( intw_rho2w(k,1)*rho_wet(i,j,k+1) +   &
                 intw_rho2w(k,2)*rho_wet(i,j,k) )*w(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Compute work1 = (<vstar>^xi2)*f3_star - (<wstar>^eta)*f2_star

IF( l_simple_coriolis ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,pdims,v,f3_comp,w,f2_comp,work1,           &
!$OMP&        intw_v2p,intw_w2rho)
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
         work1(i,j,k) = f3_comp(i,j)*                                 &
                          ( intw_v2p(j,1)*v(i,j-1,k)                  &
                           +intw_v2p(j,2)*v(i,j,k) )                  &
                       -f2_comp(i,j)*                                 &
                          ( intw_w2rho(k,1)*w(i,j,k)                  &
                           +intw_w2rho(k,2)*w(i,j,k-1) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,pdims,vstar,f3_star,wstar,f2_star,work1)
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        work1(i,j,k) = 0.5*(vstar(i,j,k)+vstar(i,j-1,k))*               &
                f3_star(i,j,k) - 0.5*(wstar(i,j,k)+wstar(i,j,k-1))*     &
                f2_star(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
END IF

CALL eg_swap_bounds( work1,wdims_s,fld_type_p,.TRUE. )

IF( l_simple_coriolis ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)         &
!$OMP& SHARED(model_levels,udims,u,r_u,beta_u_dt,work1,intw_p2u)
  DO k=1, model_levels
    DO j=udims%j_start, udims%j_end
      DO i=udims%i_start, udims%i_end
        r_u(i,j,k) = u(i,j,k) + r_u(i,j,k) + beta_u_dt *                &
                     ( intw_p2u(i,1)*work1(i,j,k)                       &
                      +intw_p2u(i,2)*work1(i+1,j,k) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,udims,u,r_u,beta_u_dt,work1,h1_xi1_u,dxi1_u)
  DO k=1, model_levels
    DO j=udims%j_start, udims%j_end
      DO i=udims%i_start, udims%i_end
        r_u(i,j,k) = u(i,j,k) + r_u(i,j,k) + beta_u_dt *                &
                     0.5*(work1(i,j,k)+work1(i+1,j,k)) /                &
                     ( h1_xi1_u(i,j,k)*dxi1_u(i) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Compute work1 = (<wstar>^eta)*f1_star - (<ustar>^xi1)*f3_star

IF( l_simple_coriolis ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,pdims,w,f1_comp,u,f3_comp,work1,           &
!$OMP&        intw_w2rho,intw_u2p)
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        work1(i,j,k) = f1_comp(i,j)*                                   &
                          ( intw_w2rho(k,1)*w(i,j,k)                   &
                           +intw_w2rho(k,2)*w(i,j,k-1) )               &
                      -f3_comp(i,j)*                                   &
                          ( intw_u2p(i,1)*u(i-1,j,k)                   &
                           +intw_u2p(i,2)*u(i,j,k) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,pdims,wstar,f1_star,ustar,f3_star,work1)    
  DO k=1, model_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        work1(i,j,k) = 0.5*(wstar(i,j,k)+wstar(i,j,k-1))*               &
                       f1_star(i,j,k) -                                 &
                       0.5*(ustar(i,j,k)+ustar(i-1,j,k))*               &
                       f3_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

CALL eg_swap_bounds( work1,wdims_s,fld_type_p,.TRUE.)

IF( l_simple_coriolis ) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,vdims,r_v,beta_v_dt,work1,v,intw_p2v)
  DO k=1, model_levels
    DO j=vdims%j_start, vdims%j_end
      DO i=vdims%i_start, vdims%i_end
        r_v(i,j,k) = v(i,j,k) + r_v(i,j,k) + beta_v_dt *                &
                     ( intw_p2v(j,1)*work1(i,j,k)                       &
                      +intw_p2v(j,2)*work1(i,j+1,k) )
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,vdims,r_v,beta_v_dt,h2_xi2_v,dxi2_v,work1,v)  
  DO k=1, model_levels
    DO j=vdims%j_start, vdims%j_end
      DO i=vdims%i_start, vdims%i_end
        r_v(i,j,k) = v(i,j,k) + r_v(i,j,k) + beta_v_dt *                &
                     0.5*(work1(i,j,k)+work1(i,j+1,k)) /                &
                     ( h2_xi2_v(i,j,k)*dxi2_v(j) )
    END DO
  END DO
END DO
!$OMP END PARALLEL DO
END IF

!----------------------------------------------------------------------
! Compute w-Coriolis acceleration term (Omega x U)w eqns:
!   (7.12) ENDGame formulation vn1.01.
!----------------------------------------------------------------------


! work = (<ustar>^xi1)*f2_star - (<vstar>^xi2)*f1_star

! approximate ustar, vstar on surface with ustar(1/2), vstar(1/2).
! Coriolis terms are also approximated with Coriolis(1/2).
! This is consistent with the approximation done for rho at the
! surface: rho(surf) ~ rho(1/2)

IF( l_simple_coriolis ) THEN
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      work1(i,j,1)= f2_comp(i,j)*                                       &
                       ( intw_u2p(i,1)*u(i-1,j,1)                       &
                        +intw_u2p(i,2)*u(i,j,1) )                       &
                   -f1_comp(i,j)*                                       &
                       ( intw_v2p(j,1)*v(i,j-1,1)                       &
                        +intw_v2p(j,2)*v(i,j,1) )
    END DO
  END DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,kp1)     &
!$OMP& SHARED(model_levels,pdims,r_w,beta_w_dt,work1,u,v,w,ih,          &
!$OMP&        f1_comp,f2_comp,intw_u2p,intw_v2p,intw_rho2w)
  DO k=1, model_levels-1
    kp1 = k+1
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        work1(i,j,kp1) = f2_comp(i,j)*                                  &
                           ( intw_u2p(i,1)*u(i-1,j,kp1)                 &
                            +intw_u2p(i,2)*u(i,j,kp1) )                 &
                        -f1_comp(i,j)*                                  &
                           ( intw_v2p(j,1)*v(i,j-1,kp1)                 &
                            +intw_v2p(j,2)*v(i,j,kp1) )
! work1(k) already computed
        r_w(i,j,k) = r_w(i,j,k) + ih*w(i,j,k) +                         &
                     beta_w_dt*( intw_rho2w(k,1)*work1(i,j,kp1)         &
                                +intw_rho2w(k,2)*work1(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      work1(i,j,1)= 0.5*(ustar(i,j,1)+ustar(i-1,j,1))*                  &
                        f2_star(i,j,1) -                                &
                    0.5*(vstar(i,j,1)+vstar(i,j-1,1))*                  &
                        f1_star(i,j,1)
    END DO
  END DO

  DO k=1, model_levels-1
    kp1 = k+1
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        work1(i,j,kp1) = 0.5*(ustar(i,j,kp1)+ustar(i-1,j,kp1))*         &
                        f2_star(i,j,kp1) -                              &
                        0.5*(vstar(i,j,kp1)+vstar(i,j-1,kp1))*          &
                        f1_star(i,j,kp1)
! work1(k) already computed
        r_w(i,j,k) = r_w(i,j,k) + ih*w(i,j,k) +                         &
                     beta_w_dt*0.5*(work1(i,j,kp1) +                    &
                     work1(i,j,k))/(h3_p_eta(i,j,k)*                    &
                     deta_xi3_theta(i,j,k)*deta_w(k))
      END DO
    END DO
  END DO
END IF 

!----------------------------------------------------------------------
! Compute components of pressure gradient in Exner-thetav form:
!   (7.2) ENDGame formulation vn2.02.
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! U Component of Pressure Gradient in R_u Eqn(7.7) of EG2.02
!----------------------------------------------------------------------

! work1: avg(thetav)^(eta) in pressure grad coeff of eqn 7.7.
! work2: avg(exner)^(xi1 eta)*(d(xi3)/d(xi1)) (last term of eqn 7.7)

! *** NB. (1) c_pd=CP=const is assumed and taken outside averaging.
!         (2) It is assumed that horizontal and vertical averages
!             commute, ie. avg(X)^(xi eta) = avg(X)^(eta xi)

! Calculate averaged thetav
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,pdims,work1,intw_w2rho,theta_wet)            
DO k=1, model_levels
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      work1(i,j,k) = intw_w2rho(k,1)*theta_wet(i,j,k) +               &
                     intw_w2rho(k,2)*theta_wet(i,j,k-1)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

CALL eg_swap_bounds( work1,wdims_s,fld_type_p,.FALSE. )

! Calculate vertical Exner gradient

!CDIR NOUNROLL
DO j=udims%j_start, udims%j_end
  DO i=udims%i_start, udims%i_end
    work2(i,j,0) = (intw_p2u(i,1)*exner_star(i,j) +                   &
                    intw_p2u(i,2)*exner_star(i+1,j)) *                &
                   (xi3_at_theta(i+1,j,0) -                           &
                    xi3_at_theta(i,j,0)) / dxi1_u(i)
  END DO
END DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,udims,work2,intw_rho2w,intw_p2u,exner,     &
!$OMP& xi3_at_theta,dxi1_u)
DO k=1, model_levels-1
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end
      work2(i,j,k) = (intw_rho2w(k,1) *                               &
                      (intw_p2u(i,1)*exner(i,j,k+1) +                 &
                       intw_p2u(i,2)*exner(i+1,j,k+1)) +              &
                      intw_rho2w(k,2) *                               &
                      (intw_p2u(i,1)*exner(i,j,k)+                    &
                       intw_p2u(i,2)*exner(i+1,j,k))) *               &
                     (xi3_at_theta(i+1,j,k) -                         &
                      xi3_at_theta(i,j,k)) / dxi1_u(i)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

k = model_levels
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end
      work2(i,j,k) = 0.0
    END DO
  END DO

! Add pressure gradient term to R_u

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,       &
!$OMP& p_grad_coeff,thetav_ave)                                       &
!$OMP& SHARED(model_levels,udims,h1_xi1_u,deta_xi3_u,intw_p2u,        &
!$OMP& work1,r_u,beta_u_dt,exner,deta_xi3,dxi1_u,work2,deta_rho)
DO k=1, model_levels
  DO j=udims%j_start, udims%j_end
    DO i=udims%i_start, udims%i_end

      p_grad_coeff = cp / ( h1_xi1_u(i,j,k)*deta_xi3_u(i,j,k) )

      thetav_ave = intw_p2u(i,1)*work1(i,j,k) +                       &
                   intw_p2u(i,2)*work1(i+1,j,k)

      r_u(i,j,k) = r_u(i,j,k) -                                       &
                   beta_u_dt * p_grad_coeff * thetav_ave *            &
                   ( (exner(i+1,j,k)*deta_xi3(i+1,j,k) -              &
                      exner(i,j,k)*deta_xi3(i,j,k)) / dxi1_u(i) -     &
                     (work2(i,j,k)-work2(i,j,k-1)) / deta_rho(k) )


    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!----------------------------------------------------------------------
! V Component of Pressure Gradient in R_v Eqn(7.8) of EG2.02
!----------------------------------------------------------------------

! work2: avg(exner)^(xi2 eta)*(d(xi3)/d(xi2)) (last term of eqn 7.8)

! Calculate vertical Exner gradient

k = 0
!CDIR NOUNROLL
  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
      work2(i,j,k) = (intw_p2v(j,1)*exner_star(i,j) +                 &
                      intw_p2v(j,2)*exner_star(i,j+1)) *              &
                     (xi3_at_theta(i,j+1,k) -                         &
                      xi3_at_theta(i,j,k)) / dxi2_v(j)
    END DO
  END DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(model_levels,vdims,work2,intw_rho2w,intw_p2v,exner,     &
!$OMP& xi3_at_theta,dxi2_v)
DO k=1, model_levels-1
  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
      work2(i,j,k) = (intw_rho2w(k,1) *                               &
                      (intw_p2v(j,1)*exner(i,j,k+1) +                 &
                       intw_p2v(j,2)*exner(i,j+1,k+1)) +              &
                      intw_rho2w(k,2) *                               &
                      (intw_p2v(j,1)*exner(i,j,k) +                   &
                       intw_p2v(j,2)*exner(i,j+1,k))) *               &
                     (xi3_at_theta(i,j+1,k) -                         &
                      xi3_at_theta(i,j,k)) / dxi2_v(j)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

k = model_levels
  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end
      work2(i,j,k) = 0.0
    END DO
  END DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,       &
!$OMP& p_grad_coeff,thetav_ave)                                       &
!$OMP& SHARED(model_levels,vdims,h2_xi2_v,deta_xi3_v,intw_p2v,        &
!$OMP& work1,r_v,beta_v_dt,exner,deta_xi3,dxi2_v,work2,deta_rho)
DO k=1, model_levels
  DO j=vdims%j_start, vdims%j_end
    DO i=vdims%i_start, vdims%i_end

      p_grad_coeff = cp / ( h2_xi2_v(i,j,k)*deta_xi3_v(i,j,k) )

      thetav_ave = intw_p2v(j,1)*work1(i,j,k) +                       &
                   intw_p2v(j,2)*work1(i,j+1,k)

      r_v(i,j,k) = r_v(i,j,k) -                                       &
                   beta_v_dt * p_grad_coeff * thetav_ave *            &
                   ( (exner(i,j+1,k)*deta_xi3(i,j+1,k) -              &
                      exner(i,j,k)*deta_xi3(i,j,k)) / dxi2_v(j) -     &
                     (work2(i,j,k)-work2(i,j,k-1)) / deta_rho(k) )

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!----------------------------------------------------------------------
! W Component of Pressure Gradient in R_w Eqn(7.9) of EG2.02
!----------------------------------------------------------------------

! vpg term at bottom, k=0

  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end

      k = 0
      r_w(i,j,k) = r_w(i,j,k) + ih*w(i,j,k) + beta_w_dt*psi_w_surf(i,j)

      k = model_levels
      r_w(i,j,k) = r_w(i,j,k) + ih*w(i,j,k) + beta_w_dt*psi_w_lid(i,j)

    END DO
  END DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,       &
!$OMP& p_grad_coeff)                                                  &
!$OMP& SHARED(model_levels,pdims,h3_p_eta,deta_xi3_theta,             &
!$OMP& r_w,beta_w_dt,theta_wet,exner,deta_w,g_theta)          
DO k=1, model_levels-1
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end

      p_grad_coeff = cp/(h3_p_eta(i,j,k)*deta_xi3_theta(i,j,k))

      r_w(i,j,k) = r_w(i,j,k) -                                       &
                   beta_w_dt * (p_grad_coeff * theta_wet(i,j,k)*      &
                   (exner(i,j,k+1)-exner(i,j,k)) / deta_w(k) +        &
                   g_theta(i,j,k))

    END DO
  END DO
END DO
!$OMP END PARALLEL DO


! Zero R_v on bloundary for channel flows
 
IF( model_domain == mt_cyclic_lam ) THEN
   IF( at_extremity(Psouth) ) THEN
      R_v(:,vdims%j_start,:) = 0.0
   END IF
   IF( at_extremity(Pnorth) ) THEN
      R_v(:,vdims%j_end,:) = 0.0
   END IF
END IF

IF (lhook) CALL dr_hook('EG_SISL_INIT_UVW',zhook_out,zhook_handle)

END SUBROUTINE eg_sisl_init_uvw
END MODULE eg_sisl_init_uvw_mod
