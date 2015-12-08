! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine BDY_IMPL3
!
!  Purpose: Calculate downward sweep of matrix for increments to
!           U, V, T and Q in the boundary layer for the
!           unconditionally stable and non-oscillatory numerical solver
!
!  Programming standard: UMDP3
!
!  Documentation: UMDP24
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE bdy_impl3 (                                                  &
! IN levels/switches
 bl_levels, l_correct,                                                  &
! IN fields
 q,qcl,qcf,q_latest,qcl_latest,qcf_latest,t,t_latest,                   &
 dtrdz_charney_grid,dtrdz_u,dtrdz_v,rhokh,rhokm_u,rhokm_v,              &
 rdz_charney_grid,rdz_u,rdz_v,gamma1,gamma2,GAMMA,                      &
 du_nt,dv_nt,                                                           &
! INOUT fields
 fqw,ftl,tau_x,tau_y,du,dv,                                             &
! OUT fields
 dqw_nt,dtl_nt,qw,tl,dqw1,dtl1,ct_ctq,ctctq1,dqw,dtl,cq_cm_u,cq_cm_v    &
 )
  
  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE level_heights_mod, ONLY:                                          &
    r_theta_levels, r_rho_levels 
  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, udims_s, vdims_s, pdims, tdims
  USE atmos_constants_mod, ONLY: cp
  USE water_constants_mod, ONLY: lc, lf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE vectlib_mod, ONLY: oneover_v 
!$ USE omp_lib
  IMPLICIT NONE

! IN arrays
  LOGICAL, INTENT(IN) ::                                                &
   l_correct

  INTEGER, INTENT(IN) ::                                                &
    bl_levels                ! IN No. of atmospheric levels for
                                 !    which boundary layer fluxes are
                                 !    calculated.

  REAL, INTENT(IN) ::                                                   &
   gamma1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
   gamma2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                 ! IN new scheme weights.
   GAMMA(bl_levels)          ! IN standard implicit scheme weights.

  REAL, INTENT(IN) ::                                                   &
   q(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                                   ! IN specific humidity
   qcl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! IN Cloud liquid water
   qcf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! IN Cloud ice (kg per kg air)
   q_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                                   ! IN specific humidity
   qcl_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              bl_levels),                                               &
                                   ! IN Cloud liquid water
   qcf_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              bl_levels),                                               &
                                   ! IN Cloud ice (kg per kg air)
   t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                                   ! IN temperature
                                   !    Latest estimates to time
                                   !    level n+1 values
   t_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                                   ! IN temperature
   dtrdz_charney_grid(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end,bl_levels),             &
                                   ! IN dz for bottom BL_LEVELS
   dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                   ! IN -g.dt/dp for model wind layers
   dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                   ! IN -g.dt/dp for model wind layers
   rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         2:bl_levels),                                                  &
                                   ! IN Exchange coeff for FTL above
                                   !    surface.
   rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
            2:bl_levels),                                               &
                                   ! IN Exchange coefficients for
                                   !    momentum, on U-grid with
                                   !    first and last rows ignored.
                                   !    for K>=2 (from KMKH).
   rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
            2:bl_levels),                                               &
                                   ! IN Exchange coefficients for
                                   !    momentum, on V-grid with
                                   !    first and last rows ignored.
                                   !    for K>=2 (from KMKH).
   rdz_charney_grid(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    bl_levels),                                         &
                                   ! IN RDZ(,1) is the reciprocal of the
                                   ! height of level 1, i.e. of the
                                   ! middle of layer 1.  For K > 1,
                                   ! RDZ(,K) is the reciprocal
                                   ! of the vertical distance
                                   ! from level K-1 to level K.
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
          2:bl_levels),                                                 &
                                   ! IN Reciprocal of the vertical
                                   !    distance from level K-1 to
                                   !    level K. (K > 1) on wind levels
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels),                                                 &
                                   ! IN Reciprocal of the vertical
                                   !    distance from level K-1 to
                                   !    level K. (K > 1) on wind levels
   du_nt(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,   &
          bl_levels),                                                   &
                                   ! IN u non-turbulent increments.
   dv_nt(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,   &
          bl_levels)
                                   ! IN v non-turbulent increments.
! INOUT arrays
  REAL, INTENT(INOUT) ::                                                &
   fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! INOUT Flux of QW (ie., for surface,
                                   !    total evaporation). Kg/sq m/s
   ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! INOUT Flux of TL (ie., for surface,
                                   !    H/Cp where H is sensible heat
                                   !    in W per sq m).
   tau_x(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
         bl_levels),                                                    &
                                   ! INOUT x-component of turbulent
                                   !    stress at levels k-1/2;
                                   !    eg. TAUX(,1) is surface stress.
                                   !    U-grid, 1st and last rows set
                                   !    to "missing data". (N/sq m)
                                   !    IN as "explicit" fluxes from
                                   !    ex_flux_uv, OUT as "implicit
   tau_y(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
         bl_levels),                                                    &
                                   ! INOUT y-component of turbulent
                                   !    stress at levels k-1/2;
                                   !    eg. TAUX(,1) is surface stress.
                                   !    V-grid, 1st and last rows set
                                   !    to "missing data". (N/sq m)
                                   !    IN as "explicit" fluxes from
                                   !    ex_flux_uv, OUT as "implicit
   du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,      &
        bl_levels),                                                     &
                                   ! INOUT BL increment to u wind field
   dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,      &
        bl_levels)
                                   ! INOUT BL increment to v wind field

! OUT arrays
  REAL, INTENT(OUT) ::                                                  &
   qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),  &
                                   ! OUT total water
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),  &
                                   ! OUT liquid water temperature
   ct_ctq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                   ! OUT Coefficient in T and q
                                   !     tri-diagonal implicit matrix
   ctctq1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
   cq_cm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                   ! OUT Coefficient in U and V
                                   !     tri-diagonal implicit matrix
   cq_cm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                   ! OUT Coefficient in U and V
                                   !     tri-diagonal implicit matrix
   dqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! OUT BL increment to q field
   dtl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! OUT BL increment to T field
   dqw_nt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                         ! NT incr to q field
   dtl_nt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                         ! NT incr to T field
   dqw1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                                   ! OUT 1 LEV BL increment to q field
   dtl1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                   ! OUT 1 LEV BL increment to T field

! LOCAL arrays
  REAL ::                                                               &
   r_theta_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,       &
             0:bl_levels),                                              &
                                               ! Vertical grids for U
   r_theta_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,       &
             0:bl_levels)! and V flux levels                            

  REAL ::                                                               &
  temp(pdims%i_end*pdims%j_end),                                        &
                            ! temp for pressure grid vector division
  temp_u((udims%i_end-udims%i_start+1)*(udims%j_end-udims%j_start+1)),  &
                                   ! temp for u grid vector division
  temp_v((vdims%i_end-vdims%i_start+1)*(vdims%j_end-vdims%j_start+1))
                                   ! temp for v grid vector division
!  Local scalars :-
  REAL ::                                                               &
   at,                                                                  &
                ! Matrix element in "T" row in eqn P244.79.
   rbt,                                                                 &
                ! Reciprocal of BT' (eqns P244.107, 110, 113).
   am,                                                                  &
                ! Matrix element in eqn P244.80.
   rbm,                                                                 &
                ! Reciprocal of BM(') (eqns P244.81, 85, 89).
   gamma1_uv,                                                           &
   gamma2_uv,                                                           &
                ! gamma1 and gamma2 shifted to u or v points
   r_sq,                                                                &
                ! square of height variables
   rr_sq    ! 1/square of height variables

  INTEGER ::                                                            &
   blm1,                                                                &
                ! BL_LEVELS minus 1.
   i,j,                                                                 &
                ! Loop counter (horizontal field index).
   k,                                                                   &
                ! Loop counter (vertical index).
   omp_block,                                                           &
                ! omp block length
   jj,                                                                  &
                ! omp block loop counter                         
   l
                ! vector counter
! Derived local parameters.

  REAL :: lcrcp,ls,lsrcp

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  PARAMETER (                                                           &
   lcrcp=lc/cp,                                                         &
                             ! Evaporation-to-dT conversion factor.
   ls=lf+lc,                                                            &
                             ! Latent heat of sublimation.
   lsrcp=ls/cp                                                          &
                             ! Sublimation-to-dT conversion factor.
    )

  IF (lhook) CALL dr_hook('BDY_IMPL3',zhook_in,zhook_handle)

!$OMP  PARALLEL DEFAULT(NONE) SHARED(l_correct,bl_levels,tdims,          &
!$OMP& dqw_nt,dtl_nt,q_latest,qcl_latest, dtrdz_v,dtrdz_u,udims, rdz_v,  &
!$OMP& gamma1,q,qcl,qcf,t_latest,t,ftl,rhokh,dtl,rdz_charney_grid,dqw,   &
!$OMP& tau_x,rhokm_u,du,rdz_u,vdims,tau_y,dv, qcf_latest,                &
!$OMP& qw,tl,r_theta_levels,r_theta_u,r_theta_v,r_rho_levels,fqw,        &
!$OMP& dtrdz_charney_grid,gamma2,ct_ctq,dqw1,dtl1,ctctq1,                &
!$OMP& gamma,cq_cm_u,cq_cm_v,du_nt,dv_nt,rhokm_v,l_vatpoles)             &
!$OMP& PRIVATE(k,j,i,r_sq,rbt,temp,temp_u,temp_v,l,                      &
!$OMP& at,blm1,am,rbm,rr_sq,jj,omp_block,gamma1_uv,gamma2_uv)  

  IF ( l_correct ) THEN

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, bl_levels
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
! Don't use QW, TL here as these are no longer at time level n
          dqw_nt(i,j,k) = q_latest(i,j,k) + qcl_latest(i,j,k)           &
                        + qcf_latest(i,j,k)                             &
                        - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k)
          dtl_nt(i,j,k) = t_latest(i,j,k)                               &
               - lcrcp * qcl_latest(i,j,k)                              &
               - lsrcp * qcf_latest(i,j,k)                              &
               - ( t(i,j,k) - lcrcp*qcl(i,j,k) - lsrcp*qcf(i,j,k) )
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

! Update explicit fluxes using predictor X* value as needed by the
! 2nd stage of the scheme. Note that: DTL=TL*-TL, DQW=QW*-QW etc

!$OMP DO SCHEDULE(STATIC)
    DO k = 2, bl_levels
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          ftl(i,j,k) = ftl(i,j,k) - rhokh(i,j,k) *                      &
            ( dtl(i,j,k) - dtl(i,j,k-1) ) * rdz_charney_grid(i,j,k)
          fqw(i,j,k) = fqw(i,j,k) - rhokh(i,j,k) *                      &
            ( dqw(i,j,k) - dqw(i,j,k-1) ) * rdz_charney_grid(i,j,k)
        END DO
      END DO

      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          tau_x(i,j,k) = tau_x(i,j,k) + rhokm_u(i,j,k) *                &
                      ( du(i,j,k) - du(i,j,k-1) ) *rdz_u(i,j,k)
        END DO
      END DO

      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          tau_y(i,j,k) = tau_y(i,j,k) + rhokm_v(i,j,k) *                &
                      ( dv(i,j,k) - dv(i,j,k-1) ) *rdz_v(i,j,k)
        END DO
      END DO

    END DO
!$OMP END DO NOWAIT

  ELSE
     
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, bl_levels
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          qw(i,j,k) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)
          tl(i,j,k) = t(i,j,k) - lcrcp*qcl(i,j,k) - lsrcp*qcf(i,j,k)
          dqw_nt(i,j,k) = q_latest(i,j,k) + qcl_latest(i,j,k)           &
                          + qcf_latest(i,j,k) - qw(i,j,k)
          dtl_nt(i,j,k) = t_latest(i,j,k)                               &
                          - lcrcp * qcl_latest(i,j,k)                   &
                          - lsrcp * qcf_latest(i,j,k)                   &
                          - tl(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  blm1 = bl_levels-1

!-----------------------------------------------------------------------
!!  1.0 Interpolate r_theta_levels to U,V columns
!-----------------------------------------------------------------------


! inlined p_to_u for OpenMP

!$OMP DO SCHEDULE(STATIC)
  DO k=0,bl_levels

    DO j=udims%j_start, udims%j_end
      DO i= udims%i_start, udims%i_end

        r_theta_u(i,j,k)= 0.5*                                          &
      ( r_theta_levels(i,j,k) + r_theta_levels(i+1,j,k) )

      END DO
    END DO

    DO j=vdims%j_start, vdims%j_end
      DO i=vdims%i_start, vdims%i_end

        r_theta_v(i,j,k)= 0.5 *                                         &
      ( r_theta_levels(i,j,k) + r_theta_levels(i,j+1,k) )

      END DO
    END DO

  END DO
!$OMP END DO      
      

!-----------------------------------------------------------------------
! 2.0 For simulations on a sphere we must use spherical geometry for
!     vertical flux-divergences.   Thus, leaving out rho for
!     simplicity, the standard cartesian flux-divergence:
!          dQ(K)/dt = -(FQ(K+1)-FQ(K))/DZ
!     becomes:
!          dQ(K)/dt = -(r_flux(K+1)^2*FQ(K+1)-r_flux(K)^2*FQ(K))
!                      / (r_full(K)^2 * DZ)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!!  3.0 Calculate matrix elements
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
! Include non-turbulent increments.
      r_sq = r_rho_levels(i,j,bl_levels)*r_rho_levels(i,j,bl_levels)
      dqw(i,j,bl_levels) = ( dtrdz_charney_grid(i,j,bl_levels) *        &
                      (r_sq * fqw(i,j,bl_levels)) +                     &
      dqw_nt(i,j,bl_levels) ) * gamma2(i,j)
      dtl(i,j,bl_levels) = ( dtrdz_charney_grid(i,j,bl_levels) *        &
                     (r_sq * ftl(i,j,bl_levels)) + dtl_nt(i,j,bl_levels)&
                           ) * gamma2(i,j)
      ct_ctq(i,j,bl_levels) = -dtrdz_charney_grid(i,j,bl_levels) *      &
             gamma1(i,j)*(rhokh(i,j,bl_levels)*r_sq)*                   &
              rdz_charney_grid(i,j,bl_levels)
      rbt = 1.0 / ( 1.0 - ct_ctq(i,j,bl_levels) )
      dqw(i,j,bl_levels) = rbt * dqw(i,j,bl_levels)
      dtl(i,j,bl_levels) = rbt * dtl(i,j,bl_levels)
      ct_ctq(i,j,bl_levels) = rbt * ct_ctq(i,j,bl_levels)
    END DO
  END DO
!$OMP END DO

  omp_block = tdims%j_end
!$ omp_block = CEILING(tdims%j_end/REAL(omp_get_num_threads()))

!$OMP DO SCHEDULE(STATIC)
  DO jj  = tdims%j_start, tdims%j_end, omp_block
    DO k = blm1, 2, -1
      l = 0
      DO j = jj, MIN(jj+omp_block-1, tdims%j_end)
        DO i = tdims%i_start,tdims%i_end
          r_sq = r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
          rr_sq = r_rho_levels(i,j,k+1)*r_rho_levels(i,j,k+1)
          dqw(i,j,k) = ( -dtrdz_charney_grid(i,j,k)*                    &
               ((rr_sq*fqw(i,j,k+1))-(r_sq*fqw(i,j,k)))+dqw_nt(i,j,k) ) &
                 *gamma2(i,j)
          dtl(i,j,k) = ( -dtrdz_charney_grid(i,j,k)*                    &
               ((rr_sq*ftl(i,j,k+1))-(r_sq*ftl(i,j,k)))+dtl_nt(i,j,k) ) &
                 *gamma2(i,j)
          at = -dtrdz_charney_grid(i,j,k) *                             &
               gamma1(i,j)*(rr_sq*rhokh(i,j,k+1))*                      &
               rdz_charney_grid(i,j,k+1)
          ct_ctq(i,j,k) = -dtrdz_charney_grid(i,j,k) *                  &
               gamma1(i,j)*(r_sq*rhokh(i,j,k))*rdz_charney_grid(i,j,k)
          l = l + 1
          temp(l) = ( 1.0 - ct_ctq(i,j,k) -                             &
               at*( 1.0 + ct_ctq(i,j,k+1) ) )
          dqw(i,j,k) = (dqw(i,j,k) - at*dqw(i,j,k+1) )
          dtl(i,j,k) = (dtl(i,j,k) - at*dtl(i,j,k+1) )
        END DO 
      END DO

      CALL ONEOVER_V(l, temp, temp)

      l = 0
      DO j = jj, MIN(jj+omp_block-1, tdims%j_end)
        DO i = tdims%i_start,tdims%i_end     
          l = l + 1
          dqw(i,j,k) = temp(l) * dqw(i,j,k)
          dtl(i,j,k) = temp(l) * dtl(i,j,k) 
          ct_ctq(i,j,k) = temp(l) * ct_ctq(i,j,k)
        END DO
      END DO
      
    END DO !blm1,2,-1
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!  Bottom model layer QW row of matrix equation.
!-----------------------------------------------------------------------

  IF ( .NOT. l_correct ) THEN

!-----------------------------------------------------------------------
! The following calculations are only done on the 1st stage (predictor).
! Their purpose is to compute the surface scalar (T, Q) increments which
! are needed by the surface scheme to compute the implicit scalar fluxes.
! The same implicit scalar fluxes are used as a boundary condition for
! discrete equations of the 2nd stage.
! Due to the dependency to the surface scalar fluxes, the 1st stage
! downward sweep remains incomplete. It is completed at the
! beginning of bdy_impl4, since there the surface scalar fluxes
! have been fully updated. This is done by sf_impl2().
! NOTE: The standard scheme solver is used for this calculation.
!       Incorporation of the new scheme for the scalar surface variables
!       would have been a more preferable choice but currently not
!       available.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
! Include non-turbulent increments.
        r_sq = r_rho_levels(i,j,bl_levels)*r_rho_levels(i,j,bl_levels)
        dqw1(i,j,bl_levels) = dtrdz_charney_grid(i,j,bl_levels)*        &
                         (r_sq*fqw(i,j,bl_levels)) +                    &
                         dqw_nt(i,j,bl_levels)
        dtl1(i,j,bl_levels) = dtrdz_charney_grid(i,j,bl_levels)*        &
                         (r_sq*ftl(i,j,bl_levels)) +                    &
                         dtl_nt(i,j,bl_levels)
        ctctq1(i,j,bl_levels) = -dtrdz_charney_grid(i,j,bl_levels)*     &
           GAMMA(bl_levels)*r_sq*rhokh(i,j,bl_levels)*                  &
           rdz_charney_grid(i,j,bl_levels)
        rbt = 1.0 / ( 1.0 - ctctq1(i,j,bl_levels) )
        dqw1(i,j,bl_levels) = rbt * dqw(i,j,bl_levels)
        dtl1(i,j,bl_levels) = rbt * dtl(i,j,bl_levels)
        ctctq1(i,j,bl_levels) = rbt * ctctq1(i,j,bl_levels)
      END DO
    END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
    DO jj  = tdims%j_start,tdims%j_end, omp_block
      DO k = blm1, 2, -1
        l = 0
        DO j = jj, MIN(jj+omp_block-1, tdims%j_end)
          DO i = tdims%i_start,tdims%i_end
            r_sq = r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
            rr_sq = r_rho_levels(i,j,k+1)*r_rho_levels(i,j,k+1)
            dqw1(i,j,k) = -dtrdz_charney_grid(i,j,k) *                  &
              ((rr_sq*fqw(i,j,k+1)) - (r_sq*fqw(i,j,k))) + dqw_nt(i,j,k)
            dtl1(i,j,k) = -dtrdz_charney_grid(i,j,k) *                  &
              ((rr_sq*ftl(i,j,k+1)) - (r_sq*ftl(i,j,k))) + dtl_nt(i,j,k)
            at = -dtrdz_charney_grid(i,j,k) *                           &
              GAMMA(k+1)*(rr_sq*rhokh(i,j,k+1))*rdz_charney_grid(i,j,k+1)
            ctctq1(i,j,k) = -dtrdz_charney_grid(i,j,k) *                &
              GAMMA(k)*(r_sq*rhokh(i,j,k))*rdz_charney_grid(i,j,k)
! pack
            l = l + 1
            temp(l) = ( 1.0 - ctctq1(i,j,k) -                           &
                 at*( 1.0 + ctctq1(i,j,k+1) ) )   
            dqw1(i,j,k) =  (dqw1(i,j,k) - at*dqw1(i,j,k+1) )
            dtl1(i,j,k) =  (dtl1(i,j,k) - at*dtl1(i,j,k+1) )
          END DO
        END DO

        CALL ONEOVER_V(l, temp, temp)
        l = 0
        DO j = jj, MIN(jj+omp_block-1, tdims%j_end)
          DO i = tdims%i_start,tdims%i_end
            l = l + 1
            dqw1(i,j,k) = temp(l) * dqw1(i,j,k)
            dtl1(i,j,k) = temp(l) * dtl1(i,j,k) 
            ctctq1(i,j,k) = temp(l) * ctctq1(i,j,k)
          END DO
        END DO
      END DO !blm1,2,-1
    END DO
!$OMP END DO
       
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        r_sq = r_rho_levels(i,j,2)*r_rho_levels(i,j,2)
        dqw1(i,j,1) = -dtrdz_charney_grid(i,j,1) * (r_sq*fqw(i,j,2)) +  &
                      dqw_nt(i,j,1)
        dtl1(i,j,1) = -dtrdz_charney_grid(i,j,1) * (r_sq*ftl(i,j,2)) +  &
                      dtl_nt(i,j,1)
        at = -dtrdz_charney_grid(i,j,1) *                               &
                   GAMMA(2)*(r_sq*rhokh(i,j,2))*rdz_charney_grid(i,j,2)
        rbt = 1.0 / ( 1.0 - at*( 1.0 + ctctq1(i,j,2) ) )
        dqw1(i,j,1) = rbt * (dqw1(i,j,1) - at*dqw1(i,j,2) )
        dtl1(i,j,1) = rbt * (dtl1(i,j,1) - at*dtl1(i,j,2) )

! Now set CT_CTQ(1) to be r^2 * BETA
        r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
        ctctq1(i,j,1) = - (r_sq * dtrdz_charney_grid(i,j,1)) * rbt
      END DO
    END DO
!$OMP END DO

  ELSE
!-----------------------------------------------------------------------

! The following calculations complete the downward sweep for  the surface
! scalar variables. They apply on the 2nd stage of the scheme.
! The equivalent calculations for the 1st stage are done at the
! beginning of bdy_impl4, since at this stage (bdy_impl3) the surface
! scalar fluxes are not fully updated. This will be done by sf_impl2().

!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        r_sq = r_rho_levels(i,j,1)*r_rho_levels(i,j,1)
        rr_sq = r_rho_levels(i,j,2)*r_rho_levels(i,j,2)
        dqw(i,j,1) = gamma2(i,j) * ( -dtrdz_charney_grid(i,j,1) *       &
            ((rr_sq*fqw(i,j,2)) - (r_sq*fqw(i,j,1))) + dqw_nt(i,j,1) )
        dtl(i,j,1) = gamma2(i,j) * ( -dtrdz_charney_grid(i,j,1) *       &
            ((rr_sq*ftl(i,j,2)) - (r_sq*ftl(i,j,1))) + dtl_nt(i,j,1) )
        at = -dtrdz_charney_grid(i,j,1) *                               &
              gamma1(i,j)*(rr_sq*rhokh(i,j,2))*rdz_charney_grid(i,j,2)
        rbt = 1.0 / ( 1.0 - at*( 1.0 + ct_ctq(i,j,2) ) )
        dqw(i,j,1) = rbt * (dqw(i,j,1) - at*dqw(i,j,2) )
        dtl(i,j,1) = rbt * (dtl(i,j,1) - at*dtl(i,j,2) )

! Now set CT_CTQ(1) to be r^2 * BETA
        r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
        ct_ctq(i,j,1) = - (r_sq * dtrdz_charney_grid(i,j,1)) * rbt
      END DO
    END DO
!$OMP END DO

  END IF

IF (l_vatpoles) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
      gamma1_uv = gamma1(i+1,j)
      gamma2_uv = gamma2(i+1,j)

      r_sq           = r_theta_u(i,j,bl_levels-1)*r_theta_u(i,j,bl_levels-1)

      du(i,j,bl_levels) = -dtrdz_u(i,j,bl_levels)*(r_sq*tau_x(i,j,bl_levels))

! addition of non-turbulent increments

      du(i,j,bl_levels) = gamma2_uv * ( du(i,j,bl_levels)               &
                           + du_nt(i,j,bl_levels) )
      cq_cm_u(i,j,bl_levels) = -dtrdz_u(i,j,bl_levels) *                &
            gamma1_uv*(r_sq*rhokm_u(i,j,bl_levels))*rdz_u(i,j,bl_levels)
      rbm = 1.0 / ( 1.0 - cq_cm_u(i,j,bl_levels) )
      du(i,j,bl_levels) = rbm * du(i,j,bl_levels)
      cq_cm_u(i,j,bl_levels) = rbm * cq_cm_u(i,j,bl_levels)
    END DO
  END DO
!$OMP END DO
ELSE
!$OMP DO SCHEDULE(STATIC)
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
      gamma1_uv = gamma1(i,j)
      gamma2_uv = gamma2(i,j)

      r_sq           = r_theta_u(i,j,bl_levels-1)*r_theta_u(i,j,bl_levels-1)

      du(i,j,bl_levels) = -dtrdz_u(i,j,bl_levels)*(r_sq*tau_x(i,j,bl_levels))

! addition of non-turbulent increments

      du(i,j,bl_levels) = gamma2_uv * ( du(i,j,bl_levels)               &
                           + du_nt(i,j,bl_levels) )
      cq_cm_u(i,j,bl_levels) = -dtrdz_u(i,j,bl_levels) *                &
            gamma1_uv*(r_sq*rhokm_u(i,j,bl_levels))*rdz_u(i,j,bl_levels)
      rbm = 1.0 / ( 1.0 - cq_cm_u(i,j,bl_levels) )
      du(i,j,bl_levels) = rbm * du(i,j,bl_levels)
      cq_cm_u(i,j,bl_levels) = rbm * cq_cm_u(i,j,bl_levels)
    END DO
  END DO
!$OMP END DO
END IF ! vatpoles

  omp_block = udims%j_end - udims%j_start + 1
!$ omp_block=CEILING((udims%j_end-udims%j_start+1)/REAL(omp_get_num_threads()))
  
IF (l_vatpoles) THEN
!$OMP DO SCHEDULE(STATIC)
  DO jj  = udims%j_start, udims%j_end, omp_block
    DO k = blm1, 2, -1
      l = 0
      DO j = jj, MIN(jj+omp_block-1, udims%j_end)
        DO i = udims%i_start, udims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
          gamma1_uv = gamma1(i+1,j)
          gamma2_uv = gamma2(i+1,j)

          r_sq           = r_theta_u(i,j,k-1)*r_theta_u(i,j,k-1)
          rr_sq           = r_theta_u(i,j,k)*r_theta_u(i,j,k)  

          du(i,j,k) = dtrdz_u(i,j,k) *                                  &
                       ( (rr_sq*tau_x(i,j,k+1)) - (r_sq*tau_x(i,j,k) ))
! addition of non-turbulent increments
          du(i,j,k) = gamma2_uv * (du(i,j,k) + du_nt(i,j,k))
          am = -dtrdz_u(i,j,k) * gamma1_uv*(rr_sq*rhokm_u(i,j,k+1))*    &
               rdz_u(i,j,k+1)
          cq_cm_u(i,j,k) = -dtrdz_u(i,j,k)*(gamma1_uv*r_sq)*            &
               rhokm_u(i,j,k)*rdz_u(i,j,k)
! pack 
          l = l + 1
          temp_u(l) = ( 1.0 - cq_cm_u(i,j,k) -                          &
               am*( 1.0 + cq_cm_u(i,j,k+1) ) )
          du(i,j,k) = ( du(i,j,k) - am*du(i,j,k+1) )
        END DO 
      END DO 

      CALL ONEOVER_V(l, temp_u, temp_u)

      l = 0 
      DO j = jj, MIN(jj+omp_block-1, udims%j_end)
        DO i = udims%i_start, udims%i_end   
          l = l + 1
          du(i,j,k) = temp_u(l) * du(i,j,k) 
          cq_cm_u(i,j,k) = temp_u(l) * cq_cm_u(i,j,k)
        END DO
      END DO ! loop over 2,BLM1

    END DO
  END DO
!$OMP END DO
ELSE
!$OMP DO SCHEDULE(STATIC)
  DO jj  = udims%j_start, udims%j_end, omp_block
    DO k = blm1, 2, -1
      l = 0
      DO j = jj, MIN(jj+omp_block-1, udims%j_end)
        DO i = udims%i_start, udims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
          gamma1_uv = gamma1(i,j)
          gamma2_uv = gamma2(i,j)

          r_sq           = r_theta_u(i,j,k-1)*r_theta_u(i,j,k-1)
          rr_sq           = r_theta_u(i,j,k)*r_theta_u(i,j,k)  

          du(i,j,k) = dtrdz_u(i,j,k) *                                  &
                       ( (rr_sq*tau_x(i,j,k+1)) - (r_sq*tau_x(i,j,k) ))
! addition of non-turbulent increments
          du(i,j,k) = gamma2_uv * (du(i,j,k) + du_nt(i,j,k))
          am = -dtrdz_u(i,j,k) * gamma1_uv*(rr_sq*rhokm_u(i,j,k+1))*    &
               rdz_u(i,j,k+1)
          cq_cm_u(i,j,k) = -dtrdz_u(i,j,k)*(gamma1_uv*r_sq)*            &
               rhokm_u(i,j,k)*rdz_u(i,j,k)
! pack 
          l = l + 1
          temp_u(l) = ( 1.0 - cq_cm_u(i,j,k) -                          &
               am*( 1.0 + cq_cm_u(i,j,k+1) ) )
          du(i,j,k) = ( du(i,j,k) - am*du(i,j,k+1) )
        END DO 
      END DO 

      CALL ONEOVER_V(l, temp_u, temp_u)

      l = 0 
      DO j = jj, MIN(jj+omp_block-1, udims%j_end)
        DO i = udims%i_start, udims%i_end   
          l = l + 1
          du(i,j,k) = temp_u(l) * du(i,j,k) 
          cq_cm_u(i,j,k) = temp_u(l) * cq_cm_u(i,j,k)
        END DO
      END DO ! loop over 2,BLM1

    END DO
  END DO
!$OMP END DO
END IF ! vatpoles
 
IF (l_vatpoles) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
      gamma1_uv = gamma1(i+1,j)
      gamma2_uv = gamma2(i+1,j)

      r_sq           = r_theta_u(i,j,1)*r_theta_u(i,j,1)

      du(i,j,1) = dtrdz_u(i,j,1) * (r_sq * tau_x(i,j,2))

! addition of non-turbulent increments

      du(i,j,1) = gamma2_uv*(du(i,j,1) + du_nt(i,j,1))
      am = -dtrdz_u(i,j,1) * gamma1_uv*(r_sq*rhokm_u(i,j,2))            &
                   *rdz_u(i,j,2)
      rbm = 1.0 / ( 1.0 - am *( 1.0 + cq_cm_u(i,j,2) ) )
      du(i,j,1) = rbm * ( du(i,j,1) - am*du(i,j,2) )

! Now set CQ_CM_U(1) to be r^2 * BETA

      r_sq = r_theta_u(i,j,0)*r_theta_u(i,j,0)
      cq_cm_u(i,j,1) = (r_sq * dtrdz_u(i,j,1)) * rbm
    END DO
  END DO
!$OMP END DO NOWAIT
ELSE
!$OMP DO SCHEDULE(STATIC)
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
      gamma1_uv = gamma1(i,j)
      gamma2_uv = gamma2(i,j)

      r_sq           = r_theta_u(i,j,1)*r_theta_u(i,j,1)

      du(i,j,1) = dtrdz_u(i,j,1) * (r_sq * tau_x(i,j,2))

! addition of non-turbulent increments

      du(i,j,1) = gamma2_uv*(du(i,j,1) + du_nt(i,j,1))
      am = -dtrdz_u(i,j,1) * gamma1_uv*(r_sq*rhokm_u(i,j,2))            &
                   *rdz_u(i,j,2)
      rbm = 1.0 / ( 1.0 - am *( 1.0 + cq_cm_u(i,j,2) ) )
      du(i,j,1) = rbm * ( du(i,j,1) - am*du(i,j,2) )

! Now set CQ_CM_U(1) to be r^2 * BETA

      r_sq = r_theta_u(i,j,0)*r_theta_u(i,j,0)
      cq_cm_u(i,j,1) = (r_sq * dtrdz_u(i,j,1)) * rbm
    END DO
  END DO
!$OMP END DO NOWAIT
END IF ! vatpoles

IF (l_vatpoles) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
      IF ( j .EQ. vdims%j_end ) THEN
           gamma1_uv = gamma1(i,tdims%j_end)
           gamma2_uv = gamma2(i,tdims%j_end)
      ELSE
           gamma1_uv = gamma1(i,j+1)
           gamma2_uv = gamma2(i,j+1)
      END IF

      r_sq = r_theta_v(i,j,bl_levels-1)*r_theta_v(i,j,bl_levels-1)

      dv(i,j,bl_levels) = -dtrdz_v(i,j,bl_levels) *                     &
                         (r_sq*tau_y(i,j,bl_levels))
! addition of non-turbulent increments
      dv(i,j,bl_levels) = gamma2_uv * ( dv(i,j,bl_levels)               &
                           + dv_nt(i,j,bl_levels) )
      cq_cm_v(i,j,bl_levels) = -dtrdz_v(i,j,bl_levels) *                &
            gamma1_uv*(r_sq*rhokm_v(i,j,bl_levels))*rdz_v(i,j,bl_levels)
      rbm = 1.0 / ( 1.0 - cq_cm_v(i,j,bl_levels) )
      dv(i,j,bl_levels) = rbm * dv(i,j,bl_levels)
      cq_cm_v(i,j,bl_levels) = rbm * cq_cm_v(i,j,bl_levels)
    END DO
  END DO
!$OMP END DO 
ELSE
!$OMP DO SCHEDULE(STATIC)
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
      gamma1_uv = gamma1(i,j)
      gamma2_uv = gamma2(i,j)

      r_sq = r_theta_v(i,j,bl_levels-1)*r_theta_v(i,j,bl_levels-1)

      dv(i,j,bl_levels) = -dtrdz_v(i,j,bl_levels) *                     &
                         (r_sq*tau_y(i,j,bl_levels))
! addition of non-turbulent increments
      dv(i,j,bl_levels) = gamma2_uv * ( dv(i,j,bl_levels)               &
                           + dv_nt(i,j,bl_levels) )
      cq_cm_v(i,j,bl_levels) = -dtrdz_v(i,j,bl_levels) *                &
            gamma1_uv*(r_sq*rhokm_v(i,j,bl_levels))*rdz_v(i,j,bl_levels)
      rbm = 1.0 / ( 1.0 - cq_cm_v(i,j,bl_levels) )
      dv(i,j,bl_levels) = rbm * dv(i,j,bl_levels)
      cq_cm_v(i,j,bl_levels) = rbm * cq_cm_v(i,j,bl_levels)
    END DO
  END DO
!$OMP END DO 
END IF ! vatpoles

  omp_block = vdims%j_end - vdims%j_start + 1
!$ omp_block=CEILING((vdims%j_end-vdims%j_start+1)/REAL(omp_get_num_threads()))

IF (l_vatpoles) THEN
!$OMP DO SCHEDULE(STATIC)
  DO jj  = vdims%j_start, vdims%j_end, omp_block
    DO k = blm1, 2, -1
      l = 0
      DO j = jj, MIN(jj+omp_block-1, vdims%j_end)
        DO i = vdims%i_start, vdims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
        IF ( j .EQ. vdims%j_end ) THEN
           gamma1_uv = gamma1(i,tdims%j_end)
           gamma2_uv = gamma2(i,tdims%j_end)
        ELSE
           gamma1_uv = gamma1(i,j+1)
           gamma2_uv = gamma2(i,j+1)
        END IF

          r_sq = r_theta_v(i,j,k-1)*r_theta_v(i,j,k-1)
          rr_sq = r_theta_v(i,j,k)*r_theta_v(i,j,k)

          dv(i,j,k) = dtrdz_v(i,j,k) *                                  &
                      ( (rr_sq*tau_y(i,j,k+1)) - (r_sq*tau_y(i,j,k)) )
! addition of non-turbulent increments
          dv(i,j,k) = gamma2_uv * (dv(i,j,k) + dv_nt(i,j,k))
          am = -dtrdz_v(i,j,k) * gamma1_uv*(rr_sq*rhokm_v(i,j,k+1))*    &
                        rdz_v(i,j,k+1)          
          cq_cm_v(i,j,k) = -dtrdz_v(i,j,k)*gamma1_uv*(r_sq*             &
              rhokm_v(i,j,k))*rdz_v(i,j,k)
! pack
          l = l + 1
          temp_v(l) = ( 1.0 - cq_cm_v(i,j,k) -                          &
              am*( 1.0 + cq_cm_v(i,j,k+1) ) )
          dv(i,j,k) = ( dv(i,j,k) - am*dv(i,j,k+1) )
        END DO
      END DO

      CALL ONEOVER_V(l, temp_v, temp_v)
      
      l = 0
      DO j = jj, MIN(jj+omp_block-1, vdims%j_end)
        DO i = vdims%i_start, vdims%i_end
          l = l + 1
          dv(i,j,k) = temp_v(l) * dv(i,j,k)
          cq_cm_v(i,j,k) = temp_v(l) * cq_cm_v(i,j,k)
        END DO
      END DO
    END DO ! loop over 2,BLM1
  END DO
!$OMP END DO
ELSE
!$OMP DO SCHEDULE(STATIC)
  DO jj  = vdims%j_start, vdims%j_end, omp_block
    DO k = blm1, 2, -1
      l = 0
      DO j = jj, MIN(jj+omp_block-1, vdims%j_end)
        DO i = vdims%i_start, vdims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
          gamma1_uv = gamma1(i,j)
          gamma2_uv = gamma2(i,j)

          r_sq = r_theta_v(i,j,k-1)*r_theta_v(i,j,k-1)
          rr_sq = r_theta_v(i,j,k)*r_theta_v(i,j,k)

          dv(i,j,k) = dtrdz_v(i,j,k) *                                  &
                      ( (rr_sq*tau_y(i,j,k+1)) - (r_sq*tau_y(i,j,k)) )
! addition of non-turbulent increments
          dv(i,j,k) = gamma2_uv * (dv(i,j,k) + dv_nt(i,j,k))
          am = -dtrdz_v(i,j,k) * gamma1_uv*(rr_sq*rhokm_v(i,j,k+1))*    &
                        rdz_v(i,j,k+1)          
          cq_cm_v(i,j,k) = -dtrdz_v(i,j,k)*gamma1_uv*(r_sq*             &
              rhokm_v(i,j,k))*rdz_v(i,j,k)
! pack
          l = l + 1
          temp_v(l) = ( 1.0 - cq_cm_v(i,j,k) -                          &
              am*( 1.0 + cq_cm_v(i,j,k+1) ) )
          dv(i,j,k) = ( dv(i,j,k) - am*dv(i,j,k+1) )
        END DO
      END DO

      CALL ONEOVER_V(l, temp_v, temp_v)
      
      l = 0
      DO j = jj, MIN(jj+omp_block-1, vdims%j_end)
        DO i = vdims%i_start, vdims%i_end
          l = l + 1
          dv(i,j,k) = temp_v(l) * dv(i,j,k)
          cq_cm_v(i,j,k) = temp_v(l) * cq_cm_v(i,j,k)
        END DO
      END DO
    END DO ! loop over 2,BLM1
  END DO
!$OMP END DO
END IF ! vatpoles

IF (l_vatpoles) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
      IF ( j .EQ. vdims%j_end ) THEN
         gamma1_uv = gamma1(i,tdims%j_end)
         gamma2_uv = gamma2(i,tdims%j_end)
      ELSE
         gamma1_uv = gamma1(i,j+1)
         gamma2_uv = gamma2(i,j+1)
      END IF

      r_sq = r_theta_v(i,j,1)*r_theta_v(i,j,1)

      dv(i,j,1) = dtrdz_v(i,j,1) * (r_sq * tau_y(i,j,2))
! addition of non-turbulent increments
      dv(i,j,1) = gamma2_uv * (dv(i,j,1) + dv_nt(i,j,1))
      am = -dtrdz_v(i,j,1) * gamma1_uv * (r_sq * rhokm_v(i,j,2))          &
                   *rdz_v(i,j,2)
      rbm = 1.0 / ( 1.0 - am *( 1.0 + cq_cm_v(i,j,2) ) )
      dv(i,j,1) = rbm * ( dv(i,j,1) - am*dv(i,j,2) )
! Now set CQ_CM_V(1) to be r^2 * BETA
      r_sq = r_theta_v(i,j,0)*r_theta_v(i,j,0)
      cq_cm_v(i,j,1) = (r_sq * dtrdz_v(i,j,1)) * rbm
    END DO
  END DO
!$OMP END DO 
ELSE
!$OMP DO SCHEDULE(STATIC)
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
      gamma1_uv = gamma1(i,j)
      gamma2_uv = gamma2(i,j)

      r_sq = r_theta_v(i,j,1)*r_theta_v(i,j,1)

      dv(i,j,1) = dtrdz_v(i,j,1) * (r_sq * tau_y(i,j,2))
! addition of non-turbulent increments
      dv(i,j,1) = gamma2_uv * (dv(i,j,1) + dv_nt(i,j,1))
      am = -dtrdz_v(i,j,1) * gamma1_uv * (r_sq * rhokm_v(i,j,2))          &
                   *rdz_v(i,j,2)
      rbm = 1.0 / ( 1.0 - am *( 1.0 + cq_cm_v(i,j,2) ) )
      dv(i,j,1) = rbm * ( dv(i,j,1) - am*dv(i,j,2) )
! Now set CQ_CM_V(1) to be r^2 * BETA
      r_sq = r_theta_v(i,j,0)*r_theta_v(i,j,0)
      cq_cm_v(i,j,1) = (r_sq * dtrdz_v(i,j,1)) * rbm
    END DO
  END DO
!$OMP END DO 
END IF ! vatpoles

!-----------------------------------------------------------------------
!! 4.0 Return fluxes back to their true value by dividing by r*r
!-----------------------------------------------------------------------

!$OMP END PARALLEL

  IF (lhook) CALL dr_hook('BDY_IMPL3',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE bdy_impl3
