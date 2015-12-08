! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
! Subroutine BDY_IMPL4
!
!  Purpose: Calculate implicit correction to boundary layer fluxes of
!           heat, moisture and momentum for the unconditionally
!           stable and non-oscillatory numerical solver.
!
!
!  Programming standard: UMDP3
!
!  Documentation: UMDP24
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE bdy_impl4 (                                                  &
! IN levels, switches
 bl_levels, l_correct, l_ftl, l_fqw, l_taux, l_tauy,                    &
! IN data :
 gamma1,gamma2,rhokm_u,rhokm_v,rdz_charney_grid,                        &
 dtrdz_charney_grid,rdz_u,rdz_v,ct_ctq,cq_cm_u,cq_cm_v,dqw_nt,dtl_nt,   &
! INOUT data :
 qw,tl,fqw,ftl,tau_x,tau_y, fqw_star,ftl_star,taux_star,tauy_star,      &
 du,dv,du_star,dv_star, dqw,dtl, rhokh,                                 &
! OUT data, NB these are really tl and qt on exit!
 t_latest,q_latest,rhokh_mix                                            &
 )

  USE level_heights_mod, ONLY: r_rho_levels
  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, udims_s, vdims_s, tdims
  USE atmos_constants_mod, ONLY: cp
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

!  Inputs :-
  INTEGER, INTENT(IN) ::                                                &
   bl_levels                   ! IN Max. no. of "boundary" levels
!                                     allowed.
  LOGICAL, INTENT(IN) ::                                                &
   l_ftl,                                                               &
   l_fqw,                                                               &
   l_taux,                                                              &
   l_tauy,                                                              &
   l_correct

  REAL, INTENT(IN) ::                                                   &
   rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
            2:bl_levels),                                               &
                                   ! IN Exchange coefficients for U
   rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
            2:bl_levels),                                               &
                                   ! IN Exchange coefficients for V
   rdz_charney_grid(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    bl_levels),                                         &
                                   ! IN RDZ(,1) is the reciprocal of the
                                   ! height of level 1, i.e. of the
                                   ! middle of layer 1.  For K > 1,
                                   ! RDZ(,K) is the reciprocal
                                   ! of the vertical distance
                                   ! from level K-1 to level K.
   dtrdz_charney_grid(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end,bl_levels),             &
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
          2:bl_levels),                                                 &
                                   ! IN  RDZ (K > 1) on U-grid.
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels),                                                 &
                                   ! IN  RDZ (K > 1) on V-grid.
   gamma1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
   gamma2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                                   ! IN new scheme weights.
   ct_ctq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                   ! IN Coefficient in T and q
                                   !       tri-diagonal implicit matrix
   cq_cm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                   ! IN Coefficient in U tri-diagonal
                                   !       implicit matrix
   cq_cm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                   ! IN Coefficient in V tri-diagonal
                                   !       implicit matrix
   dqw_nt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                        ! IN NT incr to qw
   dtl_nt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels)
                                        ! IN NT incr to TL

!  In/outs :-

  REAL, INTENT(INOUT) ::                                                &
   rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         2:bl_levels),                                                  &
                                   ! INOUT Exchange coeffs for moisture.
                                   ! shouldnt change but is scaled by
                                   ! r_sq then back
   qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                   ! INOUT Total water content, but
                                   !       replaced by specific
                                   !       humidity in LS_CLD.
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                   ! INOUT Ice/liquid water temperature,
                                   !       but replaced by T in LS_CLD.
   fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! INOUT Moisture flux between layers
                                   !       (kg per square metre per sec)
                                   !       FQW(,1) is total water flux
                                   !       from surface, 'E'.
   ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! INOUT FTL(,K) contains net
                                   !       turbulent sensible heat flux
                                   !       into layer K from below; so
                                   !       FTL(,1) is the surface
                                   !       sensible heat, H. (W/m2)
   tau_x(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
         bl_levels),                                                    &
                                   ! INOUT W'ly component of surface
                                   !       wind stress (N/sq m).(On
                                   !       UV-grid with first and last
                                   !       rows undefined or at present,
                                   !       set to  missing data
   tau_y(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
         bl_levels),                                                    &
                                   ! INOUT S'ly component of surface
                                   !       wind stress (N/sq m).  On
                                   !       UV-grid; comments as per TAUX
!                                  4 arrays below:
                                   ! INOUT temp arrays for diags
   fqw_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
   ftl_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
   taux_star(udims%i_start:udims%i_end,udims%j_start:udims%j_end,       &
              bl_levels),                                               &
   tauy_star(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,       &
              bl_levels),                                               &
   dqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! INOUT BL increment to q field
   dtl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! INOUT BL increment to T field
   du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,      &
        bl_levels),                                                     &
                                   ! INOUT BL increment to u wind field
   dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,      &
        bl_levels),                                                     &
                                   ! INOUT BL increment to v wind field
   du_star(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end, &
        bl_levels),                                                     &
                                        ! INOUT BL incr to u wind field
   dv_star(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end, &
        bl_levels)
! OUT fields
  REAL, INTENT(OUT) ::                                                  &
   q_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
        ! OUT specific humidity
        ! But at this stage it is qT = qv+qcl+qcf
   t_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
        ! OUT temperature
        ! But at this stage it is tL = t - lcrcp*qcl - lsrcp*qcf
   rhokh_mix(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                        ! OUT Exch coeffs for moisture

!-----------------------------------------------------------------------
!  Local scalars :-

  REAL :: r_sq, rbt, at ,                                               &
   gamma1_uv,                                                           &
   gamma2_uv      ! gamma1 and gamma2 shifted to u or v points

  INTEGER ::                                                            &
   i,j,                                                                 &
                  ! LOCAL Loop counter (horizontal field index).
   k          ! LOCAL Loop counter (vertical level index).

  INTEGER :: jj, j_block  ! omp blocking variables

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('BDY_IMPL4',zhook_in,zhook_handle)

  j_block = 4

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,jj,at,rbt,gamma1_uv,      &
!$OMP& gamma2_uv,r_sq)

!$OMP DO SCHEDULE(STATIC)
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      du(i,j,1) = du(i,j,1) - cq_cm_u(i,j,1)*tau_x(i,j,1)
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      dv(i,j,1) = dv(i,j,1) - cq_cm_v(i,j,1)*tau_y(i,j,1)
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO jj=udims%j_start, udims%j_end, j_block
    DO k = 2, bl_levels

      DO j = jj, MIN(jj+j_block-1,udims%j_end)
        DO i = udims%i_start, udims%i_end
          du(i,j,k) = du(i,j,k) - cq_cm_u(i,j,k)*du(i,j,k-1)
        END DO
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO jj=vdims%j_start, vdims%j_end, j_block
    DO k = 2, bl_levels

      DO j = jj, MIN(jj+j_block-1,vdims%j_end)
        DO i = vdims%i_start, vdims%i_end
          dv(i,j,k) = dv(i,j,k) - cq_cm_v(i,j,k)*dv(i,j,k-1)
        END DO
      END DO

    END DO
  END DO
!$OMP END DO


  IF ( .NOT. l_correct ) THEN
!  1st stage: predictor
!  Keep a copy of computed taux_1.

!$OMP DO SCHEDULE(STATIC)
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        taux_star(i,j,1) = tau_x(i,j,1)
      END DO
    END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        tauy_star(i,j,1) = tau_y(i,j,1)
      END DO
    END DO
!$OMP END DO 

! Save increment from 1st stage. Will be needed for the explicit flux
! at the next stage.

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, bl_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          du_star(i,j,k) = du(i,j,k)
        END DO
      END DO
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          dv_star(i,j,k) = dv(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO 
!---------------------------------------------------------------------

! Complete downward sweep of matrix for increments to TL and QW in the
! boundary layer. It needs to be done here since the scalar fluxes
! (FQW(,,1),FTL(,,1)) computed previously by the surface scheme are not
! available when the main downward sweep subroutine bdy_impl3 is first
! invoked (at the 1st stage of the new solver).

!---------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        r_sq = r_rho_levels(i,j,1)*r_rho_levels(i,j,1)
        rhokh(i,j,2) = r_sq * rhokh(i,j,2)
        fqw(i,j,1)   = r_sq * fqw(i,j,1)
        ftl(i,j,1)   = r_sq * ftl(i,j,1)
        fqw(i,j,2)   = r_sq * fqw(i,j,2)
        ftl(i,j,2)   = r_sq * ftl(i,j,2)
        dqw(i,j,1) = gamma2(i,j) * ( -dtrdz_charney_grid(i,j,1) *       &
            ( fqw(i,j,2) - fqw(i,j,1) ) + dqw_nt(i,j,1) )
        dtl(i,j,1) = gamma2(i,j) * ( -dtrdz_charney_grid(i,j,1) *       &
            ( ftl(i,j,2) - ftl(i,j,1) ) + dtl_nt(i,j,1) )
        at = -dtrdz_charney_grid(i,j,1) *                               &
                 gamma1(i,j)*rhokh(i,j,2)*rdz_charney_grid(i,j,2)
        rbt = 1.0 / ( 1.0 - at*( 1.0 + ct_ctq(i,j,2) ) )
        dqw(i,j,1) = rbt*(dqw(i,j,1) - at*dqw(i,j,2) )
        dtl(i,j,1) = rbt*(dtl(i,j,1) - at*dtl(i,j,2) )
        rhokh(i,j,2) = rhokh(i,j,2)/r_sq
        fqw(i,j,1) = fqw(i,j,1)/r_sq
        ftl(i,j,1) = ftl(i,j,1)/r_sq
        fqw(i,j,2) = fqw(i,j,2)/r_sq
        ftl(i,j,2) = ftl(i,j,2)/r_sq
      END DO
    END DO
!$OMP END DO

  ELSE ! L_CORRECT == TRUE: 2nd stage of the scheme

! Compute 2nd stage correction (total: from tn to tn+1)

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, bl_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          du(i,j,k) = du(i,j,k) + du_star(i,j,k)
        END DO
      END DO
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          dv(i,j,k) = dv(i,j,k) + dv_star(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO

! Compute total surface stress from tn to tn+1 (diagnostic):
! TAUX_star: flux from tn to t*.
! TAU_X: flux from t* to tn+1.

!$OMP DO SCHEDULE(STATIC)
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        tau_x(i,j,1) = tau_x(i,j,1)+taux_star(i,j,1)
      END DO
    END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        tau_y(i,j,1) = tau_y(i,j,1)+tauy_star(i,j,1)
      END DO
    END DO
!$OMP END DO

  END IF

! Update TL, QW and their increments
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      tl(i,j,1) = tl(i,j,1) + dtl(i,j,1)
      qw(i,j,1) = qw(i,j,1) + dqw(i,j,1)
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO jj=tdims%j_start, tdims%j_end, j_block
    DO k = 2, bl_levels
      DO j = jj, MIN(jj+j_block-1,tdims%j_end)
        DO i = tdims%i_start, tdims%i_end
          dtl(i,j,k) = dtl(i,j,k) - ct_ctq(i,j,k)*dtl(i,j,k-1)
          tl(i,j,k) = tl(i,j,k) + dtl(i,j,k)
          dqw(i,j,k) = dqw(i,j,k) - ct_ctq(i,j,k)*dqw(i,j,k-1)
          qw(i,j,k) = qw(i,j,k) + dqw(i,j,k)
        END DO
      END DO
    END DO !bl_levels
  END DO
!$OMP END DO

! Calculate stress and flux diagnostics using the new scheme equations.
! The fluxes are calculated only when requested.

  
  IF ( l_taux ) THEN
    
    IF ( .NOT. l_correct ) THEN
     IF (l_vatpoles) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
            gamma1_uv = gamma1(i+1,j)
            gamma2_uv = gamma2(i+1,j)
            taux_star(i,j,k) = gamma2_uv*tau_x(i,j,k)                 &
                         + gamma1_uv*rhokm_u(i,j,k)*rdz_u(i,j,k)      &
                         * ( du(i,j,k)-du(i,j,k-1) )
          END DO
        END DO
      END DO ! bl_levels
!$OMP END DO
     ELSE
!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
            gamma1_uv = gamma1(i,j)
            gamma2_uv = gamma2(i,j)
            taux_star(i,j,k) = gamma2_uv*tau_x(i,j,k)                 &
                         + gamma1_uv*rhokm_u(i,j,k)*rdz_u(i,j,k)      &
                         * ( du(i,j,k)-du(i,j,k-1) )
          END DO
        END DO
      END DO ! bl_levels
!$OMP END DO
     END IF ! vatpoles
    ELSE
     IF (l_vatpoles) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
! Copy gamma into an array defined on u-points
! In principle, gamma should be interpolated, but sensitivity is expected
! to be small so the adjacent p-point is used to avoid message passing
! This makes vatpoles formulation consistent with the current formulation
            gamma1_uv = gamma1(i+1,j)
            gamma2_uv = gamma2(i+1,j)
            tau_x(i,j,k)= taux_star(i,j,k)+gamma2_uv*tau_x(i,j,k)     &
                        + gamma1_uv*rhokm_u(i,j,k)*rdz_u(i,j,k)       &
                        * ( du(i,j,k)-du(i,j,k-1) )
          END DO
        END DO
      END DO ! bl_levels
!$OMP END DO
     ELSE
!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
        DO j = udims%j_start, udims%j_end
          DO i = udims%i_start, udims%i_end
            gamma1_uv = gamma1(i,j)
            gamma2_uv = gamma2(i,j)
            tau_x(i,j,k)= taux_star(i,j,k)+gamma2_uv*tau_x(i,j,k)     &
                        + gamma1_uv*rhokm_u(i,j,k)*rdz_u(i,j,k)       &
                        * ( du(i,j,k)-du(i,j,k-1) )
          END DO
        END DO
      END DO ! bl_levels
!$OMP END DO
     END IF ! vatpoles
    END IF ! L_correct
  END IF ! l_taux

  IF ( l_tauy ) THEN

    IF ( .NOT. l_correct ) THEN
     IF (l_vatpoles) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
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
            tauy_star(i,j,k) = gamma2_uv*tau_y(i,j,k)                 &
                         + gamma1_uv*rhokm_v(i,j,k)*rdz_v(i,j,k)      &
                         * (dv(i,j,k)-dv(i,j,k-1))
          END DO
        END DO
      END DO ! bl_levels
!$OMP END DO
     ELSE
!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
        DO j = vdims%j_start, vdims%j_end
          DO i = vdims%i_start, vdims%i_end
            gamma1_uv = gamma1(i,j)
            gamma2_uv = gamma2(i,j)
            tauy_star(i,j,k) = gamma2_uv*tau_y(i,j,k)                 &
                         + gamma1_uv*rhokm_v(i,j,k)*rdz_v(i,j,k)      &
                         * (dv(i,j,k)-dv(i,j,k-1))
          END DO
        END DO
      END DO ! bl_levels
!$OMP END DO
     END IF ! vatpoles
    ELSE
     IF (l_vatpoles) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
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
            tau_y(i,j,k)= tauy_star(i,j,k)+gamma2_uv*tau_y(i,j,k)     &
                         + gamma1_uv*rhokm_v(i,j,k)*rdz_v(i,j,k)      &
                         * (dv(i,j,k)-dv(i,j,k-1))
          END DO
        END DO
      END DO ! bl_levels
!$OMP END DO
     ELSE
!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
        DO j = vdims%j_start, vdims%j_end
          DO i = vdims%i_start, vdims%i_end
            gamma1_uv = gamma1(i,j)
            gamma2_uv = gamma2(i,j)
            tau_y(i,j,k)= tauy_star(i,j,k)+gamma2_uv*tau_y(i,j,k)     &
                         + gamma1_uv*rhokm_v(i,j,k)*rdz_v(i,j,k)      &
                         * (dv(i,j,k)-dv(i,j,k-1))
          END DO
        END DO
      END DO ! bl_levels
!$OMP END DO
     END IF ! vatpoles
    END IF ! L_correct
  END IF ! l_tauy

  IF ( l_ftl ) THEN

    IF ( .NOT. l_correct ) THEN

!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            ftl_star(i,j,k) = gamma2(i,j)*ftl(i,j,k)                    &
                 - gamma1(i,j)*rhokh(i,j,k)*rdz_charney_grid(i,j,k)     &
                              * (dtl(i,j,k)-dtl(i,j,k-1))
          END DO
        END DO
      END DO
!$OMP END DO
    ELSE

!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            ftl(i,j,k) = ftl_star(i,j,k)+gamma2(i,j)*ftl(i,j,k)         &
                 - gamma1(i,j)*rhokh(i,j,k)*rdz_charney_grid(i,j,k)     &
                              * (dtl(i,j,k)-dtl(i,j,k-1))
          END DO
        END DO
      END DO
!$OMP END DO

    END IF ! L_correct
  END IF ! L_ftl

  IF ( l_fqw ) THEN

    IF ( .NOT. l_correct ) THEN

!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            fqw_star(i,j,k) = gamma2(i,j)*fqw(i,j,k)                    &
                 - gamma1(i,j)*rhokh(i,j,k)*rdz_charney_grid(i,j,k)     &
                 * (dqw(i,j,k)-dqw(i,j,k-1))
          END DO
        END DO
      END DO ! bl_levels
!$OMP END DO

    ELSE

!$OMP DO SCHEDULE(STATIC)
      DO k = 2, bl_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            fqw(i,j,k) = fqw_star(i,j,k)+gamma2(i,j)*fqw(i,j,k)         &
                 - gamma1(i,j)*rhokh(i,j,k)*rdz_charney_grid(i,j,k)     &
                   * (dqw(i,j,k)-dqw(i,j,k-1))
          END DO
        END DO
      END DO ! bl_levels
!$OMP END DO

    END IF
  END IF

  IF ( l_correct ) THEN

!-----------------------------------------------------------------------
!!     Convert FTL to sensible heat flux in Watts per square metre.
!      Also, IMPL_CAL only updates FTL_TILE(*,1) and FQW_TILE(*,1)
!      over sea points, so copy this to remaining tiles
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
    DO k = 2, bl_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end 
          ftl(i,j,k) = ftl(i,j,k)*cp
        END DO
      END DO
    END DO
!$OMP END DO

!Copy T and Q from workspace to INOUT space.

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, bl_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          t_latest(i,j,k)=tl(i,j,k)
          q_latest(i,j,k)=qw(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO

!-----------------------------------------------------------------------
!    Set RHOKH, the coefficients required for tracer mixing.
!    Required 5B and after due to change in contents of RHOKH in rest
!    of routine.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
    DO k = 2, bl_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          rhokh_mix(i,j,k) = rhokh(i,j,k)*                              &
                         rdz_charney_grid(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO

  END IF ! L_CORRECT

!$OMP END PARALLEL

  IF (lhook) CALL dr_hook('BDY_IMPL4',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE bdy_impl4
