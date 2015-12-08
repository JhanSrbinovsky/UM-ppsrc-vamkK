! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE IM_BL_PT1 ----------------------------------------------
!!!
!!!  Purpose: Calculate downward sweep of matrix for increments to
!!!           T and Q in the boundary layer, using an
!!!           implicit numerical scheme.  The tridiagonal matrices are
!!!           inverted using simple Gaussian elimination.
!!!
!!!
!!!  Programming standard: UM Documentation Paper No 3
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------
!!  Arguments :-
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
SUBROUTINE im_bl_pt1 (                                                  &
 offx ,offy ,row_length,rows,n_rows,bl_levels,                          &
 global_row_length,n_proc, n_procy, proc_row_group,at_extremity,        &
 halo_i, halo_j,                                                        &
 dtrdz_charney_grid,dtrdz_u,dtrdz_v,                                    &
 rhokh,rhokm_u,rhokm_v,                                                 &
 rdz_charney_grid,rdz_u,rdz_v,GAMMA,                                    &
 dqw_nt,dtl_nt,du_nt,dv_nt,                                             &
 fqw,ftl,tau_x,tau_y,                                                   &
 ct_ctq,dqw,dtl,cq_cm_u,cq_cm_v,du,dv,                                  &
 ltimer                                                                 &
)

  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, udims_s, vdims_s, udims_l, vdims_l, tdims_l

  USE level_heights_mod, ONLY:                                          &
    r_theta_levels, r_rho_levels

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE p_to_u_mod, ONLY: p_to_u
  USE p_to_v_mod, ONLY: p_to_v
  IMPLICIT NONE

  LOGICAL :: ltimer

  INTEGER ::                                                            &
    row_length,                                                         &
                                 ! IN Number of points on a row
    rows,                                                               &
                                 ! IN Number of rows in a theta field
    n_rows,                                                             &
                                 ! IN Number of rows in a v field
    offx ,                                                              &
                                 ! IN Size of small halo in i
    offy ,                                                              &
                                 ! IN Size of small halo in j.
    halo_i,                                                             &
                                 ! IN Size of halo in i direction
    halo_j,                                                             &
                                 ! IN Size of halo in j direction
    n_proc,                                                             &
                                ! Total number of processors
    n_procy,                                                            &
                                ! Number of procs N/S
    global_row_length,                                                  &
                                ! number of points per global row
    proc_row_group,                                                     &
                                ! Group id for processors on the same row                             
                                 
    bl_levels                   ! IN No. of atmospheric levels for
                                 !    which boundary layer fluxes are
                                 !    calculated.
  LOGICAL ::                                                            &
    at_extremity(4)  ! Indicates if this processor is at north,
                         !  south,east or west of the processor grid 
 

  REAL ::                                                               &
   dtrdz_charney_grid(row_length,rows,bl_levels),                       &
                                   ! IN dz for bottom BL_LEVELS
   dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                   ! IN -g.dt/dp for model wind layers
   dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                   ! IN -g.dt/dp for model wind layers
   rhokh(row_length,rows,2:bl_levels),                                  &
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
   rdz_charney_grid(row_length,rows,bl_levels),                         &
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
   GAMMA(bl_levels),                                                    &
                                   ! IN Implicit weighting.
   dqw_nt(row_length,rows,bl_levels),                                   &
                                   ! IN Non-turbulent increment for QW.
   dtl_nt(row_length,rows,bl_levels),                                   &
                                   ! IN Non-turbulent increment for TL.
   du_nt(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,   &
          bl_levels),                                                   &
                                   ! IN u non-turbulent increments.
   dv_nt(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,   &
          bl_levels),                                                   &
                                   ! IN v non-turbulent increments.
   fqw(row_length,rows,bl_levels),                                      &
                                   ! IN Flux of QW (ie., for surface,
                                   !    total evaporation). Kg/sq m/s
   ftl(row_length,rows,bl_levels),                                      &
                                   ! IN Flux of TL (ie., for surface,
                                   !    H/Cp where H is sensible heat
                                   !    in W per sq m).
   tau_x(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
         bl_levels),                                                    &
                                   ! IN x-component of turbulent
                                   !    stress at levels k-1/2;
                                   !    eg. TAUX(,1) is surface stress.
                                   !    U-grid, 1st and last rows set
                                   !    to "missing data". (N/sq m)
                                   !    IN as "explicit" fluxes from
                                   !    ex_flux_uv, OUT as "implicit
   tau_y(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
         bl_levels)
                                   ! IN y-component of turbulent
                                   !    stress at levels k-1/2;
                                   !    eg. TAUX(,1) is surface stress.
                                   !    V-grid, 1st and last rows set
                                   !    to "missing data". (N/sq m)
                                   !    IN as "explicit" fluxes from
                                   !    ex_flux_uv, OUT as "implicit


  REAL ::                                                               &
   ct_ctq(row_length,rows,bl_levels),                                   &
                                   ! OUT Coefficient in T and q
                                   !     tri-diagonal implicit matrix
   cq_cm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                   ! OUT Coefficient in U and V
                                   !     tri-diagonal implicit matrix
   cq_cm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                   ! OUT Coefficient in U and V
                                   !     tri-diagonal implicit matrix
   dqw(row_length,rows,bl_levels),                                      &
                                   ! OUT BL increment to q field
   dtl(row_length,rows,bl_levels),                                      &
                                   ! OUT BL increment to T field
   du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,      &
        bl_levels),                                                     &
                                   ! INOUT BL increment to u wind field
   dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,      &
        bl_levels) 
                                   ! INOUT BL increment to v wind field

  REAL ::                                                               &
   r_theta_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,       &
               0:bl_levels),                                            &
                                               ! Vertical grids for U
   r_theta_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,       &
               0:bl_levels)                    ! and V flux levels

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
   r_sq,                                                                &
                ! square of height variables
   rr_sq    ! 1/square of height variables

  INTEGER ::                                                            &
   blm1,                                                                &
                ! BL_LEVELS minus 1.
   i,j,                                                                 &
                ! Loop counter (horizontal field index).
   k        ! Loop counter (vertical index).

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle






!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('IM_BL_PT1',zhook_in,zhook_handle)

  blm1 = bl_levels-1

!-----------------------------------------------------------------------
!!  1.0 Interpolate r_theta_levels to U,V columns
!-----------------------------------------------------------------------
      CALL p_to_u(r_theta_levels,                                 &
                  tdims_l%i_start,tdims_l%i_end,                  &
                  tdims_l%j_start,tdims_l%j_end,                  &
                  udims%i_start,udims%i_end,                      &
                  udims%j_start,udims%j_end,                      &
                  0,bl_levels,r_theta_u) 

      CALL p_to_v(r_theta_levels,                                 &
                  tdims_l%i_start,tdims_l%i_end,                  &
                  tdims_l%j_start,tdims_l%j_end,                  &
                  vdims%i_start,vdims%i_end,                      &
                  vdims%j_start,vdims%j_end,                      &
                  0,bl_levels,r_theta_v) 

!-----------------------------------------------------------------------
!! 2.0 For simulations on a sphere we must use spherical geometry for
!!     vertical flux-divergences.   Thus, leaving out rho for
!!     simplicity, the standard cartesian flux-divergence:
!!          dQ(K)/dt = -(FQ(K+1)-FQ(K))/DZ
!!     becomes:
!!          dQ(K)/dt = -(r_flux(K+1)^2*FQ(K+1)-r_flux(K)^2*FQ(K))
!!                      / (r_full(K)^2 * DZ)
!!     In the code below it would be clearer to include the r^2
!!     explicitly where they occur in the equations but it is
!!     computationally cheaper to multiply the fluxes by r^2 in one go
!!     at the start and then divide out the r^2 at the end.
!-----------------------------------------------------------------------
  DO k = 2, bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        r_sq = r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
        rhokh(i,j,k) = r_sq * rhokh(i,j,k)
        fqw(i,j,k)   = r_sq * fqw(i,j,k)
        ftl(i,j,k)   = r_sq * ftl(i,j,k)
      END DO
    END DO
  END DO
  DO k = 2, bl_levels
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        r_sq           = r_theta_u(i,j,k-1)*r_theta_u(i,j,k-1)
        rhokm_u(i,j,k) = r_sq * rhokm_u(i,j,k)
        tau_x(i,j,k)   = r_sq * tau_x(i,j,k)
      END DO
    END DO
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        r_sq           = r_theta_v(i,j,k-1)*r_theta_v(i,j,k-1)
        rhokm_v(i,j,k) = r_sq * rhokm_v(i,j,k)
        tau_y(i,j,k)   = r_sq * tau_y(i,j,k)
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
!!  3.0 Calculate matrix elements
!-----------------------------------------------------------------------

  DO j = 1, rows
    DO i = 1, row_length
! Include non-turbulent increments.
      dqw(i,j,bl_levels) = dtrdz_charney_grid(i,j,bl_levels) *          &
                           fqw(i,j,bl_levels) + dqw_nt(i,j,bl_levels)
      dtl(i,j,bl_levels) = dtrdz_charney_grid(i,j,bl_levels) *          &
                           ftl(i,j,bl_levels) + dtl_nt(i,j,bl_levels)

      ct_ctq(i,j,bl_levels) = -dtrdz_charney_grid(i,j,bl_levels) *      &
             GAMMA(bl_levels)*rhokh(i,j,bl_levels)*                     &
              rdz_charney_grid(i,j,bl_levels)

      rbt = 1.0 / ( 1.0 - ct_ctq(i,j,bl_levels) )

      dqw(i,j,bl_levels) = rbt * dqw(i,j,bl_levels)
      dtl(i,j,bl_levels) = rbt * dtl(i,j,bl_levels)

      ct_ctq(i,j,bl_levels) = rbt * ct_ctq(i,j,bl_levels)    ! P244.1
    END DO
  END DO


  DO k = blm1, 2, -1
    DO j = 1, rows
      DO i = 1, row_length

        dqw(i,j,k) = -dtrdz_charney_grid(i,j,k) *                       &
                    ( fqw(i,j,k+1) - fqw(i,j,k) ) + dqw_nt(i,j,k)
        dtl(i,j,k) = -dtrdz_charney_grid(i,j,k) *                       &
                    ( ftl(i,j,k+1) - ftl(i,j,k) ) + dtl_nt(i,j,k)

        at = -dtrdz_charney_grid(i,j,k) *                               &
               GAMMA(k+1)*rhokh(i,j,k+1)*rdz_charney_grid(i,j,k+1)

        ct_ctq(i,j,k) = -dtrdz_charney_grid(i,j,k) *                    &
                     GAMMA(k)*rhokh(i,j,k)*rdz_charney_grid(i,j,k)

        rbt = 1.0 / ( 1.0 - ct_ctq(i,j,k) -                             &
                               at*( 1.0 + ct_ctq(i,j,k+1) ) )

        dqw(i,j,k) = rbt * (dqw(i,j,k) - at*dqw(i,j,k+1) )
        dtl(i,j,k) = rbt * (dtl(i,j,k) - at*dtl(i,j,k+1) )

        ct_ctq(i,j,k) = rbt * ct_ctq(i,j,k)               ! P244.1
      END DO
    END DO
  END DO !blm1,2,-1

!-----------------------------------------------------------------------
!! 3.3 Bottom model layer QW row of matrix equation.
!-----------------------------------------------------------------------

  DO j = 1, rows
    DO i = 1, row_length

      dqw(i,j,1) = -dtrdz_charney_grid(i,j,1) * fqw(i,j,2) +            &
                    dqw_nt(i,j,1)
      dtl(i,j,1) = -dtrdz_charney_grid(i,j,1) * ftl(i,j,2) +            &
                    dtl_nt(i,j,1)

      at = -dtrdz_charney_grid(i,j,1) *                                 &
                 GAMMA(2)*rhokh(i,j,2)*rdz_charney_grid(i,j,2)

      rbt = 1.0 / ( 1.0 - at*( 1.0 + ct_ctq(i,j,2) ) )

      dqw(i,j,1) = rbt * (dqw(i,j,1) - at*dqw(i,j,2) )
      dtl(i,j,1) = rbt * (dtl(i,j,1) - at*dtl(i,j,2) )

! Now set CT_CTQ(1) to be r^2 * BETA
      r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
      ct_ctq(i,j,1) = - r_sq * dtrdz_charney_grid(i,j,1) * rbt

    END DO
  END DO




  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end

      du(i,j,bl_levels) = -dtrdz_u(i,j,bl_levels) *                     &
                         tau_x(i,j,bl_levels)

! addition of non-turbulent increments
      du(i,j,bl_levels) = du(i,j,bl_levels)                             &
                           + du_nt(i,j,bl_levels)

      cq_cm_u(i,j,bl_levels) = -dtrdz_u(i,j,bl_levels) *                &
            GAMMA(bl_levels)*                                           &
            rhokm_u(i,j,bl_levels)*rdz_u(i,j,bl_levels)

      rbm = 1.0 / ( 1.0 - cq_cm_u(i,j,bl_levels) )

      du(i,j,bl_levels) = rbm * du(i,j,bl_levels)

      cq_cm_u(i,j,bl_levels) = rbm * cq_cm_u(i,j,bl_levels)
    END DO
  END DO


  DO k = blm1, 2, -1

    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end

        du(i,j,k) = dtrdz_u(i,j,k) *                                    &
                       ( tau_x(i,j,k+1) - tau_x(i,j,k) )
! addition of non-turbulent increments
        du(i,j,k) = du(i,j,k) + du_nt(i,j,k)

        am = -dtrdz_u(i,j,k) * GAMMA(k+1)*rhokm_u(i,j,k+1)*             &
                        rdz_u(i,j,k+1)

        cq_cm_u(i,j,k) = -dtrdz_u(i,j,k) * GAMMA(k)*rhokm_u(i,j,k)*     &
              rdz_u(i,j,k)

        rbm = 1.0 / ( 1.0 - cq_cm_u(i,j,k) -                            &
                        am*( 1.0 + cq_cm_u(i,j,k+1) ) )

        du(i,j,k) = rbm * ( du(i,j,k) - am*du(i,j,k+1) )

        cq_cm_u(i,j,k) = rbm * cq_cm_u(i,j,k)
      END DO
    END DO
  END DO ! loop over 2,BLM1


  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end

      du(i,j,1) = dtrdz_u(i,j,1) * tau_x(i,j,2)

! addition of non-turbulent increments
      du(i,j,1) = du(i,j,1) + du_nt(i,j,1)

      am = -dtrdz_u(i,j,1) * GAMMA(2)*rhokm_u(i,j,2)                    &
                   *rdz_u(i,j,2)

      rbm = 1.0 / ( 1.0 - am *( 1.0 + cq_cm_u(i,j,2) ) )

      du(i,j,1) = rbm * ( du(i,j,1) - am*du(i,j,2) )

! Now set CQ_CM_U(1) to be r^2 * BETA
      r_sq = r_theta_u(i,j,0)*r_theta_u(i,j,0)
      cq_cm_u(i,j,1) = r_sq * dtrdz_u(i,j,1) * rbm
    END DO
  END DO




  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end

      dv(i,j,bl_levels) = -dtrdz_v(i,j,bl_levels) *                     &
                         tau_y(i,j,bl_levels)

! addition of non-turbulent increments
      dv(i,j,bl_levels) = dv(i,j,bl_levels)                             &
                           + dv_nt(i,j,bl_levels)

      cq_cm_v(i,j,bl_levels) = -dtrdz_v(i,j,bl_levels) *                &
            GAMMA(bl_levels)*                                           &
            rhokm_v(i,j,bl_levels)*rdz_v(i,j,bl_levels)

      rbm = 1.0 / ( 1.0 - cq_cm_v(i,j,bl_levels) )

      dv(i,j,bl_levels) = rbm * dv(i,j,bl_levels)

      cq_cm_v(i,j,bl_levels) = rbm * cq_cm_v(i,j,bl_levels)
    END DO
  END DO


  DO k = blm1, 2, -1

    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end

        dv(i,j,k) = dtrdz_v(i,j,k) *                                    &
                       ( tau_y(i,j,k+1) - tau_y(i,j,k) )
! addition of non-turbulent increments
        dv(i,j,k) = dv(i,j,k) + dv_nt(i,j,k)

        am = -dtrdz_v(i,j,k) * GAMMA(k+1)*rhokm_v(i,j,k+1)*             &
                        rdz_v(i,j,k+1)

        cq_cm_v(i,j,k) = -dtrdz_v(i,j,k) * GAMMA(k)*rhokm_v(i,j,k)*     &
              rdz_v(i,j,k)

        rbm = 1.0 / ( 1.0 - cq_cm_v(i,j,k) -                            &
                        am*( 1.0 + cq_cm_v(i,j,k+1) ) )

        dv(i,j,k) = rbm * ( dv(i,j,k) - am*dv(i,j,k+1) )

        cq_cm_v(i,j,k) = rbm * cq_cm_v(i,j,k)
      END DO
    END DO
  END DO ! loop over 2,BLM1


  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end

      dv(i,j,1) = dtrdz_v(i,j,1) * tau_y(i,j,2)

! addition of non-turbulent increments
      dv(i,j,1) = dv(i,j,1) + dv_nt(i,j,1)

      am = -dtrdz_v(i,j,1) * GAMMA(2)*rhokm_v(i,j,2)                    &
                   *rdz_v(i,j,2)

      rbm = 1.0 / ( 1.0 - am *( 1.0 + cq_cm_v(i,j,2) ) )

      dv(i,j,1) = rbm * ( dv(i,j,1) - am*dv(i,j,2) )

! Now set CQ_CM_V(1) to be r^2 * BETA
      r_sq = r_theta_v(i,j,0)*r_theta_v(i,j,0)
      cq_cm_v(i,j,1) = r_sq * dtrdz_v(i,j,1) * rbm
    END DO
  END DO

!-----------------------------------------------------------------------
!! 4.0 Return fluxes back to their true value by dividing by r*r
!-----------------------------------------------------------------------
  DO k = 2, bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        rr_sq = 1.0 / ( r_rho_levels(i,j,k)*r_rho_levels(i,j,k) )
        rhokh(i,j,k) = rhokh(i,j,k) * rr_sq
        fqw(i,j,k)   = fqw(i,j,k) * rr_sq
        ftl(i,j,k)   = ftl(i,j,k) * rr_sq
      END DO
    END DO
  END DO
  DO k = 2, bl_levels
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        rr_sq          = 1.0 / (r_theta_u(i,j,k-1)*r_theta_u(i,j,k-1))
        rhokm_u(i,j,k) = rr_sq * rhokm_u(i,j,k)
        tau_x(i,j,k)   = rr_sq * tau_x(i,j,k)
      END DO
    END DO
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        rr_sq          = 1.0 / (r_theta_v(i,j,k-1)*r_theta_v(i,j,k-1))
        rhokm_v(i,j,k) = rr_sq * rhokm_v(i,j,k)
        tau_y(i,j,k)   = rr_sq * tau_y(i,j,k)
      END DO
    END DO
  END DO

  IF (lhook) CALL dr_hook('IM_BL_PT1',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE im_bl_pt1
