! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE IM_BL_PT2 ----------------------------------------------
!!!
!!!  Purpose: Calculate increments for
!!!           T and Q in the boundary layer, using an
!!!           implicit numerical scheme.  The tridiagonal matrices are
!!!           inverted using simple Gaussian elimination.
!!!
!!!
!!!  Programming standard: UM Documentation Paper No 4, Version 2,
!!!                        dated 18/1/90
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
SUBROUTINE im_bl_pt2 (                                                  &
 offx ,offy ,row_length,rows,n_rows,bl_levels,                          &
 l_ftl,l_fqw,l_taux,l_tauy,                                             &
 rhokh,rhokm_u,rhokm_v,                                                 &
 rdz_charney_grid,rdz_u,rdz_v,GAMMA,                                    &
 ct_ctq,dqw,dtl,cq_cm_u,cq_cm_v,du,dv,                                  &
 fqw,ftl,tau_x,tau_y,qw,tl,                                             &
 ltimer                                                                 &
)

  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, udims_s, vdims_s, udims_l, vdims_l
  USE atmos_constants_mod, ONLY: cp

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

  LOGICAL :: ltimer

  INTEGER ::                                                            &
    row_length,                                                         &
                                 ! Local number of points on a row
    rows,                                                               &
                                 ! Local number of rows in a theta field
    n_rows,                                                             &
                                 ! Local number of rows in a v field
    offx ,                                                              &
                                 ! Size of small halo in i
    offy ,                                                              &
                                 ! Size of small halo in j.
   bl_levels                 ! IN No. of atmospheric levels for
!                                     which boundary layer fluxes are
!                                     calculated.

  LOGICAL ::                                                            &
   l_ftl,                                                               &
   l_fqw,                                                               &
   l_taux,                                                              &
   l_tauy


  REAL ::                                                               &
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
   GAMMA(bl_levels)          ! IN Implicit weighting.


  REAL ::                                                               &
   ct_ctq(row_length,rows,bl_levels),                                   &
                                 ! INOUT Coefficient in T and q
                                 !       tri-diagonal implicit matrix
   cq_cm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                 ! INOUT Coefficient in U tri-diagonal
                                 !       implicit matrix
   cq_cm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                 ! INOUT Coefficient in V tri-diagonal
                                 !       implicit matrix
   dqw(row_length,rows,bl_levels),                                      &
                                 ! INOUT BL increment to q field
   dtl(row_length,rows,bl_levels),                                      &
                                 ! INOUT BL increment to T field
   du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,      &
        bl_levels),                                                     &
                                 ! INOUT BL increment to u wind field
   dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,      &
        bl_levels),                                                     &
                                 ! INOUT BL increment to v wind field
   fqw(row_length,rows,bl_levels),                                      &
                                 ! INOUT Flux of QW (ie., for surface,
                                 !       total evaporation). Kg/sq m/s
   ftl(row_length,rows,bl_levels),                                      &
                                 ! INOUT Flux of TL (ie., for surface,
                                 !       H/Cp where H is sensible heat
                                 !       in W per sq m).
   tau_x(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
         bl_levels),                                                    &
                                 ! INOUT x-component of turbulent
                                 !       stress at levels k-1/2;
                                 !       eg. TAUX(,1) is surface stress
                                 !       on UV-grid, 1st and last rows
                                 !       set to "missing data". (N/sq m)
                                 !       IN as "explicit" fluxes from
                                 !       ex_flux_uv, OUT as "implicit
   tau_y(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
         bl_levels),                                                    &
                                 ! INOUT y-component of turbulent
                                 !       stress at levels k-1/2;
                                 !       eg. TAUX(,1) is surface stress
                                 !       on UV-grid, 1st and last rows
                                 !       set to "missing data". (N/sq m)
                                 !       IN as "explicit" fluxes from
                                 !       ex_flux_uv, OUT as "implicit
   qw(row_length,rows,bl_levels),                                       &
                                 ! INOUT Total water content (kg per
                                 !       kg air).  From P243.
   tl(row_length,rows,bl_levels)
                                 ! INOUT Liquid/frozen water
                                 !       temperature (K).  From P243.


!  Local scalars :-
  INTEGER ::                                                            &
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

  IF (lhook) CALL dr_hook('IM_BL_PT2',zhook_in,zhook_handle)

  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      du(i,j,1) = du(i,j,1) - cq_cm_u(i,j,1)*tau_x(i,j,1)
    END DO
  END DO

  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      dv(i,j,1) = dv(i,j,1) - cq_cm_v(i,j,1)*tau_y(i,j,1)
    END DO
  END DO


  DO k = 2, bl_levels

    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        du(i,j,k) = du(i,j,k) - cq_cm_u(i,j,k)*du(i,j,k-1)
      END DO
    END DO

    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        dv(i,j,k) = dv(i,j,k) - cq_cm_v(i,j,k)*dv(i,j,k-1)
      END DO
    END DO

  END DO


  DO j = 1, rows
    DO i = 1, row_length
      dtl(i,j,1) = dtl(i,j,1) - ct_ctq(i,j,1)*ftl(i,j,1)/cp
      tl(i,j,1) = tl(i,j,1) + dtl(i,j,1)
      dqw(i,j,1) = dqw(i,j,1) - ct_ctq(i,j,1)*fqw(i,j,1)
      qw(i,j,1) = qw(i,j,1) + dqw(i,j,1)
    END DO
  END DO


  DO k = 2, bl_levels
    DO j = 1, rows
      DO i = 1, row_length

        dtl(i,j,k) = dtl(i,j,k) - ct_ctq(i,j,k)*dtl(i,j,k-1)
        tl(i,j,k) = tl(i,j,k) + dtl(i,j,k)
        dqw(i,j,k) = dqw(i,j,k) - ct_ctq(i,j,k)*dqw(i,j,k-1)
        qw(i,j,k) = qw(i,j,k) + dqw(i,j,k)

      END DO
    END DO
  END DO !bl_levels


! MD
! Calculate stress and flux diagnostics. The fluxes are
! calculated only when requested.

  IF ( l_taux ) THEN
    DO k = 2, bl_levels
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          tau_x(i,j,k) = tau_x(i,j,k)                                   &
                       + GAMMA(k)*rhokm_u(i,j,k)*rdz_u(i,j,k)           &
                       * ( du(i,j,k)-du(i,j,k-1) )
        END DO
      END DO
    END DO ! bl_levels
  END IF

  IF ( l_tauy ) THEN
    DO k = 2, bl_levels
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          tau_y(i,j,k) = tau_y(i,j,k)                                   &
                       + GAMMA(k)*rhokm_v(i,j,k)*rdz_v(i,j,k)           &
                       * (dv(i,j,k)-dv(i,j,k-1))
        END DO
      END DO
    END DO ! bl_levels
  END IF

  IF ( l_ftl ) THEN
    DO k = 2, bl_levels
      DO j = 1, rows
        DO i = 1, row_length
          ftl(i,j,k) = ftl(i,j,k) - GAMMA(k)*rhokh(i,j,k)               &
                     * rdz_charney_grid(i,j,k)                          &
                     * (dtl(i,j,k)-dtl(i,j,k-1))
        END DO
      END DO
    END DO
  END IF

  IF ( l_fqw ) THEN
    DO k = 2, bl_levels
      DO j = 1, rows
        DO i = 1, row_length
          fqw(i,j,k) = fqw(i,j,k) - GAMMA(k)*rhokh(i,j,k)               &
                     * rdz_charney_grid(i,j,k)                          &
                     * (dqw(i,j,k)-dqw(i,j,k-1))
        END DO
      END DO
    END DO ! bl_levels
  END IF

  IF (lhook) CALL dr_hook('IM_BL_PT2',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE im_bl_pt2
