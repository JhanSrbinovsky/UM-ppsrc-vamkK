! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!     SUBROUTINE IMP_MIX -----------------------------------------------
!
!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: Boundary Layer
!
!    Purpose: Calculate turbulent mixing increments for a passive tracer
!             using an implicit numerical scheme.  The tridiagonal
!             matices are inverted using simple Gaussian elimination.
!
!    Programming standard: UM Documentation Paper No. 3
!
!    Documentation: UM Documentation Paper No 24.
!
!----------------------------------------------------------------------
MODULE imp_mix_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE imp_mix (                                                    &
 bl_levels,dtrdz,                                                       &
 gamma_rhokh_rdz, gamma_rhok_dep,f_field,surf_dep_flux,field            &
 )

  USE level_heights_mod, ONLY:                                          &
    r_theta_levels, r_rho_levels
  USE atm_fields_bounds_mod, ONLY:  pdims, tdims
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  !$ USE omp_lib
  IMPLICIT NONE

!  Inputs :-

  INTEGER, INTENT(IN) ::                                                &
   bl_levels                   ! IN No. of atmospheric levels for
                                   !    which boundary layer fluxes are
                                   !    calculated.
  REAL, INTENT(IN) ::                                                   &
   dtrdz(pdims%i_start:,pdims%j_start:,:)
                                   ! IN  dt/(rho*r*r*dz) for scalar
                                   !     flux divergence

!  Next 2 arrays are IN as "explicit" fluxes and OUT as "implicit"
!  fluxes.

  REAL, INTENT(INOUT) ::                                                &
   gamma_rhokh_rdz(pdims%i_start:,pdims%j_start:,2:),                   &
                                   ! INOUT Turbulent mixing coefs. above
                                   !    surface, =GAMMA(K)*RHOKH(,K)
                                   !    *RDZ(K) for K>=2 (from KMKH).
   gamma_rhok_dep(pdims%i_start:,pdims%j_start:),                       &
                                   ! INOUT Surface exchange coefficient
                                   !    for surface deposition*GAMMA(1)
   f_field(pdims%i_start:,pdims%j_start:,:),                            &
                                           ! INOUT Flux of tracer
   surf_dep_flux(pdims%i_start:,pdims%j_start:),                        &
                                           ! INOUT surface deposition
                                           !       flux
   field(pdims%i_start:,pdims%j_start:,:)
                                           ! INOUT Amount of tracer

!    Local and other symbolic constants :-
!   Workspace :-
!   4*BL_LEVELS + 4 blocks of real workspace are required.
  REAL ::                                                               &
   af(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),   &
                                           ! Elements in rows in matrix
                                   ! equation (modified during
                                   ! Gaussian elimination calculations).
   d_field(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels)
                                           ! Delta FIELD (tracer field)
                                   ! elements of vector on RHS, then
                                   ! LHS, of eqn P244.79.
!
!  Local scalars :-
  REAL ::                                                               &
   cf,                                                                  &
                ! Matrix element for local increments
   rbf,                                                                 &
                ! Reciprocal of B for local increments
   r_sq,                                                                &
                ! r*r
   rr_sq    ! 1/(r*r)

  INTEGER ::                                                            &
   blm1,                                                                &
                ! BL_LEVELS minus 1.
   i,                                                                   &
                ! Loop counter (horizontal field index).
   j,                                                                   &
                ! Offset version of I.
   k,                                                                   &
                ! Loop counter (vertical index).
   km1,                                                                 &
                ! K minus 1.
   kp1,                                                                 &
                ! K plus 1.
   jj,                                                                  &
                ! omp blocking counter                               
   omp_block                                                            
                ! omp block size 

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('IMP_MIX',zhook_in,zhook_handle)

  blm1 = bl_levels-1
  
  omp_block =  pdims%j_end -  pdims%j_start + 1

!$OMP PARALLEL DEFAULT(NONE) SHARED(f_field, dtrdz, field, pdims,       &
!$OMP& gamma_rhokh_rdz, r_rho_levels, surf_dep_flux, d_field, af,       &
!$OMP& r_theta_levels, gamma_rhok_dep, bl_levels, blm1)      &
!$OMP& PRIVATE(r_sq, cf, rbf, kp1, km1, jj, k, j, i, rr_sq, omp_block) 

!$ omp_block = (pdims%j_end -  pdims%j_start + 1)/ omp_get_num_threads()

! ----------------------------------------------------------------------
!  (A) Calculations on P-grid.
!-----------------------------------------------------------------------
! 1.0 For simulations on a sphere use spherical geometry for vertical
!     derivatives so multiply fluxes and rho*K_h by r*r
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        r_sq = r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
        gamma_rhokh_rdz(i,j,k) = r_sq * gamma_rhokh_rdz(i,j,k)
        f_field(i,j,k)         = r_sq * f_field(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO 

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
      gamma_rhok_dep(i,j) = r_sq * gamma_rhok_dep(i,j)
      f_field(i,j,1)      = r_sq * f_field(i,j,1)
      surf_dep_flux(i,j)  = r_sq * surf_dep_flux(i,j)
    END DO
  END DO
!$OMP END DO 

! ----------------------------------------------------------------------
!  4.  Calculate those matrix and vector elements on the LHS of eqn
!      which are to do with implicit solution of the tracer
!      transport problem at the surface and between all levels.
!      Begin with "upward sweep" through lower half of matrix).
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
!-----------------------------------------------------------------------
!  4.2 Lowest atmospheric layer FIELD row of matrix.
!-----------------------------------------------------------------------
          ! "Explicit" increment to FIELD(1)
      d_field(i,j,1) = -dtrdz(i,j,1) *                                  &
                     ( f_field(i,j,2) - f_field(i,j,1)                  &
                                    - surf_dep_flux(i,j) )

      cf = -dtrdz(i,j,1) * gamma_rhok_dep(i,j)
      af(i,j,1) = -dtrdz(i,j,1) * gamma_rhokh_rdz(i,j,2)
      rbf = 1.0 / ( 1.0 - af(i,j,1) - cf )
      d_field(i,j,1) = rbf * d_field(i,j,1)
      af(i,j,1) = rbf * af(i,j,1)
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!  4.3 Rows of matrix applying to FIELD transport into model layers in
!      the "middle" of the "boundary" layer, i.e. all but the bottom
!      layer and the top "boundary" layer.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO jj=pdims%j_start, pdims%j_end, omp_block
    DO k = 2, blm1
      kp1 = k+1
      km1 = k-1
      DO j = jj, MIN(jj+omp_block-1, pdims%j_end)
        DO i = pdims%i_start, pdims%i_end

!   "Explicit" flux divergence across layer giving explicit FIELD
!   increment due to mixing

          d_field(i,j,k) = -dtrdz(i,j,k) *                              &
                         (f_field(i,j,kp1) - f_field(i,j,k))
          af(i,j,k) = -dtrdz(i,j,k) * gamma_rhokh_rdz(i,j,kp1)
          cf = -dtrdz(i,j,k) * gamma_rhokh_rdz(i,j,k)
          rbf = 1.0 / ( 1.0 - af(i,j,k)                                 &
                            - cf * ( 1.0 + af(i,j,km1) ) )
          d_field(i,j,k) = rbf * ( d_field(i,j,k)                       &
                                    - cf*d_field(i,j,km1) )
          af(i,j,k) = rbf * af(i,j,k)
        END DO
      END DO
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!  4.4 Top "boundary" layer FIELD row of matrix. FIELD for this layer
!      can then be, and is, updated.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      d_field(i,j,bl_levels) = dtrdz(i,j,bl_levels) *                   &
                               f_field(i,j,bl_levels)

      cf = -dtrdz(i,j,bl_levels) * gamma_rhokh_rdz(i,j,bl_levels)
      rbf = 1.0 / ( 1.0 - cf*( 1.0 + af(i,j,blm1) ) )
      d_field(i,j,bl_levels) = rbf * ( d_field(i,j,bl_levels)           &
                                          - cf*d_field(i,j,blm1) )
      field(i,j,bl_levels) = field(i,j,bl_levels) +                     &
                             d_field(i,j,bl_levels)
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!  5.  "Downward sweep" through whole matrix.  FIELD is updated when
!      the final implicit increments have been calculated.
!-----------------------------------------------------------------------
!  5.1 Remaining FIELD rows of matrix and add implicit increments
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO jj=pdims%j_start, pdims%j_end, omp_block
    DO k = blm1, 1, -1
      DO j = jj, MIN(jj+omp_block-1,pdims%j_end)
        DO i = pdims%i_start, pdims%i_end
          d_field(i,j,k) = d_field(i,j,k) -                             &
                         af(i,j,k)*d_field(i,j,k+1)
          field(i,j,k) = field(i,j,k) + d_field(i,j,k)
        END DO
      END DO
    END DO
  END DO
!$OMP END DO

!-----------------------------------------------------------------------
!  6.  Calculate final implicit flux of tracer.
!-----------------------------------------------------------------------
!! 6.0 First divide by r*r for diagnostics.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        rr_sq = 1.0 / ( r_rho_levels(i,j,k)*r_rho_levels(i,j,k) )
        gamma_rhokh_rdz(i,j,k) = rr_sq * gamma_rhokh_rdz(i,j,k)
        f_field(i,j,k)         = rr_sq * f_field(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      rr_sq = 1.0 / ( r_theta_levels(i,j,0)*r_theta_levels(i,j,0) )
      gamma_rhok_dep(i,j) = rr_sq * gamma_rhok_dep(i,j)
      f_field(i,j,1)      = rr_sq * f_field(i,j,1)
      surf_dep_flux(i,j)  = rr_sq * surf_dep_flux(i,j)
    END DO
  END DO
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------
!  6.1 Surface fluxes.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      surf_dep_flux(i,j) = surf_dep_flux(i,j)                           &
                          - gamma_rhok_dep(i,j) * d_field(i,j,1)
    END DO
  END DO
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------
!  6.2 Fluxes at layer interfaces above the surface.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    km1 = k-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

!  Calculate and store implicit fluxes due to local mixing.

        f_field(i,j,k) = f_field(i,j,k) - gamma_rhokh_rdz(i,j,k)        &
                         * ( d_field(i,j,k) - d_field(i,j,km1) )
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL 

  IF (lhook) CALL dr_hook('IMP_MIX',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE imp_mix
END MODULE imp_mix_mod
