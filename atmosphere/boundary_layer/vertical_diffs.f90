! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! SUBROUTINE VERTICAL_DIFFS
! PURPOSE: Calculate vertical differences required by BL scheme
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!
! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 95
!   THIS CODE IS WRITTEN TO UMDP3 PROGRAMMING STANDARDS.
!
SUBROUTINE vertical_diffs (                                             &
! IN Field dimensions/logicals
  bl_levels, lq_mix_bl,                                                 &
! IN Fields.
  rho_p_rsq, rho_wet, rho_dry,                                          &
! OUT Vertical differences required by physics.
  rho_uv, rho_dry_tq,                                                   &
  dzl_charney, rdz,                                                     &
  z1_uv, z1_uv_top, z1_tq, z1_tq_top,                                   &
  rdz_charney_grid,                                                     &
  dtrdz_charney_grid,                                                   &
  dtrdz_u, dtrdz_v, rdz_u, rdz_v                                        &
  )

  USE level_heights_mod, ONLY:                                          &
    r_theta_levels, r_rho_levels
  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, pdims_s, pdims, pdims_l, tdims
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE switches, ONLY: i_modiscopt
  USE timestep_mod, ONLY: timestep
  USE bl_option_mod, ONLY: on

  USE p_to_t_mod, ONLY: p_to_t
  USE p_to_u_mod, ONLY: p_to_u
  USE p_to_v_mod, ONLY: p_to_v
  IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

! LOGICAL SWITCHES

  LOGICAL, INTENT(IN) :: lq_mix_bl ! IN switch for using mixing ratios

  INTEGER, INTENT(IN) :: bl_levels

  REAL, INTENT(IN) ::                                                   &
    rho_p_rsq(pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,bl_levels+1),               &
                        ! wet density times r^2 on rho levels (kg/m3)
    rho_wet(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels+1),                                               &
                        ! wet density on rho levels (kg/m3)
    rho_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels+1)
                        ! dry density on rho levels (kg/m3)

! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

  REAL, INTENT(OUT) ::                                                  &
    rho_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels+1),                                                &
                                ! OUT density on UV (ie. rho) levels;
                                !    used in RHOKH so dry density if
                                !    Lq_mix_bl is true
                                
    rho_dry_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               bl_levels),                                              &
                                ! OUT density on TQ (ie. theta) levels;
                                !    used in non-turb flux integration
                                !    so dry density if Lq_mix_bl is true
    dtrdz_charney_grid (tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end,bl_levels),           &
    dzl_charney(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                bl_levels),                                             &
    rdz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
           bl_levels),                                                  &
                                ! OUT  Reciprocal of distance between
    rdz_charney_grid(tdims%i_start:tdims%i_end,                         &
                     tdims%j_start:tdims%j_end,bl_levels),              &
    z1_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
    z1_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
  REAL, INTENT(OUT) :: z1_uv_top(pdims%i_start:pdims%i_end,             &
                                 pdims%j_start:pdims%j_end)
                                ! Top of lowest uv-layer
  REAL, INTENT(OUT) :: z1_tq_top(tdims%i_start:tdims%i_end,             &
                                 tdims%j_start:tdims%j_end)
                                ! Top of lowest Tq-layer

  REAL, INTENT(OUT) ::                                                  &
              ! quantities interpolated to u and v grids.
    rdz_u (udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
          2:bl_levels),                                                 &
    rdz_v (vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
          2:bl_levels),                                                 &
    dtrdz_u (udims%i_start:udims%i_end,udims%j_start:udims%j_end,       &
             bl_levels),                                                &
    dtrdz_v (vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,       &
             bl_levels)

! LOCAL VARIABLES.
  INTEGER ::                                                            &
    k,                                                                  &
    i,j

  REAL ::                                                               &
    dz_p,                                                               &
    dz_t

  REAL ::                                                               &
    work1(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end, bl_levels),                    &
    work2(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end, bl_levels+1),                  &
    work3(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end, bl_levels)

  REAL ::                                                               &
    dzl_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,          &
             bl_levels),                                                &
    dzl_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,          &
             bl_levels)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! CALCULATE rho on layer centres and layer boundaries
! ----------------------------------------------------------------------
  IF (lhook) CALL dr_hook('VERTICAL_DIFFS',zhook_in,zhook_handle)
  IF ( lq_mix_bl ) THEN
        ! conservation of mixing ratios requires use of dry density
    DO k = 1, bl_levels+1
      DO j = pdims%j_start,pdims%j_end
        DO i = pdims%i_start,pdims%i_end
          rho_uv(i,j,k) = rho_dry(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    DO k = 1, bl_levels+1
      DO j = pdims%j_start,pdims%j_end
        DO i = pdims%i_start,pdims%i_end
          rho_uv(i,j,k) = rho_wet(i,j,k)
        END DO
      END DO
    END DO
  END IF

! Interpolate RHO_uv to temperature levels for RHO_DRY_tq
!  - can't think of a better name but this will be wet or dry
!    depending on lq_mix_bl

  CALL p_to_t (                                                         &
! IN Field dimensions and pointers.
    pdims%i_end, pdims%j_end, pdims_l%halo_i, pdims_l%halo_j, 0, 0,     &
    bl_levels,                                                          &
! IN Vertical coordinate levels.
    r_theta_levels, r_rho_levels,                                       &
! IN field to be interpolated.
    rho_uv,                                                             &
! OUT Interpolated field.
    rho_dry_tq                                                          &
   )
! ----------------------------------------------------------------------
! CALCULATE DTRDZ, DZL, RDZ
! ----------------------------------------------------------------------
  k = 1
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      dz_p = r_theta_levels(i,j,k)-r_theta_levels(i,j,k-1)
      dz_t = r_rho_levels(i,j,k+1)-r_theta_levels(i,j,k-1)
      rdz(i,j,k) = 1./                                                  &
                     (r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1))
      z1_uv(i,j) = r_rho_levels(i,j,k) - r_theta_levels (i,j,k-1)
      z1_tq(i,j) = r_theta_levels(i,j,k) -                              &
                       r_theta_levels (i,j,k-1)
      rdz_charney_grid(i,j,k) = 1./ dz_p
          ! Following is dt/(rho * r^2 * Dz), used for flux-divergence
          ! equation for scalars - hence dry density if lq_mix_bl
          ! is true, wet otherwise
      dtrdz_charney_grid(i,j,k) = timestep/                             &
                ( r_theta_levels (i,j,k)*r_theta_levels (i,j,k)*        &
                                        rho_dry_tq(i,j,k)*dz_t )
! Dzl_charney( ,,1) is such that
! 0.5 * DZL_Charney = height of first temperature level
      dzl_charney(i,j,k) = 2. *                                         &
                 (r_theta_levels(i,j,1) - r_theta_levels(i,j,0))
    END DO
  END DO

!     Additionally calculate the tops of the lowest layers if using
!     fully conservative discretization.
  IF (i_modiscopt == on) THEN
    z1_uv_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
      = r_theta_levels(tdims%i_start:tdims%i_end,                       &
                       tdims%j_start:tdims%j_end,k) -                   &
        r_theta_levels(tdims%i_start:tdims%i_end,                       &
                       tdims%j_start:tdims%j_end,k-1)
    z1_tq_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)      &
      = r_rho_levels(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end,k+1) -                   &
        r_theta_levels(tdims%i_start:tdims%i_end,                       &
                       tdims%j_start:tdims%j_end,k-1)
  END IF

  DO k = 2, bl_levels
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        dz_p = r_theta_levels (i,j,k) - r_theta_levels(i,j,k-1)
        dz_t = r_rho_levels (i,j,k+1) - r_rho_levels(i,j,k)
        rdz(i,j,k) = 1./                                                &
                  (r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1))
        dzl_charney(i,j,k) = dz_t
        rdz_charney_grid(i,j,k) = 1./ dz_p
          ! Following is dt/(rho * r^2 * Dz), used for flux-divergence
          ! equation for scalars - hence dry density if lq_mix_bl
          ! is true, wet otherwise
        dtrdz_charney_grid(i,j,k) = timestep/                           &
                ( r_theta_levels (i,j,k)*r_theta_levels (i,j,k)*        &
                                        rho_dry_tq(i,j,k)*dz_t )
      END DO
    END DO
  END DO

  DO k = 1, bl_levels
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        work1(i,j,k) = r_theta_levels(i,j,k)-r_theta_levels(i,j,k-1)
      END DO
    END DO
  END DO

  DO k = 1, bl_levels
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        work3(i,j,k) = timestep/( rho_p_rsq(i,j,k) * work1(i,j,k) )
      END DO
    END DO
  END DO

  DO k = 2, bl_levels
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
            ! Note work1(K=1) is not being set here so dzl_u,v(K=1),
            ! calculated from it below, should not be used
        work1(i,j,k) = r_rho_levels(i,j,k)-r_rho_levels(i,j,k-1)
      END DO
    END DO
  END DO

! ----------------------------------------------------------------------
! horizontal interpolations to momentum points.
! ----------------------------------------------------------------------

      CALL p_to_u (work3,                                         &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,bl_levels,dtrdz_u)  


      CALL p_to_u (work1,                                         &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,bl_levels,dzl_u)  

      CALL p_to_v (work3,                                         &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   1,bl_levels,dtrdz_v) 

      CALL p_to_v (work1,                                         &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   1,bl_levels,dzl_v)

! calculate rdz_u, rdz_v from local arrays dzl_u, dzl_v.

  DO k = 2, bl_levels
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        rdz_u(i,j,k) = 1.0/dzl_u(i,j,k)
      END DO
    END DO
  END DO

  DO k = 2, bl_levels
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        rdz_v(i,j,k) = 1.0/dzl_v(i,j,k)
      END DO
    END DO
  END DO

  IF (lhook) CALL dr_hook('VERTICAL_DIFFS',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE vertical_diffs
