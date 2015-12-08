! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE BDY_EXPL1----------------------------------------------
!
!  Purpose: Calculate explicit scalling parameters and additional
!           boundary layer information required by the surface
!           exchange scheme
!
! Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 24.
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Boundary Layer
!
!---------------------------------------------------------------------
SUBROUTINE bdy_expl1 (                                                  &
! IN values defining vertical grid of model atmosphere :
 bl_levels,                                                             &
 p,p_theta_levels, rho_rsq, rho_wet, rho_dry,                           &
! IN cloud data :
 cf,q,qcf,qcl,t,                                                        &
! IN everything not covered so far :
 pstar,lq_mix_bl,                                                       &
! OUT
 dtrdz_charney_grid,rdz_charney_grid,dtrdz_u,dtrdz_v,                   &
 rdz_u,rdz_v,rho_uv,rho_dry_tq,dzl_charney,rdz,                         &
 z1_tq,z1_tq_top,z1_uv,z1_uv_top,                                       &
 p_half,deltap,qw,tl,                                                   &
 bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt                     &
 )

  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, pdims, pdims_s, tdims
  USE atmos_constants_mod, ONLY: cp
  USE water_constants_mod, ONLY: lc, lf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

!  IN data

! Defining vertical grid of model atmosphere.

  INTEGER, INTENT(IN) ::                                                &
   bl_levels
                                   ! IN Max. no. of "boundary" levels

  REAL, INTENT(IN) ::                                                   &
    p(pdims_s%i_start:pdims_s%i_end,                                    &
      pdims_s%j_start:pdims_s%j_end, bl_levels+1),                      &
    p_theta_levels(tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end, bl_levels+1),             &
    rho_rsq(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end,bl_levels+1),                 &
                                               ! IN Density * R**2
    rho_wet(pdims%i_start:pdims%i_end,                                  &
            pdims%j_start:pdims%j_end, bl_levels+1),                    &
                        ! IN wet density on rho levels (kg/m3)
    rho_dry(pdims%i_start:pdims%i_end,                                  &
            pdims%j_start:pdims%j_end, bl_levels+1)
                        ! IN dry density on rho levels (kg/m3)

! Cloud data.

  REAL, INTENT(IN) ::                                                   &
   cf(tdims%i_start:tdims%i_end,                                        &
      tdims%j_start:tdims%j_end, bl_levels),                            &
                                          ! IN Cloud fraction (decimal).
   qcf(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,bl_levels),                            &
                                     ! IN Cloud ice (kg per kg air)
   qcl(tdims%i_start:tdims%i_end,                                       &
       tdims%j_start:tdims%j_end,bl_levels),                            &
                                     ! IN Cloud liquid water
   q(tdims%i_start:tdims%i_end,                                         &
     tdims%j_start:tdims%j_end,bl_levels),                              &
                                     ! IN specific humidity
   t(tdims%i_start:tdims%i_end,                                         &
     tdims%j_start:tdims%j_end,bl_levels)       ! IN temperature

! Atmospheric + any other data not covered so far, incl control.

  REAL, INTENT(IN) ::                                                   &
   pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! IN Surface pressure (Pascals).

  LOGICAL, INTENT(IN) :: lq_mix_bl  ! IN switch for using mixing ratios

! OUT data

  REAL, INTENT(OUT) ::                                                  &
   dtrdz_charney_grid(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end,bl_levels),             &
                                ! OUT dt/(rho*r*r*dz) for scalar
                                !     flux divergence
   rdz_charney_grid(tdims%i_start:tdims%i_end,                          &
                    tdims%j_start:tdims%j_end,bl_levels),               &
                                ! OUT RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
   dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                           ! OUT dt/(rho*r*r*dz) for
   dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                           ! OUT U,V flux divergence
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
          2:bl_levels),                                                 &
                                ! OUT 1/(Z_U(K)-Z_U(K-1)) for K > 1
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels),                                                 &
                                ! OUT 1/(Z_V(K)-Z_V(K-1)) for K > 1
   rho_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels+1),                                                 &
                                ! OUT density on UV (ie. rho) levels;
                                !    used in RHOKH so dry density if
                                !    Lq_mix_bl is true
   rho_dry_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
             bl_levels),                                                &
                                ! OUT density on TQ (ie. theta) levels;
                                !    used in non-turb flux integration
                                !    so dry density if Lq_mix_bl is true
   dzl_charney(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               bl_levels),                                              &
                                ! OUT DZL(,K) is depth in m of theta
!                                 level K, i.e. distance from boundary
!                                 K-1/2 to boundary K+1/2.
   rdz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                ! OUT RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
   z1_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                                ! OUT Height of lowest theta level.
   z1_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                ! OUT Height of lowest u,v level.
   p_half(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels),                                                   &
                                ! OUT P_HALF(*,K) is pressure at half
                                ! level k-1/2. (rho levels)
   deltap(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels),                                                   &
                                ! OUT Difference in pressure between
                                ! levels
   qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                ! OUT Total water content
   tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                ! OUT Ice/liquid water temperature
   bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                ! OUT A buoyancy parameter for clear air
                                ! on p,T,q-levels (full levels).
   bq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                ! OUT A buoyancy parameter for clear air
                                ! on p,T,q-levels (full levels).
   bt_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                ! OUT A buoyancy parameter for cloudy
                                ! air on p,T,q-levels (full levels).
   bq_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                ! OUT A buoyancy parameter for cloudy
                                ! air on p,T,q-levels (full levels).
   bt_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! OUT A grid-box mean buoyancy parameter
                                ! on p,T,q-levels (full levels).
   bq_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! OUT A grid-box mean buoyancy parameter
                                ! on p,T,q-levels (full levels).
   a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                                ! OUT Saturated lapse rate factor
                                ! on p,T,q-levels (full levels).
   a_dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           bl_levels),                                                  &
                                ! OUT Saturated lapse rate factor
                                ! on p,T,q-levels (full levels).
   dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                ! OUT Derivative of q_SAT w.r.t. T
  REAL, INTENT(OUT) :: z1_uv_top(pdims%i_start:pdims%i_end,             &
                                 pdims%j_start:pdims%j_end)
                                ! Height of top of lowest uv-layer
                                ! above the surface
  REAL, INTENT(OUT) :: z1_tq_top(tdims%i_start:tdims%i_end,             &
                                 tdims%j_start:tdims%j_end)
                                ! Height of top of lowest Tq-layer
                                ! above the surface
!-----------------------------------------------------------------------
!   Symbolic constants (parameters) read in top-level routine :-

! Derived local parameters.

  REAL :: lcrcp,ls,lsrcp

  PARAMETER (                                                           &
   lcrcp=lc/cp,                                                         &
                             ! Evaporation-to-dT conversion factor.
   ls=lf+lc,                                                            &
                             ! Latent heat of sublimation.
   lsrcp=ls/cp                                                          &
                             ! Sublimation-to-dT conversion factor.
    )

!-----------------------------------------------------------------------
!  Workspace :-

!  Local scalars :-

  INTEGER ::                                                            &
   i,j,                                                                 &
                  ! LOCAL Loop counter (horizontal field index).
   k          ! LOCAL Loop counter (vertical level index).

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('BDY_EXPL1',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 1.2 Calculate layer depths and heights, and construct wind fields on
!     P-grid.
!-----------------------------------------------------------------------

! DEPENDS ON: vertical_diffs
    CALL vertical_diffs(                                                &
! IN Field dimensions/logicals
      bl_levels, lq_mix_bl,                                             &
! IN Fields.
      rho_rsq, rho_wet, rho_dry,                                        &
! OUT Vertical differences required by physics.
      rho_uv, rho_dry_tq,                                               &
      dzl_charney, rdz,                                                 &
      z1_uv, z1_uv_top, z1_tq, z1_tq_top,                               &
      rdz_charney_grid,                                                 &
      dtrdz_charney_grid,                                               &
      dtrdz_u, dtrdz_v, rdz_u, rdz_v                                    &
      )

! set pressure array.
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        p_half(i,j,1) = pstar(i,j)
        deltap(i,j,1) = p(i,j,2) - pstar(i,j)
      END DO
    END DO
    DO k = 2, bl_levels
      DO j = pdims%j_start,pdims%j_end
        DO i = pdims%i_start,pdims%i_end
          p_half(i,j,k) = p(i,j,k)
          deltap(i,j,k) = p(i,j,k+1) - p(i,j,k)
        END DO
      END DO
    END DO  ! end of loop over bl_levels

!-----------------------------------------------------------------------
! Calculate total water content, QW and Liquid water temperature, TL
!-----------------------------------------------------------------------

  DO k = 1, bl_levels
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        qw(i,j,k) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)
                                ! P243.10
        tl(i,j,k) = t(i,j,k) - lcrcp*qcl(i,j,k) - lsrcp*qcf(i,j,k)
                                ! P243.9
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
! 5.1  Calculate buoyancy parameters BT and BQ.
!-----------------------------------------------------------------------

! DEPENDS ON: buoy_tq
  CALL buoy_tq (                                                        &
! IN dimensions/logicals
   bl_levels, lq_mix_bl,                                                &
! IN fields
   p_theta_levels,t,q,qcf,qcl,cf,                                       &
! OUT fields
   bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt                   &
    )

  IF (lhook) CALL dr_hook('BDY_EXPL1',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE bdy_expl1
