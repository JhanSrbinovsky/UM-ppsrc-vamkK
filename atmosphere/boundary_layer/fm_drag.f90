! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Calculate form drag profiles for distributed
!           drag parametrization
!
!           Based on work by Wood, Brown and Hewer (2001),
!           Quart. J. Roy. Met. Soc., 127, 759--777.
!
!  Programming standard: UMDP 3
!
! Description:
! Turbulent form drag due to sub-grid orography

! Method:
! Calculates the form drag and stress profiles due to sub-grid scale
! orography. The stress is later added as an additional explicit
! stress to the boundary-layer equations in EXFXUV8A. The orographic
! stress is added to the surface stress in BDYLYR8A.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer

! Code Description:
! Language: FORTRAN 77 + common extensions.
! This code is written to UMDP3 v6 programming standards.
!---------------------------------------------------------------------
SUBROUTINE fm_drag (                                                    &
! IN levels
  land_pts, land_index, bl_levels,                                      &
! IN fields
  u_p, v_p, rho_tq, z_uv, z_tq, z0m, zh, Rib, sil_orog_land,            &
! OUT fields
  tau_fd_x, tau_fd_y                                                    &
  )

  USE atm_fields_bounds_mod, ONLY: pdims, tdims
  USE atmos_constants_mod, ONLY: vkman
  USE c_surf, ONLY: ri_crit
  USE conversions_mod, ONLY: pi
  USE bl_option_mod, ONLY: fd_stab_dep, on, orog_drag_param
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

  INTEGER, INTENT(IN):: land_pts,                                       &
                                             ! Number of land points
                        bl_levels        ! No. of levels for which
                                             ! boundary layer fluxes
                                             ! are calculated

  INTEGER, INTENT(IN):: land_index(land_pts) ! Index for compressed
                                                 ! land point array; ith
                                                 ! element holds
                                                 ! position in the FULL
                                                 ! field of the ith land
                                                 ! point to be processed

  REAL, INTENT(IN)::                                                    &
   u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                      ! Wind component in x direction
                                      ! horizontally interpolated to
                                      ! P-grid (m/s)
   v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                      ! Wind component in y direction
                                      ! horizontally interpolated to
                                      ! P-grid (m/s)
   rho_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                      ! For a vertically staggered grid
                                      ! with a u,v-level first above the
                                      ! surface, RHO_TQ(*,K) is the
                                      ! density of the k-th T,q-level
                                      ! above the surface;
                                      ! for an unstaggered grid the
                                      ! densities at the layer interface
                                      ! (half-levels) 1.5 to BL_LEVELS+0
                                      ! should be input to elements 1 to
                                      ! BL_LEVELS.
                                      ! (Value for BL_LEVELS not used
                                      ! in either case.)
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels), &
                                      ! For a vertically staggered grid
                                      ! with a u,v-level first above the
                                      ! surface, Z_UV(*,K) is the height
                                      ! of the k-th u,v-level(half level
                                      ! k-1/2) above the surface;
                                      ! for an unstaggered grid the
                                      ! heights of the half-levels
                                      ! 0.5 to BL_LEVELS-0.5 should be
                                      ! input to elements 1 to BL_LEVELS
                                      ! (1st value not used in either
                                      !  case.)
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                                      ! For a vertically staggered grid
                                      ! with a u,v-level first above the
                                      ! surface, Z_TQ(*,K) is the height
                                      ! of the k-th T,q-level (full
                                      ! level k) above the surface;
                                      ! for an unstaggered grid the
                                      ! heights of the half levels
                                      ! 1.5 to BL_LEVELS+0.5 should be
                                      ! input to elements 1 to BL_LEVELS
                                      ! (Value for BL_LEVELS not used
                                      ! in either case.)
   z0m(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                                      ! Roughness length for momentum (m
   zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
                                      ! Boundary layer height (actually
                                      ! ZH_PREV)
   Rib(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                                      ! Bulk Richardson number for
                                      ! lowest layer
   sil_orog_land(land_pts)
                                      ! Silhouette area of unresolved
                                      ! orography per unit hoz. area

  REAL, INTENT(OUT)::                                                   &
   tau_fd_x(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                                           ! X-comp of orographic stress
   tau_fd_y(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels) ! Y-comp of orographic stress
                                           ! (N/m^2)

!     Local and other symbolic constants:

!     Define local storage

!     Local arrays:

  REAL::                                                                &
    h_m(land_pts),                                                      &
                                     ! height at which form drag is
                                     ! calculated (m)
    u_hm(land_pts),                                                     &
                                     ! X-component of wind at height h_m
    v_hm(land_pts),                                                     &
                                     ! Y-component of wind at height h_m
    fp_x(land_pts),                                                     &
                                     ! X-component of pressure force
    fp_y(land_pts)               ! Y-component of pressure force

!     Local scalars:

  REAL::                                                                &
    zeta,                                                               &
                      ! Log(h_m/z0m)
    height_fac,                                                         &
                      ! Height dependency factor for form drag
    wta, wtb,                                                           &
                      ! weights for interpolation between levels
    tausx,tausy,                                                        &
                      ! Surface stress
    rib_fn        ! Richardson number function for stability
                      ! correction to drag

  INTEGER::                                                             &
    i,                                                                  &
                      ! Loop counter (horizontal field index)
    j,                                                                  &
                      ! Loop counter (offset within I loop)
    k,                                                                  &
                      ! Loop counter (vertical level index)
    l             ! Loop counter (horizontal land field index)

! Local parameters
      ! Tunable parameters in calculation of explicit orographic stresss
  REAL,PARAMETER:: alpha    = 12.0,                                     &
                                             ! Tunable parameter for form
                   beta     = 1.0,                                      &
                                             ! drag calculation.
                   fd_decay = 0.3333,                                   &
                                             ! Height scale factors for
                   max_ht_scale  = 300.0,                               &
                                             ! stress profiles
                   min_ht_scale  = 100.0
  LOGICAL,PARAMETER:: l_lowhill=.false.  ! Set to .TRUE. for Wood &
                                         ! Mason (1993) drag formula
                                         ! Set to .FALSE. for steep
                                         ! hill expression

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('FM_DRAG',zhook_in,zhook_handle)

!----------------------------------------------------------------
! 1. Calculate the height scale h_m and interpolate the wind and
!    density to this height.
!----------------------------------------------------------------

  DO l = 1, land_pts
    j = (land_index(l)-1)/pdims%i_end + 1
    i = land_index(l) - (j-1)*pdims%i_end

    h_m(l)=MIN(max_ht_scale,fd_decay*zh(i,j))
    h_m(l)=MAX(min_ht_scale,h_m(l))
    u_hm(l)=0.0
    v_hm(l)=0.0
  END DO

!     Interpolate to get U and V at z=h_m

  DO k = 2, bl_levels
    DO l = 1, land_pts

      j = (land_index(l)-1)/pdims%i_end + 1
      i = land_index(l) - (j-1)*pdims%i_end

      IF(h_m(l) <= z_uv(i,j,k).AND.h_m(l) >= z_uv(i,j,k-1)) THEN

        wta = ( h_m(l) - z_uv(i,j,k-1) )                                &
                 /( z_uv(i,j,k) - z_uv(i,j,k-1) )
        wtb = ( z_uv(i,j,k) - h_m(l) )                                  &
                 /( z_uv(i,j,k) - z_uv(i,j,k-1) )
        u_hm(l)   = wta*u_p(i,j,k) + wtb*u_p(i,j,k-1)
        v_hm(l)   = wta*v_p(i,j,k) + wtb*v_p(i,j,k-1)

      END IF

    END DO
  END DO

!-----------------------------------------------------------------------
! 2. Calculate the pressure force from the wind components and density
!    at height h_m, the frontal silhouette area and surface roughness
!    length for momentum, z_0m.
!-----------------------------------------------------------------------
  DO l = 1, land_pts

    j = (land_index(l)-1)/pdims%i_end + 1
    i = land_index(l) - (j-1)*pdims%i_end

    IF ( fd_stab_dep == on) THEN
      rib_fn=1.-Rib(i,j)/ri_crit
      IF(rib_fn >  1.)rib_fn=1.
      IF(rib_fn <  0.)rib_fn=0.
    ELSE
      rib_fn=1.
    END IF

    zeta = LOG( h_m(l)/z0m(i,j) )

    IF(l_lowhill)THEN

!          Compute Wood and Mason (1993) low-hill drag expression

      tausx=(vkman/zeta)*(vkman/zeta)*u_hm(l)*                          &
            SQRT(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
      tausy=(vkman/zeta)*(vkman/zeta)*v_hm(l)*                          &
            SQRT(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
      fp_x(l) = rho_tq(i,j,1)*alpha*beta*pi*pi                          &
              *sil_orog_land(l)*sil_orog_land(l)                        &
              *rib_fn*tausx
      fp_y(l) = rho_tq(i,j,1)*alpha*beta*pi*pi                          &
              *sil_orog_land(l)*sil_orog_land(l)                        &
              *rib_fn*tausy
    ELSE

!          Compute steep hill drag expression

      tausx=u_hm(l)*SQRT(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
      tausy=v_hm(l)*SQRT(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
      fp_x(l)=0.5*rho_tq(i,j,1)*orog_drag_param*                        &
              sil_orog_land(l)*rib_fn*tausx
      fp_y(l)=0.5*rho_tq(i,j,1)*orog_drag_param*                        &
              sil_orog_land(l)*rib_fn*tausy

    END IF

    tau_fd_x(i,j,1) = fp_x(l)
    tau_fd_y(i,j,1) = fp_y(l)

  END DO

!-----------------------------------------------------------------------
! 3. Calculate the vertical profiles of the explicit orographic stress
!-----------------------------------------------------------------------
  DO k = 2, bl_levels
    DO l = 1, land_pts

      j = (land_index(l)-1)/pdims%i_end + 1
      i = land_index(l) - (j-1)*pdims%i_end

      height_fac = EXP(z_tq(i,j,k-1)/h_m(l))
      tau_fd_x(i,j,k) = tau_fd_x(i,j,1)/height_fac
      tau_fd_y(i,j,k) = tau_fd_y(i,j,1)/height_fac

    END DO
  END DO

  IF (lhook) CALL dr_hook('FM_DRAG',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE fm_drag
