! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE conv_surf_flux_mod

IMPLICIT NONE

CONTAINS

! Calculate surface fluxes for convective diagnosis

SUBROUTINE conv_surf_flux(                                        &
  row_length, rows                                                &
, model_levels, wet_model_levels, land_points                     &

, l_mixing_ratio,l_ctile, l_flux_bc, l_spec_z0                    &
, land_mask                                                       &

, pstar, tstar_land, tstar_sea, tstar_sice, zh, flandg            &
, ice_fract, u_p, v_p, u_0_p, v_0_p                               &
, flux_e, flux_h,  z0msea, z0m_scm, z0h_scm                       &
, z_full, q, theta, exner_theta_levels                            &
! INOUT 
, tstar , fb_surf, tv1_sd, bl_vscale2                             &
 )

!-----------------------------------------------------------------------
! Purpose:
!  Calculate the surface buoyancy flux.
!  Also calculates an approximate standard deviation of virtual 
!  temperature at level 1
!
!  Part of convective diagnosis routines
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3  programming standards v8.3
!
!-----------------------------------------------------------------------

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                    &
  pdims, pdims_s, tdims, qdims

USE cv_run_mod, ONLY:                                               &
  tv1_sd_opt

USE earth_constants_mod, ONLY: g

USE atmos_constants_mod, ONLY:                                      &
  vkman, cp, kappa, r, repsilon, c_virtual

USE water_constants_mod, ONLY: lc, lf, tm

USE cv_derived_constants_mod, ONLY:   gamma_dry 


USE surf_param, ONLY : z0sice, z0h_z0m_sice
USE c_rough,    ONLY : z0hsea


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

INTEGER, INTENT(IN) ::  &
  row_length            & ! Local number of points on a row
 ,rows                  & ! Local number of rows in a theta field
 ,model_levels          & ! number of model levels
 ,wet_model_levels      & ! number of wet model levels
 ,land_points             ! number of land points

LOGICAL,INTENT(IN) ::   &
  l_mixing_ratio        & ! true moisture input as mixing ratios
                          ! false moisture input as specific humidity
 ,l_ctile               & ! true if coastal tiling
 ,l_flux_bc             & ! true if SCM using specified surface fluxes
 ,l_spec_z0               ! true if roughness length has been specIFied

LOGICAL,INTENT(IN) ::          &
  land_mask(row_length, rows)    ! T if land, F elsewhere.

REAL, INTENT(IN) ::            &
  pstar(row_length, rows)      & ! Surface pressure (Pa)
 ,tstar_land(row_length, rows) & ! Surface T on land
 ,tstar_sea(row_length, rows)  & ! Surface T on sea
 ,tstar_sice(row_length, rows) & ! Surface T on sea-ice 
 ,zh(row_length,rows)          & ! Height above surface of top
                                 !  of boundary layer (metres).
 ,ice_fract(row_length,rows)     ! fraction of sea that has ice

REAL, INTENT(IN) ::                             &
  flandg(pdims_s%i_start:pdims_s%i_end,         & ! Land fraction of gridbox
         pdims_s%j_start:pdims_s%j_end)         & ! on all points
 ,z_full(row_length,rows,1:tdims%k_end)         & ! height th lev (m)
 ,q(qdims%i_start:qdims%i_end,                  & ! water vapour (kg/kg)
    qdims%j_start:qdims%j_end,                  &
                1:qdims%k_end)                  &
 ,theta(tdims%i_start:tdims%i_end,              & ! Theta (Kelvin)
        tdims%j_start:tdims%j_end,              &
                    1:tdims%k_end)              &
 ,exner_theta_levels(tdims%i_start:tdims%i_end, & ! exner pressure theta lev
                     tdims%j_start:tdims%j_end, & !  (Pa)
                                 1:tdims%k_end)

REAL, INTENT(IN) ::            &
  u_p(row_length, rows)        & ! U(1) on P-grid.
 ,v_p(row_length, rows)        & ! V(1) on P-grid.
 ,u_0_p(row_length,rows)       & ! W'ly component of surface current
                                 !    (metres per second) on P-grid.
 ,v_0_p(row_length,rows)       & ! S'ly component of surface current
                                 !    (metres per second) on P-grid.
 ,flux_e(row_length,rows)      & ! Specified surface
                                 !    latent heat flux (W/m^2)
 ,flux_h(row_length,rows)      & ! Specified surface
                                 !    sensible heat fluxes (in W/m2)
 ,z0msea(row_length,rows)      & ! Sea roughness length for momentum (m)
 ,z0m_scm(row_length,rows)     & ! Namelist input z0m (if >0)
 ,z0h_scm(row_length,rows)       ! Namelist input z0h (if >0)

REAL, INTENT(INOUT) ::          &
  tstar(row_length,rows)        & ! Surface temperature 
                                  ! (= top soil layer temperature) (K).
 ,fb_surf(row_length, rows)     & ! Change in theta_v from surface
                                  ! to layer 1 (note diff from BL)
 ,tv1_sd( row_length*rows)      & ! Approx to standard dev of level 1
                                  ! virtual temperature (K).(unstable points)
 ,bl_vscale2(row_length*rows)     ! Velocity scale squared for 
                                  ! boundary layer eddies (m2/s2)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER ::  routinename = 'conv_surf_flux'

INTEGER ::   &
  i,j,ii     & ! Local Loop counter (horizontal field index).
 ,nunstable    ! total number of unstable points

REAL ::                         &
  qs_star(row_length, rows)     & ! Saturated sp humidity at surface
 ,qs_star_sice(row_length, rows)& ! Saturated sp humidity at sea-ice surface
 ,dqsdt(row_length,rows)        & ! d(qsat)/dT 
 ,z0(row_length,rows)           & ! roughness length (m)
 ,z0m_land(row_length,rows)     & ! roughness length for momentum over land (m)
 ,z0h_land(row_length,rows)     & ! roughness length for heat over land (m)
 ,z0m_sea(row_length,rows)      & ! roughness length for momentum over sea(m)
 ,z0h_sea(row_length,rows)        ! roughness length for heat over sea(m)


! Used in calculation to decide on unstable points 

REAL ::           &
  theta1          &  ! Potential temperature in layer 1
 ,ushear          &  ! U wind shear from level 1 to surface
 ,vshear          &  ! V wind shear from level 1 to surface
 ,wshr1           &  ! magnitude of surface wind shear
 ,wshr2           &  ! (wshr1)**2
 ,rhostar         &  ! surface air density
 ,theta_star      &  ! theta at surface
 ,wthvbar         &  ! surface buoyancy flux
 ,cd              &  ! bulk transfer coefficient for momentum
 ,ch              &  ! bulk transfer coefficient for heat
 ,rib             &  ! Ri for surface exchange
 ,ustar           &  ! surface friction velocity
 ,w_m             &  ! 
 ,w_s_cubed       &  !
 ,wstar_tmp       &  ! convective velocity scale
 ,ustar2_sea      &  ! ustar^2 over sea
 ,ustar2_land     &  ! ustar^2 over land
 ,wthvbar_sea     &  ! surface buoyancy flux over sea
 ,wthvbar_land       ! surface buoyancy flux

REAL  :: tv1_sd_temp(row_length, rows),    &
         bl_vscale2_temp(row_length, rows)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('CONV_SURF_FLUX',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Set up roughness lengths
!-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j)

IF (tv1_sd_opt == 2) THEN

  ! using grid-box mean in coastal points, hence need land and sea 
  ! separately
!$OMP DO SCHEDULE(STATIC)
  DO j=1, rows
    DO i=1, row_length
      ! Land: assume z0m=0.1m and z0h=z0m/10.
      z0m_land(i,j) = 0.1
      z0h_land(i,j) = 0.01
      ! Sea: use parametrized values
      !      z0h_sea updated later for low wind speed limit
      z0m_sea(i,j) = z0msea(i,j)
      z0h_sea(i,j) = MAX( 2.56e-9/z0msea(i,j), 7.0e-08 )
    END DO ! i
  END DO ! j
!$OMP END DO

  IF ( L_spec_z0 ) THEN
    ! Code to use z0mh_scm if namelist specifies it
!$OMP DO SCHEDULE(STATIC)
    DO j=1, rows
      DO i=1, row_length
        IF ( z0m_scm(i,j)  >   0.0 ) THEN
          z0m_sea(i,j)  = z0m_scm(i,j)
          z0m_land(i,j) = z0m_scm(i,j)
        END IF
        IF ( z0h_SCM(i,j)  >   0.0 ) THEN
          z0h_sea(i,j)  = z0h_scm(i,j)
          z0h_land(i,j) = z0h_scm(i,j)
        END IF
      END DO ! i
    END DO ! j
!$OMP END DO
  END IF
END IF

!$OMP DO SCHEDULE(STATIC)
DO j=1, rows
  DO i=1, row_length
    IF (land_mask(i,j)) THEN
      ! Approximate z0 for land as 0.1.
      z0(i,j) = 0.1
    ELSE
      z0(i,j) = z0hsea
    END IF
  END DO ! i
END DO ! j
!$OMP END DO

! Code to use z0h_scm IF Namelist specIFies it

IF ( L_spec_z0 ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j=1, rows
    DO i=1, row_length
      IF ( z0h_SCM(i,j)  >   0.0 ) THEN
        z0(i,j) = z0h_SCM(i,j)
      END IF! z0h_scm
    END DO ! i
  END DO ! j
!$OMP END DO
END IF

!$OMP END PARALLEL

!-----------------------------------------------------------------------
! 1.5 Surface buoyancy flux and calculation of unstable points
!-----------------------------------------------------------------------

IF ( .NOT. l_flux_bc) THEN        ! used by most UM runs

  IF (tv1_sd_opt == 2) THEN
  !-----------------------------------------------------------------------
  ! Calculate the surface buoyancy flux
  ! new method includes stability dependence and area mean for coastal 
  ! and sea-ice points
  !-----------------------------------------------------------------------
! DEPENDS ON: qsat_mix
    CALL qsat_mix(qs_star,tstar_sea,pstar,row_length*rows,        &
                     l_mixing_ratio)
! DEPENDS ON: qsat_mix
    CALL qsat_mix(qs_star_sice,tstar_sice,pstar,row_length*rows,  &
                     l_mixing_ratio)
!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(theta1, rhostar, Ushear, wshr1, &
!$OMP& Vshear, wthvbar_sea, rib, cd, ustar2_sea, ch, w_m, wshr2,        &
!$OMP& wthvbar_land, ustar2_land, w_s_cubed, i, j, ustar, wthvbar)

!$OMP DO SCHEDULE(STATIC)
    DO j=1,rows
      DO i=1,row_length

        theta1 = theta(i,j,1)
        rhostar = pstar(i,j) / ( R*tstar(i,j) )

        ushear = u_p(i,j) - u_0_p(i,j)
        vshear = v_p(i,j) - v_0_p(i,j)
        wshr2 = MAX (1.0E-6 , ushear*ushear + vshear*vshear)
        wshr1 = SQRT(wshr2)
        !-------------------------------------------------------------
        ! Sea 
        !-------------------------------------------------------------
        wthvbar_sea  = 0.0
        ustar2_sea   = 0.0
        IF ( flandg(i,j) < 0.99 ) THEN
          ! Include a crude stability dependence
          rib = - ( g / tstar_sea(i,j) ) *                             &
            (tstar_sea(i,j) - ( theta1*exner_theta_levels(i,j,1)       &
                    +gamma_dry*z_full(i,j,1) ) ) * z_full(i,j,1) / wshr2
          cd = ( vkman / LOG(z_full(i,j,1)/z0m_sea(i,j)) )**2
          cd = cd * (1.0 + 0.7 * (MAX(0.0, -rib))**0.33 )
          ustar2_sea = cd * wshr2
          IF ( .NOT.l_spec_z0 ) THEN
            ! include low wind speed limit (larger z0)
            z0h_sea(i,j) = MAX( 2.52e-6/(SQRT(ustar2_sea)+1.0e-05),   &
                                 z0h_sea(i,j) ) 
          END IF
          ch = vkman**2 / ( LOG(z_full(i,j,1)/z0m_sea(i,j)) *         &
                             LOG(z_full(i,j,1)/z0h_sea(i,j) ) )
          ch = ch * (1.0 + (MAX(0.0, -rib))**0.33 )
          wthvbar_sea = ch * wshr1 * ( tstar_sea(i,j) -                     &
             ( theta1*exner_theta_levels(i,j,1) + gamma_dry*z_full(i,j,1) ) &
               + 0.61*theta1*(qs_star(i,j)-q(i,j,1)) )
          !-------------------------------------------------------------
          ! Sea-ice
          !-------------------------------------------------------------
          IF ( ice_fract(i,j) > 0.01 ) THEN
            ! Include a crude stability dependence
            rib = - ( g / tstar_sice(i,j) ) *                            &
              (tstar_sice(i,j) - ( theta1 * exner_theta_levels(i,j,1)    &
                    +gamma_dry*z_full(i,j,1) ) ) * z_full(i,j,1) / wshr2
            cd = ( vkman / LOG(z_full(i,j,1)/z0sice) )**2
            cd = cd * (1.0 + 0.7 * (MAX(0.0, -rib))**0.33 )
            ustar2_sea = (1.0-ice_fract(i,j)) * ustar2_sea +          &
                               ice_fract(i,j) * cd*wshr2
            ch = vkman**2 / ( LOG(z_full(i,j,1)/z0sice) * &
                              LOG(z_full(i,j,1)/(z0sice*z0h_z0m_sice)) )
            ch = ch * (1.0 + (MAX(0.0, -rib))**0.33 )
            wthvbar_sea = (1.0-ice_fract(i,j)) * wthvbar_sea +        &
                        ice_fract(i,j) * ch*wshr1*( tstar_sice(i,j) - &
              ( theta1*exner_theta_levels(i,j,1) + gamma_dry*z_full(i,j,1) )&
                + 0.61*theta1*(qs_star_sice(i,j)-q(i,j,1)) )
          END IF
        END IF

        !-------------------------------------------------------------
        ! Land
        !-------------------------------------------------------------
        wthvbar_land = 0.0
        ustar2_land   = 0.0
        IF ( flandg(i,j) > 0.01 ) THEN
          ! Include a crude stability dependence
          rib = - ( g / tstar_land(i,j) ) * &
            (tstar_land(i,j) - ( theta1 * exner_theta_levels(i,j,1)  &
             +gamma_dry*z_full(i,j,1) ) ) * z_full(i,j,1) / wshr2
          cd = ( vkman / LOG(z_full(i,j,1)/z0m_land(i,j)) ) ** 2
          cd = cd * (1.0 + 0.7 * (MAX(0.0, -rib))**0.33 )
          ustar2_land = cd * wshr2
          ch = vkman**2 / ( LOG(z_full(i,j,1)/z0m_land(i,j)) *      &
                             LOG(z_full(i,j,1)/z0h_land(i,j) ) )
          ch = ch * (1.0 + (MAX(0.0, -rib))**0.33 )
          wthvbar_land = ch * wshr1 * &
            ( tstar_land(i,j) - ( theta1 * exner_theta_levels(i,j,1) &
             +gamma_dry*z_full(i,j,1) ) )
        END IF
        !-------------------------------------------------------------
        ! Combine to cell average values and then take sqrt for ustar
        !-------------------------------------------------------------
        IF ( flandg(i,j) < 0.01 ) THEN
          ustar = ustar2_sea
          wthvbar = wthvbar_sea
        ELSE IF ( flandg(i,j) > 0.99 ) THEN
          ustar = ustar2_land
          wthvbar = wthvbar_land
        ELSE
          ! Take area-weighted mean
          ustar = (1.0 - flandg(i,j) ) * ustar2_sea +                &
                     flandg(i,j) * ustar2_land 
          wthvbar = (1.0 - flandg(i,j) ) * wthvbar_sea +             &
                     flandg(i,j) * wthvbar_land
        END IF

        ustar = SQRT(ustar)
        fb_surf(i,j) = g * wthvbar /( rhostar * theta1*(1.0+0.61*q(i,j,1)) )

        IF (fb_surf(i,j)  >   0.0) THEN
          w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
          w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
          tv1_sd_temp(i, j) = 1.93 * wthvbar/( rhostar * w_m )   
          bl_vscale2_temp(i,j) = 2.5 * 2.52 * w_m * w_m
                  ! 2.5 from Beare (2008), 2.52 to undo the 0.25 factor
          bl_vscale2_temp(i,j) = MAX( 0.01, bl_vscale2_temp(i,j) ) ! for safety
        END IF! (fb_surf > 0.0)

      END DO ! i
    END DO ! j
!$OMP END DO

!$OMP END PARALLEL

    ii = 0
    
    DO j=1, rows
      DO i=1, row_length
        IF (fb_surf(i,j)  >   0.0) THEN
          ii= ii+1
          tv1_sd(ii) = tv1_sd_temp(i,j)
          bl_vscale2(ii) = bl_vscale2_temp(i,j)
        END IF
      END DO
    END DO

  ELSE  ! tv1_sd_opt /= 2
    ! old (neutral stability) method

    ! qsat at surface
! DEPENDS ON: qsat_mix
    CALL qsat_mix(qs_star,tstar,pstar,row_length*rows,l_mixing_ratio)

    !----------------------------------------------------------------
    ! Standard UM code for surface T boundary condition:
    ! Calculate the surface buoyancy flux using
    ! approximation for unstable cd as 1.5*neutral value (defined
    ! Garratt p54) and that ch=cd.
    !----------------------------------------------------------------
    ii=0
    DO j=1,rows
      DO i=1,row_length

        theta1 = theta(i,j,1)
        theta_star = tstar(i,j)*((100000.0/pstar(i,j))**kappa)
        rhostar = pstar(i,j) / ( r*tstar(i,j) )

        ushear = u_p(i,j) - u_0_p(i,j)
        vshear = v_p(i,j) - v_0_p(i,j)
        wshr2 = MAX (1.0e-6 , ushear*ushear + vshear*vshear)
        wshr1 = SQRT(wshr2)
        cd = 1.5 * ( vkman/ALOG(z_full(i,j,1)/z0(i,j)) )**2

        IF (land_mask(i,j)) THEN         ! land
          wthvbar = wshr1 * cd * ( theta_star - theta1 )
        ELSE                             ! sea
          wthvbar = wshr1 * cd *                                       &
                  ( theta_star - theta1 + 0.61*theta1*(qs_star(i,j)-q(i,j,1)) )
        END IF

        ustar = SQRT(cd * wshr2)
        fb_surf(i,j) = g * wthvbar / ( rhostar * theta1*(1.0+0.61*q(i,j,1)) )

        IF (fb_surf(i,j)  >   0.0) THEN
          ii= ii+1
          ! improved method uses BL depth from the previous timestep
          w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
          w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
          tv1_sd(ii) = 1.93 * wthvbar/( rhostar * w_m )           
          bl_vscale2(ii) = 2.5 * 2.52 * w_m * w_m
                 ! 2.5 from Beare (2008), 2.52 to undo the 0.25 factor
          bl_vscale2(ii) = MAX( 0.01, bl_vscale2(ii) ) ! for safety
        END IF

      END DO ! i
    END DO ! j
  END IF ! test on tv1_sd_opt /= 2

ELSE ! if l_flux_bc      (used by some SCM runs)

  !-----------------------------------------------------------
  ! Code for specified surface flux boundary condition.
  ! Assumes a saturated surface (ie. only appropriate for sea)
  !-----------------------------------------------------------

  ii=0
  DO j=1, rows
    DO i=1, row_length
      ! For taylor expansion about T0=SL(K=1)
      tstar(i,j) = theta(i,j,1)*exner_theta_levels(i,j,1)+                  &
                                 gamma_dry*z_full(i,j,1)
    END DO
  END DO

! DEPENDS ON: qsat_mix
  CALL qsat_mix(qs_star,tstar,pstar,row_length*rows,l_mixing_ratio)

  DO j=1, rows
    DO i=1, row_length
      dqsdt(i,j) = (repsilon * lc * qs_star(i,j))                          &
                      / ( r * tstar(i,j) * tstar(i,j) )
    END DO
  END DO

  ii=0
  DO j=1,rows
    DO i=1,row_length

      ushear = u_p(i,j) - u_0_p(i,j)
      vshear = v_p(i,j) - v_0_p(i,j)
      ! Need to have a higher minimum wind speed limit with
      ! specified fluxes in order not to generate huge tstar
      wshr2 = MAX (0.1, ushear*ushear + vshear*vshear)
      wshr1 = SQRT(wshr2)

      ! Calculate wthv from namelist flux_h and flux_e (in W/m2)
      wthvbar = ((flux_h(i,j)/cp)+0.61*(flux_e(i,j)/lc))             &
               * ((100000.0/pstar(i,j))**kappa)

      cd = 1.5 * ( vkman/LOG(z_full(i,j,1)/z0(i,j)) )**2

      theta1 = theta(i,j,1)
      ! Taylor expansion for qsat(T*) about SL(k=1)

      tstar(i,j) = ( theta1 + (wthvbar/(wshr1*cd))                     &
               -   0.61*theta1                                         &
               *   (qs_star(i,j)-q(i,j,1)-dqsdt(i,j)*tstar(i,j)))      &
               /   ( (100000.0/pstar(i,j))**kappa +                    &
               0.61*theta1*dqsdt(i,j) )

      rhostar = pstar(i,j) / ( r*tstar(i,j) )

      ustar = SQRT(cd * wshr2)
      fb_surf(i,j) = g * wthvbar /( rhostar * theta1*(1.0+0.61*q(i,j,1)) )

      IF (fb_surf(i,j)  >   0.0) THEN
        ii= ii+1
        IF (tv1_sd_opt == 0) THEN 
          ! old method assumes the BL depth zh=300.
          tv1_sd(ii) = 1.93 * wthvbar / ( rhostar * ( 75.0 *           &
                            fb_surf(i,j) + ustar*ustar*ustar)**(1.0/3.0) ) 
          w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
          w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
          bl_vscale2(ii) = 2.5 * 2.52 * w_m * w_m
                  ! 2.5 from Beare (2008), 2.52 to undo the 0.25 factor
        ELSE
          ! improved method uses BL depth from the previous timestep
          w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
          w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
          tv1_sd(ii) = 1.93 * wthvbar/( rhostar * w_m )           
          bl_vscale2(ii) = 2.5 * 2.52 * w_m * w_m
                  ! 2.5 from Beare (2008), 2.52 to undo the 0.25 factor
        END IF
        bl_vscale2(ii) = MAX( 0.01, bl_vscale2(ii) ) ! for safety
      END IF! (fb_surf > 0.0)

    END DO ! i
  END DO ! j

END IF!  l_flux_bc

!----------------------------------------------------------------------------

IF (lhook) CALL dr_hook('CONV_SURF_FLUX',zhook_out,zhook_handle)
RETURN
END SUBROUTINE conv_surf_flux

END MODULE conv_surf_flux_mod
