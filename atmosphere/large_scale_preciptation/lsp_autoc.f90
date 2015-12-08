! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Autoconversion of liquid to rain
! Subroutine Interface:
MODULE lsp_autoc_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_autoc(                                                   &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  qcl, qrain, t, p,                                                     &
                                          ! Water contents, temp and p
  cfliq, rhcpt,                                                         &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
!    &, cf, cfl                           ! Current cloud fractions for
!                                         ! updating
  area_liq, area_mix, area_ice,                                         &
                                          ! Cloud fraction overlaps
  rainfrac, rain_liq, rain_mix,                                         &
                                          ! Rain fractions for updating
  rain_ice, rain_clear,                                                 &
  rho, rhor, corr2,                                                     &
                                          ! Parametrization information
  lcrcp,                                                                &
                                          ! Microphysical information
  one_over_tsi,                                                         &
                                          ! Number of iterations and
                                          ! 1/(timestep*iterations)
  ptransfer, rftransfer,                                                &
                                          ! Mass and rainfrac transfer
 
!    &, cftransfer,cfltransfer            ! Cloud transfer diagnostics
  land_fract,                                                           &
                                          ! Land surface information
  n_drop_tpr, n_drop_out,                                               &
                                          ! Droplet concs.
  r_theta_levels_c, fv_cos_theta_latitude_c                             &
  )

  ! Microphysics modules
  USE lsp_autoc_consts_mod, ONLY: inhomog_rate, inhomog_lim, r_thresh,  &
                                  r_auto, n_auto, power_droplet_auto,   &
                                  power_qcl_auto, power_rho_auto,       &
                                  consts_auto, aut_pref, aut_qc, aut_nc

  USE mphys_inputs_mod,     ONLY: l_autoc_3b, l_autolim_3b, l_warm_new

  USE mphys_constants_mod,  ONLY: ec_auto, l_inhomog
 
  USE mphys_ice_mod,        ONLY: qcfmin

  USE mphys_bypass_mod,     ONLY: mp_dell, mp_delp


  ! General/atmosphere modules   
  USE atmos_constants_mod,  ONLY: r, repsilon
  USE water_constants_mod,  ONLY: lc
  USE conversions_mod,      ONLY: pi
  USE um_input_control_mod, ONLY: l_mr_physics1, l_auto_debias,         &
                                  l_use_sulphate_autoconv,              &
                                  l_use_bmass_autoconv,                 &
                                  l_use_ocff_autoconv,                  &
                                  l_use_seasalt_autoconv,               &
                                  l_use_nitrate_autoconv,               &
                                  l_use_biogenic
  USE murk_inputs_mod,      ONLY: l_murk

  ! Dr Hook modules
  USE yomhook,              ONLY: lhook, dr_hook
  USE parkind1,             ONLY: jprb, jpim

  IMPLICIT NONE

! Purpose:
!   Update cloud and rain due to autoconversion of liquid to rain

! Method:
!   Use a rate based on a power law of liquid content and droplet
!   number with a minimum liquid content for autoconversion.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

!   There are many routes through this code depending upon the
!   selected options. Vector options have been removed since
!   they are very difficult to trace through.
!   1) Calculate number of droplets based upon either a
!      prescribed number or the amount of aerosol (murk or
!      sulphate/sea-salt/biomass/fossil-fuel organic carbon/
!      nitrate) in the model.
!   2) Calculate the autoconversion rate based upon a
!      Tripoli and Cotton power law formulation.
!   3) Calculate the autoconversion limit based either on
!      the concentration of droplets above 20 um diameter or
!      on a Tripoli and Cotton type argument

! Subroutine Arguments


  INTEGER, INTENT(IN) ::                                                &
    points
                          ! Number of points to calculate

  REAL, INTENT(IN) ::                                                   &
    timestep,                                                           &
                          ! Timestep / s
    p(points),                                                          &
                          ! Air pressure / N m-2
    t(points),                                                          &
                          ! Temperature / K
    cfliq(points),                                                      &
                          ! Fraction of gridbox with liquid cloud
!    &, cf(points)        ! Bulk cloud fraction for updating
!    &, cfl(points)       ! Liquid cloud fraction for updating
    area_liq(points),                                                   &
                          ! Fraction of gridbox with liquid but no ice
    area_mix(points),                                                   &
                          ! Fraction of gridbox with liquid and ice
    area_ice(points),                                                   &
                          ! Fraction of gridbox with ice but no liquid
    rhcpt(points),                                                      &
                          ! Rhcrit on each point (no units)
    rho(points),                                                        &
                          ! Air density / kg m-3
    rhor(points),                                                       &
                          ! 1 / air density / m3 kg-1
    corr2(points),                                                      &
                          ! Temperature correction factor (no units)
    lcrcp,                                                              &
                          ! Latent heat of condensation/cP / K
    land_fract(points),                                                 &
                          ! Fraction of gridbox with land
    n_drop_tpr(points),                                                 &
                          ! Droplet number determined by
                          ! droplet taper curves
    r_theta_levels_c(points),                                           &
                          ! Distance from centre of Earth and ... 
    fv_cos_theta_latitude_c(points),                                    &
                          ! ... grid info for working out gridbox size.
    one_over_tsi
                          ! 1/(timestep*iterations)

  REAL, INTENT(INOUT) ::                                                &
    qcl(points),                                                        &
                          ! Liquid water content / kg kg-1
    qrain(points),                                                      &
                          ! Rain water content / kg kg-1
    ptransfer(points),                                                  &
                          ! Autoconversion rate / kg kg-1 s-1
!    &, cftransfer(points) ! Rate of change of bulk cloud frac / s-1
!    &, cfltransfer(points)! Rate of change of liquid cloud frac / s-1
    rftransfer(points),                                                 &
                          ! Rate of change of rain fraction / s-1
    rainfrac(points),                                                   &
                          ! Rain fraction (no units)
    rain_liq(points),                                                   &
                          ! Overlap of rain with liquid (no units)
    rain_mix(points),                                                   &
                          ! Overlap of rain with mixed phase region
    rain_ice(points),                                                   &
                          ! Overlap of rain with ice
    rain_clear(points),                                                 &
                          ! Overlap of rain with clear sky
    n_drop_out(points)    
                          ! Droplet number from autoconversion

!-----------------------------
! Local Variables

  INTEGER ::                                                            &
    kk,                                                                 &
                          ! Index for condensed points
    i
                          ! Loop counter

  REAL ::                                                               &
    n_drop(points),                                                     &
                          ! Number of droplets / m-3
    dpr(points),                                                        &
                          ! Amount of mass autoconverted / kg kg-1
    qsl_tl(points),                                                     &
                          ! Saturated humidity wrt liquid at
                          !  temperature TL / kg kg-1
    t_l(points),                                                        &
                          ! Liquid temperature / K
    r_mean(points),                                                     &
                          ! Mean radius of droplets / m
    r_mean0,                                                            &
                          ! Factor in calculation of r_mean
    a_factor(points),                                                   &
                          ! 3 r_mean / m
    b_factor(points),                                                   &
                          ! 0.5 n_drop / m-3
    n_gt_20(points),                                                    &
                          ! Concentration of droplets with
                          !  diameter gt 20 um / m-3
    autorate(points),                                                   &
                          ! Rate constants for autoc.
    autolim(points),                                                    &
                          ! Autoconversion threshold / kg kg-1
    qc(points),                                                         &
                          ! Maximum amout of liquid to remove / kg kg-1
!    &, delta_cf(points)  ! Cloud fraction transferred (no units)
!    &, delata_cfl(points)! Liquid cloud fraction transferred (no units)
    rainfracnew(points),                                                &
                          ! Updated rain fraction
    rho_1(points)
                          ! rho^(power_rho_auto - power_qcl_auto + 1)

  REAL ::                                                               &
                          ! For debiasing code
    alpha_l,                                                            &
                          ! dqsat/dT at T_L / kg kg-1 K-1
    a_l,                                                                &
                          ! 1 / (1 + L/cp alpha)
    sigma_s,                                                            &
                          ! Width of moisture PDF
    g_l,                                                                &
                          ! Factor in debiasing calculation
    gacb,                                                               &
                          ! Factor in debiasing calculation
    ac_factor,                                                          &
                          ! Multiplying factor for debiasing
    ltiny
                          ! Largest real number represented on platform
                       
  REAL ::                                                               &
    fsd,                                                                &
            ! fractional standard dev of in-cloud liquid water content
    bias,                                                               &
            ! autoconversion rate bias
    x_in_km ! grid-box size in km

!  Function calls
  REAL, EXTERNAL :: number_droplet

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

      !-----------------------------------------------
      ! Part 1. Calculate droplet number concentration
      !-----------------------------------------------

  IF (lhook) CALL dr_hook('LSP_AUTOC',zhook_in,zhook_handle)

  ltiny = LOG(TINY(1.0))      ! Initialize ltiny

  kk = 0 ! Reset index counter for condensed points

  ! Use pre-determined droplet number concentration from
  ! lsp_taper_ndrop.F90 subroutine.

  DO i = 1, points

    IF (qcl(i) > 0.0  .AND. cfliq(i) > 0.0) THEN

      kk = kk + 1

      n_drop(kk)    = n_drop_tpr(i)
      n_drop_out(i) = n_drop(kk)

    ELSE

      n_drop_out(i) = 0.0

    END IF

  END DO

      !-----------------------------------------------
      ! Part 2: Calculating the autoconversion rate and limit
      !-----------------------------------------------
  IF (l_warm_new) THEN

     !------------------------------------------------
     ! Use new autoconversion scheme
     !------------------------------------------------
    kk=0

    DO i = 1, points

      IF (qcl(i) > 0.0 .AND. cfliq(i) > 0.0) THEN

        ! Autoconversion active for this point
        kk=kk+1

        IF (n_drop(kk) > 0.0) THEN

          ! calculate the in-cloud autoconversion rate
          autorate(i) = aut_pref *                                      &
                        ((qcl(i)/cfliq(i)) ** aut_qc) *                 &
                        ((1.0E-6*n_drop(kk)) ** aut_nc)

        ELSE ! n_drop=0

          ! No autoconversion
          autorate(i) = 0.0

        END IF ! n_drop > 0

        IF (l_inhomog) THEN

          ! calculate bias factor E based on Boutle et al 2012 QJ
          x_in_km = 0.001*SQRT (   r_theta_levels_c(i) * mp_dell        &
                                 * r_theta_levels_c(i) * mp_delp        &
                                 * fv_cos_theta_latitude_c(i)     )

          IF ( cfliq(i) < 1.0 ) THEN
            fsd = (0.45-0.25*cfliq(i))*(((x_in_km*cfliq(i))**0.333)     &
                  *((0.06*x_in_km*cfliq(i))**1.5+1.0)**(-0.17))
          ELSE
            fsd = 0.11*(((x_in_km*cfliq(i))**0.333)                     &
                  *((0.06*x_in_km*cfliq(i))**1.5+1.0)**(-0.17))
          END IF
          fsd=fsd*1.414
          bias = ((1+fsd**2)**(-0.5*aut_qc))*                           &
                 ((1+fsd**2)**(0.5*aut_qc**2))

        ELSE ! no inhomog param

          bias = 1.0

        END IF !l_inhomog

        dpr(i) = autorate(i) * bias * timestep * cfliq(i)
        dpr(i) = MAX(MIN(dpr(i),qcl(i)-qcfmin),0.0)

          !-----------------------------------------------
          ! Update liquid water content and rain
          !-----------------------------------------------
        qcl(i)   = qcl(i)   - dpr(i)
        qrain(i) = qrain(i) + dpr(i)

          !-----------------------------------------------
          ! Store process rate
          !-----------------------------------------------
        ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

        IF (dpr(i) >  0.0) THEN
            !-----------------------------------------------
            ! Calculate change in rain fraction
            !-----------------------------------------------
          rainfracnew(i) = MAX(rainfrac(i),cfliq(i))
          rftransfer(i) = rftransfer(i)                                 &
          + (rainfracnew(i) - rainfrac(i)) * one_over_tsi

            !-----------------------------------------------
            ! Update rain fractions
            !-----------------------------------------------
          rainfrac(i)   = rainfracnew(i)
          rain_liq(i)   = MIN(area_liq(i),rainfrac(i))
          rain_mix(i)   = MIN(area_mix(i),rainfrac(i)-rain_liq(i))
          rain_ice(i)   =                                               &
              MIN(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
          rain_clear(i) =                                               &
              rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i)

        END IF  ! dpr gt 0

      END IF ! qcl and cfl > 0

    END DO ! points

  ELSE ! original autoconversion etc

    IF (l_autoc_3b) THEN

        !-----------------------------------------------
        ! Use the original 3B method for autoconversion
        !-----------------------------------------------
      kk=0

      DO i = 1, points

          ! The numbers 5907.24 and 4188.79 represent combinations
          ! of physical constants. Do NOT change them.
          ! Please refer to UMDP26
        IF (qcl(i) >  0.0.AND.cfliq(i) >  0.0) THEN

          kk=kk+1
          IF (n_drop(kk) > 0.0) THEN
            autorate(i) = 5907.24*ec_auto*inhomog_rate                  &
                               *n_drop(kk)**(-1.0/3.0)
            autolim(i)  = 4188.79*r_thresh**3*n_drop(kk)*inhomog_lim
          ELSE
              ! No droplets so no autoconversion
            autorate(i) = 0.0
            autolim(i)  = 0.0
          END IF  ! n_drop >0

        END IF  ! qcl > 0

      END DO  ! points

    ELSE  ! l_auto_3b
        !-----------------------------------------------
        ! Use the later 3C/3D method for autoconversion
        !-----------------------------------------------

      kk=0
        ! Calculate multiplying factor for droplet size
      r_mean0=(27.0/(80.0*pi*1000.0))**(-1.0/3.0)

      DO i = 1, points
        IF (qcl(i) >  0.0.AND.cfliq(i) >  0.0) THEN

            ! Autoconversion is active
          kk=kk+1

            ! Calculate inverse of mean droplet size
          r_mean(kk) = r_mean0 * (rho(i) * qcl(i)                       &
                     / (cfliq(i) * n_drop(kk))) ** (-1.0/3.0)

            ! Calculate numerical factors
          b_factor(kk) = 3.0 * r_mean(kk)
          a_factor(kk) = n_drop(kk) * 0.5

            !-----------------------------------------------
            ! Calculate droplet number concentration greater
            ! than threshold radius
            !-----------------------------------------------
       
          
! Only compute Exponent if it will not generate an inexact signal.
! Required for initial ENDGame runs on IBM.

          IF ( (b_factor(kk) * r_auto) < -ltiny) THEN
            
            n_gt_20(i) = (b_factor(kk)**2 * a_factor(kk)) * r_auto**2   &
                                      * EXP(-b_factor(kk) * r_auto)     &
                       + (2.0 * b_factor(kk) * a_factor(kk)) * r_auto   &
                                      * EXP(-b_factor(kk) * r_auto)     &
                       + (2.0 * a_factor(kk))                           &
                                      * EXP(-b_factor(kk) * r_auto)

          ELSE

            n_gt_20(i)=0.0
       
          END IF       

            !-----------------------------------------------
            ! Test to see if there is a sufficient concentration of
            ! droplets with diameters > threshold for autoconversion
            ! to proceed.
            !-----------------------------------------------
          IF (n_gt_20(i)  >=  n_auto) THEN
              ! Calculate autoconversion rate
            autorate(i) = consts_auto * ec_auto                         &
                      * n_drop(kk) ** power_droplet_auto
          ELSE
              ! No autoconversion
            autorate(i)=0.0
          END IF ! n_gt_20 ge n_auto

            !-----------------------------------------------
            ! Calculate the autoconversion limit
            !-----------------------------------------------
            ! Calculate value of local qcl at which the droplet
            ! concentration with radii greater than 20um will
            ! fall below a threshold (1000 m-3). This is a
            ! hardwired numerical approximation
          autolim(i)=(6.20e-31*n_drop(kk)**3)-(5.53e-22*n_drop(kk)**2)  &
                 +(4.54e-13*n_drop(kk))+(3.71e-6)-(7.59/n_drop(kk))

        END IF  ! qcl(i) >  0.0
      END DO  ! points

    END IF ! l_autoc_3b

    IF (l_autolim_3b) THEN

        !-----------------------------------------------
        ! Overwrite the autoconversion limit with 3B values
        !-----------------------------------------------
      DO i = 1, points
        IF (land_fract(i)  >=  0.5) THEN
          autolim(i)=8.621e-4
        ELSE
          autolim(i)=2.155e-4
        END IF  ! land_fract ge 0.5
      END DO

    END IF ! l_autolim3b

      !-----------------------------------------------
      ! Part 3: Optionally debias the autoconversion rate
      !-----------------------------------------------
    IF (l_auto_debias) THEN

        !-----------------------------------------------
        ! Calculate qsat(TL) in order to allow subgrid
        ! debiasing of the autoconversion rate.
        !-----------------------------------------------
      DO i = 1, points
        t_l(i) = t(i) - (lcrcp * qcl(i) )
      END DO
! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix( qsl_tl, t_l, p, points, l_mr_physics1 )

      kk=0
      DO i = 1, points
        IF (qcl(i) >  0.0.AND.cfliq(i) >  0.0) THEN

            ! Autoconversion is active
          kk=kk+1

            !-----------------------------------------------
            ! De-bias the autoconversion rate based on a
            ! triangular Smith PDF, following Wood et al
            ! (Atmos. Res., 65, 109-128, 2002).
            !-----------------------------------------------

          IF (qcl(i)  <   1.0e-15) THEN
              ! Water contents are small so set factor to 1
            ac_factor=1.0

          ELSE  ! qcl lt 1e-15
            alpha_l = repsilon * lc * qsl_tl(i) / ( r * t_l(i)**2 )
            a_l = 1.0 / (1.0+(lcrcp*alpha_l))
            sigma_s = (1.0 - rhcpt(i)) * a_l * qsl_tl(i) / SQRT(6.0)
            g_l = 1.15 * (power_qcl_auto-1.0) * sigma_s
            gacb = EXP(-1.0 * qcl(i) / g_l)

              ! Calculate autoconversion rate multiplication factor
            ac_factor = MAX(                                            &
                    (cfliq(i)**(power_qcl_auto-1.0))*1.0/(1.0-gacb),    &
                     1.0)

          END IF  ! qcl lt 1e-15

            !-----------------------------------------------
            ! Apply the debiasing factor
            !-----------------------------------------------
          autorate(i) = ac_factor * autorate(i)

        END IF ! qcl > 0

      END DO  ! points

    END IF  ! l_auto_debias

      !-----------------------------------------------
      ! Part 4. Calculate the autoconversion
      !-----------------------------------------------
    DO i = 1, points

        !-----------------------------------------------
        ! Set the dependence of autoconversion on air density
        !-----------------------------------------------
        ! power_rho_auto and power_qcl_auto are set in c_lspmic.
      rho_1(i) = rho(i) ** (power_rho_auto-power_qcl_auto+1.0)

      IF (qcl(i) >  0.0.AND.cfliq(i) >  0.0) THEN

            ! Calculate maximum amount of liquid that can be removed
            ! from the grid box
        qc(i) = MIN( autolim(i) * cfliq(i) * rhor(i) , qcl(i) )

            !-----------------------------------------------
            ! Calculate the autoconversion amount
            !-----------------------------------------------
        dpr(i) = MIN(autorate(i)                                        &
                  *(rho(i)*qcl(i)/cfliq(i))**(power_qcl_auto-1.0)       &
                  *timestep*qcl(i)*rho_1(i)/corr2(i),qcl(i)-qc(i))

          !-----------------------------------------------
          ! Update liquid water content and rain
          !-----------------------------------------------
        qcl(i)   = qcl(i)   - dpr(i)
        qrain(i) = qrain(i) + dpr(i)

          !-----------------------------------------------
          ! Store process rate
          !-----------------------------------------------
        ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

        IF (dpr(i) >  0.0) THEN
            !-----------------------------------------------
            ! Calculate change in rain fraction
            !-----------------------------------------------
          rainfracnew(i) = MAX(rainfrac(i),cfliq(i))
          rftransfer(i) = rftransfer(i)                                 &
          + (rainfracnew(i) - rainfrac(i)) * one_over_tsi

            !-----------------------------------------------
            ! Update rain fractions
            !-----------------------------------------------
          rainfrac(i)   = rainfracnew(i)
          rain_liq(i)   = MIN(area_liq(i),rainfrac(i))
          rain_mix(i)   = MIN(area_mix(i),rainfrac(i)-rain_liq(i))
          rain_ice(i)   =                                               &
              MIN(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
          rain_clear(i) =                                               &
              rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i)

        END IF  ! dpr gt 0

          !-----------------------------------------------
          ! Update cloud fractions
          !-----------------------------------------------
!            These are commented out since there is currently no
!            cloud fraction update associated with the autoconversion.
!              cf(i)  = cf(i) + delta_cf(i)
!              cfl(i) = cfl(i)+ delta_cfl(i)
!            cftransfer(i)  = cftransfer(i)  + delta_cf(i)
!     &                                      / (timestep*iterations)
!            cfltransfer(i) = cfltransfer(i) + delta_cfl(i)
!     &                                      / (timestep*iterations)

      END IF  ! qcl gt 0
    END DO  ! points

  END IF ! l_warm_new

  IF (lhook) CALL dr_hook('LSP_AUTOC',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_autoc
END MODULE lsp_autoc_mod
