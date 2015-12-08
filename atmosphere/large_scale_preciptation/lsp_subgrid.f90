! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Subgrid-scale set ups and checks
MODULE lsp_subgrid_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_subgrid(                                                 &
  points,                                                               &
                                          ! Number of points
  q, qcf_cry, qcf_agg, qcftot, t,                                       &
                                          ! Water contents and temp
  qsl, qs,                                                              &
                                          ! Saturated water contents
  q_ice, q_clear, q_ice_1, q_ice_2,                                     &
                                          ! Local vapour contents
  area_liq,area_mix,area_ice,area_clear,                                &
                                              ! Cloud frac partitions
  area_ice_1, area_ice_2,                                               &
                                          ! Subdivision of area_ice
  areamix_over_cfliq,                                                   &
                                          ! area_mix/cfliq
  rain_liq,rain_mix,rain_ice,rain_clear,                                &
                                              ! Rain overlap partitions
  cftot, cfliq, cfice, cficei,                                          &
                                          ! Cloud fractions for
  frac_ice_above,                                                       &
                                          ! partition calculations
  cf, cff, rainfrac,                                                    &
                                          ! Cloud and rain fractions
                                          ! for updating
  lsrcp,                                                                &
                                          ! Latent heat of sublim./cp
  rhcpt                                                                 &
                                          ! RH crit values
  )

  ! Microphysics modules
  USE mphys_ice_mod,     ONLY: qcfmin

  ! Cloud modules
  USE cloud_inputs_mod,  ONLY: ice_width

  ! General atmosphere modules
  USE conversions_mod,   ONLY: zerodegc

  ! Dr Hook Modules
  USE yomhook,           ONLY: lhook, dr_hook
  USE parkind1,          ONLY: jprb, jpim

  IMPLICIT NONE

! Purpose:
!   Perform the subgrid-scale setting up calculations

! Method:
!   Parametrizes the width of the vapour distribution in the part
!   of the gridbox which does not have liquid water present.
!   Calculates the overlaps within each gridbox between  the cloud
!   fraction prognostics and rainfraction diagnostic.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! The subgrid calculations are a necessary step to calculating the
! subsequent deposition and sublimation transfers and for setting up
! the partition information that is used by the subsequent transfers.

! Subroutine Arguments

  INTEGER, INTENT(IN) ::                                                &
    points            ! Number of points to calculate

  REAL, INTENT(IN) ::                                                   &
    qs(points),                                                         &
                          ! Saturated mixing ratio wrt ice
    qsl(points),                                                        &
                          ! Saturated mixing ratio wrt liquid
    cfliq(points),                                                      &
                          ! Fraction of gridbox with liquid cloud
    rainfrac(points),                                                   &
                          ! Fraction of gridbox containing rain
    frac_ice_above(points),                                             &
                               ! Ice cloud in level above this one
    lsrcp,                                                              &
                          ! Latent heat of sublimation
                          ! / heat capacity of air / K
    rhcpt(points)     ! RH crit values

  REAL, INTENT(INOUT) ::                                                &
    q(points),                                                          &
                          ! Vapour content / kg kg-1
    qcf_cry(points),                                                    &
                          ! Ice crystal content / kg kg-1
    qcf_agg(points),                                                    &
                          ! Ice aggregate content / kg kg-1
    qcftot(points),                                                     &
                          ! Total ice content before advection / kg kg-1
    t(points),                                                          &
                          ! Temperature / K
    cf(points),                                                         &
                          ! Current cloud fraction
    cff(points)       ! Current ice cloud fraction

  REAL, INTENT(OUT) ::                                                  &
    q_clear(points),                                                    &
                          ! Local vapour in clear-sky region / kg kg-1
    q_ice(points),                                                      &
                          ! Local vapour in ice-only region  / kg kg-1
    q_ice_1(points),                                                    &
                          ! Local vapour in ice-only regions that are:
    q_ice_2(points),                                                    &
                          !   1, depositing; and 2, subliming.
    cftot(points),                                                      &
                          ! Modified cloud fraction for partition calc.
    cfice(points),                                                      &
                          ! Modified ice cloud frac. for partition calc.
    cficei(points),                                                     &
                          ! 1/cfice
    area_liq(points),                                                   &
                          ! Frac of gridbox with liquid cloud but no ice
    area_mix(points),                                                   &
                          ! Frac of gridbox with liquid and ice cloud
    area_ice(points),                                                   &
                          ! Frac of gridbox with ice cloud but no liquid
    area_clear(points),                                                 &
                          ! Frac of gridbox with no cloud
    area_ice_1(points),                                                 &
                          ! Frac of gridbox where ice-only cloud is:
    area_ice_2(points),                                                 &
                          !  1, depositing; and 2, subliming.
    areamix_over_cfliq(points),                                         &
                          ! area_mix/cfliq
    rain_liq(points),                                                   &
                          ! Frac of gbox with rain and liquid but no ice
    rain_mix(points),                                                   &
                          ! Frac of gbox with rain and liquid and ice
    rain_ice(points),                                                   &
                          ! Frac of gbox with rain and ice but no liquid
    rain_clear(points)! Frac of gbox with rain but no condensate

! Local Variables

  INTEGER ::                                                            &
    i                 ! Loop counter

  REAL ::                                                               &
    tempw(points),                                                      &
                          ! Vapour content in ice and clear partitions
    temp7(points),                                                      &
                          ! Temporary in width of PDF calculation
    width(points)     ! Full width of vapour distribution in ice and
                          ! clear sky.

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('LSP_SUBGRID',zhook_in,zhook_handle)
  DO i = 1, points

        !-----------------------------------------------
        ! Check that ice cloud fraction is sensible.
        !-----------------------------------------------
        ! Difference between the way PC2 and non-PC2 code operates
        ! is kept here in order to be tracable across model versions.
        ! However, perhaps the code ought to be the same.
          ! 0.001 is to avoid divide by zero problems

      cfice(i)  = MAX( cff(i), 0.001 )
      cficei(i) = 1.0/cfice(i)
      cftot(i)  = cf(i)
      cftot(i)  = MIN( MAX( cftot(i), cfice(i) ),( cfice(i) + cfliq(i)) )

        ! -----------------------------------------------
        ! Calculate overlaps of liquid, ice and rain fractions
        ! -----------------------------------------------
    area_liq(i) = MAX(cftot(i)-cfice(i),0.0)
    area_mix(i) = MAX(cfice(i)+cfliq(i)-cftot(i),0.0)

    IF (cfliq(i) /= 0.0) THEN
      areamix_over_cfliq(i) = area_mix(i)/cfliq(i)
    END IF

    IF (cfice(i) == cftot(i)) THEN
      area_mix(i)           = cfliq(i)
      areamix_over_cfliq(i) = 1.0
    END IF

        ! Remove tiny mixed phase areas
    IF (area_mix(i) < 0.0001) THEN
      area_mix(i)           = 0.0
      areamix_over_cfliq(i) = 0.0
    END IF

    area_ice(i) = MAX(cftot(i)-cfliq(i),0.0)
    area_clear(i) = MAX(1.0-cftot(i),0.0)
    rain_liq(i) = MAX(MIN(area_liq(i),rainfrac(i)),0.0)
    rain_mix(i) = MAX(MIN(area_mix(i),rainfrac(i)-rain_liq(i)),0.0)
    rain_ice(i) =                                                       &
      MAX(MIN(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i)),0.0)
    rain_clear(i) =                                                     &
      MAX(rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i),0.0)

    IF (cfliq(i)  <   1.0) THEN

          ! -----------------------------------------------
          ! Calculate width of vapour dist. in ice and clear region
          ! -----------------------------------------------
          ! tempw is the mean vapour content in the ice only and clear
          ! sky partitions
      tempw(i) = (q(i) - cfliq(i)*qsl(i)) / (1.0 - cfliq(i))
      temp7(i) = ice_width * qsl(i)
          ! 0.001 is to avoid divide by zero problems
      width(i) = 2.0 *(1.0-rhcpt(i))*qsl(i)                             &
                     *MAX(  (1.0-0.5*qcftot(i)/temp7(i)),0.001  )
          ! The full width cannot be greater than 2q because otherwise
          ! part of the gridbox would have negative q. Also ensure that
          ! the full width is not zero (possible if rhcpt is 1).
      width(i) = MIN(width(i) , MAX(2.0*q(i),0.001*qs(i) )    )

          ! -----------------------------------------------
          ! Calculate vapour contents in ice only and clear regions
          ! -----------------------------------------------
      IF (area_ice(i)  >   0.0) THEN
        q_clear(i) = tempw(i) - 0.5*width(i) * area_ice(i)
        q_ice(i) = (q(i)-cfliq(i)*qsl(i)-area_clear(i)*q_clear(i))      &
                 / area_ice(i)
      ELSE
        q_clear(i) = tempw(i)
        q_ice(i) = 0.0               ! q_ice is a dummy value here
      END IF  ! area_ice gt 0

    ELSE ! cf_liq lt 1

          ! -----------------------------------------------
          ! Specify dummy values for q_clear and q_ice
          ! -----------------------------------------------
      width(i)   = 1.0
      q_clear(i) = 0.0
      q_ice(i)   = 0.0

    END IF ! cf_liq lt 1

        ! -------------------------------------------------
        ! Remove any small amount of ice to be tidy.
        ! -------------------------------------------------
        ! If QCF is less than QCFMIN and isn't growing by deposition
        ! (assumed to be given by RHCPT) then evaporate it.
    IF ((qcf_cry(i)+qcf_agg(i)) <  qcfmin) THEN
      IF (t(i) >  zerodegc .OR.                                         &
         (q_ice(i)  <=  qs(i) .AND. area_mix(i)  <=  0.0)               &
         .OR. (qcf_cry(i)+qcf_agg(i)) <  0.0)  THEN
        q(i) = q(i) +qcf_cry(i)+qcf_agg(i)
        t(i) = t(i) - lsrcp * (qcf_cry(i)+qcf_agg(i))
        qcf_cry(i)=0.0
        qcf_agg(i)=0.0
      END IF ! T gt 0 etc.
    END IF ! qcf_cry+qcf_agg lt qcfmin

        ! -------------------------------------------------
        ! First estimate of partition sizes for ice sublimation
        ! and deposition and vapour contents within these partitions
        ! -------------------------------------------------
    IF (q_ice(i)  >   qs(i)) THEN
          ! First estimate is to use a deposition process
      area_ice_1(i) = area_ice(i)
      area_ice_2(i) = 0.0
      q_ice_1(i)    = q_ice(i)
      q_ice_2(i)    = qs(i)       ! Dummy value

    ELSE ! q_ice gt qs
          ! First estimate is to use a sublimation process
      area_ice_1(i) = 0.0
      area_ice_2(i) = area_ice(i)
      q_ice_1(i)    = qs(i)       ! Dummy value
      q_ice_2(i)    = q_ice(i)

    END IF ! q_ice gt qs

        ! -------------------------------------------------
        ! Detailed estimate of partition sizes for ice sublimation
        ! and deposition and vapour contents within these partitions
        ! -------------------------------------------------
    IF (area_ice(i)  >   0.0) THEN
        ! Temp7 is the estimate of the proportion of the gridbox
        ! which contains ice and has local q > than qs (wrt ice)
      temp7(i) = 0.5*area_ice(i) + (q_ice(i)-qs(i)) / width(i)

      IF (temp7(i) >  0.0 .AND. temp7(i) <  area_ice(i)) THEN
            ! Calculate sizes of regions and q in each region
            ! These overwrite previous estimates
        area_ice_1(i) = temp7(i)
        area_ice_2(i) = area_ice(i) - area_ice_1(i)
        q_ice_1(i) = qs(i) + 0.5 * area_ice_1(i) * width(i)
        q_ice_2(i) = qs(i) - 0.5 * area_ice_2(i) * width(i)
      END IF ! temp7 gt 0 etc.

    END IF ! area_ice gt 0

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_SUBGRID',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_subgrid
END MODULE lsp_subgrid_mod
