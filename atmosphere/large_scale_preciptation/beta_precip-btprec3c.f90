! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE BETA_PRECIP --------------------------------------------
!     PURPOSE:
! Process fields of precipitation intensity to give scattering coefft
! in 1/metres.
! Calculated at model level (eg bottom eta level 25m)
! or level within surface layer eg screen ht ( 1.5M )
!  Programming standard: U M Doc. Paper No. 4
!
!  External documentation
!    Forecasting Research Scientific Paper NO.4
!    Diagnosis of visibility in the UK Met Office Mesoscale Model
!    and the use of a visibility analysis to constrain initial
!    conditions.  SP Ballard, BJ Wright, BW Golding    1992
!      NIMROD diagnostic:
!      Wright, B. J., 1997: Improvements to the Nimrod Visibility
!         Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!      Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!         for Nimrod. Met. Office FR Tech Rep., No. 222.
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Large Scale Precipitation
MODULE beta_precip_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE beta_precip                                                  &
           (ls_rain, ls_snow, c_rain, c_snow, qcf,                      &
                                                      !INPUT
            rho, t, pressure,                                           &
                                                      !INPUT
            lca,cca,pct,avg,                                            &
                                                      !INPUT
            p_field,points,k1stpt,                                      &
                                                      !INPUT
            beta_ls_rain, beta_ls_snow,                                 &
                                                      !OUTPUT
            beta_c_rain, beta_c_snow,error)     !OUTPUT

  ! General atmosphere modules
  USE atmos_constants_mod, ONLY: pref
  USE water_constants_mod, ONLY: rho_water
  USE conversions_mod,     ONLY: pi

  ! Microphysics modules
  USE mphys_psd_mod,       ONLY: ri, si, cr, dr, x2i, x4i, ci0, di0
  USE mphys_inputs_mod,    ONLY: x1r, x2r, ai,  bi, l_psd
  USE mphys_constants_mod, ONLY: x1i, x1ic, x4r

  ! Dr Hook Modules
  USE yomhook,             ONLY: lhook, dr_hook
  USE parkind1,            ONLY: jprb, jpim

  ! Large scale precipitation modules
  USE gammaf_mod,          ONLY: gammaf
  USE lsp_moments_mod,     ONLY: lsp_moments

  IMPLICIT NONE
!---------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
  INTEGER ::                                                            &
   p_field,                                                             &
                                        ! IN NO. points in field.
   points,                                                              &
                ! IN Number of gridpoints being processed.
   k1stpt,                                                              &
                ! IN First gridpoint processed within complete field.
   error    ! OUT Error code

  REAL ::                                                               &
   ls_rain(p_field),                                                    &
                                        ! IN Large scale Rain
   ls_snow(p_field),                                                    &
                                        ! IN Large scale Snow
   c_rain(p_field),                                                     &
                                        ! IN Convective Rain
   c_snow(p_field),                                                     &
                                        ! IN Convective Snow
   qcf(p_field),                                                        &
                                        ! IN large-scale ice / kg kg-1
   rho(p_field),                                                        &
                                        ! IN Air density / kg m-3
   t(p_field),                                                          &
                                        ! IN Temperature / K
   pressure(p_field),                                                   &
                                        ! IN Pressure
   lca(p_field),                                                        &
                                        ! IN Total Layer Cloud.
   cca(p_field)                         ! IN Convective Cloud.


  LOGICAL ::                                                            &
   pct,                                                                 &
                                        ! IN T:Cloud amounts are in %
   avg
                                        ! IN T:Precip =local*prob

!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
  REAL ::                                                               &
   beta_ls_rain(p_field),                                               &
                                        ! OUT Scattering in LS Rain.
   beta_ls_snow(p_field),                                               &
                                        ! OUT Scattering in LS Snow.
   beta_c_rain(p_field),                                                &
                                        ! OUT Scattering in Conv Rain
   beta_c_snow(p_field)                 ! OUT Scattering in Conv Snow
!---------------------------------------------------------------------

  REAL ::                                                               &
   powerr,                                                              &
   poweri,                                                              &
   factorr1,                                                            &
   factorr2,                                                            &
   factori1,                                                            &
   factori2
  REAL ::                                                               &
   inst_ls_rain(p_field),                                               &
                                          ! Local Large scale Rain
   inst_ls_snow(p_field),                                               &
                                          ! Local Large scale Snow
   inst_c_rain(p_field),                                                &
                                          ! Local Convective Rain
   inst_c_snow(p_field),                                                &
                                          ! Local Convective Snow
   lcai(p_field),                                                       &
                                          ! 1 / large-scale cloud (lca)
   m_si(p_field)
                          ! si moment of ice particle size distribution

! Local varables:-----------------------------------------------------
  REAL ::                                                               &
   pfactor(p_field),                                                    &
   gammar,                                                              &
   gammar1,                                                             &
   gammai,                                                              &
   gammai1,                                                             &
   smallvalue
  PARAMETER(smallvalue=1.0e-7)

!-----------------------------------------------------------------------
!  Define local variables ----------------------------------------------
  INTEGER :: i       ! Loop counters: I - horizontal field index;

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('BETA_PRECIP',zhook_in,zhook_handle)
  error=0
  IF((k1stpt+points-1) >  p_field)THEN
    error=1
    GO TO 9999
  END IF

  DO i = k1stpt, k1stpt+points-1
    pfactor(i)=(pref/pressure(i))**0.4
  END DO

  powerr=(3.0-x2r+x4r)/(dr+4.0-x2r+x4r)
  poweri=(1.0+si-x2i+x4i)/(bi+di0+1.0-x2i+x4i)
  factorr1=0.5*pi*x1r
  factorr2=pi/6.0*rho_water*cr*x1r
  factori1=2.0*ri*x1i
  factori2=ai*ci0*x1i

  CALL gammaf(x4r+dr+4.0,gammar)
  CALL gammaf(x4r+3.0,gammar1)
  CALL gammaf(x4i+bi+di0+1.0,gammai)
  CALL gammaf(x4i+si+1.0,gammai1)

  IF (l_psd) THEN
        ! Find inverse of layer cloud amount
    DO i = 1, p_field
      lcai(i)=1.0/MAX(lca(i),0.01)
    END DO

        ! Use the generic ice particle size distribution to
        ! calculate the si moment of the ice particle size distribution
    CALL lsp_moments(p_field,rho,t,qcf,lcai,si,m_si )

  END IF  ! l_psd

  IF (avg) THEN

    IF (pct) THEN

      DO i = k1stpt, k1stpt+points-1

        IF(lca(i)  >   smallvalue .AND. cca(i)  <   100.0) THEN
          inst_ls_rain(i) = 10000.0*ls_rain(i)/(100.0-cca(i))/lca(i)
          inst_ls_snow(i) = 10000.0*ls_snow(i)/(100.0-cca(i))/lca(i)
        ELSE
          inst_ls_rain(i) = 0.0
          inst_ls_snow(i) = 0.0
        END IF
        IF(cca(i)  >   smallvalue) THEN
          inst_c_rain(i) = 100.0*c_rain(i)/cca(i)
          inst_c_snow(i) = 100.0*c_snow(i)/cca(i)
        ELSE
          inst_c_rain(i) = 0.0
          inst_c_snow(i) = 0.0
        END IF

      END DO

    ELSE

      DO i = k1stpt, k1stpt+points-1

        IF(lca(i)  >   smallvalue .AND. cca(i)  <   1.0) THEN
          inst_ls_rain(i) = ls_rain(i)/(1.0-cca(i))/lca(i)
          inst_ls_snow(i) = ls_snow(i)/(1.0-cca(i))/lca(i)
        ELSE
          inst_ls_rain(i) = 0.0
          inst_ls_snow(i) = 0.0
        END IF
        IF(cca(i)  >   smallvalue) THEN
          inst_c_rain(i) = c_rain(i)/cca(i)
          inst_c_snow(i) = c_snow(i)/cca(i)
        ELSE
          inst_c_rain(i) = 0.0
          inst_c_snow(i) = 0.0
        END IF

      END DO

    END IF

  ELSE

    DO i = k1stpt, k1stpt+points-1
      inst_ls_rain(i) = ls_rain(i)
      inst_ls_snow(i) = ls_snow(i)
      inst_c_rain(i)  = c_rain(i)
      inst_c_snow(i)  = c_snow(i)
    END DO

  END IF

  DO i = k1stpt, k1stpt+points-1

    IF(inst_ls_rain(i)  >   smallvalue) THEN
      beta_ls_rain(i)=factorr1*gammar1*(inst_ls_rain(i)/                &
                      (pfactor(i)*factorr2*gammar))**powerr
    ELSE
      beta_ls_rain(i)=0.0
    END IF

    IF(inst_ls_snow(i)  >   smallvalue) THEN
      IF (l_psd) THEN
            ! Use the generic ice particle size distribution
        beta_ls_snow(i)=2.0*ri*m_si(i)
      ELSE
        beta_ls_snow(i)=factori1*gammai1*(inst_ls_snow(i)/              &
                        (pfactor(i)*factori2*gammai))**poweri
      END IF  ! l_psd
    ELSE
      beta_ls_snow(i)=0.0
    END IF

    IF(inst_c_rain(i)  >   smallvalue) THEN
      beta_c_rain(i) =factorr1*gammar1*(inst_c_rain(i) /                &
                      (pfactor(i)*factorr2*gammar))**powerr
    ELSE
      beta_c_rain(i) =0.0
    END IF

    IF(inst_c_snow(i)  >   smallvalue) THEN
          ! Use the 3C size distribution of snow for the convective
          ! contribution since we do not have the equivalent of qcf
          ! easily available.
      beta_c_snow(i) =factori1*gammai1*(inst_c_snow(i) /                &
                      (pfactor(i)*factori2*gammai))**poweri
    ELSE
      beta_c_snow(i) =0.0
    END IF

  END DO

  9999 CONTINUE

  IF (lhook) CALL dr_hook('BETA_PRECIP',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE beta_precip
END MODULE beta_precip_mod
