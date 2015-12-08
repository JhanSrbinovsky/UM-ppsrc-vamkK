! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates constants used in large-scale precipitation scheme.
MODULE lspcon_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lspcon( m_ci )


  ! General modules
  USE stochastic_physics_run_mod, ONLY: l_rp2
  USE conversions_mod,            ONLY: pi
  USE earth_constants_mod,        ONLY: g
  USE water_constants_mod,        ONLY: rho_water
  USE missing_data_mod,           ONLY: rmdi

!   Microphysics modules
  USE lsp_dif_mod,         ONLY: air_density0, air_viscosity0,         &
                                 air_conductivity0, air_diffusivity0,  &
                                 air_pressure0, apb1, apb2, apb3, apb4,&
                                 apb5, apb6, tw1, tw2, tw3, tw4, tw5,  &
                                 sc, vent_ice1, vent_ice2, vent_rain1, &
                                 vent_rain2, lc, lf

  USE mphys_psd_mod,       ONLY: l_calcfall, cr, dr, c1r, d1r, h1r,    &
                                 c2r, d2r, h2r, x2i, x3i, x4i, ri, si, &
                                 ci0, di0, x2ic, x3ic, x4ic, ric,      &
                                 sic, cic0, dic0, x1g, x2g, x4g,       &
                                 ag, bg, cg, dg 

  USE mphys_constants_mod, ONLY: cx, constp, m_ci_sav, rho_q_veloc,    &
                                 l_calc_mp_const, x4r, x1i, x1ic,      &
                                 lsp_ei, lsp_fi, lam_evap_enh,         &
                                 n0_murk, m0_murk

  USE mphys_inputs_mod,    ONLY: x1r, x2r, l_psd, ai, bi, aic, bic,    &
                                 lsp_eic, lsp_fic, l_rainfall_as,      &
                                 l_clark_aero, ar, arc, l_mcr_qrain

  USE lsp_autoc_consts_mod, ONLY: n0_clark, n0_haywood, m0_clark,        &
                                  m0_haywood

!   Dr hook modules
  USE yomhook,             ONLY: lhook, dr_hook
  USE parkind1,            ONLY: jprb, jpim

  USE gammaf_mod, ONLY: gammaf
  IMPLICIT NONE

!  Description:
!     Calculates constants used within the LSP_ICE routine.

!  Method:
!     Calculate powers, gamma functions and constants from combinations
!     of physical parameters selected by the comdecks.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

!  Code description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

! Declarations:

  REAL, INTENT(IN) ::  m_ci
                  ! Used to modify ice fall speed

!  Local scalars:
  INTEGER :: i
! Counter to print out the contents of CX and CONSTP
  REAL :: temp,                                                         &
! Forms input to the GAMMAF routine which calculates gamma functions.
       g1, g2, g3,                                                      &
       gb1, gb2, gb3,                                                   &
       gbc1, gbd1, gbdc1,                                               &
       gc1, gc2, gc3,                                                   &
       gd3, gdc3, gd52, gdc52,                                          &
       gdr3, gdr4, gdr52,                                               &
       gd1r3, gd2r3, gd1r4, gd2r4,                                      &
       gr2, gr4, gr5, gr6,                                              &
       gg1, ggb1, ggbd1,                                                &
       gg2, gdg52, gdg3, gg3,                                           &
! Represents the gamma function of BI+DI+1 etc.
! Fall speed of ice particles parameters
       ci,di,cic,dic

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


!- End of header
!---------------------------------------------------------

  IF (lhook) CALL dr_hook('LSPCON',zhook_in,zhook_handle)

! First initialise cx and constp to rmdi where terms are
! not being used to prevent non-initialisation
 cx(:)     = rmdi
 constp(:) = rmdi

!----------------------------------------------------------------------
! Is Generic (Field et al) ice PSD set on?
! If it is, the constants x1i and x1ic will be set to
! 0.00E+0 in the UMUI. To avoid a divide by zero error:
!----------------------------------------------------------------------

  IF (l_psd) THEN
    x1i=2.0e6
    x1ic=40.0e6
  END IF !(L_PSD)

! Do we need to calculate fall speeds?
  IF (l_calcfall) THEN
! Define fall speeds
    ci=lsp_ei*air_viscosity0**(1.0-2.0*lsp_fi)                          &
       *air_density0**(lsp_fi-1.0)                                      &
       *(2.0*g)**lsp_fi*(ai/ri)**lsp_fi
    di=lsp_fi*(bi+2.0-si)-1.0
    cic=lsp_eic*air_viscosity0**(1.0-2.0*lsp_fic)                       &
       *air_density0**(lsp_fic-1.0)                                     &
       *(2.0*g)**lsp_fic*(aic/ric)**lsp_fic
    dic=lsp_fic*(bic+2.0-sic)-1.0
! Modify fallspeeds for random parameters 2
    IF (l_rp2) THEN
      ci=ci*m_ci
      cic=cic*m_ci
    END IF
  ELSE
! Use preset parameters
    ci=ci0
    di=di0
    cic=cic0
    dic=dic0
  END IF  ! Calculation of fall speeds

! CX values. 1-20 are for the crystal population. 21-40 are for the
! aggregate population. 41-62 are for rain. 63-80 are for graupel.
! 81-99 are for the generic ice size distribution.

! Crystals
  cx(1)=(bic+1.0+x4ic-x2ic)/bic
  cx(2)=-(x4ic+1.0-x2ic)/bic
  cx(3)=dic/(bic+1.0+x4ic-x2ic)
  cx(4)=(2.0+x4ic-x2ic)/(bic+1.0+x4ic-x2ic)
  cx(5)=(5.0+dic+2.0*x4ic-2.0*x2ic)*0.5/(bic+1.0+x4ic-x2ic)
  cx(6)=(3.0+dic+x4ic-x2ic)/(bic+1.0+x4ic-x2ic)
  cx(7)=1.0/(x2ic-x4ic-1.0-bic)
  cx(8)=1.0+x4ic
  cx(9)=2.0+x4ic
  cx(10)=3.0+x4ic
  cx(11)=x2ic
  cx(12)=x3ic
  cx(13)=1.0+x4ic+bic
  cx(14)=bic

! Aggregates
  cx(23)=di/(bi+1.0+x4i-x2i)
  cx(24)=(2.0+x4i-x2i)/(bi+1.0+x4i-x2i)
  cx(25)=(5.0+di+2.0*x4i-2.0*x2i)*0.5/(bi+1.0+x4i-x2i)
  cx(26)=(3.0+di+x4i-x2i)/(bi+1.0+x4i-x2i)
  cx(27)=1.0/(x2i-x4i-1.0-bi)
  cx(28)=1.0+x4i
  cx(29)=2.0+x4i
  cx(30)=3.0+x4i
  cx(31)=x2i
  cx(32)=x3i
  cx(33)=3.0+x4i+bi
  cx(34)=2.0+x4i+bi
  cx(35)=1.0+x4i+bi
! Rain
  cx(41)=dr/(4.0+dr-x2r+x4r)
  cx(42)=1.0/(4.0+dr-x2r+x4r)
  cx(43)=6.0+x4r
  cx(44)=5.0+x4r
  cx(45)=4.0+x4r
  cx(46)=x2r
  cx(47)=2.0+x4r-x2r
  cx(48)=1.0/(4.0+x4r+dr-x2r)
  cx(49)=(dr+5.0)*0.5-x2r+x4r
  cx(50)=(3.0+dr-x2r+x4r)/(4.0+dr-x2r+x4r)
! Rain mixing ratio
  cx(51)=dr/(4.0-x2r+x4r)
  cx(52)=1.0/(4.0-x2r+x4r)
  cx(53)=(3.0+dr-x2r+x4r)/(4.0-x2r+x4r)
! Abel & Shipway Terms
  cx(56)=h1r
  cx(57)=h2r
  cx(59)=4.0+d1r+x4r
  cx(60)=4.0+d2r+x4r
  cx(61)=3.0+d1r+x4r
  cx(62)=3.0+d2r+x4r

!Note there is no space between rain and graupel.

! Graupel

  cx(63)=dg/(bg+1.0+x4g-x2g)
  cx(64)=(2.0+x4g-x2g)/(bg+1.0+x4g-x2g)
  cx(65)=(5.0+dg+2.0*x4g-2.0*x2g)*0.5/(bg+1.0+x4g-x2g)
  cx(66)=(3.0+dg+x4g-x2g)/(bg+1.0+x4g-x2g)
  cx(67)=1.0/(x2g-x4g-1.0-bg)
  cx(68)=1.0+x4g
  cx(69)=2.0+x4g
  cx(70)=3.0+x4g
  cx(71)=x2g

  cx(73)=3.0+x4g+bg
  cx(74)=2.0+x4g+bg
  cx(75)=1.0+x4g+bg


! Generic ice particle size distribution
  cx(81)=2.0+di
  cx(82)=di+bi
  cx(83)=1.0+bi
  cx(84)=1.0
  cx(85)=1.0+0.5*(di+1.0)

! Define gamma values
! Crystals
  temp=1.0+x4ic
  CALL gammaf(temp,gc1)
  temp=bic+1.0+x4ic
  CALL gammaf(temp,gbc1)
  temp=bic+dic+1.0+x4ic
  CALL gammaf(temp,gbdc1)
  temp=2.0+x4ic
  CALL gammaf(temp,gc2)
  temp=(dic+5.0+2.0*x4ic)*0.5
  CALL gammaf(temp,gdc52)
  temp=dic+3.0+x4ic
  CALL gammaf(temp,gdc3)
  temp=3.0+x4ic
  CALL gammaf(temp,gc3)
! Aggregates
  temp=1.0+x4i
  CALL gammaf(temp,g1)
  temp=bi+1.0+x4i
  CALL gammaf(temp,gb1)
  temp=bi+2.0+x4i
  CALL gammaf(temp,gb2)
  temp=bi+3.0+x4i
  CALL gammaf(temp,gb3)
  temp=bi+di+1.0+x4i
  CALL gammaf(temp,gbd1)
  temp=2.0+x4i
  CALL gammaf(temp,g2)
  temp=(di+5.0+2.0*x4i)*0.5
  CALL gammaf(temp,gd52)
  temp=di+3.0+x4i
  CALL gammaf(temp,gd3)
  temp=3.0+x4i
  CALL gammaf(temp,g3)
! Rain
  temp=dr+4.0+x4r
  CALL gammaf(temp,gdr4)
  temp=4.0+x4r
  CALL gammaf(temp,gr4)
  temp=6.0+x4r
  CALL gammaf(temp,gr6)
  temp=5.0+x4r
  CALL gammaf(temp,gr5)
  temp=2.0+x4r
  CALL gammaf(temp,gr2)
  temp=(dr+5.0+2.0*x4r)*0.5
  CALL gammaf(temp,gdr52)
  temp=dr+3.0+x4r
  CALL gammaf(temp,gdr3)
!Abel & Shipway Rain Terms
  temp= 3.0+d1r+x4r
  CALL gammaf(temp,gd1r3)
  temp= 3.0+d2r+x4r
  CALL gammaf(temp,gd2r3)
  temp= 4.0+d1r+x4r
  CALL gammaf(temp,gd1r4)
  temp= 4.0+d2r+x4r
  CALL gammaf(temp,gd2r4)

! Graupel
  temp=1.0+x4g
  CALL gammaf(temp,gg1)
  temp=bg+1.0+x4g
  CALL gammaf(temp,ggb1)
  temp=bg+dg+1.0+x4g
  CALL gammaf(temp,ggbd1)
  temp=2.0+x4g
  CALL gammaf(temp,gg2)
  temp=(dg+5.0+2.0*x4g)*0.5
  CALL gammaf(temp,gdg52)
  temp=dg+3.0+x4g
  CALL gammaf(temp,gdg3)
  temp=3.0+x4g
  CALL gammaf(temp,gg3)

! CONSTP values. 1-20 are for the crystal population. 21-40 are for the
! aggregate population. 41-60 are for rain. 61-80 are for graupel.
! 81-99 are for the generic ice size distribution

! Crystals
  constp(1)=x1ic
  constp(2)=1.0/gc1
  constp(3)=1.0/(aic*gbc1)
  constp(4)=cic*gbdc1/gbc1
  constp(5)=1.0/(aic*x1ic*gbc1)
  constp(6)=2.0*pi*x1ic
  constp(7)=vent_ice1*gc2
  constp(8)=vent_ice2*sc**(1.0/3.0)/air_viscosity0**0.5                 &
            *SQRT(cic)*gdc52
  constp(9)=pi*0.25*x1ic*cic*gdc3
  constp(10)=5.0*gc1
  constp(11)=2.0*gc2
  constp(12)=0.25*gc3
  constp(13)=pi**2*rho_water*x1ic*x1r
  constp(14)=2.0*pi*air_conductivity0/lf*x1ic
! Capacitance relative to spheres of same maximum dimension
! Formula depends on value of axial ratio
  constp(15)=arc
  IF (arc  >   1.0) THEN
! Prolate
    constp(15)=(1.0-(1.0/constp(15))**2)**0.5 /                         &
               ALOG( constp(15) + (constp(15)**2-1.0)**0.5 )
  ELSE IF (arc  ==  1.0) THEN
! Spherical
    constp(15)=1.0
  ELSE
! Oblate
    constp(15)=(1.0-constp(15)**2)**0.5                                 &
               /ASIN((1.0-constp(15)**2)**0.5)
  END IF
! Now adjust diffusional growth constants for capacitance
  constp(6)=constp(6)*constp(15)
  constp(14)=constp(14)*constp(15)
!            constp(20) is reserved for lsp_collection

! Values 16 to 23 are unused
! Aggregates

  constp(24)=ci*gbd1/gb1
  constp(25)=1.0/(ai*x1i*gb1)
  constp(26)=2.0*pi*x1i
  constp(27)=vent_ice1*g2
  constp(28)=vent_ice2*sc**(1.0/3.0)/air_viscosity0**0.5                &
             *SQRT(ci)*gd52
  constp(29)=pi*0.25*x1i*ci*gd3
  constp(30)=5.0*g1
  constp(31)=2.0*g2
  constp(32)=0.25*g3
  constp(33)=pi**2*rho_water*x1i*x1r
  constp(34)=2.0*pi*air_conductivity0/lf*x1i
! Capacitance relative to spheres of same maximum dimension
! Formula depends on value of axial ratio
  constp(35)=ar
  IF (ar  >   1.0) THEN
! Prolate
    constp(35)=(1.0-(1.0/constp(35))**2)**0.5 /                         &
               ALOG( constp(35) + (constp(35)**2-1.0)**0.5 )
  ELSE IF (ar  ==  1.0) THEN
! Spherical
    constp(35)=1.0
  ELSE
! Oblate
    constp(35)=(1.0-constp(35)**2)**0.5                                 &
               / ASIN((1.0-constp(35)**2)**0.5)
  END IF
! Now adjust diffusional growth constants for capacitance
  constp(26)=constp(26)*constp(35)
  constp(34)=constp(34)*constp(35)

  constp(36) = gc1*gb3
  constp(37) = 2.0*gc2*gb2
  constp(38) = gc3*gb1
  constp(39) = ai*x1ic*x1i*pi/4.0
!            constp(40) is reserved for lsp_collection
! Rain
  constp(41)=6.0*cr*gdr4/gr4
  constp(42)=pi*rho_water/6.0*x1r*gdr4*cr
  constp(43)=1.0/120.0*gr6
  constp(44)=1.0/24.0*gr5
  constp(45)=1.0/6.0*gr4
  constp(46)=2.0*pi*x1r
  constp(47)=vent_rain1*gr2
  constp(48)=vent_rain2*sc**(1.0/3.0)/air_viscosity0**0.5               &
             *gdr52*SQRT(cr)
  constp(49)=pi*0.25*x1r*cr*gdr3
  constp(50)=1.0/(pi*rho_water*x1r*gr4/6.0)

! Abel & Shipway Rain Section
!-----------------------------------------------------------
!First overwrite lam_evap_enh if maximum rain rate specified
!-----------------------------------------------------------

  constp(51)=c1r*x1r*pi*0.25*gd1r3
  constp(52)=c2r*x1r*pi*0.25*gd2r3
  constp(53)=gr4
  constp(54)=c1r*gd1r4
  constp(55)=c2r*gd2r4
  constp(56)=(pi/6.0)*rho_water*x1r*(lam_evap_enh**x2r)
  constp(57)=(pi/6.0)*rho_water*x1r


! Values 58 to 60 are unused

! Graupel
  constp(64) = cg*ggbd1/ggb1
  constp(65) = 1.0/(ag*x1g*ggb1)
  constp(67) = vent_rain1*gg2
  constp(68) = vent_rain2*sc**(1.0/3.0)/air_viscosity0**0.5             &
               *SQRT(cg)*gdg52
  constp(69) = pi*0.25*x1g*cg*gdg3
  constp(74) = 2.0*pi*air_conductivity0*x1g/lf
  constp(76) = gg1*gb3
  constp(77) = 2.0*gg2*gb2
  constp(78) = gg3*gb1
  constp(79) = ai*x1g*x1i*pi/4.0
  constp(80) = pi*pi/24.0*x1g*ag

! Generic ice particle size distribution
  constp(81)=pi*0.25*ci
  constp(82)=ci*ai
  constp(83)=2.0*pi*constp(35)
  constp(84)=vent_ice1
  constp(85)=vent_ice2*sc**(1.0/3.0)/air_viscosity0**0.5                &
             *SQRT(ci)
  constp(86)=pi*pi/24.0*rho_water*x1r
  constp(87)=gr6
  constp(88)=2.0*gr5
  constp(89)=gr4
  constp(90)=2.0*pi*constp(35)*air_conductivity0/lf
  constp(91)=gb3
  constp(92)=2.0*gb2
  constp(93)=gb1

!-----------------------------------------------------------------------
! Determine Abel & Shipway Prognostic Fall Velocity at the point where
! it diverges from the traditional UM parametrization
! Note: already have lambda from module
!-----------------------------------------------------------------------
! Determine rhoqv = droplet velocity multiplied by any given rho and q.
! rho_q_veloc is a constant and
! droplet velocity = rho_q_veloc * air density correction / (rho q)

  IF (l_rainfall_as .AND. .NOT. l_mcr_qrain ) THEN

    rho_q_veloc = constp(56) *  (                                       &
    ( constp(54) / ((lam_evap_enh+cx(56))**cx(59) ) ) +                 &
    ( constp(55) / ((lam_evap_enh+cx(57))**cx(60) ) )    )

  END IF

!-----------------------------------------------------------------
! Determine aerosol scheme to use for
!-----------------------------------------------------------------

  IF (l_clark_aero) THEN

    ! Use Clark et al (2008, QJRMS) constants
    n0_murk = n0_clark
    m0_murk = m0_clark

  ELSE ! l_clark_aero

    ! Use Haywood et al (2008, QJRMS) constants
    n0_murk = n0_haywood
    m0_murk = m0_haywood

  END IF ! l_clark_aero

!-----------------------------------------------------------------
! Constants are now calculated, no need to recalculate
  l_calc_mp_const = .FALSE.

!-----------------------------------------------------------------
! Update m_ci_sav to the new value defined by m_ci

  m_ci_sav = m_ci
!-----------------------------------------------------------------

  IF (lhook) CALL dr_hook('LSPCON',zhook_out,zhook_handle)
  RETURN
! End the subroutine
END SUBROUTINE lspcon

END MODULE lspcon_mod
