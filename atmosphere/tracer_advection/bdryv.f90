! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  FUNCTION BDRYV----------------------------------------------------
!LL
!LL  Purpose: Special routine to add psuedo source terms to boundary
!LL           data in limited area.
!LL  PARAMETERS ARE SPECIFIC TO UK MESOSCALE MODEL
!LL  Method:  The boundary concentrations are computed using a
!LL           simple model of transport from sources outside the
!LL           model. Analysis of the source distribution outside
!LL           the UK MES shows that it can be well represented by
!LL           a line source at constant radius from the centre of
!LL           the model, with a source distribution given by the
!LL           sum of two Gaussians. Concentrations from these are
!LL           computed assuming transport using the local windspeed u
!LL           or 1 m/s, whichever is stronger, over a distance
!LL           determined from the centroid of the source distribution, x
!LL           with a linear transformation rate k from emission to
!LL           aerosol, dry deposition at a rate determined from the
!LL           dry deposition velocity vd and mean mixed layer depth h.
!LL           Thus the max concentration is given by
!LL               Q/(uh)*k/(k+vd/h)*(1-exp(-k*x/u))
!LL           The source term is assumed to decrease with level
!LL           pressure. See forthcoming documentation for details.
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3,
!LL                        Version 7, dated 11/3/93.
!LL
!LL
!*L  Arguments:---------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Tracer Advection
      REAL FUNCTION BDRYV(                                              &
     & WDIR                                                             &
     &,WSPEED                                                           &
     &,PRESS                                                            &
     &)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      REAL                                                              &
     & WDIR                                                             &
                     ! IN Wind direction : Cartesian degrees
     &,WSPEED                                                           &
                     ! IN Wind speed m/s
     &,PRESS         ! IN Pressure

!*
!* Local, including SAVE'd, storage------------------------------------
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL                                                              &
     & ZANGLE1,WIDTH1                                                   &
                      ! Centre and width of first source Gaussian.
     &,ZANGLE2,WIDTH2                                                   &
                      ! Centre and width of second source Gaussian.
     &,CMAX,CZER                                                        & 
                 ! Max concentration and 'background'.
     &,WDIRN,RECIPROOT2PI                                               &
     &,MIXD,TRAVEL                                                      &
                      ! Average mixed layer depth and travel distance.
     &,QMAX1                                                            &
                      ! Peak height of first source Gaussian.
     &,QMAX2                                                            &
                      ! Peak height of second source Gaussian.
     &,VD                                                               &
                      ! Dry deposition velocity.
     &,K,K1                                                             &
            ! Transformation parameters.
     &,PH                                                               &
          ! Pressure height scale.
     &,KRAT,KT

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      PARAMETER(ZANGLE1=178.0)
      PARAMETER(WIDTH1=5.0)
      PARAMETER(ZANGLE2=173.0)
      PARAMETER(WIDTH2=25.0)
      PARAMETER(RECIPROOT2PI=0.3989422803)
      PARAMETER(QMAX1=4.7E5,QMAX2=9.0E4,MIXD=800.0,TRAVEL=7.7E5)
      PARAMETER(CZER=6.0)
      PARAMETER(VD=5.0E-3)
      PARAMETER(K=3.0E-6, K1=K+VD/MIXD, KRAT=K/K1,KT=-K1*TRAVEL)
      PARAMETER(PH=3.0E4)
!
!     Max concentration = Q/(uh)*k/(k+vd/h)*(1-exp(-k*x/u))
!
      IF (lhook) CALL dr_hook('BDRYV',zhook_in,zhook_handle)
      CMAX=1.0/MAX(WSPEED,1.0)/MIXD
      CMAX=CMAX*KRAT*(1-EXP(KT/MAX(WSPEED,1.0)))
      WDIRN = WDIR - ZANGLE1
      IF (WDIRN  <   -180.0) WDIRN=WDIRN+360.0
      IF (WDIRN  >    180.0) WDIRN=WDIRN-360.0
      WDIRN=WDIRN/WIDTH1
      BDRYV= QMAX1 * EXP(-WDIRN*WDIRN/2.0)
      WDIRN = WDIR - ZANGLE2
      IF (WDIRN  <   -180.0) WDIRN=WDIRN+360.0
      IF (WDIRN  >    180.0) WDIRN=WDIRN-360.0
      WDIRN=WDIRN/WIDTH2
      BDRYV= BDRYV + QMAX2 * EXP(-WDIRN*WDIRN/2.0)
!
!     Add 'background' value.
!
      BDRYV= BDRYV * CMAX + CZER
!
!     Reduce concentration with pressure altitude.
!
      BDRYV= BDRYV * EXP(-(1.E5-PRESS)/PH)
      IF (lhook) CALL dr_hook('BDRYV',zhook_out,zhook_handle)
      RETURN
      END FUNCTION BDRYV
