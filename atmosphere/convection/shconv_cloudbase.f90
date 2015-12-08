! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Cloud base fluxes
!

      SUBROUTINE shconv_cloudbase( n_sh, ntra, l_tracer, mb             &
     &,                     theta_plus,r_plus,rsat_plus,exner_plus      &
     &,                     theta_minus,r_minus,exner_minus             &
     &,                     exner_cb, wthetav_surf, dtracer_cb          &
     &,                     dthetal_cb, drt_cb,dthetav_cb               &
     &,                     wtheta_plus_cb, wthetal_cb,wq_plus_cb       &
     &,                     wqt_cb,wql_cb                               &
     &,                     wthetav_cb, wtracer_cb)

!
! Purpose:
!   Calculate the fluxes at cloud base
!
!   Called by Shallow_conv (5A version)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!

      USE atmos_constants_mod, ONLY: cp, c_virtual, rv, recip_epsilon

      USE water_constants_mod, ONLY: lc
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!

      integer, intent(in) ::                                            &
     &   n_sh                                                           &
                         ! No. of shallow convection points
     &,  ntra            ! No. of tracers

      logical, intent(in) ::                                            &
     &   l_tracer        ! true for tracers

      real, intent(in) ::                                               &
     &  mb(n_sh)                                                        &
                               ! cloud base mass flux
     &, theta_plus(n_sh)                                                &
                               ! theta at + side of cloud base (K)
     &, theta_minus(n_sh)                                               &
                               ! theta at - side of cloud base (k)
     &, r_plus(n_sh)                                                    &
                               ! water vapour at + side of cloud base
                               ! (mixing ratio kg/kg)
     &, r_minus(n_sh)                                                   &
                               ! water vapour at - side of cloud base
                               ! (mixing ratio kg/kg)
     &, rsat_plus(n_sh)                                                 &
                               ! qsat at + side of cloud base
                               ! (mixing ratio kg/kg)
     &, exner_plus(n_sh)                                                &
                               ! exner at + side of cloud base
     &, exner_minus(n_sh)                                               &
                               ! exner at - side of cloud base
     &, exner_cb(n_sh)                                                  &
                               ! exner at cloud base
     &, wthetav_surf(n_sh)                                              & 
                               ! w'thetav' at surface
     &, dtracer_cb(n_sh,ntra)  ! change in tracers across cloud base
                               ! (kg/kg)
!
! Arguments with intent INOUT:
!
!          NONE

!
! Arguments with intent OUT:
!

      real, intent(out) ::                                              &
     &  wthetal_cb(n_sh)                                                & 
                                  !  w'thetal' at cloud base
     &, wqt_cb(n_sh)                                                    & 
                                  !  w'qt'  at cloud base
     &, wql_cb(n_sh)                                                    & 
                                  !  w'ql' at cloud base
     &, wq_plus_cb(n_sh)                                                & 
                                  !  w'q' at cloud base
     &, wtheta_plus_cb(n_sh)                                            & 
                                  !  w'theta' at cloud base
     &, wthetav_cb(n_sh)                                                & 
                                  !  w'thetav' at cloud base
     &, dthetal_cb(n_sh)                                                &
                               ! dthetal across cloud base
     &, dthetav_cb(n_sh)                                                &
                               ! dthetav across cloud base
     &, drt_cb(n_sh)                                                    &
                               ! drt across cloud base
     &, wtracer_cb(n_sh,ntra)  ! w'tracer' at cloud base

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
      real                                                              &
     & lc_o_cp                 ! lc/cp


      real                                                              &
     &  T_plus(n_sh)                                                    &
                               ! T at + side of cloud base
     &, T_minus(n_sh)                                                   &
                               ! T at - side of cloud base
     &, thetav_plus(n_sh)                                               &
                               ! thetav at + side of cloud base
     &, thetav_minus(n_sh)     ! thetav at - side of cloud base

! temporary variables
      real                                                              &
        term_a                                                          &
      , term_c                                                          &
      , term_d                                                          &
      , r_rsat_term                                                     &
                         ! r-rsat   term
      , drsatdt_topt     ! drsat/dT at the top of the transition region

!
! Loop counters
!

      integer :: i, ktra

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! 1.0 Initialise arrays
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SHCONV_CLOUDBASE',zhook_in,zhook_handle)

! set up constants to be used in routine

      lc_o_cp = lc/cp

!-----------------------------------------------------------------------
! calculate terms independent of level
!-----------------------------------------------------------------------

      Do i = 1,n_sh

        t_plus(i) = theta_plus(i)*exner_plus(i)
        t_minus(i) = theta_minus(i)*exner_minus(i)

! Calculation for mixing ratio

        thetav_plus(i)  = theta_plus(i) *                               &
                         (1.+r_plus(i)*recip_epsilon)/(1.+r_plus(i))
        thetav_minus(i) = theta_minus(i)*                               &
                         (1.+r_minus(i)*recip_epsilon)/(1.+r_minus(i))
! values across cloud base

        dthetav_cb(i) = thetav_plus(i) - thetav_minus(i)


! ----------------------------------------------------------------------
! New bit trying to interpolate to actual location of top of transition
! zone.
! Assume minus theta level is the bottom of the transition zone
!
! (a) assume gradient in theta & q roughly constant through transition
!      zone
! (b) May try later d/dpstar method of interpolation.
! ----------------------------------------------------------------------

!
! The convective parametrisation scheme assumes that the enviromental
! air has qcl=0 and qcf=0. i.e thetal=theta and rt=r (this may not be
! the case in the model).

        dthetal_cb(i) = theta_plus(i)-theta_minus(i)

        drt_cb(i)     = r_plus(i) - r_minus(i)

!
! expression for dqsat/dT
!

        drsatdt_topt = Lc *rsat_plus(i)                                 &
     &                                 /(Rv * t_plus(i)*t_plus(i))

        r_rsat_term = r_plus(i) - rsat_plus(i)

!
! assume w'ql'+ is value at cloud base
!

        wql_cb(i) = 0.4*wthetav_surf(i)/(lc_o_cp/exner_plus(i) -        &
     &                        (1.+c_virtual)*theta_plus(i))


        term_a = -0.7*mb(i)*dthetav_cb(i)

        term_c = wql_cb(i)*(1.+lc_o_cp*drsatdt_topt)/0.4

        term_d = mb(i)*r_rsat_term

        wqt_cb(i) = (term_c - term_d                                    &
     &           + exner_plus(i)*drsatdt_topt*term_a )/                 &
     &     (1.+c_virtual*exner_plus(i)*theta_minus(i)*drsatdt_topt)

        wthetal_cb(i) = term_a                                          &
     &                     - c_virtual*theta_minus(i)*wqt_cb(i)

        wthetav_cb(i) = term_a

! this then implies as far as in cloud in concerned

        wq_plus_cb(i) = wqt_cb(i) - wql_cb(i)

        wtheta_plus_cb(i) = wthetal_cb(i) + lc_o_cp*wql_cb(i)           &
     &                                                /exner_cb(i)


      End do


! Note expect wthetal & wtheta to be negative at cloud base and wqt
! wq and wql to be positive.

!-----------------------------------------------------------------------
! tracers - assumed to behave as other conserved variables

      If (l_tracer) then
        Do ktra = 1,ntra
          Do i = 1,n_sh
            wtracer_cb(i,ktra) = -0.7*mb(i)*dtracer_cb(i,ktra)
          End do
        End do
      End If


!-----------------------------------------------------------------------
!    End Subroutine
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SHCONV_CLOUDBASE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE shconv_cloudbase
