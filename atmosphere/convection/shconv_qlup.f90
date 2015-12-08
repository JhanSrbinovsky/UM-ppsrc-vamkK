! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Liquid water content of the cumulus ensemble, ql_up
!

      SUBROUTINE shconv_qlup( n_sh, max_cldlev                          &
     &,                     mb,  mb_o_wsc, wstar_up, zcld               &
     &,                     wthetal_cb,wqt_cb                           &
     &,                     dthetal_cb,  dqt_cb                         &
     &,                     dthetal_cld, dqt_cld                        &
     &,                     theta_cld, r_cld, exner_cld                 &
     &,                     ftheta_func,f0_func,f1_func                 &
     &,                     ql_up)

!
! Purpose:
!   Calculate the cumulus ensemble liquid water mixing ratio, ql_up
!   after first calculating the cumuls ensemble thetal, thetav and
!   qT. Assumes all moisture variables are in mixing ratio not
!   specific humidity.
!
!   Also using relationship between h'up calculate theta_up & rv_up.
!
!   Called by Shallow_conv (5A version)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!

      USE atmos_constants_mod, ONLY: repsilon, c_virtual, cp, rv
      USE earth_constants_mod, ONLY: g
      
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
     &,  max_cldlev      ! Maximum number of convective cloud levels

      real, intent(in) ::                                               &
     &  mb(n_sh)                                                        &
                                  ! cloud base mass flux
     &, mb_o_wsc(n_sh)                                                  &
                                  ! mb/wstar_up
     &, wstar_up(n_sh)                                                  &
                                  ! wstar_up
     &, zcld(n_sh)                                                      &
                                  ! cloud depth
     &, dthetal_cb(n_sh)                                                &
                                  ! dthetal across cloud base
     &, dqt_cb(n_sh)                                                    &
                                  ! dqt across cloud base
     &, dthetal_cld(n_sh)                                               &
                                  ! dthetal across cloud
     &, dqt_cld(n_sh)                                                   &
                                  ! dqt across cloud
     &, wthetal_cb(n_sh)                                                & 
                                  ! w'thetal' at cloud base
     &, wqt_cb(n_sh)                                                    & 
                                  ! w'qt'  at cloud base

! fields on cloud levels
     &, theta_cld(n_sh,max_cldlev)                                      &
                                     ! incloud theta values (K)
     &, r_cld(n_sh,max_cldlev)                                          &
                                     ! incloud q values (K)
     &, exner_cld(n_sh,max_cldlev)                                      &
                                     ! incloud exner values
     &, ftheta_func(n_sh,max_cldlev)                                    &
                                      ! ftheta function from LES
     &, f0_func(n_sh,max_cldlev)                                        &
                                      ! f0 function from LES
     &, f1_func(n_sh,max_cldlev)      ! f1 function from LES

!
! Arguments with intent OUT:
!

      real, intent(out) ::                                              &
     &  ql_up(n_sh,max_cldlev)    ! cumulus ensemble ql up    (kg/kg)


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

      real ::                                                           &
     &  rv_up(n_sh,max_cldlev)                                          &
                                  ! cumulus ensemble rv up     (kg/kg)
     &, theta_up(n_sh,max_cldlev)                                       &
                                  ! cumulus ensemble theta up     (K)
     &, thetal_up(n_sh,max_cldlev)                                      &
                                       ! thetal  up
     &, thetav_up(n_sh,max_cldlev)                                      &
                                       ! thetav  up
     &, qt_up(n_sh,max_cldlev)         ! qt  up

      real ::                                                           &
     &  thetal_star2(n_sh)                                              &
                               ! thetal star (squared)
     &, qt_star2(n_sh)                                                  &
                               ! qt star (squared)
     &, thetal_cb_up2(n_sh)                                             &
                               ! thetal at cloud base (squared)
     &, qt_cb_up2(n_sh)        ! qt at cloud base (squared)

      real ::                                                           &
     &  thetav_cld(n_sh,max_cldlev)    ! thetav on cloud levels


! temporary variables

      real ::                                                           &
     &  qt_up2                                                          &
     &, thetal_up2                                                      &
     &, virt_term(n_sh)                                                 &
                                     ! virtual_term
     &, d_term                 ! denominator term
      real                                                              &
     & lc_o_cp                 ! lc/cp



! Loop counters
!

      integer :: i,k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
! 1.0 Initialise arrays
!-----------------------------------------------------------------------


      IF (lhook) CALL dr_hook('SHCONV_QLUP',zhook_in,zhook_handle)
      Do k = 1,max_cldlev
        Do i = 1,n_sh

          ql_up(i,k) = 0.0
          rv_up(i,k) = 0.0
          theta_up(i,k) = 0.0

        End do
      End do

! set up constants to be used in routine

      lc_o_cp = lc/cp

!-----------------------------------------------------------------------
! 2.0 calculate terms independent of level
!-----------------------------------------------------------------------

      Do i = 1,n_sh
!
! thetal_star & qt_star  squared
!
        thetal_star2(i) = mb_o_wsc(i)*dthetal_cld(i)*dthetal_cld(i)

        qt_star2(i)     = mb_o_wsc(i)*dqt_cld(i)*dqt_cld(i)


! upward cloud base fluxes

        qt_cb_up2(i) = -0.175 * wqt_cb(i) * dqt_cb(i) / mb(i)

        thetal_cb_up2(i) =-0.175 * wthetal_cb(i) * dthetal_cb(i)/ mb(i)

!
! virtual term
!
        virt_term(i) = sqrt(wstar_up(i)/mb(i))*wstar_up(i)*wstar_up(i)  &
     &                                     /(g*zcld(i))
      End do


!-----------------------------------------------------------------------
!   Loop over levels - calculating various terms requried for ql
!-----------------------------------------------------------------------
! Note at present doing caluations over a lot of points not in cloud
!-----------------------------------------------------------------------

      Do k = 1,max_cldlev
        Do i = 1,n_sh
!
! Calculate virtual potential temperature in cloud
!
          thetav_cld(i,k) = theta_cld(i,k)*                             &
     &               (1.+r_cld(i,k)/repsilon)/(1.+r_cld(i,k))

!
! Ensemble thetal up
!
          thetal_up2     = thetal_cb_up2(i)*f0_func(i,k)                &
     &                      + thetal_star2(i)*f1_func(i,k)

! Take negative root as cumulus updraughts have lower thetal than the
! environment as cloud liquid water forms

! thetal_up2 should be positive but as it is the result of combining
! 2 terms one of which can be negative the result can be negative.

          if (thetal_up2 <  0.0) then
            thetal_up(i,k) = 0.0
          else
            thetal_up(i,k) = -1.*sqrt(thetal_up2)
          endif
!
! Ensemble qt up
!
          qt_up2     = qt_cb_up2(i)*f0_func(i,k)                        &
     &                      + qt_star2(i)*f1_func(i,k)

! Take positive root as cumulus updraughts are greater than mean
! environment
! qt_up2 should be positive but as it is the result of combining
! 2 terms one of which can be negative the result can be negative.

          if (qt_up2 <  0.0) then
            qt_up(i,k) = 0.0
          else
            qt_up(i,k) = sqrt(qt_up2)
          endif

!
! Ensemble virtual potential temperature up
!
          thetav_up(i,k) = thetav_cld(i,k)*virt_term(i)                 &
     &                                    *ftheta_func(i,k)

!
! original calculation of ql_up using thetav_up & qt_up
!

          d_term = ( lc_o_cp/exner_cld(i,k) - theta_cld(i,k)/repsilon )

          ql_up(i,k)=(thetav_up(i,k)-thetal_up(i,k)                     &
     &                   -c_virtual*theta_cld(i,k)*qt_up(i,k))          &
     &                           /d_term

        End do
      End do

!-----------------------------------------------------------------------
!     End Subroutine
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SHCONV_QLUP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE shconv_qlup
