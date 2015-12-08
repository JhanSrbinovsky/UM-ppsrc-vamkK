! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Shallow convection scheme - wql & dqsat/dT in cloud
!

      SUBROUTINE SHALLOW_WQL(n_sh, nlev, max_cldlev                     &
     &,                            ncld_thlev                           &
     &,                            zcld, wstar_up, root_mb_o_wsc        &
     &,                            wql_cb, wql_inv                      &
       ! fields on cloud levels
     &,                           qsat_cld, q_cld, theta_cld            &
     &,                           t_cld, fql_func                       &
       ! In/Output fields
     &,                           wql_cld                               &
       ! Output fields
     &,                       dqsatdt_cld)



!
! Purpose:
!   Shallow convection scheme - routine to calculate the liquid water
!   budget.
!   Make use of the mass flux and qlup calculated earlier.
!
!   Called by SHALLOW_CONV
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
      USE earth_constants_mod, ONLY: g

      USE atmos_constants_mod, ONLY: repsilon, rv

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
     &  n_sh                                                            &
                         ! No. of shallow convection points
     &, nlev                                                            &
                         ! No. of model layers
     &, max_cldlev       ! maximum noumber of in cloud levels

      integer, intent(in) ::                                            &
     &  ncld_thlev(n_sh) ! number of theta cloud levels

      real, intent(in)    ::                                            &
     &  zcld(n_sh)                                                      &
                               ! Depth of cloud layer (m)
     &, wstar_up(n_sh)                                                  &
                               ! w* convective velocity scale (wsc)
     &, wql_cb(n_sh)                                                    &
                               ! flux of wql across cloud base
     &, wql_inv(n_sh)                                                   &
                               ! flux of wql across inversion
     &, root_mb_o_wsc(n_sh)     ! sqrt of (mb/wsc)


      real, intent(in) ::                                               &
     &  theta_cld(n_sh,nlev)                                            &
                              ! theta
     &, q_cld(n_sh,nlev)                                                &
                              ! q - mixing ratio
     &, qsat_cld(n_sh,nlev)                                             &
                              ! qsat in cloud
     &, t_cld(n_sh,nlev)                                                &
                              ! temperature in cloud
     &, fql_func(n_sh,nlev)   ! fql function


!
! Arguments with intent INOUT:
!
      real, intent(inout) ::                                            &
     &  wql_cld(n_sh,nlev)     ! wql in cloud (includes cloud base)


!
! Arguments with intent OUT:
!
      real, intent(out) ::                                              &
     &  dqsatdt_cld(n_sh,max_cldlev+1)    ! dqsat/dt

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

      real ::                                                           &
     &  thetav_cld(n_sh,max_cldlev+1)                                   &
     &, temp(n_sh)

!
! Loop counters
!

      integer :: i,k,itop

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! 1.0  Calculate dqsat/dT  in cloud & store ratio w*/zcld
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SHALLOW_WQL',zhook_in,zhook_handle)

!
!   dqsatdt on cloud levels  (expression for mixing ratio)
!
      Do k = 1,max_cldlev+1
        Do i = 1,n_sh
          dqsatdt_cld(i,k) = qsat_cld(i,k)*lc                           &
     &                     /(Rv*t_cld(i,k)*t_cld(i,k))
          thetav_cld(i,k)=theta_cld(i,k)*(1.0+q_cld(i,k)/repsilon)      &
     &                 /(1.+q_cld(i,k))
        End do
      End do
!-----------------------------------------------------------------------
! 2.0 Calculate  in cloud values of wql on uv levels using
!      similarity expression
!-----------------------------------------------------------------------
! values for new water flux parametrisation

        k=1
        Do i=1,n_sh
          wql_cld(i,k) = wql_cb(i)
          temp(i) =root_mb_o_wsc(i)*(wstar_up(i)**3)/(zcld(i)*g)
        End do

! taking thetav/theta=1.

        Do k=2,max_cldlev+1
          Do i=1,n_sh

            wql_cld(i,k) = 0.5*temp(i)*fql_func(i,k)*theta_cld(i,k)     &
     &                                        /thetav_cld(i,k)

          End do
        End do

! value at inversion

        Do i=1,n_sh
          itop=ncld_thlev(i)+1
          wql_cld(i,itop) = wql_inv(i)
        End do

!-----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('SHALLOW_WQL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SHALLOW_WQL
