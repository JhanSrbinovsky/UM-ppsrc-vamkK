! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  New Shallow convection scheme - Calculations of LES functions
!

      SUBROUTINE les_shall_func_thlev(n_sh,max_cldlev,nlev              &
     &,     eta,wsc_o_mb                                                &
     &,     int_type                                                    &
     &,     fw_func, g_func,k_func, f0_func, f1_func                    &
     &,     ftheta_func)

!
! Purpose:
!   Calculate various functions of height based on fits to LES data
!   for use in the new turbulence base shallow convection scheme.
!   Done for all cloud levels, rest of array set to zero.
!
!
!   Called by SHALLOW_CONV5A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3  programming standards
!


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
     &, max_cldlev                                                      &
                      ! Maximum number of convective cloud levels
     &, nlev                                                            &
                      ! Maximum number of levels
     &, int_type      ! method to use for integrating
                      ! 1 - trapezium rule (only available method).

      real, intent(in) ::                                               &
     &  eta(n_sh,nlev)                                                  &
                       ! non-dimension height of the cloud levels
     &, wsc_o_mb(n_sh)     ! wstar_up over cloud base mass flux

!
! Arguments with intent INOUT:
!
!     None
!
! Arguments with intent OUT:
!
      real, intent(out) ::                                              &
     &  g_func(n_sh,nlev)                                               &
                             ! Used in turbulent transport
     &, k_func(n_sh,nlev)                                               &
                            ! Used in eddy visocisty for CMT
     &, fw_func(n_sh,nlev)                                              &
                             ! Used in ensemble vertical velocity
     &, f0_func(n_sh,nlev)                                              &
                             ! Used in thermodynamic properties
     &, f1_func(n_sh,nlev)                                              &
                             ! Used in thermodynamic properties
     &, ftheta_func(n_sh,nlev) ! used in virtual potential temp

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
      integer :: rlev              ! Number of levels for integration
                                   ! needs to be higher than cloud
                                   ! levels and fixed given resolution
                                   ! of model
      parameter (rlev = 30)    ! ?


      integer :: i,j,k                                                  &
                                   ! loop counters
     &,  ibelow                                                         &
     &, k_minus, k_centre, k_plus


! variables used to precalculate some of the constants in the
! functions

      real ::                                                           &
     & gg1, gg2, cc1, dfw_at_eta, inc_eta

! variables used to store exponentails calculated
      real ::                                                           &
     &  wexp6(n_sh,max_cldlev)                                          &
                                 ! exp(-6eta)
     &, wexp10(n_sh,max_cldlev)  ! exp(-10eta)

! variables used in integral calculation

      real ::                                                           &
     &   temp                                                           &
                          ! temporary store
     &,  eta_reg(0:rlev)                                                &
                          ! eta on regular grid
     &,  fw(0:2*rlev)                                                   &
                          ! fw on regular grid
     &,  dfw_deta(0:rlev)                                               &
                          ! dfw/deta on regular grid
     &,  fob(n_sh)                                                      &
     &,  fob_fot(n_sh)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!
!-----------------------------------------------------------------------
! 1.0 Initialise arrays - set all functions to zero
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('LES_SHALL_FUNC_THLEV',zhook_in,zhook_handle)

!
! initialise fob and fot
!
      Do i = 1,n_sh
        fob(i) = 0.1*wsc_o_mb(i)
!       fot(i) = 2. - fob(i)        ! for information only
        fob_fot(i) = 0.2*wsc_o_mb(i) - 2.
      End Do

!-----------------------------------------------------------------------
! 2.0 Functions not defined as integrals i.e. can be evaluated
!     directly
!-----------------------------------------------------------------------
!
!  K   = 0.115( 1+eta)(1-exp(-10eta))
!
!  Fng = 1 -0.3(1+1.75eta)(1-exp(-10eta))
!
!  B   = 1.4(1+eta-0.6(eta**2))(1-exp(-10eta))
!
!  G   = 1.45(eta - 0.665*eta*eta+ (4.67-1.33*6eta)exp(-6eta)/36.) -C
!    (c = 1.45*4.67/36.)
!
!  f0  = 1. -[0.473(fob-fot)+fob]eta
!                             + (0.473/0.64)*(fob-fot)(1-exp(-0.64eta))
!
!{ df0/deta = 0.473(fob -fot)exp(-0.64eta)-0.473(fob -fot)-fob}
!
!
!  f1  = 9*(eta**1.5)
!
!  ftheta = 2.25 -1.75*exp(-4eta)
!
!-----------------------------------------------------------------------
!
      gg1 = 4.67/36.
      gg2 = 1.33/6.
      cc1 = 0.473/0.64

      Do k =1,max_cldlev
        Do i = 1,n_sh
          If (eta(i,k) <= 1.0) then    ! in cloud
            wexp10(i,k) = 1. - exp(-10.*eta(i,k))
            wexp6(i,k)  = exp(-6.*eta(i,k))

            g_func(i,k)   = 1.45*(eta(i,k) -0.665*eta(i,k)*eta(i,k)     &
     &                    + (gg1 - gg2*eta(i,k))*wexp6(i,k) -gg1)

             k_func(i,k)   = 0.115*(1. + eta(i,k))*wexp10(i,k)

            f0_func(i,k)  = 1. -(0.473*fob_fot(i)+fob(i))*eta(i,k)      &
     &                   + cc1*fob_fot(i)*(1.-exp(-0.64*eta(i,k)))

            f1_func(i,k)  = 9.*(eta(i,k)**1.5)

            ftheta_func(i,k)  = 2.25 - 1.75*exp(-4.*eta(i,k))

          Else        ! non cloud levels

            k_func(i,k)   = 0.0
            g_func(i,k)   = 0.0
            f0_func(i,k)  = 0.0
            f1_func(i,k)  = 0.0
            ftheta_func(i,k)  = 0.0

          End if
        End Do
      End Do

!-----------------------------------------------------------------------
! 3.0 Functions defined as differentials with no analytic integral
!-----------------------------------------------------------------------
!
! dfw/deta = 20.(1-exp(-5eta))/(1+2eta)**2  - 1
!
!
!  How many levels for rlev ?
!  Quick look at function in wave indicates that there may be problems
!  near the lowest end of the function where eta is ~ 0.0
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Note fw(0) = 0.0
      dfw_deta(0) = -1.
      fw(0) = 0.0
      eta_reg(0) = 0.0
!-----------------------------------------------------------------------
! calculate function at n regularly spaced points


! (a) Trapezium rule method  - checked

        inc_eta=1.0/(float(rlev))
        Do k = 1,rlev
          eta_reg(k) = k*inc_eta
          temp = (1.+2.*eta_reg(k))*(1.+2.*eta_reg(k))
          dfw_deta(k) = 20.*(1.-exp(-5.*eta_reg(k)))/temp  -1.
        End do
        Do k = 1,rlev
          fw(k) = fw(k-1)+ 0.5*inc_eta*(dfw_deta(k) + dfw_deta(k-1))
        End do

!-----------------------------------------------------------------------
! Calculation at actual points (uses trapezium assumption for residual
!  part after dfinding nearest point
!  (i) find nearest regular point below required eta
!  (ii) work out value of differential at eta
!  (iii) calculate value at point using this info

      Do k =1,max_cldlev
        Do i = 1,n_sh
          If (eta(i,k) <= 1.0) then  ! in cloud

            ibelow = int(eta(i,k)/inc_eta)
            temp = (1+2.*eta(i,k))*(1+2.*eta(i,k))
            dfw_at_eta = 20.*(1.-exp(-5.*eta(i,k)))/temp  -1.
            fw_func(i,k) = fw(ibelow) + 0.5*(eta(i,k)-eta_reg(ibelow))  &
     &                        *(dfw_at_eta + dfw_deta(ibelow))

          Else                  ! non cloud level
            fw_func(i,k) = 0.0
          End if           ! test on eta <1.0
        End Do
      End Do


!-----------------------------------------------------------------------
!   End Subroutine
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('LES_SHALL_FUNC_THLEV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE les_shall_func_thlev
