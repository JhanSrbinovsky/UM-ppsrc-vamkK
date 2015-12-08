! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  New Shallow convection scheme - Calculations of LES functions
!

      SUBROUTINE les_shall_func_rhlev(n_sh,max_cldlev,nlev              &
     &,       eta_rho, int_type                                         &
     &,     k_func, fng_func, b_func                                    &
     &,     fw_func, g_func,fql_func, gql_func )

!
! Purpose:
!   Calculate various functions of height based on fits to LES data
!   for use in the new turbulence base shallow convection scheme.
!   Values on eta rho levels
!   Done for all cloud levels, rest of array set to zero.
!
!
!   Called by SHALLOW_CONV5A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
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
     &, int_type      ! Method to use for integration
                      ! 1 - trapezium (only available method).

      real, intent(in) ::                                               &
     &  eta_rho(n_sh,nlev) ! non-dimension height of the cloud levels

!
! Arguments with intent INOUT:
!
!     None
!
! Arguments with intent OUT:
!
      real, intent(out) ::                                              &
     &  k_func(n_sh,nlev)                                               &
                             ! Used in turbulent transport
     &, fng_func(n_sh,nlev)                                             &
                             ! Used in turbulent transport
     &, b_func(n_sh,nlev)                                               &
                             ! Used in turbulent transport
     &, g_func(n_sh,nlev)                                               &
                             ! Used in turbulent transport
     &, fw_func(n_sh,nlev)                                              &
                             ! Used in ensemble vertical velocity
     &, fql_func(n_sh,nlev)                                             &
                             ! used in wql calculation
     &, gql_func(n_sh,nlev)  ! used in liquid water flux parametrization

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
      integer :: rlev              ! Number of levels for integration
                                   ! needs to be higher than cloud
                                   ! levels and fixed given resolution
                                   ! of model
      parameter (rlev = 30)        ! ?

      integer :: i,k                                                    &
                                 ! loop counters
     &,  ibelow                                                         &
     &, k_minus, k_centre, k_plus

! variables used to precalculate some of the constants in the
! functions

      real ::                                                           &
     & gg1, gg2, cc1,  dfw_at_eta, inc_eta

! variables used to store exponentials calculated
      real ::                                                           &
     &  wexp6(n_sh,max_cldlev+1)                                        &
                                   ! exp(-6eta)
     &, wexp10(n_sh,max_cldlev+1)  ! exp(-10eta)

! variables used in integral calculation

      real ::                                                           &
     &   temp                                                           &
                   ! temporary store
     &,  eta_reg(0:rlev)                                                &
                          ! eta on regular grid
     &,  fw(0:2*rlev)                                                   &
                            ! fw on regular grid
     &,  dfw_deta(0:rlev) ! dfw/deta on regular grid

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



!-----------------------------------------------------------------------
! 1.0 Initialise arrays - Not required as loop code sets all of array
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 1.0 Functions not defined as integrals i.e. can be evaluated
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
!
!   where c =1.45*4.67/36.
!
!
!  fql = (1+eta-0.6*eta**2)(1-0.6exp(-5eta)
!
!  gql = 1-0.6exp(-4.5eta)
!-----------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('LES_SHALL_FUNC_RHLEV',zhook_in,zhook_handle)

      gg1 = 4.67/36.
      gg2 = 1.33/6.
      cc1 = 0.473/0.64

      Do k =1,max_cldlev+1
        Do i = 1,n_sh
          If (eta_rho(i,k) <  1.0) then      ! in cloud
            wexp10(i,k) = 1. - exp(-10.*eta_rho(i,k))
            wexp6(i,k)  = exp(-6.*eta_rho(i,k))

            k_func(i,k)   = 0.115*(1. + eta_rho(i,k))*wexp10(i,k)
            fng_func(i,k) = 1. - 0.3*(1.+1.75*eta_rho(i,k))*wexp10(i,k)
            b_func(i,k)   = 1.4*(1. + eta_rho(i,k)                      &
     &                       -0.6*eta_rho(i,k)*eta_rho(i,k))            &
     &                                                  *wexp10(i,k)

            g_func(i,k)   = 1.45*(eta_rho(i,k)                          &
     &                         -0.665*eta_rho(i,k)*eta_rho(i,k)         &
     &                    + (gg1 - gg2*eta_rho(i,k))*wexp6(i,k) - gg1)


          Else     ! not in cloud set function to zero

            k_func(i,k)   = 0.0
            fng_func(i,k) = 0.0
            b_func(i,k)   = 0.0
            g_func(i,k)   = 0.0

          End if
         If (eta_rho(i,k) <= 1.0) then      ! in cloud
            fql_func(i,k) = (1.+eta_rho(i,k)                            &
     &                             -0.6*eta_rho(i,k)*eta_rho(i,k))      &
     &                               *( 1.-0.6*exp(-5.*eta_rho(i,k)) )
            gql_func(i,k) = 1. -0.6*exp(-4.5*eta_rho(i,k))

          Else     ! not in cloud set function to zero
            fql_func(i,k)    = 0.0
            gql_func(i,k)    = 0.0
          End if
        End Do
      End Do

!-----------------------------------------------------------------------
! 2.0 Functions defined as differentials with no analytic integral
!-----------------------------------------------------------------------
!
! dfw/deta = 20.(1-exp(-5eta))/(1+2eta)**2  - 1
!
!
!  How many levels for rlev ?
!  Quick look at function in wave indicates that there may be problems
!  near the lowest end of the function where eta is ~ 0.0
!
! dg/deta = 1.45(1-1.33*eta)*(1-exp(-6eta))
!
!-----------------------------------------------------------------------
! Note fw(0) = 0.0

      dfw_deta(0) = -1.
      fw(0) = 0.0
      eta_reg(0) = 0.0

!-----------------------------------------------------------------------
! calculate function at n regularly spaced points
! (Checked code now correct- Note was incorrect at start)


! (a) Trapezium rule method - checked this code

        inc_eta=1.0/(float(rlev))
        Do k = 1,rlev
          eta_reg(k) = k*inc_eta
          temp = (1+2.*eta_reg(k))*(1+2.*eta_reg(k))
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

      Do k =1,max_cldlev+1
        Do i = 1,n_sh
          If (eta_rho(i,k) <= 1.0) then  ! in cloud
            ibelow = int(eta_rho(i,k)/inc_eta)
            temp = (1+2.*eta_rho(i,k))*(1+2.*eta_rho(i,k))
            dfw_at_eta = 20.*(1.-exp(-5.*eta_rho(i,k)))/temp  -1.
            fw_func(i,k)=fw(ibelow) + 0.5*(eta_rho(i,k)-eta_reg(ibelow))&
     &                        *(dfw_at_eta + dfw_deta(ibelow))
          Else                  ! non cloud level
            fw_func(i,k) = 0.0
          End if           ! test on eta <1.0
        End Do
      End Do

!-----------------------------------------------------------------------
!   End Subroutine
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('LES_SHALL_FUNC_RHLEV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE les_shall_func_rhlev
