! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Switches and parameters controlling the deep TCS scheme.


MODULE dts_cntl_mod

  ! Description:
  !   Switches and parameters controlling the deep TCS scheme.
  !
  !   Settings as used in the last 10 year deep turbulence climate run testing 
  !   the deep TCS scheme. These were considered the best values at the time 
  !   but the deep TCS scheme did not give acceptable performance relative to 
  !   the deep mass flux scheme. 
  !
  !   These switches are hard-wired and are not available in the umui, the
  !   only way to alter them is via a code branch.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.
  !
  ! Declarations:

  IMPLICIT NONE
  SAVE

!------------------------------------------------------------------------------
! Deep turbulence scheme sensitivity testing
!------------------------------------------------------------------------------
  INTEGER :: idts_rhfac = 3     ! sensitivity to rhfac (UM8.1 default was 0)
  INTEGER :: idts_zsc   = 0     ! sensitivity to zsc
  INTEGER :: idts_dampfac = 3   ! sensitivity to idampw2 (UM8.1 default was 0)
  INTEGER :: idts_gamma_evp = 0 ! sensitivity to expression for gamma_evp
  INTEGER :: idts_mb = 0        ! sensitivity to expression for mb cal

  REAL :: dts_dif_pow = 1.0 ! Sensitivity to thetav difference between parcel 
                            ! and environ  (UM8.1 default was 1.5)

  REAL :: dts_mb = 0.04     ! Sensitivity to factor used in expression for mb

  REAL :: dts_qfac = 1.0    ! Sensitivity to qfac (UM8.1 default was 1.2)

  REAL :: dts_gamma_evp = 0.0  ! Sensitivity to value of gamma_evp

  REAL :: dts_gamma_fac = 0.30 ! Sensitivity to value of gamma_fac
                               !   (UM8.1 default was 0.125)

END MODULE dts_cntl_mod
