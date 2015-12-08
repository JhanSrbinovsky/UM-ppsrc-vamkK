! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate initial downdraught mass flux
!
SUBROUTINE flx_init_4a5a(npnts,kct,iccb,icct,flx,flx_dd_k,bddi    &
                   ,flx_strt)

USE cv_run_mod,                                                   &
    ONLY: dd_opt

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description : 
!    Calculate initial downdraught mass flux
!
!   Method    : See Unified Model documentation paper 27.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.3.
!-----------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                & ! Vector length
 ,kct                    ! Convective cloud top

INTEGER, INTENT(IN) :: &
  iccb(npnts)          & ! Convective cloud base level 
 ,icct(npnts)            ! Convective cloud top level 


LOGICAL, INTENT(IN) :: &
  bddi(npnts)            ! Mask for points where downdraught may initiate

REAL, INTENT(IN) ::    &
  flx(npnts,kct+1)       ! updraught mass flux (Pa/s)

REAL, INTENT(OUT) ::    &
  flx_dd_k(npnts)       & ! Downdraught mass flux of layer k (Pa/s)
 ,flx_strt(npnts)         ! Updraught mass flux at level downdraught starts
                          ! (Pa/s)

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER ::        & 
  i               & ! loop counter
 ,kddref            ! reference level for downdraught massflux

REAL :: flxscale    ! The scaling factor for the initial downdraught massflux


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Calculate downdraught massflux based of a reference level which is
! 3/4 cloud depth
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('FLX_INIT_4A5A',zhook_in,zhook_handle)

IF (dd_opt == 1) THEN
  flxscale=0.1
ELSE
  flxscale=0.05
END IF

DO i=1,npnts
  IF (bddi(i)) THEN
    kddref = iccb(i) + 0.75*(icct(i) - iccb(i))
    IF (kddref  >=  icct(i)-1) THEN
      kddref=icct(i)-1
    END IF
    flx_strt(i) = flx(i,kddref)
    flx_dd_k(i) = flx_strt(i) * flxscale
  END IF
END DO

IF (lhook) CALL dr_hook('FLX_INIT_4A5A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE flx_init_4a5a
