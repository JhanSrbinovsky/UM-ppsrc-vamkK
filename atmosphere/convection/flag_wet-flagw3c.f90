! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate a mask for when condensation is liquid
!
! Subroutine Interface:
SUBROUTINE flag_wet (np_field, npnts, nlev,                                  &
                     th, exner_layer_centres,                                &
                     bwater )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE cv_run_mod, ONLY:                                                        &
    tice

IMPLICIT NONE

!
! Description : Calculates a mask for when condensation is liquid
!
! Method:
!   If 0.5 * (TK + TK+1) > TICE  then any condensation  in layer k+1 is liquid
!
!   If 0.5 * (TK + TK+1) < TICE  then any condensation  in layer k+1 is ice
!
!  Returns bwater - logical true if liquid rather than ice
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

INTEGER, INTENT(IN) :: &
  np_field             &  ! Full vector length
 ,npnts                &  ! Vector length
 ,nlev                    ! Number of model layers

REAL, INTENT(IN) ::                    &
  th(np_field,nlev)                    & ! Potential temperature (K)
 ,exner_layer_centres(np_field,0:nlev)   ! Exner ratio at layer centres
                                         ! (starting with the surface).
LOGICAL, INTENT(OUT) ::  &
  bwater(npnts,2:nlev)     ! mask for those points at which condensate
                           ! is liquid

!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER ::       &
  i,k               ! Loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
!  Calculate mask
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('FLAG_WET',zhook_in,zhook_handle)

DO k=1,nlev-1
  DO i=1,npnts

    bwater(i,k+1) = 0.5*(th(i,k)*exner_layer_centres(i,k) +               &
                         th(i,k+1)*exner_layer_centres(i,k+1))  >  tice

  END DO  ! npnts
END DO    !  nlev -1

IF (lhook) CALL dr_hook('FLAG_WET',zhook_out,zhook_handle)
RETURN
END SUBROUTINE flag_wet
