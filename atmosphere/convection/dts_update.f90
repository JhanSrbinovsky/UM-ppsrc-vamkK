! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! applies increment dpdt to oldprof
!
! Subroutine Interface:
SUBROUTINE dts_update(n_dp,nlev,dts_ntpar,ntparmax,                     &
                      timestep,oldprof,dpdt,newprof)

! Modules none used


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Applies the increment
! new profile = old profile + timestep  x increment
!
! Inputs:
! -------
! oldprof, dpdt
!
! Outputs:
! --------
! newprof
!
! Weaknesses: 
! -----------
! Relies on the old profile being on theta levels, and the increments
! also on theta levels
!
! Other information: 
! ------------------ 
! called from FLUX CALCULATION section of deep_turb_conv
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  n_dp                 & ! No. of deep convection points
 ,nlev                 & ! No. of model layers
 ,dts_ntpar(n_dp)      & ! Top level of initial parcel ascent
 ,ntparmax               ! Maximum of ntpar across all conv points

REAL, INTENT(IN) ::    &
  timestep             & ! timestep
 ,oldprof(n_dp,nlev)   & ! profile before increment on theta levels
 ,dpdt(n_dp,nlev)        ! on rho levels
        
REAL, INTENT(OUT) ::  &
  newprof(n_dp,nlev)    ! updated profile

! ------------------------------------------------------------------------------
! Local loop variables
INTEGER ::            &
  i_dp,k                ! loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ------------------------------------------------------------------------------
! First set the new profile to be the old profile

IF (lhook) CALL dr_hook('DTS_UPDATE',zhook_in,zhook_handle)
  newprof(:,:) = oldprof(:,:) 

  DO k=1,ntparmax
    DO i_dp=1,n_dp
    ! Add to this the increment
      IF(k <= dts_ntpar(i_dp)) THEN
        newprof(i_dp,k)=newprof(i_dp,k)+timestep*dpdt(i_dp,k)
      END IF
    END DO
  END DO
  IF (lhook) CALL dr_hook('DTS_UPDATE',zhook_out,zhook_handle)
  RETURN

END subroutine dts_update
! ------------------------------------------------------------------------------
