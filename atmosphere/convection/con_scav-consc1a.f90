! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Purpose: Scavenge aerosol by convective precipitation.
!
SUBROUTINE con_scav(                                                      &
                     timestep,                                            &
                     cblevel,ctlevel,                                     &
                     rain,snow,aerosol                                    &
                      )

USE atm_fields_bounds_mod, ONLY:                                          &
   tdims

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!  Description: 
!      Scavenge aerosol by convective precipitation.
!  Method
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection (murk aerosol routine)
!
! Code description:
!   Language: Fortran 90. 
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER, INTENT(IN) ::            &
 cblevel(tdims%i_end,tdims%j_end) & ! Convective cloud base level.
,ctlevel(tdims%i_end,tdims%j_end)   ! Convective cloud top level.

REAL, INTENT(IN) ::               &
 timestep                         & ! Timestep (s).
,rain(tdims%i_end,tdims%j_end)    & ! Rate of rainfall in this layer from
                                    !  above (kg per sq m per s).
,snow(tdims%i_end,tdims%j_end)      ! Rate of snowfall in this layer from
                                    !  above (kg per sq m per s). 

REAL, INTENT(INOUT) ::                   &
 aerosol(tdims%i_start:tdims%i_end,      & ! Aerosol mixing ratio
         tdims%j_start:tdims%j_end,      &
                     1:tdims%k_end)

!-----------------------------------------------------------------------
!  Define local scalars.
!-----------------------------------------------------------------------
!  (a) Reals effectively expanded to workspace by the Cray (using
!      vector registers).
REAL, PARAMETER ::                                                &
 krain = 1.0e-4                                                   &
,ksnow = 1.0e-4

REAL ::     &
 rrain,rsnow         ! Real workspace.

!  (b) Others.
INTEGER :: i,j,k    ! Loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
IF (lhook) CALL dr_hook('CON_SCAV',zhook_in,zhook_handle)

! Overall rate = KRAIN*(R) where R is in mm/hr=kg/m2/s*3600.0
rrain=krain*timestep*3600.0
rsnow=ksnow*timestep*3600.0

DO j = 1, tdims%j_end
  DO i = 1,tdims%i_end
    IF (ctlevel(i,j) > 0) THEN
      DO k=1,MIN(ctlevel(i,j),tdims%k_end)
        aerosol(i,j,k)=aerosol(i,j,k)/(1.0+rrain*rain(i,j)+rsnow*snow(i,j))
      END DO
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook('CON_SCAV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE con_scav
