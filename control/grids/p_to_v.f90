! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     subroutine P_TO_V  for calculating variables held at p points
!     at v points, on both the new dynamics and endgame grids
!
!     This routine does interior points of array not halos,
!     but requires halo information to be set.
!
!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: Grids
!
!     global code has E-W wrap around
!     LAM code has most east u point set to zero
MODULE p_to_v_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE p_to_v(array_on_p_points,                              &
                  ini_start,ini_end,                              &
                  inj_start,inj_end,                              &
                  outi_start,outi_end,                            &
                  outj_start,outj_end,                            &
                  outk_start,outk_end,                            &
                  array_on_v_points)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN)  ::   ini_start,ini_end
INTEGER, INTENT(IN)  ::   inj_start,inj_end
INTEGER, INTENT(IN)  ::   outi_start,outi_end
INTEGER, INTENT(IN)  ::   outj_start,outj_end
INTEGER, INTENT(IN)  ::   outk_start,outk_end

REAL, INTENT(IN) ::  array_on_p_points(ini_start:ini_end,         &
                                       inj_start:inj_end,         &
                                       outk_start:outk_end)

REAL, INTENT(OUT)::  array_on_v_points(outi_start:outi_end,       &
                                       outj_start:outj_end,       &
                                       outk_start:outk_end)

! local variables

INTEGER  ::   i 
INTEGER  ::   j
INTEGER  ::   k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('P_TO_V',zhook_in,zhook_handle)

DO k=outk_start,outk_end
  DO j=outj_start,outj_end
    DO i=outi_start,outi_end


      array_on_v_points(i,j,k)= 0.5 *                             &
    ( array_on_p_points(i,j,k) + array_on_p_points(i,j+1,k) )




    END DO

  END DO
END DO

IF (lhook) CALL dr_hook('P_TO_V',zhook_out,zhook_handle)
RETURN
END SUBROUTINE p_to_v
END MODULE p_to_v_mod
