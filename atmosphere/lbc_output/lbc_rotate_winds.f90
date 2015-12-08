! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Rotates LBC winds to rotated grid

! Subroutine Interface:

SUBROUTINE lbc_rotate_winds (                                     &
           u_p,                                                   &
           v_p,                                                   &
           lbc_size_p,                                            &
           lbc_levels,                                            &
           coeff1,                                                &
           coeff2           )


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Description:
!   Rotate LBC wind components for a rotated grid.

! Method:
!   Rotation Coefficients calculated in LBC_Interp_Coeffs are used to
!   transform u and v to rotated u and v.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.1 programming standards.

! Declarations:

! Subroutine arguments

INTEGER :: lbc_size_p                      !  size of lbc grid
INTEGER :: lbc_levels                      !  no of lbc levels

REAL    :: u_p ( lbc_size_p, lbc_levels )  ! u lbcs on p grid
REAL    :: v_p ( lbc_size_p, lbc_levels )  ! v lbcs on p grid
REAL    :: coeff1 (lbc_size_p)             ! \ rotation
REAL    :: coeff2 (lbc_size_p)             ! / coefficients

! Local scalars:

INTEGER :: i      ! loop index for points
INTEGER :: level  ! loop index for levels
REAL    :: u, v   ! u and v component

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header

IF (lhook) CALL dr_hook('LBC_ROTATE_WINDS',zhook_in,zhook_handle)

! ---------------------------------
! Perform the rotation of the winds
! ---------------------------------

!     w_lltoeq not used here.
!     use loop from w_lltoeq and modify to reduce arrays required.

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP&         SHARED(lbc_levels, lbc_size_p, u_p, v_p, coeff1, coeff2) &
!$OMP&         PRIVATE(i, u, v, level)

DO level = 1, lbc_levels

  DO i = 1, lbc_size_p

    u = u_p(i,level)
    v = v_p(i,level)

    u_p(i,level) = coeff1(i) * u - coeff2(i) * v
    v_p(i,level) = coeff1(i) * v + coeff2(i) * u

  END DO

END DO
!$OMP END PARALLEL DO

!     u_p and v_p now contain winds for the rotated grid

IF (lhook) CALL dr_hook('LBC_ROTATE_WINDS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE lbc_rotate_winds
