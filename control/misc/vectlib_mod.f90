! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

MODULE vectlib_mod

! Description:
! This routine acts as an interface to vector versions 
! of intrinsics functions on a platform.
!
! Supported libraries:
! IBM's VMASS (compile using VMASS def)
! Intel's MKL (compile with MKL def)
!
! Default compiles to equivalent do loop over array 
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Misc

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

CONTAINS

SUBROUTINE exp_v(n,x,y)
  USE um_types
  IMPLICIT NONE

! Sets y(i) to the exponential function of x(i), for i=1,..,n

  REAL (KIND=real64) :: y(n), x(n)
  INTEGER :: n
  INTEGER (KIND=integer32) :: l_n
  INTEGER i

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle






  IF (lhook) CALL dr_hook('VECTLIB_MOD:EXP_V',zhook_in,zhook_handle)
  l_n=n








  DO i=1, n
    y(i) = EXP(x(i))
  END DO


  IF (lhook) CALL dr_hook('VECTLIB_MOD:EXP_V',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE exp_v

!-----------------------------------------------------------

SUBROUTINE powr_v(n, x, power, z)
  USE um_types
  IMPLICIT NONE

! Sets z(i) to x(i) raised to the power y(i), for i=1,..,n

  REAL (KIND=real64) :: z(n), x(n), y(n), power
  INTEGER :: n, i
  INTEGER (KIND=integer32) :: l_n

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle






  IF (lhook) CALL dr_hook('VECTLIB_MOD:POWR_V',zhook_in,zhook_handle)
  l_n=n

  DO i=1, n
    z(i) = x(i)**power
  END DO

  IF (lhook) CALL dr_hook('VECTLIB_MOD:POWR_V',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE powr_v

!-----------------------------------------------------------

SUBROUTINE rtor_v(n, x, y, z)
  USE um_types
  IMPLICIT NONE

! Sets z(i) to x(i) raised to the power y(i), for i=1,..,n

  REAL (KIND=real64) :: z(n), x(n), y(n)
  INTEGER :: n
  INTEGER (KIND=integer32) :: l_n
  INTEGER i

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('VECTLIB_MOD:RTOR_V',zhook_in,zhook_handle)
  l_n=n

  DO i=1, n
    z(i) = x(i)**y(i)
  END DO

  IF (lhook) CALL dr_hook('VECTLIB_MOD:RTOR_V',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE rtor_v

!-----------------------------------------------------------

SUBROUTINE sqrt_v(n, x, y)
  USE um_types
  IMPLICIT NONE

! Sets y(i) to the square root of x(i), for i=1,..,n

  REAL (KIND=real64) :: x(n), y(n)
  INTEGER :: n
  INTEGER (KIND=integer32) :: l_n
  INTEGER i

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('VECTLIB_MOD:SQRT_V',zhook_in,zhook_handle)
  l_n=n


  DO i=1, n 
    y(i) = SQRT(x(i))
  END DO

  IF (lhook) CALL dr_hook('VECTLIB_MOD:SQRT_V',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE sqrt_v

!-----------------------------------------------------------

SUBROUTINE oneover_v(n, x, y)
  USE um_types
  IMPLICIT NONE

! Sets y(i) to the reciprocal of x(i), for i=1,..,n

  REAL (KIND=real64) :: x(n), y(n)
  INTEGER :: n
  INTEGER (KIND=integer32) :: l_n
  INTEGER i 

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('VECTLIB_MOD:ONEOVER_V',zhook_in,zhook_handle)
  l_n=n
  
  DO i=1, n
    y(i) = 1/x(i) 
  END DO
  IF (lhook) CALL dr_hook('VECTLIB_MOD:ONEOVER_V',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE oneover_v

!-----------------------------------------------------------

SUBROUTINE log_v (n, x, y)
  USE um_types
  IMPLICIT NONE

! Sets y(i) to the natural logarithm of x(i), for i=1,..,n

  REAL (KIND=real64) :: x(n), y(n)
  INTEGER :: n
  INTEGER (KIND=integer32) :: l_n
  INTEGER i 

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('VECTLIB_MOD:LOG_V',zhook_in,zhook_handle)
  l_n=n

  DO i=1, n 
    y(i) = LOG(x(i))
  END DO

  IF (lhook) CALL dr_hook('VECTLIB_MOD:LOG_V',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE log_v

!-----------------------------------------------------------

SUBROUTINE sin_v(n,x,y)
  USE um_types
  IMPLICIT NONE

! Sets y(i) to the sin function of x(i), for i=1,..,n

  REAL (KIND=real64) :: y(n), x(n)
  INTEGER :: n
  INTEGER (KIND=integer32) :: l_n
  INTEGER i 

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('VECTLIB_MOD:SIN_V',zhook_in,zhook_handle)
  l_n=n

  DO i=1, n
    y(i) = SIN(x(i)) 
  END DO
  
  IF (lhook) CALL dr_hook('VECTLIB_MOD:SIN_V',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE sin_v

!-----------------------------------------------------------

SUBROUTINE cos_v(n,x,y)
  USE um_types
  IMPLICIT NONE

! Sets y(i) to the cos function of x(i), for i=1,..,n

  REAL (KIND=real64) :: y(n), x(n)
  INTEGER :: n
  INTEGER (KIND=integer32) :: l_n
  INTEGER i 

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('VECTLIB_MOD:COS_V',zhook_in,zhook_handle)
  l_n=n

  DO i=1, n 
    y(i) = COS(x(i)) 
  END DO

  IF (lhook) CALL dr_hook('VECTLIB_MOD:COS_V',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE cos_v

!-----------------------------------------------------------
SUBROUTINE acos_v(n,x,y)
  USE um_types
  IMPLICIT NONE

! Sets y(i) to the cos function of x(i), for i=1,..,n

  REAL (KIND=real64) :: y(n), x(n)
  INTEGER :: n
  INTEGER (KIND=integer32) :: l_n
  INTEGER i 

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('VECTLIB_MOD:ACOS_V',zhook_in,zhook_handle)
  l_n=n

  DO i=1, n 
    y(i) = ACOS(x(i)) 
  END DO

  IF (lhook) CALL dr_hook('VECTLIB_MOD:ACOS_V',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE acos_v

!-----------------------------------------------------------

SUBROUTINE asin_v(n,x,y)
  USE um_types
  IMPLICIT NONE

! Sets y(i) to the asin function of x(i), for i=1,..,n

  REAL (KIND=real64) :: y(n), x(n)
  INTEGER :: n
  INTEGER (KIND=integer32) :: l_n
  INTEGER i 

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('VECTLIB_MOD:ASIN_V',zhook_in,zhook_handle)
  l_n=n

  DO i=1, n
    y(i) = ASIN(x(i)) 
  END DO

  IF (lhook) CALL dr_hook('VECTLIB_MOD:ASIN_V',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE asin_v

!-----------------------------------------------------------

END MODULE vectlib_mod
