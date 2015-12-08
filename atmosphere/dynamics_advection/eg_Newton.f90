! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
SUBROUTINE eg_newton(p,theta,rho,t0_ref,g,z,                          &
                     intw_rho2w, intw_w2rho, p_0,n,prof_type,u_term)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod, ONLY : ereport
IMPLICIT NONE
!
! Description:
!  
! The code solves the 1D initialization problem using
! Newtons method
!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


INTEGER, PARAMETER                      :: itmax = 1000
REAL,    PARAMETER                      :: tol   = 1.0e-11
INTEGER,                     INTENT(IN) :: n, prof_type
REAL,                        INTENT(IN) :: intw_rho2w(n,2)
REAL,                        INTENT(IN) :: intw_w2rho(n,2)
REAL,                        INTENT(IN) :: p_0
REAL,                        INTENT(IN) :: u_term(1:n)

REAL                                    :: t0_ref(0:n)
REAL                                    :: g(0:n)
REAL                                    :: theta(0:n)   
REAL                                    :: z(0:n)   
REAL                                    :: p(0:n+1)  
REAL                                    :: rho(n)
INTEGER                                 :: i, k
REAL                                    :: err
REAL                                    :: x(3*n+3)
REAL                                    :: f(3*n+3)
REAL                                    :: dx(3*n+3)
INTEGER                                 :: it


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_NEWTON',zhook_in,zhook_handle)


! Set initial guess (interleave the variables to
! reduce the bandwidth of the Jacobian matrix)

k = 1
DO i = 0, n
   IF( i == 0 ) THEN
      x(k) = rho(1)
   ELSE
      x(k) = rho(i)
   END IF
   k    = k + 1
   x(k) = p(i)
   k    = k + 1
   x(k) = theta(i)
   k    = k + 1
END DO

! Perform Newton iterations until converged


dx = 0.0
DO it = 1, itmax

! DEPENDS ON: eg_1d_equations
   CALL eg_1d_equations(f,x,t0_ref,g,z,intw_rho2w,intw_w2rho,         &
                        p_0,n,prof_type, u_term)

   err = MAXVAL(ABS(f))+MAXVAL(ABS(dx))
   IF( err < tol ) EXIT

! DEPENDS ON: eg_invert_jacob
   CALL eg_invert_jacob(dx,f,x,t0_ref,g,z,                            &
                        intw_rho2w, intw_w2rho,n,prof_type)!, u_term)

   x = ABS(x - dx)
END DO

IF(it >= itmax .AND. err > tol ) THEN

     Call ereport("eg_newton", 1,                                     &
                  "Newton iteration failed to converge" )
END IF
! Copy solution back (discarding the surface rho value)

k = 1
DO i = 0, n
   IF( i /= 0 ) rho(i)   = x(k)
   k        = k + 1
   p(i)     = x(k)
   k        = k + 1
   theta(i) = x(k)
   k        = k + 1
END DO

IF (lhook) CALL dr_hook('EG_NEWTON',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_newton
