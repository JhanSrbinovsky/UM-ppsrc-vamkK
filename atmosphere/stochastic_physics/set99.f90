! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE set99_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE set99(trigs,ifax,n)

!  computes factors of n & trigonometric
!  functions required by fourier

USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport

USE c_skeb2_mod, ONLY: nblock, nfacts
USE ereport_mod, ONLY: ereport
  IMPLICIT NONE

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


! nfacts and nblock values

REAL    :: trigs(n), del, angle
INTEGER :: ifax(nfacts), jfax(nfacts), n, nhl, k, nu, nfax, ifac, i
INTEGER :: lfax(1:7) = (/6,8,5,4,3,2,1/)

INTEGER                       ::  icode
CHARACTER (Len=80)            ::  cmessage
CHARACTER (Len=* ), PARAMETER ::  RoutineName='SET99'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SET99',zhook_in,zhook_handle)
del=4.0*ASIN(1.0)/FLOAT(n)
nhl=(n/2)-1
DO k=0,nhl
   angle=FLOAT(k)*del
   trigs(2*k+1)=COS(angle)
   trigs(2*k+2)=SIN(angle)
   END DO

!     find factors of n (8,6,5,4,3,2; only one 8 allowed)
!     look for the eight first and store factors in descending order

nu=n
nfax=1
IF( mod(nu,8) == 0 ) THEN
   jfax(1)=8
   nu=nu/8
   nfax=2
   END IF

ifac=6
DO WHILE ( nu /= 1 .AND. nfax < nfacts)
  IF( MOD(nu,ifac) == 0 ) THEN
    jfax(nfax)=ifac
    nu=nu/ifac
    nfax=nfax+1
  ELSE
    ifac=ifac-1
    IF( ifac == 1 ) THEN
      WRITE(6,'(A4,I4,A27)')'1n =',n,' - contains illegal factors'
      cmessage = 'Too may factors'
      icode    = 1
      CALL EReport(RoutineName,icode,cmessage)  
    END IF
  END IF
END DO
IF( nfax > nfacts-1 ) THEN
  WRITE(6,'(A4,I4,A27)')'1n =',n,' - contains illegal factors'
  cmessage = 'Too may factors'
  icode    = 1
  CALL EReport(RoutineName,icode,cmessage)
END IF

!     now reverse order of factors

nfax=nfax-1
ifax(1)=nfax
DO i=1,nfax
   ifax(nfax+2-i)=jfax(i)
   END DO
ifax(nfacts)=n
IF (lhook) CALL dr_hook('SET99',zhook_out,zhook_handle)
RETURN

END SUBROUTINE set99

END MODULE set99_mod
