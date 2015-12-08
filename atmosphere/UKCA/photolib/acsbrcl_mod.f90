! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acsbrcl_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate BrCl temperature dependent cross sections,

! Method:
!     Parameterisation from Maric et al [1994]

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA

!  Code Description:
!    Language:  Fortran 95
!    This code is written to UMDP3 standards.

CONTAINS

SUBROUTINE acsbrcl(temp,jpwav,wavenm,abrcl)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
INTEGER, INTENT(IN) :: jpwav
REAL, INTENT(IN)    :: temp
REAL, INTENT(IN)    :: wavenm(jpwav) ! Contains the wavelength intervals.
! abrcl contains the cross sections at model wavelengths. 
! Contains the result on exit.
REAL, INTENT(INOUT) :: abrcl(jpwav)

! Local variables
REAL, PARAMETER :: a1=7.34e-20
REAL, PARAMETER :: a2=4.35e-19
REAL, PARAMETER :: a3=1.12e-19
REAL, PARAMETER :: b1= 68.6
REAL, PARAMETER :: b2=123.6
REAL, PARAMETER :: b3= 84.8
REAL, PARAMETER :: c1=227.6
REAL, PARAMETER :: c2=372.5
REAL, PARAMETER :: c3=442.4
REAL, PARAMETER :: we=443.1
REAL, PARAMETER :: h=6.63e-34
REAL, PARAMETER :: c=3.0e10
REAL, PARAMETER :: boltz=1.38e-23

REAL :: lam
REAL :: tant
REAL :: tc

INTEGER :: jw

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



!     Check that temperature is in range.
IF (lhook) CALL dr_hook('ACSBRCL',zhook_in,zhook_handle)
tc = temp
tc=MIN(tc, 300.0)
tc=MAX(tc, 200.0)

tant=TANH(h*c*we/(2.0*boltz*tc))

!     Wavelength loop 200nm - 516nm (estimated upper limit for Br-Cl bond)
DO jw = 60,136
  lam =wavenm(jw)
  abrcl(jw)=                                                                   &
    a1*SQRT(tant)*EXP(-b1*tant*(LOG(c1/lam))**2.0)                             &
    + a2*SQRT(tant)*EXP(-b2*tant*(LOG(c2/lam))**2.0)                           &
    + a3*SQRT(tant)*EXP(-b3*tant*(LOG(c3/lam))**2.0)
END DO
IF (lhook) CALL dr_hook('ACSBRCL',zhook_out,zhook_handle)
RETURN

END SUBROUTINE acsbrcl
END MODULE acsbrcl_mod
