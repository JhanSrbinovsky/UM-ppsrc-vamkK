! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
MODULE acscs2_mod

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


! Description:
!     Calculate T-dependent CS2 cross sections, following JPL (2002)

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

SUBROUTINE acscs2(t,wavenm,acs2)
USE ukca_parpho_mod, ONLY: jpwav
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Subroutine interface
REAL, INTENT(IN)  :: t               ! Temperature in kelvin
REAL, INTENT(IN)  :: wavenm(jpwav)   ! Wavelength of each interval in nm
! Absorption cross-section of CS2 in cm^2 for temperature T for each
! interval wavelength.
REAL, INTENT(OUT) :: acs2(jpwav)

! Local variables

INTEGER :: i

REAL :: arg
REAL :: tc
!     CS2: JPL 2002 T=298K values.
REAL, PARAMETER :: acs2d(jpwav) = (/                                           &
  (0.0,i=1,86),                                                                &
  0.02,  0.05,  0.12,  0.30,  0.64,                                            &
  1.36,  2.95,  4.17,  8.53,  9.44,  4.52,  8.63,  3.80,                       &
  1.38,  0.49,  0.35,  0.24,  0.25,  0.12,  0.01,  0.02,                       &
  (0.0,i=1,96)/)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('ACSCS2',zhook_in,zhook_handle)

! Linear interpolation between the cross sections given.
acs2 = 1.0e-20 * acs2d

IF (lhook) CALL dr_hook('ACSCS2',zhook_out,zhook_handle)

END SUBROUTINE acscs2
END MODULE acscs2_mod
