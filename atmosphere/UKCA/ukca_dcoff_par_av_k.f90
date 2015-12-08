! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Calculates the average (over kth moment) diffusion coefficient
!    for a log-normally distributed aerosol population with
!    geometric mean diameter DG, geometric mean standard
!    deviation SSIGMA.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_DCOFF_PAR_AV_K(NBOX,K,DP,SIGMA,T,DVISC,MFPA,      &
                                     DCOFF_PAR_AV_K)
!--------------------------------------------------------------------
!
! Purpose
! -------
! Calculates the average (over kth moment) diffusion coefficient
! for a log-normally distributed aerosol population with
! geometric mean diameter DG, geometric mean standard
! deviation SIGMA following method in Regional Particulate
! Model as described by Binkowski & Shankar (1995).
! Also follows expression in Seinfeld & Pandis pg 474
! (Stokes-Einstein relation with slip-flow correction)
!
! Parameters
! ----------
! None
!
! Inputs
! ------
! NBOX           : Number of boxes in domain
! K              : Index of moment for calculation
! DP             : Geometric mean particle diameter for mode (m)
! SIGMA          : Geometric standard deviation for mode
! T              : Temperature of air (K)
! DVISC          : Dynamic viscosity of air (kg m-1 s-1)
! MFPA           : Mean free path of air (m)
!
! Outputs
! -------
! DCOFF_PAR_AV_K : Avg. ptcl diffusion coefficient (m^2 s-1)
!
! Local variables
! ---------------
! KNG            : Knudsen number for geo. mean sized particle
! LNSQSG         : ln(SIGMA)*ln(SIGMA)
! PREF           : Prefactor term to expression
! TERM1,TERM2    : Terms in average diff coeff expression
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! ZBOLTZ         : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
! PPI            : 3.1415927.....
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS,   ONLY: zboltz, ppi
      USE yomhook,          ONLY: lhook, dr_hook
      USE parkind1,         ONLY: jprb, jpim
      IMPLICIT NONE
!
! Subroutine interface
      INTEGER, INTENT(IN) :: K
      INTEGER, INTENT(IN) :: NBOX
      REAL, INTENT(IN)    :: DP(NBOX)
      REAL, INTENT(IN)    :: SIGMA
      REAL, INTENT(IN)    :: T(NBOX)
      REAL, INTENT(IN)    :: DVISC(NBOX)
      REAL, INTENT(IN)    :: MFPA(NBOX)
      REAL, INTENT(OUT)   :: DCOFF_PAR_AV_K(NBOX)

! Local variables
      REAL    :: LNSQSG
      REAL    :: TERM1
      REAL    :: TERM2
      REAL    :: TERM3
      REAL    :: TERM4
      REAL    :: KNG(NBOX)
      REAL    :: PREF(NBOX)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_DCOFF_PAR_AV_K',zhook_in,zhook_handle)

      LNSQSG=LOG(SIGMA)*LOG(SIGMA)
      TERM1=(-2.0*FLOAT(K)+1.0)/2.0
      TERM2=(-4.0*FLOAT(K)+4.0)/2.0
      TERM3=EXP(TERM1*LNSQSG)
      TERM4=1.246*EXP(TERM2*LNSQSG)

      KNG(:)=2.0*MFPA(:)/DP(:)
      PREF(:)=ZBOLTZ*T(:)/(3.0*PPI*DVISC(:)*DP(:))
      DCOFF_PAR_AV_K(:)=PREF(:)*(TERM3+TERM4*KNG(:))

      IF (lhook) CALL dr_hook('UKCA_DCOFF_PAR_AV_K',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_DCOFF_PAR_AV_K
