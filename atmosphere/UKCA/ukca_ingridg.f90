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
!    Generate sectional particle mass grid in molecules per particle.
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
      SUBROUTINE UKCA_INGRIDG(NBINS,MBSMALL,MBLARGE,MBMID,MBLO,MBHI)
!-----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Generate sectional particle mass grid (in molecules per particle)
!
!     Inputs
!     ------
!     NBINS   : Number of aerosol ptcl size bins
!     MBSMALL : Number of molecules for ptcl at smallest bin mid-pt
!     MBLARGE : Number of molecules for ptcl at largest  bin mid-pt
!
!     Outputs
!     -------
!     MBLO, MBMID, MBHI : Mass of lower, mid and upper bin edge
!
!     Local variables
!     ---------------
!     LGMRANGE : Diff. in logarithms of masses of largest/smallest bins
!     LGMGRID  : Diff. in logarithms of masses of limits of each bin
!-----------------------------------------------------------------------

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN)  :: NBINS
      REAL, INTENT(IN)     :: MBSMALL
      REAL, INTENT(IN)     :: MBLARGE
      REAL, INTENT(OUT)    :: MBMID(NBINS)
      REAL, INTENT(OUT)    :: MBLO(NBINS)
      REAL, INTENT(OUT)    :: MBHI(NBINS)

! Local variables
      INTEGER :: JV
      REAL    :: LGMRANGE
      REAL    :: LGMGRID
      
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
     
      IF (lhook) CALL dr_hook('UKCA_INGRIDG',zhook_in,zhook_handle)

      LGMRANGE=LOG(MBLARGE)-LOG(MBSMALL)
      DO JV=1,NBINS+1
       LGMGRID=LOG(MBSMALL)+LGMRANGE*FLOAT(JV-1)/FLOAT(NBINS)
       IF(JV < (NBINS+1)) MBLO(JV)=EXP(LGMGRID)
       IF(JV > 1) MBHI(JV-1)=EXP(LGMGRID)
      ENDDO
      DO JV=1,NBINS
       MBMID(JV)=EXP(0.5*LOG(MBLO(JV)*MBHI(JV))) ! geometric mean
      ENDDO

      IF (lhook) CALL dr_hook('UKCA_INGRIDG',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_INGRIDG
