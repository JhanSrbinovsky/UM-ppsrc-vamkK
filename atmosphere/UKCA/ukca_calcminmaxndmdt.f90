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
!    Used to calculate and print out max,min,mean of
!    number conc. (ND) & total mass/ptcl (MDT) for each aerosol mode
!    Used for error checking in box model.
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
      SUBROUTINE UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
!----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Used to calculate and print out max,min,mean of
!     number conc. (ND) & total mass/ptcl (MDT) for each aerosol mode
!     Used for error checking in box model
!
!     Inputs
!     ------
!     NBOX       : Number of grid boxes
!     ND         : Aerosol ptcl number density for mode (cm^-3)
!     MDT        : Avg tot mass of aerosol ptcl in mode (particle^-1)
!     myproc     : Processor number
!
!     Outputs
!     -------
!     None
!
!     Local variables
!     ---------------
!     XMINN,XMAXN,XMEANN : min,max,mean of particle number conc.
!     XMINM,XMAXM,XMEANM : min,max,mean of avg total mass per particle
!
!--------------------------------------------------------------------
      USE UKCA_MODE_SETUP,   ONLY: mode, nmodes
      USE parkind1,          ONLY: jprb, jpim
      USE yomhook,           ONLY: lhook, dr_hook
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: myproc
      REAL, INTENT(IN)    :: ND(NBOX,NMODES)
      REAL, INTENT(IN)    :: MDT(NBOX,NMODES)

! Local variables
      INTEGER :: IMODE
      INTEGER :: JL
      REAL    :: XMINN
      REAL    :: XMAXN
      REAL    :: XMEANN
      REAL    :: XMINM
      REAL    :: XMAXM
      REAL    :: XMEANM

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_CALCMINMAXNDMDT',zhook_in,zhook_handle)
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        XMINN=1.0e16
        XMAXN=-1.0e16
        XMEANN=0.0e0
        XMINM=1.0e16
        XMAXM=-1.0e16
        XMEANM=0.
        DO JL=1,NBOX
         XMINN=MIN(ND(JL,IMODE),XMINN)
         XMAXN=MAX(ND(JL,IMODE),XMAXN)
         XMEANN=XMEANN+ND(JL,IMODE)/FLOAT(NBOX)
         XMINM=MIN(MDT(JL,IMODE),XMINM)
         XMAXM=MAX(MDT(JL,IMODE),XMAXM)
         XMEANM=XMEANM+MDT(JL,IMODE)/FLOAT(NBOX)
        END DO
        WRITE(6,'(1i4,6e15.6)') IMODE,XMINN,XMAXN,                      &
                                XMEANN,XMINM,XMAXM,XMEANM
       END IF
      END DO

      IF (lhook) CALL dr_hook('UKCA_CALCMINMAXNDMDT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CALCMINMAXNDMDT
