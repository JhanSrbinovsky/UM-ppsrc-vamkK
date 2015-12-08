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
!    each condensable gas phase concentration for
!    error checking in box model.
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
      SUBROUTINE UKCA_CALCMINMAXGC(STR_AT,NBOX,GC,myproc)
!----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Used to calculate and print out max,min,mean of
!     each condensable gas phase concentration for
!     error checking in box model
!
!     Inputs
!     ------
!     STR_AT     : String indicating which process the code is up to
!     NBOX       : Number of grid boxes
!     GC         : Condensable cpt number density (molecules cm-3)
!     myproc     : Processor number
!
!     Outputs
!     -------
!     None
!
!     Local variables
!     ---------------
!     GCMIN,GCMAX,GCMEAN : min,max,mean of gas phase condensable conc.
!
!--------------------------------------------------------------------
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      USE PrintStatus_mod
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      IMPLICIT NONE

! Arguments
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: myproc
      REAL, INTENT(IN)    :: GC(NBOX,NCHEMG)
      CHARACTER(len=30), INTENT(IN) :: STR_AT

! Local variables
      INTEGER :: JL
      INTEGER :: JV
      REAL    :: GCMIN
      REAL    :: GCMAX
      REAL    :: GCMEAN

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_CALCMINMAXGC',zhook_in,zhook_handle)

      DO JV=1,NCHEMG
       IF(CONDENSABLE(JV)) THEN
        IF (PrintStatus >= PrStatus_Diag) THEN
         WRITE(6,'(A)') '**********************************'
         WRITE(6,'(A)') STR_AT
        END IF
        GCMIN=1.0e9
        GCMAX=-1.0e9
        GCMEAN=0.0e0
        DO JL=1,NBOX
         GCMIN=MIN(GC(JL,JV),GCMIN)
         GCMAX=MAX(GC(JL,JV),GCMAX)
         GCMEAN=GCMEAN+GC(JL,JV)/FLOAT(NBOX)
        ENDDO
        IF (PrintStatus >= PrStatus_Diag) THEN
         WRITE(6,'(A20,I6,3E12.3)') 'GC:JV,Min,max,mean=',JV,GCMIN,     &
                                    GCMAX,GCMEAN
         WRITE(6,'(A)') '**********************************'
        END IF
       END IF
      END DO

      IF (lhook) CALL dr_hook('UKCA_CALCMINMAXGC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CALCMINMAXGC
