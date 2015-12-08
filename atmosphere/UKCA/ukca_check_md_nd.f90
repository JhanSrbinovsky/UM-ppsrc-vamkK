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
!    Checks component and total average masses MD, MDT [per ptcl] and
!    number concentrations ND for each mode for bad values.
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
      SUBROUTINE UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
!-----------------------------------------------------------
!
!   Purpose
!   -------
!   Checks cpt and total average masses MD, MDT [per ptcl] and
!   number concentrations ND for each mode for bad values.
!
!   Bad is taken to be:
!
!   MD : greater than 10^20 or less than zero.
!   MDT: greater than 10^20 or less than or equal to zero.
!   ND : greater than 10^9 per cc or less than zero.
!
!   Inputs
!   ------
!   NBOX      : Number of grid boxes
!   PROCESS   : Character string with process just completed
!   ND        : Aerosol ptcl number density for mode (cm^-3)
!   MD        : Avg cpt mass of aerosol ptcl in size mode (particle^-1)
!   MDT       : Avg tot mass of aerosol ptcl in size mode (particle^-1)
!   myproc    : Processor number
!
!   Outputs
!   -------
!   None
!
!   Inputted by module UKCA_MODE_SETUP
!   ----------------------------------
!   NMODES    : Number of possible aerosol modes
!   NCP       : Number of possible aerosol components
!   MODE      : Logical variable defining which modes are set.
!   COMPONENT : Logical variable defining which cpt are in which dsts
!   MLO       : Lo-interf masses for initial radius grid
!-----------------------------------------------------------

      USE UKCA_MODE_SETUP,   ONLY: nmodes, ncp, mode, component, mlo
      USE parkind1,          ONLY: jprb, jpim
      USE yomhook,           ONLY: lhook, dr_hook
      USE ereport_mod,       ONLY: ereport
      USE PrintStatus_mod,   ONLY: PrintStatus, PrStatus_Oper
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: myproc
      REAL, INTENT(IN)    :: ND(NBOX,NMODES)
      REAL, INTENT(IN)    :: MDT(NBOX,NMODES)
      REAL, INTENT(IN)    :: MD(NBOX,NMODES,NCP)
      CHARACTER(len=30), INTENT(IN) :: PROCESS

! Local variables
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: TOTERR
      INTEGER :: ERRCODE
      LOGICAL :: MASK1(NBOX)
      LOGICAL :: MASK2(NBOX)
      LOGICAL :: MASK12(NBOX)
      LOGICAL :: MASK3(NBOX)
      LOGICAL :: MASK4(NBOX)
      LOGICAL :: MASK5(NBOX)
      REAL    :: SUMMD(NBOX,NMODES)
      REAL    :: MDTMIN
      REAL    :: MDTMAX
      REAL    :: NDMAX
      CHARACTER(LEN=72) :: cmessage
      
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
! .. also checks if the sum of MD over defined cpts is zero
      SUMMD(:,:)=0.0
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          SUMMD(:,IMODE)=SUMMD(:,IMODE)+MD(:,IMODE,ICP)
         END IF
        END DO
       END IF
      END DO
!
      IF (lhook) CALL dr_hook('UKCA_CHECK_MD_ND',zhook_in,zhook_handle)
!
      TOTERR=0
      MDTMAX=1.0E20
      NDMAX =1.0E9
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
!
! ..    set minimum allowed value of MDT
        MDTMIN=MLO(IMODE)*0.001 ! MDTMIN equiv. to DPLIM0*0.1
!
        MASK1(:)=((ND (:,IMODE) <  NDMAX).AND.                          &
                  (ND (:,IMODE) >= 0.0E0))
        MASK2(:)=((MDT(:,IMODE) <  MDTMAX).AND.                         &
                  (MDT(:,IMODE) >= MDTMIN*0.001)) ! 0.01*DPLIM0 here
        MASK12(:)=((.NOT.MASK1(:)).OR.(.NOT.MASK2(:)))
        MASK4(:)=(SUMMD(:,IMODE) <= MDTMIN*0.001) ! 0.01*DPLIM0 here
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          MASK3(:)=((MD(:,IMODE,ICP) < MDTMAX).AND.                     &
                    (MD(:,IMODE,ICP) >= 0.0))
          MASK5(:)=((MASK12(:).OR.(.NOT.MASK3(:))).OR.MASK4(:))
          DO JL=1,NBOX
           IF(MASK5(JL)) THEN
             IF (PrintStatus >= PrStatus_Oper) THEN
               write(6,'(A)') PROCESS
               WRITE(6,'(3i6,4e13.5)') JL,IMODE,ICP,                    &
                 ND(JL,IMODE),MDT(JL,IMODE),MD(JL,IMODE,ICP),           &
                 (ND(JL,IMODE)*MD(JL,IMODE,ICP))
             END IF
            TOTERR=TOTERR+1
           END IF ! if bad value for either ND,MDT or MD
          END DO ! loop over boxes
         END IF ! COMPONENT(IMODE,ICP)
        END DO
       END IF ! MODE(IMODE)
      END DO ! loop over modes

      IF(TOTERR > 0) THEN
       cmessage='Extreme values of ND, MDT, or MD found '
       IF (PrintStatus >= PrStatus_Oper)                                & 
         WRITE(6,'(A40,A10,I10)') cmessage,' Toterr: ',toterr
       errcode=-1
       CALL ereport('UKCA_CHECK_MD_ND',errcode,cmessage)
      END IF

      IF (lhook) CALL dr_hook('UKCA_CHECK_MD_ND',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CHECK_MD_ND
