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
!    Cloud processing of aerosol.
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
      SUBROUTINE UKCA_CLOUDPROC(NBOX,ND,MD,MDT,DRYDP,                   &
       LOWCLOUD,VFAC,ACT,VERBOSE,IACTMETHOD,BUD_AER_MAS)
!----------------------------------------------------------------------
!
!  Purpose
!  -------
!  If gridbox is in-cloud then accounts for the fact that
!  particles in the Aitken soluble mode will be activated
!  to form cloud droplets.
!  Calculates the fraction of the number & mass of Aitsol
!  particles which have dry radii larger than the activation
!  radius "ACT". Then transfers this mass and number
!  from Aitsol to accsol mode (updates number and mass of each)
!  so that Aitken soluble mode particles which were large
!  enough will then receive in-cloud produced S(VI) as
!  they will now be in the accsol mode.
!
!  Currently uses constant value for activation radius as in
!  first versions of GLOMAP (as published in Spracklen et al (2005)
!  This is IACTMETHOD=1.
!
!  Later want to use input updraft velocity to diagnose
!  maximum supersaturation and hence activation radius using
!  Nenes & Seinfeld (2003) parameterisation (IACTMETHOD=2).
!
!  Inputs
!  ------
!  NBOX        : Number of grid boxes
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  DRYDP       : Geometric mean dry diameter of particles in each mode (m)
!  LOWCLOUD    : Horizontal low cloud fraction
!  VFAC        : Vertical low cloud fraction
!  ACT         : Particle dry radius above which activation is assumed
!  VERBOSE     : If =1 prints min/max (ND,MDT etc) after each process
!                for 1st grid box (for box model tests)
!  IACTMETHOD  : Switch for activation method (0=off,1=fixed ract,2=NSO3 scheme)
!
!  Outputs
!  -------
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  BUD_AER_MAS : Aerosol mass budgets
!
!  Local variables
!  ---------------
!  DP       : Number median diameter of mode (m)
!  LNRATN   : Log of ratio of threshold diameter to no. median diameter
!  ERFNUM   : LNRATN/(sqrt(2)*log(sigmag))
!  FRAC_N   : Fraction of ptcl number which is within bounds
!  DELN     : Number concentration to transfer due to cloud processing
!  LOG2SG   : log(sigmag)*log(sigmag)
!  DP2      : Volume median diameter of mode (m)
!  LNRATM   : Log of ratio of threshold diameter to mass median diameter
!  ERFMAS   : LNRATM/(sqrt(2)*log(sigmag))
!  FRAC_M   : Fraction of ptcl mass which is within bounds
!  DM       : Cpt mass concentration to transfer due to cloud processing
!  NEWN     : New particle number conc in mode IMODE
!  NEWNP1   : New particle number conc in mode IMODE+1
!
!  Inputted by module UKCA_MODE_SETUP
!  ----------------------------------
!  NMODES      : Number of possible aerosol modes
!  NCP         : Number of possible aerosol components
!  MODE        : Logical variable defining which modes are set.
!  COMPONENT   : Logical variable defining which cpt are in which dsts
!  DDPMID      : Mid-point of size mode = exp(0.5*(lndp0+lndp1)) (m)
!  SIGMAG      : Geometric standard deviation for each mode
!  MFRAC_0     : Initial mass fraction to set when no particles.
!  MMID        : Mass of particle with dp=dpmed_g=exp(0.5*(lndp0+lndp1)) (ptcl^-1)
!  NUM_EPS     : Value of NEWN below which do not recalculate MD (per cc)
!              : or carry out process
!  CP_SU       : Component where sulfate is stored
!  CP_BC       : Component where black carbon is stored
!  CP_OC       : Component where organic carbon is stored
!  CP_SO       : Component where condensible organic species is stored
!
!  Inputted by module UKCA_SETUP_INDICES
!  -------------------------------------
!  Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_MODE_SETUP,    ONLY: nmodes, ncp, mode, component,       &
                                    ddpmid, sigmag, mfrac_0, mmid,      &
                                    num_eps, cp_su, cp_bc, cp_oc, cp_so
      USE UKCA_SETUP_INDICES

      USE yomhook,            ONLY: lhook, dr_hook
      USE parkind1,           ONLY: jprb, jpim
      IMPLICIT NONE

! .. Input/output variables
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: VERBOSE
      INTEGER, INTENT(IN) :: IACTMETHOD
      REAL, INTENT(IN)    :: DRYDP(NBOX,NMODES)
      REAL, INTENT(IN)    :: LOWCLOUD(NBOX),VFAC(NBOX)
      REAL, INTENT(IN)    :: ACT
      REAL, INTENT(INOUT) :: ND(NBOX,NMODES)
      REAL, INTENT(INOUT) :: MDT(NBOX,NMODES)
      REAL, INTENT(INOUT) :: MD(NBOX,NMODES,NCP)
      REAL, INTENT(INOUT) :: BUD_AER_MAS(NBOX,0:NBUDAER)

! .. Local variables
      INTEGER :: JL
      INTEGER :: IMODE
      INTEGER :: ICP
      REAL    :: F
      REAL    :: DP
      REAL    :: LNRATN
      REAL    :: ERFNUM
      REAL    :: FRAC_N
      REAL    :: DELN
      REAL    :: LOG2SG
      REAL    :: DP2
      REAL    :: LNRATM
      REAL    :: ERFMAS
      REAL    :: FRAC_M
      REAL    :: DM(NCP)
      REAL    :: NEWN
      REAL    :: NEWNP1
      REAL    :: ERF       ! function
      REAL    :: DERF      ! function

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_CLOUDPROC',zhook_in,zhook_handle)

      IF(IACTMETHOD == 1) THEN
       DO JL=1,NBOX
        F=LOWCLOUD(JL)*VFAC(JL)
        IF(F > 0.0) THEN ! if in cloud
         IMODE=2 ! apply for Aitsol mode (transfer to accsol)
! .. calculate fraction of Aitsol mode no & mass with r>ACT
         DP=DRYDP(JL,IMODE)
         IF(ND(JL,IMODE) > NUM_EPS(IMODE)) THEN
          LNRATN=LOG(ACT*2.0/DP) ! use 2nd threshold to work out fraction
          ERFNUM=LNRATN/SQRT(2.0)/LOG(SIGMAG(IMODE))
          FRAC_N=0.5*(1.0+DERF(dble(ERFNUM))) ! fraction remaining in mode
          IF(FRAC_N < 0.5) FRAC_N=0.5 ! limit DELN to be max = half # of ptcls
          IF(FRAC_N < 0.99) THEN ! if more than 1% activated
           DELN=ND(JL,IMODE)*(1.0-FRAC_N)
           DELN=DELN*F ! modify change in number taking into account gridbox
!                        cloud fraction
           LOG2SG=LOG(SIGMAG(IMODE))*LOG(SIGMAG(IMODE))
           DP2=EXP(LOG(DP)+3.0*LOG2SG) ! volume median diameter
           LNRATM=LOG(ACT*2.0/DP2)
           ERFMAS=LNRATM/SQRT(2.0)/LOG(SIGMAG(IMODE))
           FRAC_M=0.5*(1.0+DERF(dble(ERFMAS)))
           NEWN=ND(JL,IMODE)-DELN
           NEWNP1=ND(JL,IMODE+1)+DELN
           IF(NEWN > NUM_EPS(IMODE)) THEN
            DO ICP=1,NCP
             IF(COMPONENT(IMODE,ICP)) THEN
! .. calculate amount of each component to transfer to next mode up
!         (use old no/mass)
              DM(ICP)=MD(JL,IMODE,ICP)*ND(JL,IMODE)*(1.0-FRAC_M)
              DM(ICP)=DM(ICP)*F ! modify change in mass taking into account
!                                 gridbox cloud fraction
             ELSE
              DM(ICP)=0.0
             END IF
             IF((ICP == CP_SU).AND.(NMASPROCSUINTR23 > 0))              &
                          BUD_AER_MAS(JL,NMASPROCSUINTR23)=             &
                          BUD_AER_MAS(JL,NMASPROCSUINTR23)+DM(ICP)
             IF((ICP == CP_BC).AND.(NMASPROCBCINTR23 > 0))              &
                          BUD_AER_MAS(JL,NMASPROCBCINTR23)=             &
                          BUD_AER_MAS(JL,NMASPROCBCINTR23)+DM(ICP)
             IF((ICP == CP_OC).AND.(NMASPROCOCINTR23 > 0))              &
                          BUD_AER_MAS(JL,NMASPROCOCINTR23)=             &
                          BUD_AER_MAS(JL,NMASPROCOCINTR23)+DM(ICP)
             IF((ICP == CP_SO).AND.(NMASPROCSOINTR23 > 0))              &
                          BUD_AER_MAS(JL,NMASPROCSOINTR23)=             &
                          BUD_AER_MAS(JL,NMASPROCSOINTR23)+DM(ICP)
            END DO
!
! .. first remove mass to be transferred from mode IMODE
            MDT(JL,IMODE)=0.0
            DO ICP=1,NCP
             IF(COMPONENT(IMODE,ICP)) THEN
              MD(JL,IMODE,ICP)=                                         &
                  (ND(JL,IMODE)*MD(JL,IMODE,ICP)-DM(ICP))/NEWN
              MDT(JL,IMODE)=MDT(JL,IMODE)+MD(JL,IMODE,ICP)
             ELSE
              MD(JL,IMODE,ICP)=0.0
             END IF ! COMPONENT(IMODE,ICP)
            END DO
! .. now set new number to mode IMODE
            ND(JL,IMODE)=NEWN ! set particle number to new value
!
! .. then add mass to be transferred to mode IMODE+1
            MDT(JL,IMODE+1)=0.0
            DO ICP=1,NCP
             IF(COMPONENT(IMODE+1,ICP)) THEN
              MD(JL,IMODE+1,ICP)=                                       &
          (ND(JL,IMODE+1)*MD(JL,IMODE+1,ICP)+DM(ICP))/NEWNP1
              MDT(JL,IMODE+1)=MDT(JL,IMODE+1)+MD(JL,IMODE+1,ICP)
             ELSE
              MD(JL,IMODE+1,ICP)=0.0
             END IF ! COMPONENT(IMODE+1,ICP)
            END DO
! .. now set new number to mode IMODE+1
            ND(JL,IMODE+1)=NEWNP1
           END IF ! IF NEWN>0
!
          END IF ! if FRAC_N<0.99 (if more than 1% activated)
         ELSE
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            MD(JL,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
           END IF
          END DO
          MDT(JL,IMODE)=MMID(IMODE)
         END IF ! if significant number of particles in lower mode
        END IF ! if F>0 (if in cloud)
       END DO ! end loop over boxes
      END IF ! if iactmethod=1

      IF (lhook) CALL dr_hook('UKCA_CLOUDPROC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CLOUDPROC
