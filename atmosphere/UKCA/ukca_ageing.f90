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
!   Carries out ageing of particles in insoluble mode.
!   Calculates the number of particles which will be coated
!   by the total soluble material and transfers number and
!   mass accordingly.
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
      SUBROUTINE UKCA_AGEING(NBOX,ND,MD,MDT,                            &
       AGETERM1,AGETERM2,WETDP,VERBOSE,BUD_AER_MAS)
!----------------------------------------------------------
!
! Purpose
! ------
! Carries out ageing of particles in insoluble mode.
! Calculates the number of particles which will be coated
! by the total soluble material and transfers number and
! mass accordingly.
! Number of ptcls given by assuming 1 particle ages with
! a 1 molecule layer thickness of particle of geometric
! mean size of mode.  i.e. 1 ptcl requires r_g^2/r_molec^2
! = 10^4,10^6,10^8 molecules for Aitken, accum, coarse modes
!
! Inputs
! ------
! NBOX      : Number of grid boxes
! ND        : Aerosol ptcl no. concentration (ptcls per cc)
! MD        : Component median aerosol mass (molecules per ptcl)
! MDT       : Total median aerosol mass (molecules per ptcl)
! AGETERM1  : Depletion rate of each component (molecules cpt/cc/DTZ)
!             from condensation onto the 3 insoluble modes.
! AGETERM2  : Rate of accomodation of material to each insoluble mode
!             as a result of coagulation with smaller soluble modes
!             (in molecules cpt /cm3/DTZ)
! WETDP     : Wet diameter corresponding to mean sized particle
! VERBOSE   : Switch for whether to do various test print statements
!
! Outputs
! -------
! ND        : Updated number concentration [each mode] (ptcls/cc)
! MD        : Updated avg cpt   mass conc. [each mode] (molecules/ptcl)
! MDT       : Updated avg total mass conc. [each mode] (molecules/ptcl)
! BUD_AER_MAS : Aerosol mass budgets
!
! Local variables
! ---------------
! AGE1PTCL  : Mass of soluble material needed to age 1 ptcl (molecules)
! NAGED     : # of insoluble ptcls aged to soluble mode [total] (/cc)
! NAGED_JV  : # of insoluble ptcls aged to soluble mode [by gas](/cc)
! TOTAGE_JV : Ageing flux to cpt by particular gas (molecules gas/cc)
! TOTAGE1   : Total ageing flux to cpt by conden.  (molecules gas/cc)
! TOTAGE2   : Total ageing flux to cpt by coaguln. (molecules gas/cc)
! TOTAGE    : Total ageing flux to cpt (coag+cond) (molecules gas/cc)
! NDINSNEW  : # in ins mode after reduction due to ageing (/cc)
! NDSOLNEW  : # in corresp. sol mode after reduction due to ageing (/cc)
! F_MM      : Ratio of molar masses of condensable gas to aerosol cpt
! CP_COAG_ADDED : Switch for whether added on ageing flux by
!                 coagulation to that cpt already
!                 (loop over jv --- need to make sure only count once)
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! CONC_EPS  : Likewise, threshold for soluble material (molecules per cc)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES    : Number of possible aerosol modes
! NCP       : Number of possible aerosol components
! DIMEN     : Molecular diameters of condensable components
! MODE      : Which modes are being carried
! SOLUBLE   : Logical variable defining which cpts are soluble
! COMPONENT : Which components are in each of modes
! MMID      : Mid-point masses for initial radius grid
! MFRAC_0   : Initial mass fraction to set when no particles.
! MM        : Molar masses of components (kg per mole)
! NUM_EPS   : Value of NEWN below which do not recalculate MD (per cc)
!             or carry out process
! CP_SU     : Component where sulfate is stored
! CP_BC     : Component in which black carbon is stored
! CP_OC     : Component in which primary organic carbon is stored
! CP_DU     : Component where dust is stored
! CP_SO     : Component where secondary organic carbon is stored
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! NCHEMG    : Number of gas phase tracers for gas phase chemistry scheme
! MM_GAS    : Array of molar masses for gas phase species (kg/mol)
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
!
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! .. Subroutine Interface
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: VERBOSE
      REAL, INTENT(IN)    :: AGETERM1(NBOX,3,NCHEMG)
      REAL, INTENT(IN)    :: AGETERM2(NBOX,4,3,NCP)
      REAL, INTENT(IN)    :: WETDP(NBOX,NMODES)
      REAL, INTENT(INOUT) :: ND(NBOX,NMODES)
      REAL, INTENT(INOUT) :: MD(NBOX,NMODES,NCP)
      REAL, INTENT(INOUT) :: MDT(NBOX,NMODES)
      REAL, INTENT(INOUT) :: BUD_AER_MAS(NBOX,0:NBUDAER)
!
! .. Local variables
      INTEGER :: JL
      INTEGER :: JV
      INTEGER :: IMODE
      INTEGER :: JMODE
      INTEGER :: ICP
      INTEGER :: CP_COAG_ADDED(NCP)
      REAL    :: TOTAGE(NCP)
      REAL    :: TOTAGE_JV
      REAL    :: TOTAGE1(NCP)
      REAL    :: TOTAGE2(NCP)
      REAL    :: AGE1PTCL
      REAL    :: NAGED
      REAL    :: NAGED_JV(NCHEMG)
      REAL    :: NDINSNEW
      REAL    :: NDSOLNEW
      REAL    :: F_MM(NCP)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
      IF (lhook) CALL dr_hook('UKCA_AGEING',zhook_in,zhook_handle)
      DO IMODE=5,7 ! loop over insoluble modes
       IF(MODE(IMODE)) THEN
        DO JL=1,NBOX
         IF(ND(JL,IMODE) > NUM_EPS(IMODE)) THEN
          CP_COAG_ADDED(:)=0
          TOTAGE(:)=0.0
          TOTAGE1(:)=0.0 ! from condensation
          TOTAGE2(:)=0.0 ! from coagulation
          F_MM(:)=0.0 ! need to initialise F_MM to 0 so it is defined
                      ! for components which don't condense
          NAGED=0.0
          DO JV=1,NCHEMG
           TOTAGE_JV=0.0
           NAGED_JV(JV)=0.0
           IF(CONDENSABLE(JV)) THEN

! .. Below add on amount of soluble material taken up by insoluble modes
! .. as result of condensation of this gas onto insoluble modes
! .. This is stored in AGETERM1(:,3,NCP) for the 3 insoluble modes
            ICP=CONDENSABLE_CHOICE(JV)
            F_MM(ICP)=MM_GAS(JV)/MM(ICP)
            TOTAGE_JV=AGETERM1(JL,IMODE-4,JV)/F_MM(ICP)
            TOTAGE (ICP)=TOTAGE (ICP)+TOTAGE_JV
            TOTAGE1(ICP)=TOTAGE1(ICP)+TOTAGE_JV
! .. condensable material taken up by aerosol assumed 100% soluble.
! .. AGETERM1 is from condensation onto insoluble modes
! .. Divide by F_MM because AGETERM1 is in "molecules" of cpt, whereas
! .. TOTAGE needs to be in molecules of condensible gas (H2SO4/SEC_ORG)
! .. since AGE1PTCL is in these units.

! .. Below add on amount of soluble material taken up by insoluble modes
! .. as result of coagulation with soluble modes
!
! .. IF(CP_COAG_ADDED(ICP) == 0) checks if already added on this
! .. aerosol component's coagulated material (make sure don't dblecount)
!
            IF(CP_COAG_ADDED(ICP) == 0) THEN
             DO JMODE=1,4
              TOTAGE_JV=TOTAGE_JV+                                      &
                 AGETERM2(JL,JMODE,IMODE-4,ICP)/F_MM(ICP)
! .. add on amount coagulated
              TOTAGE (ICP)=TOTAGE (ICP)+                                &
                 AGETERM2(JL,JMODE,IMODE-4,ICP)/F_MM(ICP)
              TOTAGE2(ICP)=TOTAGE2(ICP)+                                &
                 AGETERM2(JL,JMODE,IMODE-4,ICP)/F_MM(ICP)
             END DO
             CP_COAG_ADDED(ICP)=1
            ENDIF
! .. condensable material taken up by aerosol assumed 100% soluble.
! .. AGETERM2 is from coagn of soluble modes with larger insoluble modes
! .. Divide by F_MM because AGETERM2 is in "molecules" of cpt, whereas
! .. TOTAGE needs to be in molecules of condensible gas (H2SO4/SEC_ORG)
! .. since AGE1PTCL is in these units.
!
            AGE1PTCL=WETDP(JL,IMODE)*WETDP(JL,IMODE)/DIMEN(JV)/DIMEN(JV)
! above is number of molecules of condensable gas to age 1 particle
            NAGED_JV(JV)=TOTAGE_JV/AGE1PTCL/10.0
!!! divide by 10 is to make particle age by 10 monolayers instead of 1
!!          NAGED_JV(JV)=TOTAGE_JV/AGE1PTCL
! here keep as single monolayer ageing.
            NAGED=NAGED+NAGED_JV(JV)

           END IF ! if gas phase species is condensable
          END DO ! loop over gas phase species
!
          IF(NAGED > NUM_EPS(IMODE)) THEN ! if significant ptcls to age
!
           IF(NAGED > ND(JL,IMODE)) THEN
            NAGED=ND(JL,IMODE) ! limit so no -ves
            TOTAGE(:)=TOTAGE(:)*ND(JL,IMODE)/NAGED
! above reduces ageing if limited by insoluble particles
           END IF
!
           NDINSNEW=ND(JL,IMODE)-NAGED
! set new insoluble mode no. (but don't update ND yet)
           NDSOLNEW=ND(JL,IMODE-3)+NAGED
! set new   soluble mode no. (but don't update ND yet)
!

           IF(SUM(TOTAGE) > 0.0) THEN
            DO ICP=1,NCP
! below transfers aged cpt masses from ins modes (doesn't include SU)
             IF(COMPONENT(IMODE-3,ICP)) THEN
! above if statement checks whether cpt is in corresponding soluble mode
              IF(IMODE == 5) THEN
               IF((ICP == CP_SU).AND.(NMASAGEDSUINTR52 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDSUINTR52)=                          &
             BUD_AER_MAS(JL,NMASAGEDSUINTR52)+TOTAGE(ICP)*F_MM(ICP)
               IF((ICP == CP_OC).AND.(NMASAGEDOCINTR52 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDOCINTR52)=                          &
             BUD_AER_MAS(JL,NMASAGEDOCINTR52)+TOTAGE(ICP)*F_MM(ICP)
               IF((ICP == CP_SO).AND.(NMASAGEDSOINTR52 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDSOINTR52)=                          &
             BUD_AER_MAS(JL,NMASAGEDSOINTR52)+TOTAGE(ICP)*F_MM(ICP)
!! above accounts for transfer of mass of condensed/coagulated material
!! multiply above by F_MM 'cos TOTAGE is in molecules of gas phase
!! species (H2SO4/SEC_ORG) whereas needs to be in "molecules" of
!! aerosol component (CP_SU/CP_OC/CP_SU)
               IF((ICP == CP_BC).AND.(NMASAGEDBCINTR52 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDBCINTR52)=                          &
             BUD_AER_MAS(JL,NMASAGEDBCINTR52)+NAGED*MD(JL,IMODE,ICP)
               IF((ICP == CP_OC).AND.(NMASAGEDOCINTR52 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDOCINTR52)=                          &
             BUD_AER_MAS(JL,NMASAGEDOCINTR52)+NAGED*MD(JL,IMODE,ICP)
!! .. above 2 lines account for transfer of mass due to aged aerosol
              ENDIF ! if mode is Aitken-insoluble
              IF(IMODE == 6) THEN
               IF((ICP == CP_SU).AND.(NMASAGEDSUINTR63 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDSUINTR63)=                          &
             BUD_AER_MAS(JL,NMASAGEDSUINTR63)+TOTAGE(ICP)*F_MM(ICP)
               IF((ICP == CP_OC).AND.(NMASAGEDOCINTR63 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDOCINTR63)=                          &
             BUD_AER_MAS(JL,NMASAGEDOCINTR63)+TOTAGE(ICP)*F_MM(ICP)
               IF((ICP == CP_SO).AND.(NMASAGEDSOINTR63 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDSOINTR63)=                          &
             BUD_AER_MAS(JL,NMASAGEDSOINTR63)+TOTAGE(ICP)*F_MM(ICP)
!! above accounts for transfer of mass of condensed/coagulated material
!! multiply above by F_MM 'cos TOTAGE is in molecules of gas phase
!! species (H2SO4/SEC_ORG) whereas needs to be in "molecules" of
!! aerosol component (CP_SU/CP_OC/CP_SU)
               IF((ICP == CP_DU).AND.(NMASAGEDDUINTR63 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDDUINTR63)=                          &
             BUD_AER_MAS(JL,NMASAGEDDUINTR63)+NAGED*MD(JL,IMODE,ICP)
!! .. above 1 line accounts for transfer of mass due to aged aerosol
              ENDIF ! if mode is accum.-insoluble
              IF(IMODE == 7) THEN
               IF((ICP == CP_SU).AND.(NMASAGEDSUINTR74 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDSUINTR74)=                          &
             BUD_AER_MAS(JL,NMASAGEDSUINTR74)+TOTAGE(ICP)*F_MM(ICP)
               IF((ICP == CP_OC).AND.(NMASAGEDOCINTR74 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDOCINTR74)=                          &
             BUD_AER_MAS(JL,NMASAGEDOCINTR74)+TOTAGE(ICP)*F_MM(ICP)
               IF((ICP == CP_SO).AND.(NMASAGEDSOINTR74 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDSOINTR74)=                          &
             BUD_AER_MAS(JL,NMASAGEDSOINTR74)+TOTAGE(ICP)*F_MM(ICP)
!! above accounts for transfer of mass of condensed/coagulated material
!! multiply above by F_MM 'cos TOTAGE is in molecules of gas phase
!! species (H2SO4/SEC_ORG) whereas needs to be in "molecules" of
!! aerosol component (CP_SU/CP_OC/CP_SU)
               IF((ICP == CP_DU).AND.(NMASAGEDDUINTR74 > 0))            &
             BUD_AER_MAS(JL,NMASAGEDDUINTR74)=                          &
             BUD_AER_MAS(JL,NMASAGEDDUINTR74)+NAGED*MD(JL,IMODE,ICP)
!! .. above 1 line accounts for transfer of mass due to aged aerosol
              ENDIF ! if mode is coarse-insoluble

!! .. below calculates new cpt total masses in soluble modes due to
!! .. transfer of mass from Ait-ins/acc-ins/cor-ins
!! .. (n.b. insoluble avg. masses unchanged)
!! .. (n.b. soluble mode mass not updated due to condensation
!! ..  onto insoluble mode yet)

              IF(COMPONENT(IMODE,ICP)) THEN ! if in sol mode & ins mode
! if sol. mode cpt is in insoluble mode then update corresponding
! soluble mode cpt mass due to trans from ins. mode & cond onto ins
! (n.b. coag already transferred in coagulation routine using COAG_MODE)
               MD(JL,IMODE-3,ICP)=(ND(JL,IMODE-3)*MD(JL,IMODE-3,ICP)    &
                +NAGED*MD(JL,IMODE,ICP)+TOTAGE1(ICP)*F_MM(ICP))/NDSOLNEW
              ELSE ! if in sol mode but not in insoluble mode
! if sol. mode component is not in insoluble mode then just update
! soluble mode cpt mass due to cond onto ins (& change in number)
! (n.b. coag already transferred in coagulation routine using COAG_MODE)
               MD(JL,IMODE-3,ICP)=(ND(JL,IMODE-3)*MD(JL,IMODE-3,ICP)    &
                                       +TOTAGE1(ICP)*F_MM(ICP))/NDSOLNEW
              END IF

             END IF ! if component is in corresponding soluble mode
            END DO ! loop over components
           END IF ! if total amount of accomodated material > 0
!
!! update ND and MDT for ins mode
           ND(JL,IMODE  )=NDINSNEW
           MDT(JL,IMODE)=0.0
           DO ICP=1,NCP
            IF(COMPONENT(IMODE,ICP)) THEN
             MDT(JL,IMODE)=MDT(JL,IMODE)+MD(JL,IMODE,ICP)
            END IF
           END DO ! end loop over cpts
!
!! update ND and MDT for sol mode
           ND(JL,IMODE-3)=NDSOLNEW
           MDT(JL,IMODE-3)=0.0
           DO ICP=1,NCP
            IF(COMPONENT(IMODE-3,ICP)) THEN
             MDT(JL,IMODE-3)=MDT(JL,IMODE-3)+MD(JL,IMODE-3,ICP)
            END IF
           END DO ! end loop over cpts
!
          END IF ! if number of aged particles > NUM_EPS(IMODE)
         END IF ! if some particles in insoluble modes (ND>epsilon)
        END DO ! end loop over boxes
       END IF ! if insoluble mode is present (IMODE=5,7)
      END DO ! Loop IMODE=5,7
!
      IF (lhook) CALL dr_hook('UKCA_AGEING',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_AGEING
