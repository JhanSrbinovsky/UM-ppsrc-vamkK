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
!    Calculates condensation of condensable cpt vapours onto pre-existing
!    aerosol particles. Includes switch for using either Fuchs (1964) or
!    modified Fuchs and Sutugin (1971) calculation of CC.
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
      SUBROUTINE UKCA_CONDEN(NBOX,GC,ND,MD,MDT,                         &
       DTZ,DRYDP,WETDP,TSQRT,RHOA,AIRDM3,DELGC_COND,                    &
       IFUCHS,AGETERM1,BUD_AER_MAS,S_COND_S,PMID,T,S,AIRD,              &
       IDCMFP,ICONDIAM)
!----------------------------------------------------------------------
!
! Purpose
! -------
! Calculates condensation of condensable cpt vapours onto pre-existing
! aerosol particles. Includes switch for using either Fuchs (1964) or
! modified Fuchs and Sutugin (1971) calculation of CC.
!
! Parameters
! ----------
! SE_SOL : Sticking efficiency for   soluble modes [set to 1.0 as in M7]
! SE_INS : Sticking efficiency for insoluble modes [set to 0.3 as in M7]
!
! Inputs
! ------
! NBOX     : Number of grid boxes
! GC       : Condensable cpt number density (molecules cm-3)
! ND       : Aerosol ptcl no. concentration (ptcls per cc)
! MD       : Component median aerosol mass (molecules per ptcl)
! MDT      : Total median aerosol mass (molecules per ptcl)
! DTZ      : Time Step for nucl/cond competition (s)
! DRYDP    : Avg dry diameter for each aerosol mode (m)
! WETDP    : Avg wet diameter for each aerosol mode (m)
! TSQRT    : Square-root of mid-level temperature (K)
! RHOA     : Air density (kg/m3)
! AIRDM3   : Number density of air (per m3)
! IFUCHS   : Switch for Fuchs (1964) or Fuchs-Sutugin (1971) for CC
! BUD_AER_MAS : Aerosol mass fluxes (molecules/cc/DTC)
! S_COND_S :  Condensation sink
! PMID        : Centre level pressure (Pa)
! T           : Centre level temperature (K)
! S        : Specific humidity (kg/kg)
! AIRD     : Number density of air (per cm3)
! IDCMFP      : Switch : diffusion/mfp  (1=as Gbin v1, 2=as Gbin v1_1)
! ICONDIAM : Switch : wet diam in UKCA_CONDEN (1=g.mean,2=condiam.)
!
! Outputs:
! -------
! MD    : Avg aerosol ptcl mass in mode (molecules per ptcl)
! GC    : Condensable cpt number density (molecules per cc)
! AGETERM1: stores mass of soluble material which has condensed onto
!           each of the insoluble modes for use in UKCA_AGEING
!           (molecules per cc)
! MDT      : Total median aerosol mass (molecules per ptcl)
! DELGC_COND : Change in vapour conc due to cond (molecules cpt/cm3)
! BUD_AER_MAS : Aerosol mass fluxes (molecules/cc/DTC)
!
! Local variables:
! ---------------
! DMOL    : Molecular diameter of condensable cpt (m)
! CC      : Conden. coeff. for condensable cpt onto particle (m^3s^-1)
! RP      : Radius of aerosol particle (m)
! SE      : Sticking efficiency (accomodation coeff)
! SE_SOL  : Sticking efficiency (accomodation coeff) for soluble mode
! SE_INS  : Sticking efficiency (accomodation coeff) for insoluble mode
! MMCG    : Molar mass of condensing gas (kg/mole)
! NC      : Product of number conc and condensation coefficient
! SUMNC   : Sum of NC over all modes
! DELTAMS : Mass of condensing gas taken up by this   soluble mode
! DELTAMI : Mass of condensing gas taken up by this insoluble mode
! DELTAM  : Mass of condensing gas taken up by both modes (-->soluble)
!    n.b. DELTAMS,DELTAMI,DELTAM all in molecules per cc.
! MASK1-3 : Logical array to define regions of domain to work on
! VPKEL    : vapour pressure of H2SO4 includes kelvin effect (Pa)
! GC_EQ    : equilibrium vapour pressure (mol cm-3)
! PTCLCONC : Product of MDH2SO4 and ND - aerosol conc in mol cm-3
! SURFTEN : Surface tension of H2SO4 aerosol = (0.0728 J m-2)
!
! References
! ----------
! Gong et al, JGR, 108(D1), 4007, doi:10.1029/2001JD002002, 2003.
! Raes et al, J. Aerosol Sci., 23 (7), pp. 759--771, 1992.
! Fuchs & Sutugin, Topics in aerosol research, 1971.
! Fuchs, "Mechanics of aerosols", Pergamon Press, 1964.
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! RA       : Dry air gas constant = 287.05 Jkg^-1 K^-1
! RR       : Universal gas constant = 8.314 Jmol^-1 K^-1
! PPI      : 3.1415927...........
! AVC      : Avogadros constant (mol-1)
! CONC_EPS : Threshold for condensable conc (molecules per cc)
! ZBOLTZ   : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES   : Number of possible aerosol modes
! NCP      : Number of possible aerosol components
! MODE     : Defines which modes are set
! COMPONENT: Defines which cpts are allowed in each mode
! CONDENSABLE : Logical variable defining which cpts are condensable
! MODESOL  : Defines whether the mode is soluble or not (=1 or 0)
! MM       : Molar masses of components (kg/mole)
! DIMEN    : Molecular diamters of condensable components (m)
! NUM_EPS  : Value of NEWN below which do not recalculate MD (per cc)
!                                             or carry out process
! CP_SU    : Index of component in which sulfate is stored
! CP_OC    : Index of component in which 1st OC cpt is stored
! CP_SO    : Index of component in which 2nd OC cpt is stored
! SIGMAG   : Geometric standard deviation of mode
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! MH2SO4   : Index of MM_GAS, WTRATC and S0G for H2SO4
! NCHEMG   : Number of gas phase tracers for gas phase chemistry scheme
! MM_GAS   : Array of molar masses for gas phase species (kg/mol)
! DIMEN    : Molecular diamters of condensable components (m)
! Various indices for budget terms in BUD_AER_MAS
!
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS,       ONLY: ra, rr, ppi, avc, conc_eps, zboltz
      USE UKCA_MODE_SETUP,      ONLY: nmodes, ncp, mode, component,     &
                                      sigmag, modesol, mm, num_eps,     &
                                      cp_su, cp_oc, cp_so
      USE UKCA_SETUP_INDICES
      USE yomhook,              ONLY: lhook, dr_hook
      USE parkind1,             ONLY: jprb, jpim
      IMPLICIT NONE
!
! Subroutine interface:
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: IFUCHS
      INTEGER, INTENT(IN) :: IDCMFP
      INTEGER, INTENT(IN) :: ICONDIAM
      REAL, INTENT(IN)    :: ND(NBOX,NMODES)
      REAL, INTENT(IN)    :: TSQRT(NBOX)
      REAL, INTENT(IN)    :: RHOA(NBOX)
      REAL, INTENT(IN)    :: AIRDM3(NBOX)
      REAL, INTENT(IN)    :: DTZ
      REAL, INTENT(IN)    :: WETDP(NBOX,NMODES)
      REAL, INTENT(IN)    :: DRYDP(NBOX,NMODES)
      REAL, INTENT(IN)    :: PMID(NBOX)
      REAL, INTENT(IN)    :: T(NBOX)
      REAL, INTENT(IN)    :: S(NBOX)
      REAL, INTENT(IN)    :: AIRD(NBOX)
      REAL, INTENT(INOUT) :: MD(NBOX,NMODES,NCP)
      REAL, INTENT(INOUT) :: MDT(NBOX,NMODES)
      REAL, INTENT(INOUT) :: GC(NBOX,NCHEMG)
      REAL, INTENT(INOUT) :: BUD_AER_MAS(NBOX,0:NBUDAER)
      REAL, INTENT(OUT)   :: DELGC_COND(NBOX,NCHEMG)
      REAL, INTENT(OUT)   :: AGETERM1(NBOX,3,NCHEMG)
      REAL, INTENT(OUT)   :: S_COND_S(NBOX)
!
! Local variables
      INTEGER :: ICP
      INTEGER :: IMODE
      INTEGER :: JV
      INTEGER :: JL
      LOGICAL :: MASK1(NBOX)
      LOGICAL :: MASK2(NBOX)
      LOGICAL :: MASK3(NBOX)
      LOGICAL :: MASK3I(NBOX)
      LOGICAL :: MASK4(NBOX)
      LOGICAL :: MASK2E(NBOX)
      LOGICAL :: MASK3E(NBOX)
      REAL    :: DMOL
      REAL    :: MMCG
      REAL    :: CC(NBOX)
      REAL    :: RP(NBOX)
      REAL    :: SUMNC(NBOX)
      REAL    :: NC(NBOX,NMODES)
      REAL    :: DELTAM(NBOX)
      REAL    :: DELTAMS(NBOX)
      REAL    :: DELTAMI(NBOX)
      REAL    :: SE
      REAL    :: Y2
      REAL    :: AA
      REAL, PARAMETER :: SE_SOL=1.0
!!      REAL, PARAMETER :: SE_INS=0.3
      REAL, PARAMETER :: SE_INS=1.0
!
      REAL, PARAMETER :: SURFTEN=0.0728
!
      REAL :: SINKARR(NBOX)
      REAL :: DIFVOL
!
!! added below for stratosphere
      REAL :: VPKEL(NBOX)
      REAL :: DIFKEL(NBOX)
      REAL :: PTCLCONC(NBOX)
      REAL :: DELGC_EVAP(NBOX,NCHEMG)
      REAL :: DELGC_EVAP_MODE(NBOX)
      REAL :: FRAC_EVAP_MODE(NBOX,NMODES)
      REAL :: WTS(NBOX)
      REAL :: RHOSOL_STRAT(NBOX)
      REAL :: GC_EQ(NBOX)
      REAL :: DIFF_GC_MODE(NBOX)
      REAL :: MAXEVAP(NBOX)
!

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_CONDEN',zhook_in,zhook_handle)

      AGETERM1(:,:,:)=0.0
      S_COND_S(:)=0.0      
      SINKARR(:)=0.0
!
      DO JV=1,NCHEMG
       IF(CONDENSABLE(JV)) THEN

! .. Set component into which component will condense
        ICP=CONDENSABLE_CHOICE(JV)

        DMOL=DIMEN(JV)
        MMCG=MM_GAS(JV)
        DELGC_COND(:,JV)=0.0

        MASK1(:) = GC(:,JV) > CONC_EPS

        SUMNC(:)=0.0
        DO IMODE=1,NMODES
         IF(MODE(IMODE)) THEN

          IF(ICONDIAM == 1) AA=0.0 ! as v1_gm4c, take g.m. number radius
          IF(ICONDIAM == 2) THEN
           IF(IMODE == 1) AA=2.0 ! continuum regime -- 2nd radial moment
           IF(IMODE == 2) AA=1.9
           IF(IMODE == 3) AA=1.5
           IF(IMODE == 4) AA=1.1 ! molecular regime -- 1st radial moment
           IF(IMODE == 5) AA=1.9
           IF(IMODE == 6) AA=1.5
           IF(IMODE == 7) AA=1.1
! .. these values are taken from Figure 1 Lehtinen et al (2003)
          END IF

!!          Y2=EXP(2.0*LOG(SIGMAG(IMODE))*LOG(SIGMAG(IMODE)))
!! original runs tried just setting equivalent to AA=2.0 (above)
          Y2=EXP(0.5*AA*AA*LOG(SIGMAG(IMODE))*LOG(SIGMAG(IMODE)))
! now use values of AA for each mode from Lehtinen et al (2003)

          NC(:,IMODE)=0.0
          MASK2(:) = MASK1(:) .AND.(ND(:,IMODE) > NUM_EPS(IMODE))

          RP(:)=WETDP(:,IMODE)*0.5*Y2 ! radius to use for condensation

          IF(MODESOL(IMODE) == 1) SE=SE_SOL
          IF(MODESOL(IMODE) == 0) SE=SE_INS

!!          IF(JV == MH2SO4  ) DIFVOL=DH2SO4 ! as H2SO4+H2O hydrate
!!          IF(JV == MSEC_ORG) DIFVOL=DSECOR ! as OH-a-pinene radical
          IF(JV == MH2SO4  ) DIFVOL=51.96  ! values from dist_data (BIN)
          IF(JV == MSEC_ORG) DIFVOL=204.14 ! values from dist_data (BIN)

! ..  Calculate change in condensable cpt conc (molecules cm^-3)
! DEPENDS ON: ukca_cond_coff_v
          CALL UKCA_COND_COFF_V(NBOX,MASK2,RP,TSQRT,AIRDM3,RHOA,        &
                 MMCG,SE,DMOL,IFUCHS,CC,SINKARR,PMID,T,DIFVOL,IDCMFP)
          WHERE(MASK2(:))
           NC(:,IMODE)=ND(:,IMODE)*CC(:)
           SUMNC(:)=SUMNC(:)+NC(:,IMODE)
          ENDWHERE
!!          IF(CONDENSABLE_CHOICE(JV).EQ.1) THEN ! if H2SO4
          IF(JV.EQ.MH2SO4) THEN ! if H2SO4
            S_COND_S(:)=S_COND_S(:)+ND(:,IMODE)*SINKARR(:)
          END IF

         END IF ! if mode is present
        END DO ! Over modes
!
        WHERE(MASK1(:))                                                 &
         DELGC_COND(:,JV)=GC(:,JV)*(1.0-EXP(-SUMNC(:)*DTZ))

! .. Update condensable cpt concentration (molecules cm^-3)
        MASK2(:) = MASK1(:) .AND. (DELGC_COND(:,JV) > CONC_EPS)
        WHERE(MASK2(:) .AND. DELGC_COND(:,JV) > GC(:,JV))
         DELGC_COND(:,JV)=DELGC_COND(:,JV)/GC(:,JV) ! make sure no -ves
        ENDWHERE
        WHERE(MASK2(:)) GC(:,JV)=GC(:,JV)-DELGC_COND(:,JV)

        DO IMODE=1,4 ! loop over sol modes (do cond sol -> ins here too)
         IF(MODE(IMODE)) THEN

!         Calculate increase in total & cpt masses in each soluble mode
          DELTAMS(:)=0.0
          DELTAMI(:)=0.0

          MASK3 (:) = MASK2(:) .AND.(ND(:,IMODE  ) > NUM_EPS(IMODE))
          MASK3I(:) = MASK2(:) .AND.(ND(:,IMODE+3) > NUM_EPS(IMODE))

          IF(IMODE == 1) THEN

           IF((ICP == CP_SU).AND.(NMASCONDSUNUCSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSUNUCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSUNUCSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_OC).AND.(NMASCONDOCNUCSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDOCNUCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDOCNUCSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_SO).AND.(NMASCONDSONUCSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSONUCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSONUCSOL)+DELTAMS(:)

            ENDWHERE
           END IF

          END IF ! if IMODE=1

          IF(IMODE == 2) THEN

           IF((ICP == CP_SU).AND.(NMASCONDSUAITSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSUAITSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSUAITSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_SU).AND.(NMASCONDSUAITINS > 0)) THEN
            WHERE(MASK3I(:))

             DELTAMI(:)=DELGC_COND(:,JV)*NC(:,IMODE+3)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSUAITINS)=                           &
             BUD_AER_MAS(:,NMASCONDSUAITINS)+DELTAMI(:)

             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)

            ENDWHERE
           END IF

           IF((ICP == CP_OC).AND.(NMASCONDOCAITSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDOCAITSOL)=                           &
             BUD_AER_MAS(:,NMASCONDOCAITSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_OC).AND.(NMASCONDOCAITINS > 0)) THEN
            WHERE(MASK3I(:))

             DELTAMI(:)=DELGC_COND(:,JV)*NC(:,IMODE+3)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDOCAITINS)=                           &
             BUD_AER_MAS(:,NMASCONDOCAITINS)+DELTAMI(:)

             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)

            ENDWHERE
           END IF

           IF((ICP == CP_SO).AND.(NMASCONDSOAITSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSOAITSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSOAITSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_SO).AND.(NMASCONDSOAITINS > 0)) THEN
            WHERE(MASK3I(:))

             DELTAMI(:)=DELGC_COND(:,JV)*NC(:,IMODE+3)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSOAITINS)=                           &
             BUD_AER_MAS(:,NMASCONDSOAITINS)+DELTAMI(:)

             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)

            ENDWHERE
           END IF

          END IF ! if IMODE=2

          IF(IMODE == 3) THEN

           IF((ICP == CP_SU).AND.(NMASCONDSUACCSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSUACCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSUACCSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_SU).AND.(NMASCONDSUACCINS > 0)) THEN
            WHERE(MASK3I(:))

             DELTAMI(:)=DELGC_COND(:,JV)*NC(:,IMODE+3)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSUACCINS)=                           &
             BUD_AER_MAS(:,NMASCONDSUACCINS)+DELTAMI(:)

             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)

            ENDWHERE
           END IF

           IF((ICP == CP_OC).AND.(NMASCONDOCACCSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDOCACCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDOCACCSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_OC).AND.(NMASCONDOCACCINS > 0)) THEN
            WHERE(MASK3I(:))

             DELTAMI(:)=DELGC_COND(:,JV)*NC(:,IMODE+3)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDOCACCINS)=                           &
             BUD_AER_MAS(:,NMASCONDOCACCINS)+DELTAMI(:)

             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)

            ENDWHERE
           END IF

           IF((ICP == CP_SO).AND.(NMASCONDSOACCSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSOACCSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSOACCSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_SO).AND.(NMASCONDSOACCINS > 0)) THEN
            WHERE(MASK3I(:))

             DELTAMI(:)=DELGC_COND(:,JV)*NC(:,IMODE+3)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSOACCINS)=                           &
             BUD_AER_MAS(:,NMASCONDSOACCINS)+DELTAMI(:)

             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)

            ENDWHERE
           END IF

          END IF ! if IMODE=3

          IF(IMODE == 4) THEN

           IF((ICP == CP_SU).AND.(NMASCONDSUCORSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSUCORSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSUCORSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_SU).AND.(NMASCONDSUCORINS > 0)) THEN
            WHERE(MASK3I(:))

             DELTAMI(:)=DELGC_COND(:,JV)*NC(:,IMODE+3)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSUCORINS)=                           &
             BUD_AER_MAS(:,NMASCONDSUCORINS)+DELTAMI(:)

             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)

            ENDWHERE
           END IF

           IF((ICP == CP_OC).AND.(NMASCONDOCCORSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDOCCORSOL)=                           &
             BUD_AER_MAS(:,NMASCONDOCCORSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_OC).AND.(NMASCONDOCCORINS > 0)) THEN
            WHERE(MASK3I(:))

             DELTAMI(:)=DELGC_COND(:,JV)*NC(:,IMODE+3)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDOCCORINS)=                           &
             BUD_AER_MAS(:,NMASCONDOCCORINS)+DELTAMI(:)

             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)

            ENDWHERE
           END IF

           IF((ICP == CP_SO).AND.(NMASCONDSOCORSOL > 0)) THEN
            WHERE(MASK3(:))

             DELTAMS(:)=DELGC_COND(:,JV)*NC(:,IMODE)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSOCORSOL)=                           &
             BUD_AER_MAS(:,NMASCONDSOCORSOL)+DELTAMS(:)

            ENDWHERE
           END IF

           IF((ICP == CP_SO).AND.(NMASCONDSOCORINS > 0)) THEN
            WHERE(MASK3I(:))

             DELTAMI(:)=DELGC_COND(:,JV)*NC(:,IMODE+3)/SUMNC(:)

             BUD_AER_MAS(:,NMASCONDSOCORINS)=                           &
             BUD_AER_MAS(:,NMASCONDSOCORINS)+DELTAMI(:)

             AGETERM1(:,IMODE-1,JV)=DELTAMI(:)

            ENDWHERE
           END IF

          END IF ! if IMODE=4

!         All mass condensed onto sol. & ins. goes to sol. mode
          DELTAM(:)=DELTAMS(:)+DELTAMI(:)

          WHERE(MASK3(:))
           MD(:,IMODE,ICP)=                                             &
            (MD(:,IMODE,ICP)*ND(:,IMODE)+DELTAMS(:))/ND(:,IMODE)
           MDT(:,IMODE)=                                                &
            (MDT(:,IMODE)*ND(:,IMODE)+DELTAMS(:))/ND(:,IMODE)
          ENDWHERE

!!          WHERE(MASK3I(:))
!!           MD(:,IMODE,ICP)=                                             &
!!            (MD(:,IMODE,ICP)*ND(:,IMODE)+DELTAMI(:))/ND(:,IMODE)
!!           MDT(:,IMODE)=                                                &
!!            (MDT(:,IMODE)*ND(:,IMODE)+DELTAMI(:))/ND(:,IMODE)
!!          ENDWHERE
!!
!! **** HERE HAVE COMMENTED OUT THE UPDATING OF THE SOLUBLE MODE
!! **** MD,MDT FOR THE CONDENSATION ONTO INSOLUBLE MODES BECAUSE IN
!! **** THE CASE WHERE THERE ARE NO PARTICLES IN THE AITKEN SOLUBLE
!! **** MODE BUT THERE ARE IN THE AITKEN INSOLUBLE MODE, THE UPDATING
!! **** OF THE SOLUBLE MODE MD,MDT WILL NOT BE POSSIBLE UNTIL THE
!! **** SOLUBLE MODE ND HAS BEEN UPDATED DUE TO THE TRANSFER OF AGED
!! **** PARTICLES. THIS IS DONE IN THE UKCA_AGEING ROUTINE -- SO THIS IS
!! **** WHERE THE UPDATE OF THE SOLUBLE MODE MD,MDT FOR THE CONDENSATION
!! **** ONTO THE SHOULD BE DONE ALSO.
!!
         END IF ! if mode is present

        END DO ! IMODE=1,4 (soluble modes)

       END IF ! IF CONDENSABLE(JV)

      END DO ! DO JV=1,NCHEMG

      IF (lhook) CALL dr_hook('UKCA_CONDEN',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CONDEN
