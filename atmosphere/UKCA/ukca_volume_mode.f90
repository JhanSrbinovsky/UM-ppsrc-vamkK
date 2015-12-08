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
!    Calculate (wet) volume corresponding to mid-pt particles in each
!    aerosol mode.
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
      SUBROUTINE UKCA_VOLUME_MODE(NBOX,ND,MD,MDT,                       &
       RH,WVOL,WETDP,RHOPAR,IWVOLMETHOD,                                &
       DVOL,DRYDP,MDWAT,PVOL,PVOL_WAT,VERBOSE,                          &
       T,PMID,S)
!----------------------------------------------------------------------
!
! Calculate (wet) volume corresponding to mid-pt particles in each mode.
! Two methods can be followed (as specified by IWVOLMETHOD:
!
! 1) Quick method --- soluble particle mass assumed to behave as SO4
!                     and wet volume calculated from Kohler equation
!                     ln(rh)=A/Dp - B/Dp^3 --- in sub-saturated
!                     conditions can be approximated as
!                     ln(rh)=-B/Dp^3, [Dp is wet diameter of particle]
!                     so  wetvol=(pi/6)Dp^3=-(pi/6)B/ln(rh)
!                     Seinfeld & Pandis have B=3.44e13.nu.(ms/Ms)
!                     which we set as B=3.44e13.nu.(MDSOL/AVC) where
!                     MDSOL is all soluble mass (per particle),
!                     # of dissociating ions, nu=3 and then
!                     wetvol=-8.974e-29*MDSOL/ln(rh)z
!                     Also add on dry volume of insoluble
!                     molecules to give overall (wet) ptcl volume.
!
! 2) Better method -- calculate the water-content and density of the
!                     aerosol using ZSR and water activity coeffs for
!                     H+,SO42-,Na+,Cl- from Jacobsen (pg 610).
!                     Complete dissociation of solute ions is assumed
!                     to give the electrolyte concentrations and
!                     associated water content associated with each.
!                     The dry volume of insoluble molecules is then
!                     added to give overall (wet) ptcl volume.
!                     Note that OC is assumed water-insoluble in the
!                     insoluble mode but is assumed to have aged
!                     chemically in the aerosol to become hygroscopic
!                     once transferred to the soluble distribution.
!                     To respresent this, in the ZSR calculation,
!                     the concentration of SO4 ions is incremented
!                     by FHYG_AOM*MD(:,:,CP_OC)/F_AO -- i.e. the aged
!                     OC is assumed to take up water at a fraction
!                     FHYG_AOM (set at 0.65) of SO4.
!
! In each case sphericity is assumed to give the ptcl wet radius.
!
! Purpose
! -------
! Calculate avg wet volume WVOL & wet diameter WETDP for each aerosol
! particle size mode from the relative humidity and avg. number of
! molecules per particle (MD) of each cpt in each mode.
!
! (N.b. Local rel. hum. values are corrected to lie between 0.1 & 0.9)
!
! Inputs
! ------
! NBOX     : Number of grid boxes
! ND       : Aerosol ptcl no. concentration (ptcls per cc)
! MD       : Component median aerosol mass (molecules per ptcl)
! MDT      : Total median aerosol mass (molecules per ptcl)
! RH       : Relative humidity (corrected to lie within range 0.1-0.9)
! IWVOLMETHOD: Chosen method for calculation of wet volume
! DVOL     : Dry volume of particle (m^3)
! DRYDP    : Dry diameter of particle (m)
! VERBOSE  : Switch for level of verbosity
!
! Outputs
! ------
! WVOL     : Avg wet volume of size mode (m3)
! WETDP    : Avg wet diameter of size mode (m)
! MDWAT    : Molecular concentration of water (molecules per particle)
! RHOPAR   : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
! PVOL     : Partial volumes of each component in each mode (m3)
! PVOL_WAT : Partial volume of water in each mode (m3)
!
! Local Variables
! ---------------
! CORRH    : Locally corrected RH (to lie between 0.1 and 0.9)
! FHYG_AOM : Hygroscopicity of aged organic (fraction relative to SO4)
! F_AO     : No. of "moles" of POM in 1 mole of aged organic species
! MM_AGE_ORG: Molar mass of aged organic species (kg/mol)
! MM_POM   : Molar mass of particulate organic matter [for CP_OC,CP_SO]
! IONS     : Logical indicating presence of ions
! CL       : Ion concentrations (moles/cm3 of air)
! MASK     : Mask where in domain to calculate values
! MDSOL    : Mass per particle (total over soluble cpts) (mlcls/ptcl)
! RHOSOL   : Density of particle solution [excl. insoluble cpts] (kg/m3)
! B        : Factor in solute term in Kohler equation =
!              BCONST*(no. of ions)*(solute mass)/(solute molar mass)
! DENOM    : Temporary variable calculating denominator of expression
! DENOM2   : Temporary variable calculating denominator of expression
! RHOTMP   : Temporary variable in calculation of particle density
! RHOTMP2  : Temporary variable in calculation of particle density
! WVOL_SOL : Contribution to wet volume from soluble components (m^3)
! WC       : Water content for aerosol (moles/cm3 of air)
! WDPCUB   : Cube of particle wet diameter (m^3)
! SIXOVRPIX: (6.0/pi)*{ 1.0/EXP((9/2)*LOG^2(SIGMA_G)) }
!
! Inputted by module UKCA_CONSTANTS
! ---------------------------------
! PPI      : 3.1415927....
! AVC      : Avogadro's constant (molecules per mole)
! CONC_EPS : Value of soluble material mass conc. below which
!            assume no soluble mass
! RHOW     : Density of water (=1000.0 kgm^-3)
! MMW      : Molecular mass of water (=0.018 kg/mol)
! BCONST   : Value of constant in B term of Kohler equation
! NU_H2SO4 : Number of dissociated ions for H2SO4
! CONVERT  : Conversion from micron^3 to m^3
! RHOSUL   : Mass density of a pure H2SO4 aerosol (kg per m^3)
!
! Inputted by module UKCA_MODE_SETUP
! ----------------------------------
! NMODES   : Number of possible aerosol modes
! NCP      : Number of possible aerosol components
! MODE     : Logical variable denoting where mode is defined
! COMPONENT: Logical variable denoting where cpt is defined
! SOLUBLE  : Logical variable defining which cpts are soluble
! MM       : Molar masses of components (kg per mole)
! NO_IONS  : Number of dissociating ions in solute (H2SO4=3,NaCl=2)
! MODESOL  : Defines which modes are soluble (integer)
! RHOCOMP  : Densities (dry) of each component (kg/m^3)
! MMID     : Avg mass of mode when rmed_g=exp(0.5*(lnr0+lnr1)) (ptcl^-1)
! X        : EXP((9/2)*LOG^2(SIGMA_G))
! NUM_EPS  : Value of NEWN below which don't recalculate MD (per cc)
!                                            or carry out process
! CP_SU    : Index of component containing SO4    component
! CP_OC    : Index of component containing 1st OC component
! CP_CL    : Index of component containing NaCl   component
! CP_SO    : Index of component containing 2nd OC component
!
! Inputted by module UKCA_SETUP_INDICES
! -------------------------------------
! NCHEMG   : Number of gas phase tracers for gas phase chemistry scheme
! MM_GAS   : Array of molar masses for gas phase species (kg/mol)
!
!--------------------------------------------------------------------

      USE UKCA_CONSTANTS,       ONLY: ppi, avc, conc_eps, rhow, mmw,    &
                                      bconst, nu_h2so4, convert, rhosul
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      USE yomhook,              ONLY: lhook, dr_hook
      USE parkind1,             ONLY: jprb, jpim
      USE ereport_mod,          ONLY: ereport
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: IWVOLMETHOD
      INTEGER, INTENT(IN) :: VERBOSE
      REAL, INTENT(IN)    :: ND(NBOX,NMODES)
      REAL, INTENT(IN)    :: MD(NBOX,NMODES,NCP)
      REAL, INTENT(IN)    :: MDT(NBOX,NMODES)
      REAL, INTENT(IN)    :: RH(NBOX)
      REAL, INTENT(IN)    :: DVOL(NBOX,NMODES)
      REAL, INTENT(IN)    :: DRYDP(NBOX,NMODES)
      REAL, INTENT(IN)    :: T(NBOX)
      REAL, INTENT(IN)    :: PMID(NBOX)
      REAL, INTENT(IN)    :: S(NBOX)
      REAL, INTENT(OUT)   :: MDWAT(NBOX,NMODES)
      REAL, INTENT(OUT)   :: WVOL(NBOX,NMODES)
      REAL, INTENT(OUT)   :: WETDP(NBOX,NMODES)
      REAL, INTENT(OUT)   :: RHOPAR(NBOX,NMODES)
      REAL, INTENT(OUT)   :: PVOL(NBOX,NMODES,NCP)
      REAL, INTENT(OUT)   :: PVOL_WAT(NBOX,NMODES)

! Local variables
      INTEGER :: errcode                ! Variable passed to ereport
      INTEGER :: I
      INTEGER :: IMODE
      INTEGER :: ICP
      INTEGER :: JJL(1)
      INTEGER :: IIMODE
      INTEGER :: IICP
      INTEGER :: IERR(NBOX)
      LOGICAL :: MASK(NBOX)
      LOGICAL :: MASK_SOL(NBOX)
      LOGICAL :: MASK_NOSOL(NBOX)
      LOGICAL :: MASK_ERROR(NBOX)
      LOGICAL :: IONS(NBOX,-NANION:NCATION) !ION PRESENCE SWITCHES
      REAL    :: CORRH(NBOX)
      REAL    :: B(NBOX)
      REAL    :: MDSOL(NBOX)
      REAL    :: RHOSOL(NBOX)
      REAL    :: WC(NBOX)
      REAL    :: WVOL_SOL(NBOX)
      REAL    :: RHOTMP(NBOX)
      REAL    :: RHOTMP2(NBOX)
      REAL    :: DENOM(NBOX)
      REAL    :: DENOM2(NBOX)
      REAL    :: CBRT
      REAL    :: WDPCUB(NBOX)
      REAL    :: PIOVRSIX
      REAL    :: SIXOVRPIX(NMODES)
      REAL    :: MDMM(NBOX,NCP)
      REAL    :: MM_OVRAVC(NCP)
      REAL    :: MMWOVRAVC
      REAL    :: MM_OVRAVCRHOCP(NCP)
      REAL    :: MMWOVRAVCRHOW
      REAL    :: MM_RHOCP(NCP)
      REAL    :: MMWRHOW
      REAL    :: F_AO
      REAL    :: CL(NBOX,-NANION:NCATION) !ION CONCS (MOL/CC OF AIR)
!
! extra local vars for stratospheric calculation
      REAL    :: RP(NBOX)
      REAL    :: WTS(NBOX)
      REAL    :: RHOSOL_STRAT(NBOX)
      REAL    :: VPKEL(NBOX)
      REAL    :: MASSH2SO4KG(NBOX)
      REAL    :: MASSWATERKG(NBOX)

      CHARACTER(LEN=72) :: cmessage
      REAL, PARAMETER :: FHYG_AOM=0.65
      REAL, PARAMETER :: MM_AGE_ORG=0.150
      REAL, PARAMETER :: MM_POM=0.0168
      REAL, PARAMETER :: PUTLS=1.5e4 ! set max pressure for UTLS region to 150hPa
!

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_VOLUME_MODE',zhook_in,zhook_handle)

      IF((IWVOLMETHOD /= 1).AND.(IWVOLMETHOD /= 2)) THEN
       cmessage=' IWVOLMETHOD NOT 1 OR 2'
       errcode=1
       CALL EREPORT('UKCA_VOLUME_MODE',errcode,cmessage)
      END IF

!at this point in the code, the value of RP does not matter
      RP(:)=100.0E-9 ! dummy value
! DEPENDS ON: ukca_vapour
      CALL UKCA_VAPOUR(NBOX,T,PMID,S,RP,VPKEL,WTS,RHOSOL_STRAT)
! Note that here we only want to get WTS and RHOSOL_STRAT which
!  are independent of particle size and composition -- so we only
!  need to call this once here -0- and VPKEL is not used.

      PIOVRSIX=PPI/6.0
      SIXOVRPIX(:)=1.0/(X(:)*PIOVRSIX)
      MM_OVRAVC(:)=MM(:)/AVC
      MMWOVRAVC   =MMW/AVC
      MM_OVRAVCRHOCP(:)=MM_OVRAVC(:)/RHOCOMP(:)
      MMWOVRAVCRHOW    =MMWOVRAVC/RHOW
      MM_RHOCP(:)=MM(:)*RHOCOMP(:)
      MMWRHOW    =MMW*RHOW

      DENOM(:)=0.0         ! define on all points to avoid error
      DENOM2(:)=0.0

! Correct relative humidities to lie within the range of 10-90%
      CORRH(:)=RH(:)
      WHERE (CORRH(:) > 0.9) CORRH(:)=0.9
      WHERE (CORRH(:) < 0.1) CORRH(:)=0.1

      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN

        MASK(:) =(ND(:,IMODE) > NUM_EPS(IMODE))

        IF(MODESOL(IMODE) == 1) THEN

! .. first calculate the total mass per particle over all soluble cpts
         MDSOL(:)=0.0
         DO ICP=1,NCP
          IF(COMPONENT(IMODE,ICP)) THEN
           IF(SOLUBLE(ICP)) THEN
! .. calculate total mass (in molecules per particle) of all solutes
            WHERE(MASK(:)) MDSOL(:)=MDSOL(:)+MD(:,IMODE,ICP)
           END IF
          END IF
         END DO

! .. make mask "MASK_SOL"   for where ND>NUM_EPS and some soluble mass
! .. make mask "MASK_NOSOL" for where ND>NUM_EPS and no   soluble mass
         MASK_SOL  (:) = (MASK(:) .AND. (MDSOL(:) > 0.0))
         MASK_NOSOL(:) = (MASK(:) .AND. (MDSOL(:) == 0.0))

         IF(IWVOLMETHOD == 1) THEN
! .. assume all soluble components are H2SO4 & uses Kohler theory

!     Kohler equation is lnS=A/Dp - B/Dp^3
!     In sub-saturated environment, B/Dp^3 > A/Dp
!     Reasonable assumption B/Dp^3 >> A/Dp for all but smallest ptcls
!     Then,   Dp=(-B/ln(rh))^(1/3)
!       or  volp=-(pi/6)*B/ln(rh)
!
!     (CONVERT converts from micron^3 to m^3 in B term, see S&P pg 787).

          WHERE(MASK(:))
           B(:)=BCONST*NU_H2SO4*MDSOL(:)/AVC
           WVOL_SOL(:)=-PIOVRSIX*B(:)/LOG(CORRH(:))*CONVERT
! .. Calculate the volume of the solution WVOL_SOL
          ENDWHERE

          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            IF(.NOT.SOLUBLE(ICP)) THEN
             WHERE(MASK(:))
              PVOL(:,IMODE,ICP)=MD(:,IMODE,ICP)*MM_OVRAVCRHOCP(ICP)
! .. Calculate partial volumes of the insoluble components
              WVOL(:,IMODE)=WVOL_SOL(:)+PVOL(:,IMODE,ICP)
! .. Wet volume set as solution volume plus insoluble cpt partial volumes
             ENDWHERE
            END IF
           END IF
          END DO

          WHERE(MASK(:))
           MDWAT(:,IMODE)=(WVOL(:,IMODE)-DVOL(:,IMODE))/MMWOVRAVCRHOW
! .. Water content calculated from difference between dry & wet volumes
           DENOM (:)=MDWAT(:,IMODE)*MMW
           RHOTMP(:)=MDWAT(:,IMODE)*MMWRHOW
! .. Set initial value of DENOM & RHOTMP for density calculation
! .. note only need to define DENOM,RHOTMP where MASK(:)
          ELSEWHERE
! below is for gridboxes with ND<NUM_EPS
! -- define DENOM,RHOTMP,DENOM2,RHOTMP2=0 (they will not be used)
           MDWAT(:,IMODE)=0.0
           DENOM(:)=0.0
           RHOTMP(:)=0.0
           DENOM2(:)=0.0
           RHOTMP2(:)=0.0
          ENDWHERE
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            WHERE(MASK(:))
             MDMM(:,ICP)=MD(:,IMODE,ICP)*MM(ICP)
             RHOTMP(:)=RHOTMP(:)+MD(:,IMODE,ICP)*MM_RHOCP(ICP)
             DENOM (:)=DENOM (:)+MDMM(:,ICP)
! .. increment RHOTMP & DENOM by each component for density calculation
            ENDWHERE
           END IF
          END DO
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            IF(SOLUBLE(ICP)) THEN
             WHERE(MASK_SOL(:))
              PVOL(:,IMODE,ICP)=MDMM(:,ICP)*WVOL_SOL(:)/DENOM(:)
! .. Calculate partial volume of soluble components as fraction of
! .. solution volume according to fraction of total mass (including water)
             ELSEWHERE
              PVOL(:,IMODE,ICP)=DVOL(:,IMODE)*MFRAC_0(IMODE,ICP)
! where ND<NUM_EPS, over-write PVOL with default value
             ENDWHERE
             WHERE(MASK_NOSOL(:))
              PVOL(:,IMODE,ICP)=0.0
! .. In case where soluble mode has some insoluble cpt mass but no soluble
! .. component mass, DENOM=0.0 so just set PVOL for soluble cpt =0.0 here
             ENDWHERE
            END IF
           END IF
          END DO
          WHERE(MASK_SOL(:))
           PVOL_WAT(:,IMODE)=MDWAT(:,IMODE)*MMW*WVOL_SOL(:)/DENOM(:)
! .. Calculate partial volume of water as fraction of
! .. solution volume according to fraction of total mass (including water)
           RHOPAR(:,IMODE)=RHOTMP(:)/DENOM(:)
! .. RHOPAR total particle density [including H2O & insoluble cpts] (kgm^-3)
! .. Calculated according to mass-weighted average of cpt & water densities
          ELSEWHERE
           PVOL_WAT(:,IMODE)=0.0
           WVOL(:,IMODE)=DVOL(:,IMODE)
           RHOPAR(:,IMODE)=RHOSUL
! where ND<NUM_EPS, set PVOL_WAT=0.0, set WVOL=DVOL and RHOPAR=RHOSUL
          ENDWHERE

         END IF ! if IWVOLMETHOD = 1

         F_AO=MM_AGE_ORG/MM_POM

         IF(IWVOLMETHOD == 2) THEN

! .. Use composition information to calculate water uptake by
! .. each component according to ZSR using water activity data from
! .. Jacobsen page 610 (Table B.10) for binary electrolyte molalities

!***********************************************************************
!**  Liquid Phase Species:
!**
!**  **Cations**           **Anions**             **Neutrals**
!**  1: H                  -1: HSO4               0: H2O
!**  2: NH4                -2: SO4
!**  3: Na                 -3: NO3
!**                        -4: Cl
!***********************************************************************

          MASK(:) =(ND(:,IMODE) > NUM_EPS(IMODE))

          DO I=-NANION,NCATION
           CL(:,I)=0.0 ! set all concentrations to zero initially
          END DO

          IF(COMPONENT(IMODE,CP_SU)) THEN ! assume all H2SO4 --> SO4
           CL(:,-2)=MD(:,IMODE,CP_SU)/AVC   ! [SO4] in moles/cc (air)
          END IF

          IF(COMPONENT(IMODE,CP_SO)) THEN
           CL(:,-2)=CL(:,-2)+(FHYG_AOM/AVC)*(MD(:,IMODE,CP_SO)/F_AO)
! .. Increment concentration of SO4 ions to represent the
! .. presence of hygroscopic aged organic aerosol mass in CP_SO.
! .. Assume it has uptake behaviour at fraction FHYG_AOM of SO4.
! .. Need to divide by F_AO because MD of CP_SO is in "moles" of POM
! .. whereas CL needs to be in moles of aged organic (MM=0.15kg/mol).
          END IF

          IF(COMPONENT(IMODE,CP_OC)) THEN
           CL(:,-2)=CL(:,-2)+(FHYG_AOM/AVC)*(MD(:,IMODE,CP_OC)/F_AO)
! .. Increment concentration of SO4 ions to represent the
! .. presence of hygroscopic aged organic aerosol mass in CP_OC.
! .. Assume it has uptake behaviour at fraction FHYG_AOM of SO4.
! .. Need to divide by F_AO because MD of CP_OC is in "moles" of POM
! .. whereas CL needs to be in moles of aged organic (MM=0.15kg/mol).
!
! .. This effectively says that by the time the primary carbonaceous
! .. aerosol has been microphysically aged to the soluble mode, the
! .. organic component has been chemically aged to become hygroscopic.
! .. (it is assumed to be water-insoluble in the insoluble mode)
!
          END IF

          IF(COMPONENT(IMODE,CP_CL)) THEN ! assume complete dissociation
           CL(:,3)=MD(:,IMODE,CP_CL)/AVC  ! [Na] in moles per cc (air)
           CL(:,-4)=MD(:,IMODE,CP_CL)/AVC ! [Cl] in moles per cc (air)
          END IF

! SET H+ FOR CHARGE BALANCE  -- CL(1) is [H] in moles per cc(air)
          CL(:,1)=MAX((2.0*CL(:,-2)+CL(:,-1)+CL(:,-3)+CL(:,-4)          &
                       -CL(:,2)-CL(:,3)),0.0)

          DO I=-NANION,NCATION
           IONS(:,I)=(CL(:,I) > 0.)
          END DO

! DEPENDS ON: ukca_water_content_v
          CALL UKCA_WATER_CONTENT_V(NBOX,MASK,CL,CORRH,IONS,WC)

          WHERE(MASK(:)) MDWAT(:,IMODE)=WC(:)*AVC
!
! over-write MDWAT in strat with value from WTS retured from UKCA_VAPOUR
          WHERE(MASK(:).AND.(PMID(:).LT.PUTLS)) ! P<PUTLS
           MASSH2SO4KG(:)=MD(:,IMODE,CP_SU)*MM(CP_SU)/AVC ! in kg/particle
           MASSWATERKG(:)=(100.0/WTS(:)-1.0)*MASSH2SO4KG(:) ! kg/ptcl
           MDWAT(:,IMODE)=MASSWATERKG(:)/MMWOVRAVC ! in molecules per particle
          ENDWHERE
!
          WHERE(MASK(:))
! .. calculate solution density (avg over each cpt mass contribution)
           RHOTMP (:)=MDWAT(:,IMODE)*MMWRHOW
           RHOTMP2(:)=RHOTMP(:)
           DENOM (:)=MDWAT(:,IMODE)*MMW
           DENOM2(:)=DENOM(:)
! .. note only need to define DENOM,DENOM2,RHOTMP,RHOTMP2 where MASK(:)
          ELSEWHERE
           MDWAT(:,IMODE)=0.0
! .. set MDWAT to be zero where ND<NUM_EPS
           RHOTMP(:)=0.0
           DENOM (:)=0.0
           RHOTMP2(:)=0.0
           DENOM2 (:)=0.0
! .. also set RHOTMP,DENOM,RHOTMP2,DENOM2=0 (they will not be used)
          ENDWHERE

          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            IF(SOLUBLE(ICP)) THEN
             WHERE(MASK(:))
              RHOTMP(:)=RHOTMP(:)+MD(:,IMODE,ICP)*MM_RHOCP(ICP)
              DENOM (:)=DENOM (:)+MD(:,IMODE,ICP)*MM(ICP)
             ENDWHERE
            END IF
            WHERE(MASK(:))
             RHOTMP2(:)=RHOTMP2(:)+MD(:,IMODE,ICP)*MM_RHOCP(ICP)
             DENOM2 (:)=DENOM2 (:)+MD(:,IMODE,ICP)*MM(ICP)
            ENDWHERE
           END IF
          END DO

! .. RHOSOL is density of ptcl solution [exclud. insoluble cpts] (kgm/3)
!
!  .. section below added to check for DENOM <= 0.0 or DENOM2 <= 0.0

          MASK_ERROR(:)= (MASK_SOL(:) .AND. (DENOM(:) <= 0.0))
          IERR(:)=0
          WHERE(MASK_ERROR(:))
           IERR(:)=1
          ENDWHERE
          IF(SUM(IERR(:)) > 0) THEN ! error (DENOM<=0 when some soluble cpt)
           DO I=1,NBOX
            IF((MASK(I)).AND.(DENOM(I) <= 0.0)) THEN
             cmessage = 'DENOM(I) <= 0.0'
             WRITE(6,*) cmessage
             WRITE(6,*) 'I,RHOTMP(I),DENOM(I),MASK(I) =',               &
                                       I,RHOTMP(I),DENOM(I),MASK(I)
             WRITE(6,*) 'IMODE,MODESOL(IMODE),ND(I,IMODE)',             &
                        'NUM_EPS(IMODE)=',                              &
                        IMODE,MODESOL(IMODE),ND(I,IMODE),NUM_EPS(IMODE)
             DO ICP=1,NCP
              IF(COMPONENT(IMODE,ICP)) THEN
               WRITE(6,*)'ICP,MM(ICP),MM_RHOCP(ICP),MD(I,IMODE,ICP)=',  &
                       ICP,MM(ICP),MM_RHOCP(ICP),MD(I,IMODE,ICP)
              ENDIF
             ENDDO
             WRITE(6,*)'MDWAT(I,IMODE),MMW=',MDWAT(I,IMODE),MMW
             errcode=1

             CALL EREPORT('UKCA_VOLUME_MODE',errcode,cmessage)
            END IF ! if DENOM(I) <= 0.0

            IF((MASK(I)).AND.(DENOM2(I) <= 0.0)) THEN
             cmessage = 'DENOM2(I) <= 0.0'
             WRITE(6,*) 'I,RHOTMP2(I),DENOM2(I),MASK(I)=',              &
                         I,RHOTMP2(I),DENOM2(I),MASK(I)
             WRITE(6,*) 'MODE,MODESOL(IMODE)=',IMODE,MODESOL(IMODE)
             DO ICP=1,NCP
              IF(COMPONENT(IMODE,ICP)) THEN
               WRITE(6,*)'ICP,MM(ICP),MM_RHOCP(ICP),MD(I,IMODE,ICP)=',  &
                       ICP,MM(ICP),MM_RHOCP(ICP),MD(I,IMODE,ICP)
              END IF
             END DO
             WRITE(6,*)'MDWAT(I,IMODE),MMW=',MDWAT(I,IMODE),MMW
             errcode=1

             CALL EREPORT('UKCA_VOLUME_MODE',errcode,cmessage)
            END IF ! if DENOM2(I) <= 0.0
           END DO ! loop over NBOX (I)
          END IF

! .. where ND>NUM_EPS and there is some soluble material
          WHERE(MASK_SOL(:))
           RHOSOL(:)=RHOTMP(:)/DENOM(:)
! .. only need to define RHOSOL where ND>NUM_EPS and MDSOL>0 (MASK_SOL)
          ENDWHERE
!
!  over-write RHOSOL in strat with RHOSOL_STRAT from UKCA_VAPOUR
          WHERE(MASK_SOL(:).AND.(PMID(:).LT.PUTLS)) ! P<PUTLS
           RHOSOL(:)=RHOSOL_STRAT(:)
          ENDWHERE
!
          WVOL(:,IMODE)=0.0        ! initialise wet volume to zero
          DO ICP=1,NCP
           IF(COMPONENT(IMODE,ICP)) THEN
            IF(SOLUBLE(ICP)) THEN
             WHERE(MASK_SOL(:))
              PVOL(:,IMODE,ICP)=MD(:,IMODE,ICP)*MM_OVRAVC(ICP)/RHOSOL(:)
! .. for soluble cpts set PVOL according to cpt mass and solution density
              WVOL(:,IMODE)=WVOL(:,IMODE)+PVOL(:,IMODE,ICP)
             ELSEWHERE
              PVOL(:,IMODE,ICP)=DVOL(:,IMODE)*MFRAC_0(IMODE,ICP)
! where ND<NUM_EPS, over-write PVOL with default value
             ENDWHERE
             WHERE(MASK_NOSOL(:))
              PVOL(:,IMODE,ICP)=0.0
! .. In case where soluble mode has some insoluble cpt mass but no soluble
! .. component mass, DENOM=0.0 so just set PVOL for soluble cpt =0.0 here
             ENDWHERE
            ELSE
             WHERE(MASK(:))
              PVOL(:,IMODE,ICP)=MD(:,IMODE,ICP)*MM_OVRAVCRHOCP(ICP)
! .. for insoluble cpts set PVOL according to cpt mass and cpt density
              WVOL(:,IMODE)=WVOL(:,IMODE)+PVOL(:,IMODE,ICP)
             ELSEWHERE
              PVOL(:,IMODE,ICP)=DVOL(:,IMODE)*MFRAC_0(IMODE,ICP)
! where ND<NUM_EPS, over-write PVOL with default value
             ENDWHERE
            END IF
! .. initially set WVOL to sum of partial volumes of each component
           END IF
          END DO

          WHERE(MASK_SOL(:))
           PVOL_WAT(:,IMODE)=MDWAT(:,IMODE)*MMWOVRAVC/RHOSOL(:)
! .. for water, set PVOL according to water mass and solution density
           WVOL(:,IMODE)=WVOL(:,IMODE)+PVOL_WAT(:,IMODE)
! .. add on partial volume of water to give total wet volume
           RHOPAR(:,IMODE)=RHOTMP2(:)/DENOM2(:)
! .. RHOPAR is total particle density [incl H2O & insoluble cpts] (kgm^-3)
! .. calculated as mass-weighted mean of component & water densities
          ELSEWHERE
           PVOL_WAT(:,IMODE)=0.0
           WVOL(:,IMODE)=DVOL(:,IMODE)
           RHOPAR(:,IMODE)=RHOSUL
! where ND<NUM_EPS, set PVOL_WAT=0.0, set WVOL=DVOL and RHOPAR=RHOSUL
          ENDWHERE

         END IF ! if IWVOLMETHOD=2

        ELSE  ! if mode not soluble

! where mode not soluble, set PVOL_WAT=0.0, set WVOL=DVOL and MDWAT=0.0
         PVOL_WAT(:,IMODE)=0.0
         WVOL(:,IMODE)=DVOL(:,IMODE)
         MDWAT(:,IMODE)=0.0

         DO ICP=1,NCP
          IF(COMPONENT(IMODE,ICP)) THEN
           PVOL(:,IMODE,ICP)=MD(:,IMODE,ICP)*MM_OVRAVCRHOCP(ICP)
! .. all cpts insoluble in insoluble modes,
! .. so set PVOL according to cpt mass & cpt density
          END IF
         END DO

         RHOTMP2(:)=0.0
         DENOM2 (:)=0.0
         DO ICP=1,NCP
          IF(COMPONENT(IMODE,ICP)) THEN
           WHERE(MASK(:))
            RHOTMP2(:)=RHOTMP2(:)+MD(:,IMODE,ICP)*MM_RHOCP(ICP)
            DENOM2 (:)=DENOM2 (:)+MD(:,IMODE,ICP)*MM(ICP)
           ENDWHERE
          END IF
         END DO

         WHERE(MASK(:))
          RHOPAR(:,IMODE)=RHOTMP2(:)/DENOM2(:)
! .. RHOPAR is total particle density [here only insoluble cpts] (kgm^-3)
! .. calculated as mass-weighted mean of component densities
         ELSEWHERE
          RHOPAR(:,IMODE)=RHOSUL
! .. if insignificant # of particles, set density to default value
         ENDWHERE

        END IF ! if mode is soluble

       ELSE

! .. set PVOL_WAT=MDWAT=0.0 and WVOL=DVOL if mode not present
        PVOL_WAT(:,IMODE)=0.0
        WVOL(:,IMODE)=DVOL(:,IMODE)
        MDWAT(:,IMODE)=0.0

! .. set density to default value if mode not present
        RHOPAR(:,IMODE)=RHOSUL
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          PVOL(:,IMODE,ICP)=DVOL(:,IMODE)*MFRAC_0(IMODE,ICP)
! .. set partial volume to default value if mode not present
         END IF
        END DO
       END IF ! if(mode(imode))
      END DO

      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        WDPCUB(:)=SIXOVRPIX(IMODE)*WVOL(:,IMODE)
        WETDP(:,IMODE)=WDPCUB(:)**(1.0/3.0)
       ELSE
        WETDP(:,IMODE)=DRYDP(:,IMODE)
       END IF  ! if(mode(imode))
      END DO ! loop over modes

      DO IMODE=1,NMODES
       IF((minval(wetdp (:,IMODE)) <= 0.0).OR.                          &
          (minval(drydp (:,IMODE)) <= 0.0).OR.                          &
          (minval(wvol  (:,IMODE)) <= 0.0).OR.                          &
          (minval(dvol  (:,IMODE)) <= 0.0).OR.                          &
          (minval(rhopar(:,IMODE)) <= 0.0)) THEN
        cmessage='Minval <= 0 - See output'
        WRITE(6,*) 'Error in VOLUME_MODE'
        WRITE(6,*) 'wetdp/drydp/wvol/dvol/rhopar <= 0'
        WRITE(6,*) 'In Volume_mode: MDWAT (min,max,sum), IMODE=',IMODE
        WRITE(6,*) minval(mdwat(:,IMODE)),maxval(mdwat(:,IMODE)),       &
                   sum(mdwat(:,IMODE))
        WRITE(6,*) 'Location of min: ',minloc(mdwat(:,IMODE))
        WRITE(6,*) 'Location of max: ',maxloc(mdwat(:,IMODE))
        WRITE(6,*) 'In Volume_mode: wetdp (min,max,sum)'
        WRITE(6,*) minval(wetdp(:,IMODE)),maxval(wetdp(:,IMODE)),       &
                   sum(wetdp(:,IMODE))
        WRITE(6,*) 'Location of min: ',minloc(wetdp(:,IMODE))
        WRITE(6,*) 'Location of max: ',maxloc(wetdp(:,IMODE))
        WRITE(6,*) 'In Volume_mode: drydp (min,max,sum)'
        WRITE(6,*) minval(drydp(:,IMODE)),maxval(drydp(:,IMODE)),       &
                   sum(drydp(:,IMODE))
        WRITE(6,*) 'Location of min: ',minloc(drydp(:,IMODE))
        WRITE(6,*) 'Location of max: ',maxloc(drydp(:,IMODE))
        WRITE(6,*) 'In Volume_mode: rhopar(min,max,sum)'
        WRITE(6,*) minval(rhopar(:,IMODE)),maxval(rhopar(:,IMODE)),     &
                   sum(rhopar(:,IMODE))
        WRITE(6,*) 'Location of min: ',minloc(rhopar(:,IMODE))
        WRITE(6,*) 'Location of max: ',maxloc(rhopar(:,IMODE))
        WRITE(6,*) 'In Volume_mode: dvol(min,max,sum)'
        WRITE(6,*) minval(dvol (:,IMODE)),maxval(dvol (:,IMODE)),       &
                   sum(dvol (:,IMODE))
        WRITE(6,*) 'Location of min: ',minloc(dvol (:,IMODE))
        WRITE(6,*) 'Location of max: ',maxloc(dvol (:,IMODE))
        WRITE(6,*) 'In Volume_mode: wvol(min,max,sum)'
        WRITE(6,*) minval(wvol (:,IMODE)),maxval(wvol (:,IMODE)),       &
                   sum(wvol (:,IMODE))
        WRITE(6,*) 'Location of min: ',minloc(wvol (:,IMODE))
        WRITE(6,*) 'Location of max: ',maxloc(wvol (:,IMODE))

        JJL=minloc(wvol (:,IMODE))
        DO IIMODE=1,NMODES
         IF(MODE(IIMODE)) THEN
          WRITE(6,'(1a14,1i5,2e15.6)') 'IIMODE,ND,MDT=',                &
                             IIMODE,ND(JJL(1),IIMODE),MDT(JJL(1),IIMODE)
          DO IICP=1,NCP
           IF(COMPONENT(IIMODE,IICP)) THEN
            WRITE(6,'(1a17,2i5,1e15.6)') 'IIMODE,IICP,MD()=',           &
                           IIMODE,IICP,MD(JJL(1),IIMODE,IICP)
           END IF ! if component present in mode (IICP)
          END DO ! loop over cpts (IICP)
         END IF ! if mode present (IIMODE)
        END DO ! loop over modes (IIMODE)
        errcode=1
        CALL EREPORT('UKCA_VOLUME_MODE',errcode,cmessage)

       END IF ! if wetdp/drydp/wvol/dvol/rhopar <= 0
      END DO ! loop over modes

      IF (lhook) CALL dr_hook('UKCA_VOLUME_MODE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_VOLUME_MODE
