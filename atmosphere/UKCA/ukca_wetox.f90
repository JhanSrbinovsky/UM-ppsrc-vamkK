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
!    Calculates conversion of dissolved SO2 to sulfate by reaction with
!    H2O2 and O3. Uses Henry's law constants & oxidation rate constants
!    from Seinfeld & Pandis "Atmospheric Chemistry & Physics"
!    Note that this code is for use in the CTM and box model framework
!    or where the SO2 wet oxidation is to be done in the aerosol
!    code (WETOX_IN_AER=1).
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
      SUBROUTINE UKCA_WETOX(NBOX,ND,DRYDP,DELSO2,DELSO2_2,              &
       DELH2O2,LOWCLOUD,VFAC,SO2VMR,H2O2,ZH2O2,ZO3,ZHO2,PMID,           &
       T,AIRD,S,DTC,LDAY,LWC,ISO2WETOXBYO3)
!----------------------------------------------------------------------
!
!     Purpose
!     -------
!     Calculates conversion of dissolved SO2 to sulfate by reaction with
!     H2O2 and O3. Uses Henry's law constants & oxidation rate constants
!     from Seinfeld & Pandis "Atmospheric Chemistry & Physics"
!     (see pgs 337-367).
!     Assumes cloud LWC=0.2 gm^-3 and pH cloud = 4/5 SO2 dependant.
!
!     Also replenishes H2O2 concentration by HO2 self-reaction according
!     to expression in Jones et al (2001), JGR Atmos. vol. 106, no.D17,
!     pp 20,293--20,310.
!
!     Henry's law constant values at 298K from Seinfeld & Pandis pg 341.
!     Heats of dissolution values for Henry's law constants at other
!     temperatures from Seinfeld & Pandis page 342.
!
!     Note that this code is for use in the CTM and box model framework
!     or where the SO2 wet oxidation is to be done in the aerosol
!     code (WETOX_IN_AER=1).
!
!     For usual UKCA runs in the UM, this code is not called.
!
!     Note also, that this code currently has an assumption that the
!     clouds indicated by LOWCLOUD and VFAC are low clouds --- in fact
!     in the current version a liquid water content of 0.2 g/m3 is set
!     for all such clouds. It is intended to be used when climatological
!     monthly-mean cloud fields from ISCCP satellite database are used.
!
!     Inputs
!     ------
!     NBOX       : Number of grid boxes
!     ND         : Initial no. concentration of aerosol mode (ptcls/cc)
!     DRYDP      : Geometric mean dry diameter for each mode (m)
!     LOWCLOUD   : Horizontal low cloud fraction
!     VFAC       : Vertical low cloud fraction
!     SO2VMR     : Volume mixing ratio of SO2(g)
!     H2O2       : Number density of H2O2(g) [molecules per cc]
!     ZO3        : Background vmr of O3
!     ZH2O2      : Background concentration of H2O2 [molecules per cc]
!     ZHO2       : Background concentration of HO2  [molecules per cc]
!     PMID       : Air pressure (Pa)
!     T          : Air temperature (K)
!     AIRD       : Number density of air molecules
!     S          : Air specific humidity (kg/kg)
!     DTC        : Time step (s)
!     LDAY       : 1-Day, 0-Night
!     LWC        : Cloud liquid water content [kg/m3]
!     ISO2WETOXBYO3: Switch (1/0) for whether to include in-cloud oxdtn
!                    by O3 (pH value prescribed according to SO2 concn)
!
!     Outputs
!     -------
!     DELSO2     : Mass of S(IV) --> S(VI) by H2O2 (molecules/cc/DTC)
!     DELSO2_2   : Mass of S(IV) --> S(VI) by O3 (molecules/cc/DTC)
!     DELH2O2    : Total change in H2O2 mass due to reaction with S(IV)
!                : and replenishment via HO2 (molecules/cc/DTC)
!
!     Local Variables
!     ---------------
!     F          : cloud fraction
!     H2O2VMR    : volume mixing ratio of H2O2(g)
!     PHCLOUD    : pH of cloud (assumed to be 4/5 if SO2>or<0.5ppbv)
!     HSO2_298   : Henry's law constant for SO2 dissol. at 298K (M/atm)
!     DELHSO2_KCAL: Heat of dissolution for SO2 (kCal per mole)
!     DELHSO2_R  : Heat of dissolution for HSO2 / R (K)
!     HH2O2_298  : Henry's law constant for H2O2 dissol. at 298K (M/atm)
!     DELHH2O2_KCAL: Heat of dissolution for H2O2 (kCal per mole)
!     DELHH2O2_R : Heat of dissolution for HH2O2 / R (K)
!     K_S1       : constant in dissociation of SO2.H2O to H+ + HSO3-
!     K_S2       : constant in dissociation of HSO3 to H+ + SO3 2-
!     K1         : constant in equation for S(IV) oxidation
!     K2         : constant in equation for S(IV) oxidation
!     JPERKCAL   : number of Joules per kilocalorie (J kcal^-1)
!     HSO2       : Henry's law constant for SO2 dissol. at T (M/atm)
!     HH2O2      : Henry's law constant for H2O2 dissol. at T (M/atm)
!     HPLUSM     : Concentration of H+ ions in solution (M)
!     HSO3M      : Concentration of HSO3- ions in solution (M)
!     H2O2M      : Concentration of H2O2 in solution (M)
!     DSIVMDT    : Rate of oxidation of S(IV) to S(VI) (M s^-1)
!     PSO2ATM    : Partial pressure of SO2(g) in atmospheres
!     PH2O2ATM   : Partial pressure of H2O2(g) in atmospheres
!     TERM       : Temporary variable in Henry's law constant calcn
!     SO2        : Number concentration of SO2 (molecules per cc)
!     MAXDSL     : Maximum change in sulfate ( =min{SO2,H2O2} )
!     NWP        : No. concentration of water vapour (molecules per cc)
!     DH2O2TMP   : Rate of re-formation of H2O2 via HO2
!     OLDH2O2    : Temporary variable to store H2O2 concentration
!
!     Inputted by module UKCA_CONSTANTS
!     ---------------------------------
!     AVC        : Avogadro's constant
!     RR         : Universal gas constant (J/mol/K)
!     ZBOLTZ     : Boltzman Constant (kg m2 s-2 K-1 molec-1)
!     PPI        : 3.1415927
!
!     Inputted by module UKCA_MODE_SETUP
!     ----------------------------------
!     NMODES     : Number of possible aerosol modes
!     MODE       : Defines which modes are set
!     NUM_EPS    : Value of NEWN below which don't apply process
!
!     Inputted by module UKCA_SETUP_INDICES
!     -------------------------------------
!     MH2O2F     : Index for semi-prognostic H2O2
!--------------------------------------------------------------------
      USE UKCA_CONSTANTS,       ONLY: avc, rr, zboltz, ppi
      USE UKCA_MODE_SETUP,      ONLY: nmodes, mode, num_eps
      USE UKCA_SETUP_INDICES,   ONLY: mh2o2f
      USE ukca_option_mod,      ONLY: L_ukca_chem
      USE Control_Max_Sizes
      USE yomhook,              ONLY: lhook, dr_hook
      USE parkind1,             ONLY: jprb, jpim
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: LDAY(NBOX)
      INTEGER, INTENT(IN) :: ISO2WETOXBYO3
      REAL, INTENT(IN)    :: ND(NBOX,NMODES)
      REAL, INTENT(IN)    :: DRYDP(NBOX,NMODES)
      REAL, INTENT(IN)    :: LOWCLOUD(NBOX)
      REAL, INTENT(IN)    :: VFAC(NBOX)
      REAL, INTENT(IN)    :: SO2VMR(NBOX)
      REAL, INTENT(IN)    :: ZO3(NBOX)
      REAL, INTENT(IN)    :: ZH2O2(NBOX)
      REAL, INTENT(IN)    :: ZHO2(NBOX)
      REAL, INTENT(IN)    :: PMID(NBOX)
      REAL, INTENT(IN)    :: T(NBOX)
      REAL, INTENT(IN)    :: AIRD(NBOX)
      REAL, INTENT(IN)    :: S(NBOX)
      REAL, INTENT(IN)    :: DTC
      REAL, INTENT(IN)    :: LWC(NBOX)
      REAL, INTENT(INOUT) :: H2O2(NBOX)
      REAL, INTENT(OUT)   :: DELSO2(NBOX)
      REAL, INTENT(OUT)   :: DELSO2_2(NBOX)
      REAL, INTENT(OUT)   :: DELH2O2(NBOX)
!
!     Local variables
      INTEGER :: JL
      REAL    :: F
      REAL    :: PHCLOUD
      REAL    :: DELHSO2_R
      REAL    :: DELHO3_R
      REAL    :: DELHH2O2_R
      REAL    :: HSO2
      REAL    :: HO3
      REAL    :: HH2O2
      REAL    :: HPLUSM
      REAL    :: HSO3M
      REAL    :: O3M
      REAL    :: H2O2M
      REAL    :: DSIVMDT
      REAL    :: DSIVMDT2
      REAL    :: H2O2VMR
      REAL    :: O3VMR
      REAL    :: PSO2ATM
      REAL    :: PO3ATM
      REAL    :: PH2O2ATM
      REAL    :: TERM
      REAL    :: SO2
      REAL    :: MAXDSL
      REAL    :: NWP
      REAL    :: DH2O2TMP
      REAL    :: OLDH2O2
      REAL    :: SO2VMR2
      REAL    :: PSO2ATM2
      REAL    :: SO2H2O
      REAL    :: SO32M
      REAL    :: LWC_GM3(NBOX)
! .. below are new local variables added
      REAL    :: DC
      REAL    :: FAC_C
      REAL    :: DROPR_3
      REAL    :: C_3
      REAL    :: R_3
      REAL    :: DROPR_4
      REAL    :: C_4
      REAL    :: R_4
      REAL, PARAMETER :: HSO2_298     =  1.23
      REAL, PARAMETER :: DELHSO2_KCAL = -6.23
      REAL, PARAMETER :: HO3_298      =  1.13e-2
      REAL, PARAMETER :: DELHO3_KCAL  = -5.04
      REAL, PARAMETER :: HH2O2_298    =  7.45e4
      REAL, PARAMETER :: DELHH2O2_KCAL=-14.5
      REAL, PARAMETER :: K_S1         =  1.3e-2
      REAL, PARAMETER :: K_S2         =  6.6e-8
      REAL, PARAMETER :: K1_H2O2      =  7.5e7
      REAL, PARAMETER :: K2_H2O2      = 13.0
      REAL, PARAMETER :: K0_O3        =  2.4e4
      REAL, PARAMETER :: K1_O3        =  3.7e5
      REAL, PARAMETER :: K2_O3        =  1.5e9
      REAL, PARAMETER :: JPERKCAL     =  4.187e3

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_WETOX',zhook_in,zhook_handle)

!      LWC(:)=0.2e-3 ! set to constant at 0.2 g/m3 [of air] at present ** NOW INTENT(IN)
      LWC_GM3(:)=LWC(:)*1.0e3
! .. LWC is in kg/m3 [of air], LWC_GM3 is in g/m3 [of air]

      DO JL=1,NBOX
!
       OLDH2O2=H2O2(JL)
       DELSO2(JL)=0.0
       DELSO2_2(JL)=0.0
!      Calculate in-cloud (wet) oxidation of SO2
       F=LOWCLOUD(JL)*VFAC(JL)
! .. F is the cloud fraction in the gridbox
       IF(F > 0.0) THEN
        IF(MODE(3).OR.MODE(4)) THEN ! if accum-sol or coarse-sol mode
! .. below is for if some particles in either mode in this gridbox.
         IF((ND(JL,3) > NUM_EPS(3)).OR.(ND(JL,4) > NUM_EPS(4))) THEN
          IF (SO2VMR(JL) > 5e-10) PHCLOUD=4.0
          IF (SO2VMR(JL) < 5e-10) PHCLOUD=5.0
! .. assumes pH of cloud is either 4 or 5 depending on if SO2>/<0.5ppb
          H2O2VMR=H2O2(JL)/AIRD(JL)
          O3VMR=ZO3(JL) ! ZO3 is already vmr

! .. below calculates partial pressures of SO2,O3,H2O2 in atmospheres
          PSO2ATM=SO2VMR(JL)*PMID(JL)/1.0e5
          PO3ATM=O3VMR*PMID(JL)/1.0e5
          PH2O2ATM=H2O2VMR*PMID(JL)/1.0e5

! .. below calculates heats of dissolution for Henry's law coeffs
! .. (values taken from page 342 Seinfeld & Pandis).
          DELHSO2_R=DELHSO2_KCAL*JPERKCAL/RR
          DELHO3_R=DELHO3_KCAL*JPERKCAL/RR
          DELHH2O2_R=DELHH2O2_KCAL*JPERKCAL/RR

! .. below stores temperature correction for Henry's law coeffs
          TERM=(1.0/T(JL))-(1.0/298.0)

! .. below calculates Henry's law coefficients for SO2,O3,H2O2
          HSO2=HSO2_298*EXP(-DELHSO2_R*TERM)
          HO3=HO3_298*EXP(-DELHO3_R*TERM)
          HH2O2=HH2O2_298*EXP(-DELHH2O2_R*TERM)

! .. below stores dissolved ion concentrations for H+,HSO3-
          HPLUSM=10.0**(-PHCLOUD)
          HSO3M=K_S1*HSO2*PSO2ATM/HPLUSM       ! [6.37 S&P pg 349]

! .. below calculates dissolved H2O2 and O3 concentration
          H2O2M=HH2O2*PH2O2ATM
          O3M=HO3*PO3ATM

!---------------------------------------------------------------------
! This section below calculates in-cloud oxidation of SO2 by H2O2
!
! .. below calculates oxidation rate of S(IV) by H2O2 (M/s)
! .. (see Seinfeld & Pandis, page 366, equation 6.84)
          DSIVMDT=K1_H2O2*HPLUSM*H2O2M*HSO3M/(1.0+K2_H2O2*HPLUSM)

! .. below then converts this rate to a change in concentration
! .. in SO2 in the gas phase over the chemistry time step DTC.
!! comment out this line to calculate following GLOMAP-bin method below
!!          DELSO2(JL)=DSIVMDT*F*(LWC_GM3(JL)*1.0e-6)*AVC*DTC*1.0e-3
! .. 0.001 per litre = 1 per cc
! .. 1.0e-6 converts from g/m3 [of air] to g/cm3 [of air]
! .. 1.0e-3 is 1/rho_water where rho_water = 1000 g/dm3 (=1000kg/m3)

! set SO2 molecular concentration (per cc)
          SO2=AIRD(JL)*SO2VMR(JL)

          DC=1.0E-5
          FAC_C=4.0*PPI*DC*1.0e6
! 1.0E6 converts ND cm^-3 ---> m^-3 ---> ensures R_3 is dimensionless
!        n.b. DROPR,DC,ND,DTC are in m,m^2s^-1,cm^-3,s

! below is for acc-sol mode
          DROPR_3=40.0*(0.5*DRYDP(JL,3)) ! assume cloud radius = 40 * g.m.radius
          C_3=FAC_C*DROPR_3
          R_3=-ND(JL,3)*C_3*DTC

! below is for cor-sol mode
          DROPR_4=40.0*(0.5*DRYDP(JL,4)) ! assume cloud radius = 40 * g.m.radius
          C_4=FAC_C*DROPR_4
          R_4=-ND(JL,4)*C_4*DTC

! below calculates as in GLOMAP-bin -- here include F as cloud fraction
          DELSO2(JL)=F*(1.0-EXP(R_3+R_4))*SO2

          IF(SO2 <= H2O2(JL)) MAXDSL=SO2
          IF(SO2 > H2O2(JL)) MAXDSL=H2O2(JL)
          IF(DELSO2(JL) > MAXDSL) DELSO2(JL)=MAXDSL
! .. above limits DELSO2 to be max(SO2,H2O2)

! .. below reduces SO2, H2O2 gas phase according to above rate
          H2O2(JL)=H2O2(JL)-DELSO2(JL)
          SO2VMR2=SO2VMR(JL)-DELSO2(JL)/AIRD(JL)

          IF(ISO2WETOXBYO3.EQ.1) THEN ! switch for SO2 wetox by O3
!---------------------------------------------------------------------
! This section below then calculates in-cloud oxidation of SO2 by O3

! .. below calculates new partial pressure of SO2 in atmospheres
           PSO2ATM2=SO2VMR2*PMID(JL)/1.0e5

! .. below calculates new dissolved concentration of SO2
           SO2H2O=HSO2*PSO2ATM2

! .. below then calculates new dissolved ions conc of HSO3^- & SO3^{2-}
           HSO3M=K_S1*SO2H2O/HPLUSM
           SO32M=K_S2*HSO3M/HPLUSM

! .. below calculates oxidation rate of S(IV) by O3 (M/s)
! .. (see Seinfeld & Pandis, page 363, equation 6.80)
           DSIVMDT2=(K0_O3*SO2H2O+K1_O3*HSO3M+K2_O3*SO32M)*O3M

! .. below then converts this rate to a change in concentration
! .. in SO2 in the gas phase over the chemistry time step DTC.
           DELSO2_2(JL)=DSIVMDT2*F*(LWC_GM3(JL)*1.0e-6)*AVC*DTC*1.0e-3
! .. 0.001 per liter = 1 per cc
! .. 1.0e-6 converts from g/m3 [of air] to g/cm3 [of air]
! .. 1.0e-3 is 1/rho_water where rho_water = 1000 g/dm3 (=1000kg/m3)

           IF(SO2VMR2 <= ZO3(JL)) MAXDSL=SO2VMR2*AIRD(JL)
           IF(SO2VMR2 > ZO3(JL)) MAXDSL=ZO3(JL)*AIRD(JL)
           IF(DELSO2_2(JL) > MAXDSL) DELSO2_2(JL)=MAXDSL
! .. above limits DELSO2_2 to available SO2,O3

! .. Note that SO2 does not need to be updated again here as
! .. we are not updating the main SO2 concentrations in this subroutine.
! .. This is done using the returned arrays DELSO2 and DELSO2_2 in the
! .. UKCA_AERO_STEP routine.

          END IF ! if switch for SO2 in-cloud oxidation by O3 = 1

         END IF ! if some particles in either soluble accum or coarse
        END IF ! if either soluble accum or coarse are defined
       END IF ! if grid box is in cloud

!---------------------------------------------------------------------
! This section calculates replenishment of H2O2 by HO2 self-reaction.
! (only does this for off-line version with semi-prognostic H2O2)

       IF(MH2O2F > 0 .AND. (.NOT. L_ukca_chem)) THEN ! if using semi-prognostic H2O2
        IF (LDAY(JL) == 1) THEN ! if daytime then H2O2 is reformed

!        Below calculates water vapour concentration (cm-3) NWP
!            (1.609*S(JL)) is vmr of H2O (g)
         NWP = AIRD(JL)*(1.609*S(JL))

!        Calculate rate of H2O2 replenishment according to background
!        HO2 concentration ZHO2 (follows expression in Jones et al 2001)
         DH2O2TMP = (2.3E-13*EXP(600.0/T(JL)) +                         &
             1.9E-33*AIRD(JL)*EXP(890.0/T(JL))) *                       &
             (1.0 + 1.4E-21*NWP*EXP(2200.0/T(JL)))*                     &
             ZHO2(JL)*ZHO2(JL)

         H2O2(JL)=H2O2(JL)+DH2O2TMP*DTC
! .. above replenishes H2O2 due to production by HO2 self-reaction
         IF (H2O2(JL) > ZH2O2(JL)) H2O2(JL)=ZH2O2(JL)
! .. above limits H2O2 to be at most ZH2O2 (background concentration)
        END IF ! if day-time
       END IF ! if MH2O2F>0 (using semi-prognostic H2O2)

!! n.b. the line below has been moved here cf v1_gm4_coupled (bug-fix).
       DELH2O2(JL)=H2O2(JL)-OLDH2O2 ! calculate net change in H2O2

      END DO

      IF (lhook) CALL dr_hook('UKCA_WETOX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_WETOX
