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
!  Does one chemistry time step of the UKCA-MODE aerosol model.
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
      SUBROUTINE UKCA_AERO_STEP(NBOX,                                   &
       ND,MDT,MD,MDWAT,S0G,DRYDP,WETDP,RHOPAR,DVOL,WVOL,SM,             &
       AIRD,AIRDM3,RHOA,MFPA,DVISC,T,TSQRT,RH,S,PMID,PUPPER,PLOWER,     &
       EMC,EMCBM,ZO3,ZHO2,ZH2O2,USTR,US10M,ZNOT,DELTA_Z,                &
       SURTP,LAND_FRAC,SURF,SEAICE,                                     &
       CRAIN,DRAIN,CRAIN_UP,DRAIN_UP,FCONV_CONV,LOWCLOUD,VFAC,CLF,      &
       EMANSO2,EMVOLCONSO2,EMVOLEXPSO2,EMBIOMSO2,ISO2EMS,               &
       SILT,CLAY,SAND,SWC,SNOWICE,                                      &
       PSOURCE,LAI,TEX,MODE2_RADIUS,MODE3_RADIUS,                       &
       MODE2_NUMBER,MODE3_NUMBER,NDUSTEMINBOX,                          &
       DTC,DTM,DTZ,NMTS,NZTS,LDAY,RMOIS,RJOUR,ACT,BUD_AER_MAS,          &
       RAINOUT_ON,                                                      &
       IMSCAV_ON,WETOX_ON,DDEPAER_ON,SEDI_ON,ISO2WETOXBYO3,             &
       DRYOX_IN_AER,WETOX_IN_AER,DELSO2,DELSO2_2,                       &
       COND_ON,NUCL_ON,COAG_ON,BLN_ON,ICOAG,IMERGE,IFUCHS,IWVOLMETHOD,  &
       IDCMFP,ICONDIAM,IBLN,I_NUC_METHOD,                               &
       IACTMETHOD,IDDEPAER,INUCSCAV,VERBOSE,CHECKMD_ND,INTRAOFF,        &
       INTEROFF,IDUSTEMS,S0G_DOT_CONDENSABLE,LWC,CLWC,PVOL,PVOL_WAT,    &
       JLABOVE,ILSCAT,N_MERGE_1D,DLONARR,DLATARR,HEIGHT,HTPBLG,myproc)

!-----------------------------------------------------------------------
!  Inputs
!  ------
!  NBOX        : Number of grid boxes
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  MDWAT       : Molecular concentration of water (molecules per ptcl)
!  S0G         : Partial masses of gas phase species (kg per gridbox)
!  DRYDP       : Geometric mean dry diameter for each mode (m)
!  WETDP       : Geometric mean wet diameter for each mode (m)
!  RHOPAR      : Particle density [incl. H2O & insoluble cpts] (kgm^-3)
!  DVOL        : Geometric mean dry volume for each mode (m^3)
!  WVOL        : Geometric mean wet volume for each mode (m^3)
!  SM          : Grid box mass of air (kg)
!  AIRD        : Number density of air (per cm3)
!  AIRDM3      : Number density of air (per m3)
!  RHOA        : Air density (kg/m3)
!  MFPA        : Mean free path of air (m)
!  DVISC       : Dynamic viscosity of air (kg m^-1 s^-1)
!  T           : Centre level temperature (K)
!  TSQRT       : Square-root of centre level temperature (K)
!  RH          : Relative humidity (dimensionless 0-1)
!  S           : Specific humidity (kg/kg)
!  PMID        : Centre level pressure (Pa)
!  PUPPER      : Upper interface pressure (Pa)
!  PLOWER      : Lower interface pressure (Pa)
!  EMC         : BC/OC ems rates from bio- & fossil-fuels (kgC/box/s)
!  EMCBM       : BC/OC ems rates from biomass burning (kgC/box/s)
!  ZO3         : Backgrnd vmr of O3
!  ZHO2        : Backgrnd conc. of HO2 (molecules per cc)
!  ZH2O2       : Backgrnd conc. of H2O2 (molecules per cc)
!  USTR        : Surface friction velocity (m/s)
!  US10M       : Scalar wind at 10m (ms-1)
!  ZNOT        : Roughness length (m)
!  SURTP       : Surface type: 0=seasurf,1=landsurf,2=oversea,3=overland
!  LAND_FRAC   : Fraction of horizontal gridbox area covered by land
!  SURF        : Surface area of box (horizontal) (m^2)
!  SEAICE      : Fraction of horizontal gridbox area containing seaice
!  CRAIN       : Rain rate for conv precip. in box (kgm^-2s^-1)
!  DRAIN       : Rain rate for dyn. precip. in box (kgm^-2s^-1)
!  CRAIN_UP    : Rain rate for conv precip. in box above (kgm^-2s^-1)
!  DRAIN_UP    : Rain rate for dyn. precip. in box above (kgm^-2s^-1)
!  FCONV_CONV  : Fraction of box condensate --> rain in 6 hours (conv)
!  LOWCLOUD    : Horizontal low cloud fraction
!  VFAC        : Vertical low cloud fraction
!  EMANSO2     : Anthrop. SO2 ems rates, low sources (kgSO2/box/s)
!  EMVOLCONSO2(NBOX) Volcanic SO2 ems rates (cont. src) (kgSO2/box/s)
!  EMVOLEXPSO2(NBOX) Volcanic SO2 ems rates (expl. src) (kgSO2/box/s)
!  EMBIOMSO2(NBOX) Biomass SO2 ems rates (kgSO2/box/s)
!  ISO2EMS     : Switch for scheme for primary H2SO4 aerosol emissions
!  SAND        : Fraction of grid box soil which is sand
!  SILT        : Fraction of grid box soil which is silt
!  CLAY        : Fraction of grid box soil which is clay
!  SWC         : Soil water content (m3water/m3soil)
!  SNOWICE     : Fraction of gridbox that is *not* covered by snow/ice
!  PSOURCE     : Preferential source area grid from Tegen et al (2002)
!  LAI         : Leaf area index grid
!  MODE2_RADIUS: Geometric mean radius for dust mode 2 (AEROCOM daily)
!  MODE3_RADIUS: Geometric mean radius for dust mode 3 (AEROCOM daily)
!  MODE2_NUMBER: Particle number flux for dust mode 2 (AEROCOM daily)
!  MODE3_NUMBER: Particle number flux for dust mode 3 (AEROCOM daily)
!  NDUSTEMINBOX: Number of non-zero 1x1 emissions fluxes in gridbox
!  DTC         : Chemistry time step (s)
!  DTM         : Microphysics time step (s)
!  DTZ         : Competition (cond/nucl) time step (s)
!  NMTS        : Number of microphysics timesteps per DTC
!  NZTS        : Number of competition timesteps per DTM
!  LDAY        : 1-Day, 0-Night
!  RMOIS       : Month of year (1.0-12.0)
!  RJOUR       : Day of month (1.0-31.0)
!  MODESOL     : Specifies whether mode is soluble/insoluble (=1/0)
!  ACT         : Particle dry radius above which activation is assumed
!  RAINOUT_ON  : Switch : rainout (nucl. scav.) is on/off (1/0)
!  IMSCAV_ON   : Switch : impaction scavenging is on/off (1/0)
!  WETOX_ON    : Switch : wet oxidation [+ cloudproc] is on/off (1/0)
!  DDEPAER_ON  : Switch : aerosol dry deposition is on/off (1/0)
!  SEDI_ON     : Switch : aerosol sedimentation is on/off (1/0)
!  ISO2WETOXBYO3:Switch : so2in-cloud ox by o3 (assumed pH) on/off (1/0)
!  DRYOX_IN_AER: Switch : update gas phase condensible concentrations
!                         on competition tstep in aerosol code? (1/0)
!  WETOX_IN_AER: Switch : calc. wet ox. of SO2 in aerosol code? (1/0)
!  COND_ON     : Switch : vapour condensation is  on/off (1/0)
!  NUCL_ON     : Switch : binary nucleation is on/off (1/0)
!  COAG_ON     : Switch : coagulation is on/off (1/0)
!  BLN_ON      : Switch : Boundary layer nucleation is on/off (1/0)
!  ICOAG       : KIJ method (1:GLOMAP, 2: M7, 3: UMorig, 4:UMorig MFPP)
!  IMERGE      : Switch to use mid-pts (=1), edges (2) or dynamic (=3)
!  IFUCHS      : Switch : Fuchs (1964) or Fuchs-Sutugin (1971) for CC
!  IWVOLMETHOD : Switch : wet volume method (1=as-H2SO4,2=multi-cpt)
!  IACTMETHOD  : Switch : activation method (0=off,1=fixed ract,2=NSO3)
!  IDDEPAER    : Switch : dry dep method (1=as in Spr05, 2=incl. sedi)
!  IDCMFP      : Switch : diffusion/mfp  (1=as Gbin v1, 2=as Gbin v1_1)
!  ICONDIAM    : Switch : wet diam in UKCA_CONDEN (1=g.mean,2=condiam.)
!  IBLN        : Switch : BLN method (1=activation,2=kinetic,3=PNAS)
!  I_NUC_METHOD: Switch for nucleation (how to combine BHN and BLN)
! (1=initial Pandis94 approach (no BLN even if switched on) -- Do not use!!
! (2=binary homogeneous nucleation applying BLN to BL only if switched on)
!   note there is an additional switch i_bhn_method (local to CALCNUCRATE)
!   to switch between using Kulmala98 or Vehkamakki02 for BHN rate
! (3=use same rate at all levels either activation(IBLN=1), kinetic(IBLN=2),
!  PNAS(IBLN=3),eucaari-kinetic(IBLN=4),eucaari-org1(IBLN=5),eucaari-org2(IBLN=6)
!  note that if I_NUC_METHOD=3 and IBLN=3 then also add on BHN rate as in PNAS.
!  I_NUC_METHOD: Switch : 1=use Pandis94 (no BLN), 2=use BHN with BLN if on 
!                       : 3=apply BLN parameterization throughout column 
!  INUCSCAV    : Switch : scheme for removal by nucl scav
!                (1=as GLOMAP [accsol,corsol], 2=scav params as Stier05]
!  VERBOSE     : If =1 prints min/max (ND,MDT etc) after each process
!                for 1st grid box (for box model tests)
!  INTRAOFF    : Switch to turn off intra-modal coagulation
!  INTEROFF    : Switch to turn off intra-modal coagulation
!  CHECKMD_ND  : Switch : check for values of MD, ND out of range? (1/0)
!  IDUSTEMS    : Switch : dust emissions scheme
!              :          (1=Pringle scheme, 2=AEROCOM00 daily)
!  S0G_DOT_CONDENSABLE : Gas phase chemistry tendencies for change in
!                        condensable gas phase species (vmr per s)
!                        due to chem, dry and wet deposition.
!  LWC         : Cloud liquid water content [kg/m3]
!  CLWC        : Cloud liquid water content [kg/kg]
!  JLABOVE     : Index of box directly above this grid box
!  ILSCAT      : Land surface category (based on 9 UM landsurf types)
!  N_MERGE_1D  : Count number of mode-merges in each box, each mode
!  DLONARR     : Longitude on NBOX array
!  DLATARR     : Latitude  on NBOX array
!  DELSO2      : S(IV) --> S(VI) by H2O2 (molecules per cc)
!                [input from gas phase chem. scheme if WETOX_IN_AER=0]
!  DELSO2_2    : S(IV) --> S(VI) by O3u  (molecules per cc)
!                [input from gas phase chem. scheme if WETOX_IN_AER=0]
!  HEIGHT      : Mid-level height of gridbox
!  HTPBLG      : Height of boundary-layer in gridbox vertical-column
!  myproc      : Processor number
!
!  Outputs
!  -------
!  ND          : Aerosol ptcl number density for mode (cm^-3)
!  MDT         : Avg tot mass of aerosol ptcl in mode (particle^-1)
!  MD          : Avg cpt mass of aerosol ptcl in mode (particle^-1)
!  MDWAT       : molecular concentration of water (molecules per ptcl)
!  S0G         : Partial masses of gas phase species (kg per gridbox)
!  BUD_AER_MAS: Aerosol budget terms for mass
!  DELSO2      : S(IV) oxidised to S(VI) by H2O2 (molecules per cc)
!                [input as zero if WETOX_IN_AER=1]
!  DELSO2_2    : S(IV) oxidised to S(VI) by O3 (molecules per cc)
!                [input as zero if WETOX_IN_AER=1]
!  PVOL        : Partial volumes of each component in each mode (m3)
!  PVOL_WAT    : Partial volume of water in each mode (m3)
!
!  Local Variables
!  ---------------
!  SO2         : Sulfur dioxide conc (molecules per cc)
!  SO2_VMR     : Sulfur dioxode vmr
!  H2SO4       : Sulfuric acid vapour conc (moleculesH2SO4/cc)
!  SEC_ORG     : 1st stage oxidation product from terpene (molecules/cc)
!  H2O2        : Hygrogen peroxide conc (per cc)
!  GC          : Condensable vapour conc (molecules cpt per cc)
!  GCOLD       : Condensable vapour conc b4 process (molecules cpt/cm3)
!  DELH2O2     : Change in hygrogen peroxide conc due to process (/cm3)
!  FRAC_AQ_ACC : Fraction of wet oxidised SO2 -> soluble accum. mode
!  FRAC_AQ_COR : Fraction of wet oxidised SO2 -> soluble coarse mode
!  TOTWETOX    : Total wet oxidation rate of SO2 (molecules/cm3)
!  DELGC_COND: Change in vapour conc due to cond (molecules cpt/cm3)
!  DELGC_NUCL: Change in vapour conc due to nucl (molecules cpt/cm3)
!  DELH2SO4_NUCL: Change in H2SO4 conc due to nucl (molecules/cm3)
!  DELTAS0G    : Overall change in partial masses of condensable gases
!                after all NZTS competition steps (kg/box).
!  S0G_TO_GC   : Molar mass ratio : condensing gas to aerosol phase cpt
!                e.g. MM_GAS(MSEC_ORG)=0.150 kg/mol but it may condense
!                into organic carbon component with MM(CP_OC)=0.0168.
!                Need to apply conversion when going from gas-->aerosol
!  JRATE       : H2SO4 depletion rate by BH nucleation (molecules/cm3/s)
!  AGETERM1    : Depletion rate of each component (molecules cpt/cc/DTZ)
!                from condensation onto the 3 insoluble modes.
!                Used for calculation of ageing rate in UKCA_AGEING.
!  AGETERM2    : Rate of accomodation of material to each insoluble mode
!                as a result of coagulation with smaller soluble modes
!                (in molecules cpt /cm3/DTZ)
!                Used for calculation of ageing rate in UKCA_AGEING.
!  PROCESS     : Character string to store process currently being done
!                for error checking in box model
!  KII_ARR     : Coag coeff for intra-modal coag (IMODE-IMODE) (cm^3/s)
!  KIJ_ARR     : Coag coeff for inter-modal coag (IMODE-JMODE) (cm^3/s)
!
!  Inputted by module UKCA_CONSTANTS
!  ---------------------------------
!  PPI         : 3.1415927...........
!  AVC         : Avogadros constant (mol-1)
!  ZBOLTZ      : Boltzmann's constant (kg m2 s-2 K-1 molec-1)
!  VKARMN      : Von Karman's constant = 0.4
!  RA          : Dry air gas constant = 287.05 Jkg^-1 K^-1
!  RR          : Universal gas constant = 8.314 J/mol/K
!  GG          : Gravitational acceleration = 9.80665 ms^-2
!  NMOL        : Number of molecules per particle at nucleation
!  CONC_EPS    : Threshold for molecular conc. (molecules per cc)
!  EMS_EPS     : Threshold for emissions fluxes (kg/gridbox/s)
!  DN_EPS      : Value of DELN below which do not carry out process
!
!  Inputted by module UKCA_MODE_SETUP
!  ----------------------------------
!  NMODES      : Number of possible aerosol modes
!  NCP         : Number of possible aerosol components
!  MODE        : Logical variable defining which modes are set.
!  COMPONENT   : Logical variable defining which cpt are in which dsts
!  SOLUBLE     : Logical variable defining which cpts are soluble
!  MM          : Molar masses of components (kg/mole)
!  RHOCOMP     : Densities (dry) of each component (kg/m^3)
!  NO_IONS     : Number of dissociating ions in soluble components
!  DDPLIM0     : Lower limit for dry diameter in mode (m)
!  DDPLIM1     : Upper limit for dry diameter in mode (m)
!  DDPMID      : Mid-point of size mode = exp(0.5*(lndp0+lndp1)) (m)
!  MFRAC_0     : Initial mass fraction to set when no particles.
!  SIGMA       : Geometric standard deviation for each mode
!  X           : EXP((9/2)*LOG^2(SIGMA_G))
!  NUM_EPS     : Value of NEWN below which do not carry out process
!  MMID        : Mass of particle with dry diameter
!                dp=dpmed_g=exp(0.5*(lndp0+lndp1)) (ptcl^-1)
!  COAG_MODE   : Switch to defines which mode an IMODE-JMODE
!                coagulation goes into
!  COLLEFF4    : Array of aerosol-raindrop collision efficiencies
!                for impaction scavenging (l-u table) [NCOLL,NROW]
!  RADDROP     : Raindrop radii (microns) at mid-points in raindrop grid
!                of dimension NROW -- for impaction scavenging
!  NCOLL       : # of columns in COLLEFF4 (=20) corresp to aerosol grid
!  NROW        : # of rows in COLLEFF4 (=19) corresp to raindrop grid
!  CP_SU       : Component where sulfate is stored
!  CP_BC       : Component where black carbon is stored
!  CP_OC       : Component where organic carbon is stored
!  CP_CL       : Component where NaCl is stored
!  CP_DU       : Component where dust is stored
!  CP_SO       : Component where condensible organic species is stored
!
!  Inputted by module UKCA_SETUP_INDICES
!  -------------------------------------
!  ICHEM       : Switch for gas phase chemistry scheme (0=off,1=on)
!  NADVG       : # of advected gas phase tracers
!  NCHEMG      : # of gas phase tracers for gas phase chemistry scheme
!  MM_GAS      : Molar masses for gas phase species (kg/mol)
!  MSOTWO      : Index of MM_GAS, WTRATC and S0G for SO2
!  MH2SO4      : Index of MM_GAS, WTRATC and S0G for H2SO4
!  MH2O2       : Index of MM_GAS, WTRATC and S0G for H2O2
!  MH2O2F      : Index of MM_GAS, WTRATC and S0G for H2O2F
!  CONDENSABLE : Logical variable defining which cpts are condensable
!  DIMEN       : Molecular diamters of condensable components (m)
!  Various indices for budget terms in BUD_AER_MAS
!
!-----------------------------------------------------------------------
      USE UKCA_CONSTANTS
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE
!
! .. Subroutine Interface
      INTEGER, INTENT(IN) :: NBOX
      INTEGER, INTENT(IN) :: IDUSTEMS
      INTEGER, INTENT(IN) :: ISO2EMS
      INTEGER, INTENT(IN) :: NMTS
      INTEGER, INTENT(IN) :: NZTS
      INTEGER, INTENT(IN) :: LDAY(NBOX)
      INTEGER, INTENT(IN) :: RAINOUT_ON
      INTEGER, INTENT(IN) :: IMSCAV_ON
      INTEGER, INTENT(IN) :: WETOX_ON
      INTEGER, INTENT(IN) :: DDEPAER_ON
      INTEGER, INTENT(IN) :: SEDI_ON
      INTEGER, INTENT(IN) :: ISO2WETOXBYO3
      INTEGER, INTENT(IN) :: DRYOX_IN_AER
      INTEGER, INTENT(IN) :: WETOX_IN_AER
      INTEGER, INTENT(IN) :: COND_ON
      INTEGER, INTENT(IN) :: NUCL_ON
      INTEGER, INTENT(IN) :: COAG_ON
      INTEGER, INTENT(IN) :: BLN_ON
      INTEGER, INTENT(IN) :: ICOAG
      INTEGER, INTENT(IN) :: IMERGE
      INTEGER, INTENT(IN) :: IFUCHS
      INTEGER, INTENT(IN) :: IWVOLMETHOD
      INTEGER, INTENT(IN) :: IACTMETHOD
      INTEGER, INTENT(IN) :: IDDEPAER
      INTEGER, INTENT(IN) :: IDCMFP
      INTEGER, INTENT(IN) :: ICONDIAM
      INTEGER, INTENT(IN) :: IBLN
      INTEGER, INTENT(IN) :: I_NUC_METHOD
      INTEGER, INTENT(IN) :: INUCSCAV
      INTEGER, INTENT(IN) :: VERBOSE
      INTEGER, INTENT(IN) :: CHECKMD_ND
      INTEGER, INTENT(IN) :: INTRAOFF
      INTEGER, INTENT(IN) :: INTEROFF
      INTEGER, INTENT(IN) :: JLABOVE(NBOX)
      INTEGER, INTENT(IN) :: ILSCAT(NBOX)
      INTEGER, INTENT(IN) :: N_MERGE_1D(NBOX,NMODES)
      INTEGER, INTENT(IN) :: NDUSTEMINBOX(NBOX)
      INTEGER, INTENT(IN) :: myproc

      REAL, INTENT(IN) :: RMOIS
      REAL, INTENT(IN) :: RJOUR
      REAL, INTENT(IN) :: DRYDP(NBOX,NMODES)
      REAL, INTENT(IN) :: WETDP(NBOX,NMODES)
      REAL, INTENT(IN) :: RHOPAR(NBOX,NMODES)
      REAL, INTENT(IN) :: DVOL(NBOX,NMODES)
      REAL, INTENT(IN) :: WVOL(NBOX,NMODES)
      REAL, INTENT(IN) :: SM(NBOX)
      REAL, INTENT(IN) :: AIRD(NBOX)
      REAL, INTENT(IN) :: AIRDM3(NBOX)
      REAL, INTENT(IN) :: RHOA(NBOX)
      REAL, INTENT(IN) :: MFPA(NBOX)
      REAL, INTENT(IN) :: DVISC(NBOX)
      REAL, INTENT(IN) :: T(NBOX)
      REAL, INTENT(IN) :: TSQRT(NBOX)
      REAL, INTENT(IN) :: RH(NBOX)
      REAL, INTENT(IN) :: S(NBOX)
      REAL, INTENT(IN) :: PMID(NBOX)
      REAL, INTENT(IN) :: PUPPER(NBOX)
      REAL, INTENT(IN) :: PLOWER(NBOX)
      REAL, INTENT(IN) :: ZO3(NBOX)
      REAL, INTENT(IN) :: ZHO2(NBOX)
      REAL, INTENT(IN) :: ZH2O2(NBOX)
      REAL, INTENT(IN) :: USTR(NBOX)
      REAL, INTENT(IN) :: US10M(NBOX)
      REAL, INTENT(IN) :: ZNOT(NBOX)
      REAL, INTENT(IN) :: DELTA_Z(NBOX)
      REAL, INTENT(IN) :: SURTP(NBOX)
      REAL, INTENT(IN) :: LAND_FRAC(NBOX)
      REAL, INTENT(IN) :: SURF(NBOX)
      REAL, INTENT(IN) :: SEAICE(NBOX)
      REAL, INTENT(IN) :: CRAIN(NBOX)
      REAL, INTENT(IN) :: DRAIN(NBOX)
      REAL, INTENT(IN) :: CRAIN_UP(NBOX)
      REAL, INTENT(IN) :: DRAIN_UP(NBOX)
      REAL, INTENT(IN) :: FCONV_CONV(NBOX)
      REAL, INTENT(IN) :: LOWCLOUD(NBOX)
      REAL, INTENT(IN) :: VFAC(NBOX)
      REAL, INTENT(IN) :: CLF(NBOX)
      REAL, INTENT(IN) :: EMANSO2(NBOX,6)
      REAL, INTENT(IN) :: EMVOLCONSO2(NBOX)
      REAL, INTENT(IN) :: EMVOLEXPSO2(NBOX)
      REAL, INTENT(IN) :: EMBIOMSO2(NBOX)
      REAL, INTENT(IN) :: EMC(NBOX,4)
      REAL, INTENT(IN) :: EMCBM(NBOX,2)
      REAL, INTENT(IN) :: SILT(NBOX)
      REAL, INTENT(IN) :: SAND(NBOX)
      REAL, INTENT(IN) :: CLAY(NBOX)
      REAL, INTENT(IN) :: SWC(NBOX)
      REAL, INTENT(IN) :: SNOWICE(NBOX)
      REAL, INTENT(IN) :: PSOURCE(NBOX)
      REAL, INTENT(IN) :: LAI(NBOX)
      REAL, INTENT(IN) :: TEX(NBOX)
      REAL, INTENT(IN) :: MODE2_RADIUS(NBOX,10)
      REAL, INTENT(IN) :: MODE3_RADIUS(NBOX,10)
      REAL, INTENT(IN) :: MODE2_NUMBER(NBOX,10)
      REAL, INTENT(IN) :: MODE3_NUMBER(NBOX,10)
      REAL, INTENT(IN) :: DTC
      REAL, INTENT(IN) :: DTM
      REAL, INTENT(IN) :: DTZ
      REAL, INTENT(IN) :: ACT
      REAL, INTENT(IN) :: S0G_DOT_CONDENSABLE(NBOX,NCHEMG)
      REAL, INTENT(IN) :: LWC(NBOX)
      REAL, INTENT(IN) :: CLWC(NBOX)
      REAL, INTENT(IN) :: DLONARR(NBOX)
      REAL, INTENT(IN) :: DLATARR(NBOX)
      REAL, INTENT(IN) :: HTPBLG(NBOX)
      REAL, INTENT(IN) :: HEIGHT(NBOX)
!
! .. Outputs
      REAL, INTENT(INOUT) :: ND(NBOX,NMODES)
      REAL, INTENT(INOUT) :: MDT(NBOX,NMODES)
      REAL, INTENT(INOUT) :: MD(NBOX,NMODES,NCP)
      REAL, INTENT(INOUT) :: MDWAT(NBOX,NMODES)
      REAL, INTENT(INOUT) :: S0G(NBOX,NADVG)
      REAL, INTENT(INOUT) :: DELSO2(NBOX)
      REAL, INTENT(INOUT) :: DELSO2_2(NBOX)
      REAL, INTENT(INOUT) :: PVOL(NBOX,NMODES,NCP)
      REAL, INTENT(INOUT) :: PVOL_WAT(NBOX,NMODES)
      REAL, INTENT(INOUT) :: BUD_AER_MAS(NBOX,0:NBUDAER)
!
!     Local variables

      INTEGER :: errcode                ! Variable passed to ereport

      INTEGER :: JV
      INTEGER :: ICP
      INTEGER :: CP_SO4
      INTEGER :: IMTS
      INTEGER :: IZTS
      REAL :: SO2(NBOX)
      REAL :: H2O2(NBOX)
      REAL :: H2SO4(NBOX)
      REAL :: SEC_ORG(NBOX)
      REAL :: DELH2O2(NBOX)
      REAL :: SO2_VMR(NBOX)
      REAL :: GC(NBOX,NCHEMG)
      REAL :: GCOLD(NBOX,NCHEMG)
      REAL :: DELGC_COND(NBOX,NCHEMG)
      REAL :: DELGC_NUCL(NBOX,NCHEMG)
      REAL :: DELTAS0G(NBOX)
      REAL :: DELH2SO4_NUCL(NBOX)
      REAL :: JRATE(NBOX)
      REAL :: AGETERM1(NBOX,3,NCHEMG)
      REAL :: AGETERM2(NBOX,4,3,NCP)
      REAL :: S0G_TO_GC
      REAL :: FRAC_AQ_ACC(NBOX)
      REAL :: FRAC_AQ_COR(NBOX)
      REAL :: TOTWETOX(NBOX)
      REAL :: KII_ARR(NBOX,NMODES)
      REAL :: KIJ_ARR(NBOX,NMODES,NMODES)
      CHARACTER(LEN=30) :: PROCESS
      CHARACTER(LEN=2)  :: STRIZTS
      LOGICAL :: MASK(NBOX)
      LOGICAL :: MASK1(NBOX)
      LOGICAL :: MASK2(NBOX)
      LOGICAL :: MASK3(NBOX)
      REAL :: S_COND_S(NBOX)

      CHARACTER(LEN=72) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!     CHECKMD_ND set to 1/0 for whether/not to check for
!                values of ND, MD and MDT out of range.
      IF (lhook) CALL dr_hook('UKCA_AERO_STEP',zhook_in,zhook_handle)
      IF(CHECKMD_ND == 1) THEN
        PROCESS='At start of UKCA_AERO_STEP    '
! DEPENDS ON: ukca_check_md_nd
        CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
      ENDIF
!
!     VERBOSE included so that max,min,mean of ND/MD can be
!             written out and monitored after each process
!             when error checking in box model.
      IF(VERBOSE >= 2) THEN
       write(6,'(A40)') 'At start of UKCA_AERO_STEP :   ND,MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
       CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
      ENDIF
!
!     Calculate dry diameter & volume for initial remode check
! DEPENDS ON: ukca_calc_drydiam
      CALL UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)
!
!     Initial remode call (check if advection taken DRYDP out of bounds)
! DEPENDS ON: ukca_remode
      CALL UKCA_REMODE(NBOX,ND,MD,MDT,DRYDP,WETDP,VERBOSE,              &
        IMERGE,BUD_AER_MAS,N_MERGE_1D,PMID)
!
!     Recalculate dry and wet diameter and volume after re-moding
! DEPENDS ON: ukca_calc_drydiam
      CALL UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)
!
! DEPENDS ON: ukca_volume_mode
      CALL UKCA_VOLUME_MODE(NBOX,ND,MD,MDT,                             &
        RH,WVOL,WETDP,RHOPAR,IWVOLMETHOD,                               &
        DVOL,DRYDP,MDWAT,PVOL,PVOL_WAT,VERBOSE,T,PMID,S)
!
      IF(CHECKMD_ND == 1) THEN
       PROCESS='Done UKCA_REMODE (1st call)   '
! DEPENDS ON: ukca_check_md_nd
       CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
      ENDIF
!
      IF(VERBOSE >= 2) THEN
       write(6,'(A40)') 'After UKCA_REMODE1 has updated ND,MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
       CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
      ENDIF
!
!     Recalculate dry and wet diameter and volume
! DEPENDS ON: ukca_calc_drydiam
      CALL UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)
!
! DEPENDS ON: ukca_volume_mode
      CALL UKCA_VOLUME_MODE(NBOX,ND,MD,MDT,                             &
        RH,WVOL,WETDP,RHOPAR,IWVOLMETHOD,                               &
        DVOL,DRYDP,MDWAT,PVOL,PVOL_WAT,VERBOSE,T,PMID,S)
!
!     Calculate impaction scavenging of aerosol (washout)
      IF(IMSCAV_ON == 1) THEN
!
! DEPENDS ON: ukca_impc_scav
        CALL UKCA_IMPC_SCAV(NBOX,ND,MD,                                 &
       CRAIN,DRAIN,WETDP,DTC,BUD_AER_MAS)
!
        IF(CHECKMD_ND == 1) THEN
         PROCESS='Done impaction scavenging     '
! DEPENDS ON: ukca_check_md_nd
         CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
        ENDIF
!
        IF(VERBOSE >= 2) THEN
         write(6,'(A35)') 'After UKCA_IMPC_SCAV hasupdatedND'
! DEPENDS ON: ukca_calcminmaxndmdt
         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
        ENDIF
!
      ENDIF
!
!     Calculate in cloud, Nucleation Scavenging
      IF(RAINOUT_ON == 1) THEN
!
! DEPENDS ON: ukca_rainout
         CALL UKCA_RAINOUT(NBOX,ND,MD,MDT,FCONV_CONV,DELTA_Z,RHOA,      &
       CRAIN,DRAIN,CRAIN_UP,DRAIN_UP,CLWC,CLF,T,DTC,                    &
       BUD_AER_MAS,INUCSCAV,DRYDP,WETDP)
!
        IF(CHECKMD_ND == 1) THEN
          PROCESS='Done nucleation scavenging    '
! DEPENDS ON: ukca_check_md_nd
          CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
        ENDIF
!
        IF(VERBOSE >= 2) THEN
         write(6,'(A35)') 'After UKCA_RAINOUT has updated ND'
! DEPENDS ON: ukca_calcminmaxndmdt
         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
        ENDIF
!
      ENDIF
!
      IF(ICHEM == 1) THEN
!
!      Calculate aqueous chemistry
       IF(WETOX_ON == 1) THEN
!
!       Initialise variables for aqueous chemistry
        SO2     (:)=0.0
        H2O2    (:)=0.0
        DELH2O2 (:)=0.0
!
        IF(WETOX_IN_AER == 1) THEN
!
         IF(MSOTWO > 0) THEN
          MASK(:)=(S0G(:,MSOTWO) > 0.0)
          WHERE(MASK(:))
           SO2(:)=S0G(:,MSOTWO)*AIRD(:)/SM(:)
          ENDWHERE
          SO2_VMR(:)=SO2(:)/AIRD(:)
         ENDIF
!
         IF(MH2O2 > 0) THEN

! .. if using coupled H2O2
          MASK(:)=(S0G(:,MH2O2) > 0.0)
          WHERE(MASK(:))
           H2O2(:)=S0G(:,MH2O2)*AIRD(:)/SM(:)
          ENDWHERE
!
         ELSE
!
          IF(MH2O2F > 0) THEN
!
! .. if using semi-prognostic H2O2
           MASK(:)=(S0G(:,MH2O2F) > 0.0)
           WHERE(MASK(:))
            H2O2(:)=S0G(:,MH2O2F)*AIRD(:)/SM(:)
           ENDWHERE
!
          ENDIF
!
         ENDIF
!
! .. below calculates DELSO2, DELSO2_2, DELH2O2 in the
! .. aerosol module (WETOX_IN_AER=1)
!
! DEPENDS ON: ukca_wetox
         CALL UKCA_WETOX(NBOX,ND,DRYDP,DELSO2,DELSO2_2,                 &
            DELH2O2,LOWCLOUD,VFAC,SO2_VMR,H2O2,ZH2O2,ZO3,               &
            ZHO2,PMID,T,AIRD,S,DTC,LDAY,LWC,ISO2WETOXBYO3)
!
         IF(MSOTWO > 0) THEN
!
          MASK(:)=(S0G(:,MSOTWO) > 0.0)
          WHERE(MASK(:))
!
!          Update SO2 partial mass after reaction with H2O2.
           S0G(:,MSOTWO)=S0G(:,MSOTWO)-DELSO2  (:)*SM(:)/AIRD(:)
!          Update SO2 partial mass after reaction with O3.
           S0G(:,MSOTWO)=S0G(:,MSOTWO)-DELSO2_2(:)*SM(:)/AIRD(:)
!
          ENDWHERE
!
          MASK(:)=(S0G(:,MSOTWO) < 0.0)
          WHERE(MASK(:)) S0G(:,MSOTWO)=0.0
!
         ENDIF
!
!        Update H2O2 partial mass after SO2 wetox & replenishment
!                                      (both included in DELH2O2)
         IF(MH2O2 > 0) THEN
! .. if using fully coupled H2O2
          MASK(:)=(S0G(:,MH2O2) > 0.0)
          WHERE(MASK(:))
           S0G(:,MH2O2)=S0G(:,MH2O2)+DELH2O2(:)*SM(:)/AIRD(:)
          ENDWHERE
         ENDIF
!
         IF(MH2O2F > 0) THEN
! .. if using semi-prognostic H2O2
          MASK(:)=(S0G(:,MH2O2F) > 0.0)
          WHERE(MASK(:))
           S0G(:,MH2O2F)=S0G(:,MH2O2F)+DELH2O2(:)*SM(:)/AIRD(:)
          ENDWHERE
         ENDIF
!
         IF( (MH2O2 <= 0).AND.(MH2O2F <= 0) ) THEN
          cmessage = ' AQUEOUS on but MH2O2<=0 and MH2O2F <=0'
          errcode = 1
          CALL ereport('UKCA_AERO_CTL',errcode,cmessage)
         ENDIF
!
        ENDIF ! IF WETOX_IN_AER=1
! n.b. if WETOX_IN_AER.NE.1, then DELSO2,DELSO2_2 passed in to AERO_STEP
!
       ELSE
        DELSO2  (:)=0.0
        DELSO2_2(:)=0.0
       ENDIF ! if WETOX_ON=1
!
       TOTWETOX(:)=DELSO2(:)+DELSO2_2(:) ! this is zero if WETOX_ON=0
!
       IF(IACTMETHOD > 0) THEN ! if cloud processing on
!
! .. below cloud-processes those aerosol in the Aitken soluble
! .. mode which are larger than ACT to accumulation mode
! .. representing the fact that they have activated so
! .. making minimum between Aitsol and accsol respond to activation
!
! DEPENDS ON: ukca_cloudproc
         CALL UKCA_CLOUDPROC(NBOX,ND,MD,MDT,DRYDP,                      &
       LOWCLOUD,VFAC,ACT,VERBOSE,IACTMETHOD,BUD_AER_MAS)
!
       ENDIF ! IF IACTMETHOD > 0
!
       MASK1(:)=((ND(:,3)+ND(:,4)) > NUM_EPS(4))
       MASK2(:)=(MASK1(:).AND.(ND(:,3) > NUM_EPS(3)))
       MASK3(:)=(MASK1(:).AND.(ND(:,4) > NUM_EPS(4)))
!
       WHERE(MASK1(:))
!
!       Calculate number fractions to partition between acc and cor
        FRAC_AQ_ACC(:)=ND(:,3)/(ND(:,3)+ND(:,4))
        FRAC_AQ_COR(:)=ND(:,4)/(ND(:,3)+ND(:,4))
!
       ENDWHERE
!
       WHERE(MASK2(:))
!       Update soluble accum mode H2SO4 mass due to aqueous chem.
        MD(:,3,CP_SU)=(MD(:,3,CP_SU)*ND(:,3)+                           &
                       FRAC_AQ_ACC(:)*TOTWETOX(:))/ND(:,3)
        MDT(:,3)=(MDT(:,3)*ND(:,3)+FRAC_AQ_ACC(:)*TOTWETOX(:))/ND(:,3)
       ENDWHERE ! if some particles in soluble accum mode
!
       IF(NMASCLPRSUACCSOL1 > 0) THEN
        WHERE(MASK2(:))
         BUD_AER_MAS(:,NMASCLPRSUACCSOL1)=                              &
           FRAC_AQ_ACC(:)*DELSO2  (:)*SM(:)/AIRD(:)
        ENDWHERE
       ENDIF
!
       IF(NMASCLPRSUACCSOL2 > 0) THEN
        WHERE(MASK2(:))
         BUD_AER_MAS(:,NMASCLPRSUACCSOL2)=                              &
           FRAC_AQ_ACC(:)*DELSO2_2(:)*SM(:)/AIRD(:)
        ENDWHERE
       ENDIF
!
       WHERE(MASK3(:))
!       Update soluble coarse mode H2SO4 mass due to aqueous chem.
        MD(:,4,CP_SU)=(MD(:,4,CP_SU)*ND(:,4)+                           &
                       FRAC_AQ_COR(:)*TOTWETOX(:))/ND(:,4)
        MDT(:,4)=(MDT(:,4)*ND(:,4)+FRAC_AQ_COR(:)*TOTWETOX(:))/ND(:,4)
       ENDWHERE ! if some particles in soluble coarse mode
!
       IF(NMASCLPRSUCORSOL1 > 0) THEN
        WHERE(MASK3(:))
         BUD_AER_MAS(:,NMASCLPRSUCORSOL1)=                              &
           FRAC_AQ_COR(:)*DELSO2  (:)*SM(:)/AIRD(:)
        ENDWHERE
       ENDIF
!
       IF(NMASCLPRSUCORSOL2 > 0) THEN
        WHERE(MASK3(:))
         BUD_AER_MAS(:,NMASCLPRSUCORSOL2)=                              &
           FRAC_AQ_COR(:)*DELSO2_2(:)*SM(:)/AIRD(:)
        ENDWHERE
       ENDIF
!
       IF(CHECKMD_ND == 1) THEN
        PROCESS='Done aqueous phase chemistry  '
! DEPENDS ON: ukca_check_md_nd
        CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
       ENDIF
!
       IF(VERBOSE >= 2) THEN
        write(6,'(A35)') 'After UKCA_WETOX has updated MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
        CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
       ENDIF
!
      ENDIF ! ICHEM=1
!
!     Calculate aerosol dry deposition
      IF(DDEPAER_ON == 1) THEN
       IF(IDDEPAER == 1) THEN
! DEPENDS ON: ukca_ddepaer
        CALL UKCA_DDEPAER(NBOX,ND,MD,MDT,                               &
          RHOPAR,ZNOT,SEAICE,                                           &
          DTC,WETDP,USTR,PMID,PUPPER,PLOWER,T,SURTP,                    &
          RHOA,MFPA,DVISC,BUD_AER_MAS,ILSCAT)
       ENDIF
       IF(IDDEPAER == 2) THEN
! DEPENDS ON: ukca_ddepaer_incl_sedi
        CALL UKCA_DDEPAER_INCL_SEDI(NBOX,ND,MD,MDT,RHOPAR,ZNOT,         &
          DTC,WETDP,USTR,PMID,PUPPER,PLOWER,T,SURTP,SEAICE,             &
          RHOA,MFPA,DVISC,BUD_AER_MAS,JLABOVE,ILSCAT,SEDI_ON,SM)
       ENDIF
!
       IF(CHECKMD_ND == 1) THEN
        PROCESS='Done aerosol dry deposition   '
! DEPENDS ON: ukca_check_md_nd
        CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
       ENDIF
!
       IF(VERBOSE >= 2) THEN
        write(6,'(A30)') 'After DDEPAER has updated ND'
! DEPENDS ON: ukca_calcminmaxndmdt
        CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
       ENDIF
!
      ENDIF ! DDEPAER_ON=1
!
!     Loop over number of microphysics timesteps NMTS
      DO IMTS=1,NMTS
!
!     Calculate rate of condensation of gases onto existing particles
!     and nucleation rate
!
       DO JV=1,NCHEMG ! loop over gas phase components
        IF(CONDENSABLE(JV)) THEN
!
!        Set index of cpt into which condensable gas to be stored
         ICP=CONDENSABLE_CHOICE(JV)
!
!        Calculate ratio of gas phase to aerosol component molar masses
         S0G_TO_GC=MM_GAS(JV)/MM(ICP)
!
         MASK(:)=(S0G(:,JV) > 0.0)
!
         WHERE(MASK(:))
!
!         GC is gas phase molecular concentration of condensables
!         but molecules refer to molar mass of aerosol cpt
          GC(:,JV)=S0G_TO_GC*S0G(:,JV)*AIRD(:)/SM(:)
!
         ENDWHERE
!
         WHERE(.NOT.MASK(:))
!
          GC(:,JV)=0.0
!
         ENDWHERE
!
         GCOLD(:,JV)=GC(:,JV)
!
        ENDIF
       ENDDO
!
!      Split DTM substep into NZTS subsubsteps to allow for competition
!      between nucleation and condensation (and gas phase production
!      if DRYOX_IN_AER=1) to compete on short timesteps.
!
! DEPENDS ON: ukca_calc_coag_kernel
       CALL UKCA_CALC_COAG_KERNEL(NBOX,KII_ARR,KIJ_ARR,                 &
        DRYDP,DVOL,WETDP,WVOL,RHOPAR,MFPA,DVISC,T,PMID,                 &
        COAG_ON,ICOAG)
!
       DO IZTS=1,NZTS
!
        IF(VERBOSE >= 2) THEN
         WRITE(STRIZTS,'(1i2.2)') IZTS
         PROCESS='Start of loop, IZTS='//STRIZTS//' (GC)   '
! DEPENDS ON: ukca_calcminmaxgc
         CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GC   ,myproc)
         PROCESS='Start of loop, IZTS='//STRIZTS//' (GCOLD)'
! DEPENDS ON: ukca_calcminmaxgc
         CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GCOLD,myproc)
        ENDIF
!
        IF(ICHEM == 1) THEN
!
         IF(DRYOX_IN_AER > 0) THEN
!
          DO JV=1,NCHEMG
           IF(CONDENSABLE(JV)) THEN
!
!           Set index of cpt into which condensable gas will be stored
            ICP=CONDENSABLE_CHOICE(JV)
!
!           Calculate ratio of gas phase to aerosol cpt molar masses
            S0G_TO_GC=MM_GAS(JV)/MM(ICP)
!
            GC(:,JV)=GC(:,JV)+DTZ*S0G_DOT_CONDENSABLE(:,JV)             &
                             *AIRD(:)*S0G_TO_GC
!
! .. update condensable gas phase species by chemical tendencies on
! .. competition tstep
!
! .. n.b. S0G_DOT_CONDENSABLE is in vmr of S0G(:,JV) per s,
! .. so need to * by S0G_TO_GC
!
           ENDIF ! if CONDENSABLE(JV)
          ENDDO ! loop JV=1,NCHEMG
!
          IF(VERBOSE >= 2) THEN
           WRITE(STRIZTS,'(1i2.2)') IZTS
           PROCESS='Done chem upd, IZTS='//STRIZTS//' (GC)   '
! DEPENDS ON: ukca_calcminmaxgc
           CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GC   ,myproc)
           PROCESS='Done chem upd, IZTS='//STRIZTS//' (GCOLD)'
! DEPENDS ON: ukca_calcminmaxgc
           CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GCOLD,myproc)
          ENDIF
!
         ENDIF ! DRYOX_IN_AER
!
!        Carry out uptake of condensable gas phase species onto aerosol
         IF(COND_ON == 1) THEN
! DEPENDS ON: ukca_conden
          CALL UKCA_CONDEN(NBOX,GC,ND,MD,MDT,                           &
       DTZ,DRYDP,WETDP,TSQRT,RHOA,AIRDM3,DELGC_COND,                    &
       IFUCHS,AGETERM1,BUD_AER_MAS,S_COND_S,PMID,T,S,AIRD,              &
       IDCMFP,ICONDIAM)
!
          IF(VERBOSE >= 2) THEN
           WRITE(STRIZTS,'(1i2.2)') IZTS
           PROCESS='Done cond upd, IZTS='//STRIZTS//' (GC)   '
! DEPENDS ON: ukca_calcminmaxgc
           CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GC   ,myproc)
           PROCESS='Done cond upd, IZTS='//STRIZTS//' (GCOLD)'
! DEPENDS ON: ukca_calcminmaxgc
           CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GCOLD,myproc)
          ENDIF
!
          IF(CHECKMD_ND == 1) THEN
           PROCESS='Done conden of H2SO4 & SEC_ORG'
! DEPENDS ON: ukca_check_md_nd
           CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
          ENDIF
!
          IF(VERBOSE >= 2) THEN
           WRITE(6,'(A40)') 'After UKCA_CONDEN has updated MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
           CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
          ENDIF
!
         ELSE
          AGETERM1(:,:,:)=0.0
         ENDIF ! if COND_ON=1
!
!        Calculate rate of binary H2SO4-H2O nucleation
         IF(NUCL_ON == 1) THEN
!
          IF(MH2SO4 > 0) THEN
!
!          Set index of aerosol cpt which condensed H2SO4 to be stored
           CP_SO4=CONDENSABLE_CHOICE(MH2SO4)
!
           IF(CP_SO4 > 0) THEN
!
!           Calculate ratio of gas phase to aerosol cpt molar mass (H2SO4)
            S0G_TO_GC=MM_GAS(MH2SO4)/MM(CP_SO4)
!
            MASK(:)=(GC(:,MH2SO4) > 0.0)
!
            WHERE(MASK(:))
             H2SO4(:)=GC(:,MH2SO4)/S0G_TO_GC
            ENDWHERE
! .. Set (competition-step-updated) H2SO4 molecular concentration
! .. (using MM_GAS) for calculation of BHN rate in UKCA_CALCNUCRATE
!
            WHERE(.NOT.MASK(:))
             H2SO4(:)=0.0
            ENDWHERE
!
            IF(MSEC_ORG > 0) THEN
             SEC_ORG(:)=GC(:,MSEC_ORG)
            ELSE
             SEC_ORG(:)=0.0
            ENDIF

! DEPENDS ON: ukca_calcnucrate
            CALL UKCA_CALCNUCRATE(NBOX,DTZ,T,S,RH,AIRD,H2SO4,           &
                          DELH2SO4_NUCL,SEC_ORG,JRATE,BLN_ON,IBLN,      &
                          I_NUC_METHOD,                                 &
                          HEIGHT,HTPBLG,S_COND_S)
!
            MASK(:)=(H2SO4(:) > 0.0)
!
            WHERE(MASK(:))
             GC(:,MH2SO4)=H2SO4(:)*S0G_TO_GC
            ENDWHERE
! .. Update GC due to updated value of H2SO4 (bug fixed in gm5c)
!
            WHERE(.NOT.MASK(:))
             GC(:,MH2SO4)=0.0
            ENDWHERE
!
            DELGC_NUCL(:,MH2SO4)=DELH2SO4_NUCL(:)*S0G_TO_GC
!
            IF(NMASNUCLSUNUCSOL > 0)                                    &
             BUD_AER_MAS(:,NMASNUCLSUNUCSOL)=                           &
             BUD_AER_MAS(:,NMASNUCLSUNUCSOL)+DELGC_NUCL(:,MH2SO4)
!
            IF(VERBOSE >= 2) THEN
             WRITE(STRIZTS,'(1i2.2)') IZTS
             PROCESS='Done nucl upd, IZTS='//STRIZTS//' (GC)   '
!
! DEPENDS ON: ukca_calcminmaxgc
             CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GC   ,myproc)
             PROCESS='Done nucl upd, IZTS='//STRIZTS//' (GCOLD)'
! DEPENDS ON: ukca_calcminmaxgc
             CALL UKCA_CALCMINMAXGC(PROCESS,NBOX,GCOLD,myproc)
            ENDIF
           ELSE
            errcode=1
            cmessage='CP_SO4 <= 0'
            WRITE(6,'(A25,I6)') cmessage,'CP_SO4 = ',CP_SO4

            CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
           ENDIF
!
          ELSE
           errcode=1
           cmessage='MH2SO4 <= 0'
           WRITE(6,'(A25,E12.3)') cmessage,'MH2SO4 = ',MH2SO4

           CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
          ENDIF ! if MH2SO4 <=0
!
         ELSE
          IF(MH2SO4 > 0) DELGC_NUCL(:,MH2SO4)=0.0
          JRATE(:)=0.0
         ENDIF
        ELSE
         IF(MH2SO4 > 0) DELGC_NUCL(:,MH2SO4)=0.0
         JRATE(:)=0.0
        ENDIF
!
        IF(VERBOSE >= 2) THEN
         write(6,'(A35)') 'After DRYDIAM & VOL has updated ND'
! DEPENDS ON: ukca_calcminmaxndmdt
         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
        ENDIF
!
        IF((COAG_ON == 1).OR.(NUCL_ON == 1)) THEN
         IF(VERBOSE >= 2) THEN
          write(6,'(A35)') 'About to call UKCA_COAGWITHNUCL'
! DEPENDS ON: ukca_calcminmaxndmdt
          CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
         ENDIF
!
!        Update ND & MD due to combined coagulation-nucleation
! DEPENDS ON: ukca_coagwithnucl
         CALL UKCA_COAGWITHNUCL(NBOX,ND,MD,MDT,DELGC_NUCL,DTZ,JRATE,    &
       AGETERM2,INTRAOFF,INTEROFF,VERBOSE,BUD_AER_MAS,KII_ARR,KIJ_ARR)
!
         IF(CHECKMD_ND == 1) THEN
          PROCESS='Done combined coag./nucleation'
! DEPENDS ON: ukca_check_md_nd
          CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
         ENDIF
!
         IF(VERBOSE >= 2) THEN
          IF((COAG_ON == 1).AND.(NUCL_ON == 1)) THEN
           write(6,'(A40)') 'After COAG & NUCL have updatedND,MD,MDT'
          ENDIF
          IF((COAG_ON == 0).AND.(NUCL_ON == 1)) THEN
           write(6,'(A35)') 'After NUCL has updated ND,MD,MDT'
          ENDIF
          IF((COAG_ON == 1).AND.(NUCL_ON == 0)) THEN
           write(6,'(A35)') 'After COAG has updated ND,MD,MDT'
          ENDIF
! DEPENDS ON: ukca_calcminmaxndmdt
          CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
         ENDIF
        ELSE
         AGETERM2(:,:,:,:)=0.0
        ENDIF ! if COAG_ON = 1 or NUCL_ON = 1
!
!       Apply ageing -- transfer ND,MD from insol. to sol. modes
! DEPENDS ON: ukca_ageing
        CALL UKCA_AGEING(NBOX,ND,MD,MDT,                                &
          AGETERM1,AGETERM2,WETDP,VERBOSE,BUD_AER_MAS)
!
        IF(CHECKMD_ND == 1) THEN
         PROCESS='Done ageing of insoluble modes'
! DEPENDS ON: ukca_check_md_nd
         CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
        ENDIF
!
        IF(VERBOSE >= 2) THEN
         write(6,'(A40)') 'After UKCA_AGEING has updated ND,MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
         CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
        ENDIF
!
       ENDDO ! end loop over competition subsubtimesteps
!
!      Recalculate dry and wet diameter and volume
! DEPENDS ON: ukca_calc_drydiam
       CALL UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)
!
! DEPENDS ON: ukca_volume_mode
       CALL UKCA_VOLUME_MODE(NBOX,ND,MD,MDT,                            &
         RH,WVOL,WETDP,RHOPAR,IWVOLMETHOD,                              &
         DVOL,DRYDP,MDWAT,PVOL,PVOL_WAT,VERBOSE,T,PMID,S)
!
!      Apply mode-merging where necessary
! DEPENDS ON: ukca_remode
       CALL UKCA_REMODE(NBOX,ND,MD,MDT,DRYDP,WETDP,VERBOSE,             &
        IMERGE,BUD_AER_MAS,N_MERGE_1D,PMID)
!
       IF(CHECKMD_ND == 1) THEN
        PROCESS='Done UKCA_REMODE              '
! DEPENDS ON: ukca_check_md_nd
        CALL UKCA_CHECK_MD_ND(NBOX,PROCESS,ND,MD,MDT,myproc)
       ENDIF
!
       IF(VERBOSE >= 2) THEN
        write(6,'(A40)') 'After UKCA_REMODE has updated ND,MD,MDT'
! DEPENDS ON: ukca_calcminmaxndmdt
        CALL UKCA_CALCMINMAXNDMDT(NBOX,ND,MDT,myproc)
       ENDIF
!
       IF(ICHEM == 1) THEN
!
        DO JV=1,NCHEMG
         IF(CONDENSABLE(JV)) THEN
!
          ICP=CONDENSABLE_CHOICE(JV)
!
!         Calculate ratio of gas phase to aerosol cpt molar masses
          S0G_TO_GC=MM_GAS(JV)/MM(ICP)
!
          DELTAS0G(:)=(GC(:,JV)-GCOLD(:,JV))*(SM(:)/AIRD(:))/S0G_TO_GC
!
          MASK(:)=(DELTAS0G(:) < -S0G(:,JV))
! above limits deltaS0G to be -S0G (stop -ves)
          WHERE(MASK(:))
           DELTAS0G(:)=-S0G(:,JV)
          ENDWHERE
!
          S0G(:,JV)=S0G(:,JV)+DELTAS0G(:)
!
         ENDIF ! IF CONDENSABLE(JV)
        ENDDO ! loop JV=1,NCHEMG
!
       ENDIF ! if ICHEM = 1
!
      ENDDO ! end loop over NMTS

      IF (lhook) CALL dr_hook('UKCA_AERO_STEP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_AERO_STEP
