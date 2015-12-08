! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
MODULE UKCA_MODE_EMS_UM_MOD

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


!  Description:
!   Provide an aerosol primary emission array for GLOMAP-mode and store
!   emission diagnostics in the STASHwork array.
!
!  Method:
!   Emission arrays for sulphate and carbonaceous aerosols are assembled
!   from the input arrays depending on the model (chosen by i_mode_setup). 
!   A call to UKCA_MODE_EMS returns mass and number emission arrays for
!   sulphate, sea-salt, organic carbon, black carbon and dust as required.
!   The number and mass fluxes are then assembled for each tracer in the
!   array em_field_mode. Diagnostics for emitted component mass are
!   stored in the STASHwork array through calls to COPYDIAG_3D. 
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
!  University of Leeds, University of Oxford, and The Met Office.
!  See:  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90

CONTAINS

! Subroutine Interface:
      SUBROUTINE UKCA_MODE_EMS_UM(imonth, iday,                         &
                      row_length, rows, model_levels,                   &
                      n_mode_tracers, ndiv,                             &
                      n_use_emissions, n_chem_emissions,                &
                      area,                                             &
                      mass,                                             &
                      p_theta_levels,                                   &
                      t_theta_levels,                                   &
                      seaice_frac,                                      &
                      rough_length,                                     &
                      u_scalar_10m,                                     &
                      land_fraction,                                    &
                      all_emissions,                                    &
                      em_index,                                         &
                      SO2_volc_3D,                                      &
                      BC_biom_3D,                                       &
                      OC_biom_3D,                                       &
                      dust_flux,                                        &
                      em_field_mode,                                    &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                      len_stashwork38,STASHwork38)

      USE UKCA_CONSTANTS,    ONLY: avc, zboltz, mm_da
      USE UKCA_MODE_SETUP,   ONLY: ii_nd, ii_md, mode_choice,           &
                                   cp_cl, cp_bc, cp_oc, cp_su,          &
                                   cp_du, mm, mode, nmodes, ncp,        &
                                   component
      USE UKCA_OPTION_MOD,   ONLY: i_mode_setup, l_ukca_primsu,         &
                                   l_ukca_primss, l_ukca_primbcoc,      &
                                   l_ukca_primdu, mode_parfrac,         &
                                   i_mode_ss_scheme,                    &
                                   l_bcoc_ff, l_bcoc_bf, l_bcoc_bm
                                 
      USE UKCA_D1_defs,      ONLY: imode_first, ukcad1codes,            &
                                   mode_diag_sect, nmax_mode_diags
      USE ukca_mode_ems_mod, ONLY: ukca_mode_ems
      USE run_aerosol_mod,   ONLY: so2_high_level
      USE ereport_mod,       ONLY: ereport
      USE printstatus_mod
      USE um_parvars,        ONLY: at_extremity

      USE yomhook,           ONLY: lhook, dr_hook
      USE parkind1,          ONLY: jprb, jpim
! version_mod items required by cstash.h
      USE version_mod,       ONLY: nproftp, nprofdp, nprofup,           &
                                   ndiagpm, ntimep, NTimSerP,           &
                                   nlevp, npslevp, npslistp,            &
                                   outfile_s, outfile_e
      USE Submodel_Mod

      IMPLICIT NONE

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
! COMDECK CSTASH
! Description:
!  Contains variables and arrays involved in STASH request specification
!  and STASH processing in the UM. Includes namelist STASH.
!
!
! Declarations:
! Imported global variables:
!    None, but see note above.

! Global parameters:

! Global scalars:
      INTEGER      NDIAG   ! No. of diagnostics
      INTEGER      NTPROF  ! No. of time profiles
      INTEGER      NSERIES ! No. of stash time series
      INTEGER      NDPROF  ! No. of domain profiles
      INTEGER      NUPROF  ! No. of useage profiles

! Global dynamic arrays:

!   STASH specification table (JSTASH file):
!   NDIAGPM set in VERSION_MOD module
      INTEGER      MODL_B(NDIAGPM)  ! Internal model no.
      INTEGER      ISEC_B(NDIAGPM)  ! Section
      INTEGER      ITEM_B(NDIAGPM)  ! Item
      INTEGER      ITIM_B(NDIAGPM)  ! Time profile number
      INTEGER      IDOM_B(NDIAGPM)  ! Domain profile number
      INTEGER      IUSE_B(NDIAGPM)  ! Useage profile number

!   Time profile information:

      CHARACTER(LEN=8)  TIMPRO(NPROFTP)         ! Name of profile
      INTEGER      ITYP_T(NPROFTP)         ! Type of profile
      INTEGER      INTV_T(NPROFTP)         ! Time Interval
      CHARACTER(LEN=2)  UNT1_T(NPROFTP)         ! Units for time interval
      INTEGER      ISAM_T(NPROFTP)         ! Sampling period
      CHARACTER(LEN=2)  UNT2_T(NPROFTP)         ! Units for sampling period
      INTEGER      IOPT_T(NPROFTP)         ! Output option
      INTEGER      ISTR_T(NPROFTP)         ! Output Start time
      INTEGER      IEND_T(NPROFTP)         ! Output End time
      INTEGER      ISDT_T(6, NPROFTP)      ! Output Start date
      INTEGER      IEDT_T(6, NPROFTP)      ! Output End date
      INTEGER      IFRE_T(NPROFTP)         ! Output frequency
      INTEGER      IOFF_T(NPROFTP)         ! Offset for sampling
      CHARACTER(LEN=2)  UNT3_T(NPROFTP)         ! Units for output times
      INTEGER      ITIM_T(NPROFTP)         ! No. of times in times table
      INTEGER      ISER_T(NTIMEP ,NPROFTP) ! Times table (with units)
      INTEGER      MODL_T(NPROFTP)         ! Indicates internal model
                                           !  for each times table

!   Domain profile information:

      CHARACTER(LEN=8) DOMPRO  (NPROFDP)           ! Name of domain profile
      INTEGER     IOPL_D  (NPROFDP)           ! Levels option
      INTEGER     LEVB_D  (NPROFDP)           ! Bottom level
      INTEGER     LEVT_D  (NPROFDP)           ! Top level
      INTEGER     IOPA_D  (NPROFDP)           ! Area option
      INTEGER     INTH_D  (NPROFDP)           ! North boundary
      INTEGER     ISTH_D  (NPROFDP)           ! South boundary
      INTEGER     IEST_D  (NPROFDP)           ! East boundary
      INTEGER     IWST_D  (NPROFDP)           ! West boundary
      INTEGER     IMSK_D  (NPROFDP)           ! Mask type
      INTEGER     IMN_D   (NPROFDP)           ! Meaning option
      INTEGER     IWT_D   (NPROFDP)           ! Weighting option
      CHARACTER(LEN=1) TS_D    (NPROFDP)           ! Time series profile
      INTEGER     IG_TS
      INTEGER     I1_TS
      INTEGER     I51_TS
      INTEGER     BLIM_TS (NTimSerP)
      INTEGER     TLIM_TS (NTimSerP)
      REAL        BLIMR_TS(NTimSerP)
      REAL        TLIMR_TS(NTimSerP)
      INTEGER     NLIM_TS (NTimSerP)
      INTEGER     SLIM_TS (NTimSerP)
      INTEGER     ELIM_TS (NTimSerP)
      INTEGER     WLIM_TS (NTimSerP)
      INTEGER     ILEV_D  (NPROFDP)           ! Output levels code
      INTEGER     LEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      REAL       RLEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      INTEGER     PLT_D   (NPROFDP)
      INTEGER     PLLEN_D (NPROFDP)
      INTEGER     PLPOS_D (NPROFDP)
      INTEGER     PSLIST_D(NPSLEVP ,NPSLISTP)
      INTEGER     NPSLISTS
      EQUIVALENCE        (RLEVLST_D,LEVLST_D)

! Useage information:

      CHARACTER(LEN=8) USEPRO(NPROFUP)   ! Name of useage profile
      INTEGER     LOCN_U(NPROFUP)   ! Storage location of profile
      INTEGER     IUNT_U(NPROFUP)   ! Unit no.

! Information from ppxref file:

      INTEGER      MODEL_ST       ! Internal model number
      INTEGER      ISPACE         ! Space code
      INTEGER      ITIMA          ! Time availability code
      INTEGER      IGP            ! Grid of data code
      INTEGER      ILEV           ! Level type code
      INTEGER      IBOT           ! First level code
      INTEGER      ITOP           ! Last level code
      INTEGER      IFLAG          ! Level compression flag
      INTEGER      IOPN(6)        ! Sectional option code
      INTEGER      VMSK           ! Integer equiv of bin vers mask
      INTEGER      IPSEUDO        ! Pseudo dimension type
      INTEGER      IPFIRST        ! First pseudo dim code
      INTEGER      IPLAST         ! Last pseudo dim code
      INTEGER      PTR_PROG       ! Section zero point back
      INTEGER      HALO_TYPE      ! Type of halo the field has

! PP output file units
      INTEGER      PPlen2LkUp(OUTFILE_S:OUTFILE_E)
      CHARACTER(LEN=1)  FTOutUnit (OUTFILE_S:OUTFILE_E)

! COMMON blocks:
      COMMON/STCHA/ TIMPRO,UNT1_T,UNT2_T,UNT3_T,DOMPRO,TS_D,            &
     &  USEPRO,FTOutUnit

      COMMON/STSH/                                                      &
     &  NDIAG   ,MODL_B  ,ISEC_B ,ITEM_B  ,ITIM_B  ,IDOM_B  ,IUSE_B,    &
     &  NTPROF  ,ITYP_T  ,INTV_T ,ISAM_T  ,ITIM_T  ,                    &
     &  IOPT_T  ,ISTR_T  ,IEND_T ,IFRE_T  ,IOFF_T, ISER_T  ,MODL_T  ,   &
     &  NDPROF  ,IOPL_D  ,LEVB_D ,ISDT_T  ,IEDT_T  ,                    &
     &  IOPA_D  ,INTH_D  ,ISTH_D ,IEST_D  ,IWST_D  ,                    &
     &  IMSK_D  ,IMN_D   ,IWT_D  ,                                      &
     &  LEVT_D  ,LEVLST_D,                                              &
     &  PLT_D   ,PLLEN_D ,PLPOS_D,PSLIST_D,NPSLISTS,                    &
     &  BLIM_TS ,TLIM_TS ,BLIMR_TS,TLIMR_TS,IG_TS   ,I1_TS   ,          &
     &  NLIM_TS ,SLIM_TS ,ELIM_TS ,WLIM_TS ,I51_TS  ,NSERIES ,          &
     &  NUPROF  ,LOCN_U  ,IUNT_U ,                                      &
     &  MODEL_ST,ISPACE  ,ITIMA  ,IGP     ,                             &
     &  ILEV    ,IBOT    ,ITOP   ,IFLAG   ,IOPN    ,VMSK    ,           &
     &  IPSEUDO ,IPFIRST ,IPLAST ,PTR_PROG, HALO_TYPE,                  &
     & PPlen2LkUp

! CSTASH end
! TYPSTS starts
! submodel_mod must be included before this file
!Applicable to all configurations
!STASH related variables for describing output requests and space
!management.
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS


! TYPSTSZ end

! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! Whether a calculation is needed for SF above
      LOGICAL :: SF_CALC(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
! TYPSTS end

! Inputs
      INTEGER, INTENT(IN) :: imonth
                                   ! Month
      INTEGER, INTENT(IN) :: iday                                     
! Day
      INTEGER, INTENT(IN) :: row_length                               
! No of pts in a row
      INTEGER, INTENT(IN) :: rows                                     
! No of rows
      INTEGER, INTENT(IN) :: model_levels                             
! No of model levels
      INTEGER, INTENT(IN) :: n_mode_tracers                           
! No of mode tracers
      INTEGER, INTENT(IN) :: ndiv                                     
! No of dust divisions
      INTEGER, INTENT(IN) :: n_use_emissions                          
! No of emissions used
      INTEGER, INTENT(IN) :: n_chem_emissions                         
! No of chem emissions
      INTEGER, INTENT(IN) :: em_index(n_chem_emissions)               
! Index to item no
      INTEGER, INTENT(IN) :: len_stashwork38                          
! Length of diag. array

      REAL, INTENT(IN) :: area(row_length,rows,model_levels)          
! area m^2
      REAL, INTENT(IN) :: mass(row_length,rows, model_levels)         
! mass of each cell
      REAL, INTENT(IN) :: p_theta_levels(row_length,rows,model_levels) 
! pressure
      REAL, INTENT(IN) :: t_theta_levels(row_length,rows,model_levels) 
! temperature
      REAL, INTENT(IN) :: all_emissions(row_length,rows,n_chem_emissions) 
! 2D emissions fields
      REAL, INTENT(IN) :: so2_volc_3d(row_length,rows,model_levels)   
! 3D volcanic SO2 emiss
      REAL, INTENT(IN) :: bc_biom_3d(row_length,rows,model_levels)    
! Biomass burn BC emiss
      REAL, INTENT(IN) :: oc_biom_3d(row_length,rows,model_levels)    
! Biomass burn OC emiss
      REAL, INTENT(IN) :: dust_flux(row_length,rows,ndiv)             
! dust emissions (kg/m2/s)
      REAL, INTENT(IN) :: seaice_frac(row_length, rows)               
! sea ice
      REAL, INTENT(IN) :: u_scalar_10m(row_length, rows)              
! wind at 10m
      REAL, INTENT(IN) :: rough_length(row_length, rows)              
! roughness length
      REAL, INTENT(IN) :: land_fraction(row_length,rows)              
! land_fraction

      REAL, INTENT(OUT) :: em_field_mode(row_length,rows,               &
                                         model_levels,n_mode_tracers) 
! emissions
      REAL, INTENT(INOUT) :: STASHwork38(len_stashwork38)             
! Diagnostic array

! Local variables
      INTEGER :: primsu_on
! Switch for whether primary sulfate particle emissions are on/off
      INTEGER :: primbcoc_on
! Switch for whether primary carbonaceous particle emissions are on/off
      INTEGER :: primdu_on
! Switch for whether primary dust particle emissions are on/off
      INTEGER :: primss_on
! Switch for whether primary sea-salt ptcl emissions are on/off
      INTEGER :: verbose
! Switch to determine level of debug o/p (0=none, 1, 2)
! Now set via printstatus from UMUI Output Choices panel, see below
      INTEGER, PARAMETER :: idustems=0
! Switch for using Pringle scheme (=1) or AEROCOMdaily (=2)

      REAL, PARAMETER  :: emfactor=1.4
! To convert emissions of OC and BC from kg[OM] to kg[C]

      REAL     :: emanso2(row_length,rows,6)
! Anthrop. SO2 ems rates, low sources (kgSO2/m2/s)
      REAL     :: emvolconso2(row_length,rows,model_levels)
! Volcanic SO2 ems rates (cont. src) (kgSO2/m2/s)
      REAL     :: emvolexpso2(row_length,rows,model_levels)
! Volcanic SO2 ems rates (expl. src) (kgSO2/m2/s)
      REAL     :: embiomso2(row_length,rows,model_levels)
! Biomass SO2 ems rates (kgSO2/box/s)
      REAL     :: emc(row_length,rows,4)
! BC/OC emission rates from bio- & fossil-fuels (kgC/m2/s)
      REAL     :: emcbm(row_length,rows,model_levels,2)
! BC/OC emission rates from biomass burning (kgC/m2/s)
      REAL     :: emcdu(row_length,rows,ndiv)
! Dust emission rates (kg m-2 s-1)
      REAL     :: aird(row_length,rows,model_levels)
! Number density of air (per cm3)
      REAL     :: rhoa(row_length,rows,model_levels)
! Air density (kg/m3)

! Now use ISO2EMS=1 when chosen only sulphate and sea-salt in 4 modes
!                   (no BC/OC) -- i_mode_setup=1 -- small sizes
!                   make up for lack of primary BC/OC emissions
! 
! For all other options use ISO2EMS=3 as standard as in GLOMAP.
! N.B. parfrac variable now determines fraction of SO2 emissions as SO4.
!
! ISO2EMS == 1 : 3.0% of SO2 ems --> primary SO4 ems --> 15%/85% to
!                10/70 nm g.m.diam. modes as in Spracklen (2005)
!                and Binkowski & Shankar (1995).
!
! ISO2EMS == 2 : 2.5% of SO2 ems --> primary SO4 ems
!                road/off-road/domestic      all -->   30nm gm.diam mode
!                industrial/power-plant/ship all --> 1000nm gm.diam mode
!                as for original AEROCOM size recommendations.
!
! ISO2EMS >= 3 : 2.5% of SO2 ems --> primary SO4 ems --> 50%/50% to
!                150/1500nm  g.m.diam. modes as for Stier et al (2005)
!                modified AEROCOM sizdis recommendations.
!
! Note also that currently, ISO2EMS also controls size assumptions for
! primary carbonaceous aerosol (this needs to be changed):
!
! ISO2EMS /= 3 : biofuel & biomass BC/OC emissions -->  80nm g.m.diam.
!                fossil-fuel BC/OC emissions --> 30nm g.m.diam.
!
! ISO2EMS == 3 : biofuel & biomass BC/OC emissions --> 150nm g.m.diam.
!                fossil-fuel BC/OC emissions --> 60nm g.m.diam.
! See ISO2EMS settings below...

      REAL    :: parfrac
! PARFRAC is fraction of mass of SO2 emitted as primary SO4

      INTEGER, SAVE :: lbcff            ! Index for BC fossil fuel emissions
      INTEGER, SAVE :: lbcbf            ! Index for BC biofuel emissions
      INTEGER, SAVE :: locff            ! Index for OC fossil fuel emissions
      INTEGER, SAVE :: locbf            ! Index for OC biofuel emissions
      INTEGER, SAVE :: lso2emlo         ! Index for SO2-low emissions
      INTEGER, SAVE :: lso2emhi         ! Index for SO2-low emissions
      INTEGER :: itra                   ! loop counter
      INTEGER :: i                      ! loop counter
      INTEGER :: icp                    ! loop counter for components
      INTEGER :: item1                  ! stash item number
      INTEGER :: section                ! stash section number
      INTEGER :: imode                  ! loop counter for modes
      INTEGER :: im_index               ! internal modal index
      INTEGER :: iso2ems                ! determines emission assumptions
                                        !  (see below)
      INTEGER :: n_reqd_tracers         ! No of required tracers

      CHARACTER(LEN=72) :: cmessage     ! Error message
      INTEGER :: errcode
      INTEGER :: ierr
      LOGICAL, SAVE :: firstcall=.TRUE. ! Indicates first call

! Emission fluxes for mass and number in kg/m^2/s (mass) or
!  equivalent kg/m^2/s (no)
      REAL :: aer_mas_primsu(row_length,rows,model_levels,nmodes)
      REAL :: aer_mas_primbc(row_length,rows,model_levels,nmodes)
      REAL :: aer_mas_primoc(row_length,rows,model_levels,nmodes)
      REAL :: aer_mas_primss(row_length,rows,model_levels,nmodes)
      REAL :: aer_mas_primdu(row_length,rows,model_levels,nmodes)
      REAL :: aer_num_primsu(row_length,rows,model_levels,nmodes)
      REAL :: aer_num_primcar(row_length,rows,model_levels,nmodes)
      REAL :: aer_num_primss(row_length,rows,model_levels,nmodes)
      REAL :: aer_num_primdu(row_length,rows,model_levels,nmodes)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!-------------------------------------

      IF (lhook) CALL dr_hook('UKCA_MODE_EMS_UM',zhook_in,zhook_handle)

      im_index = internal_model_index(atmos_im)

! Set VERBOSE from printstatus as defined in UMUI
      verbose = MIN(printstatus-1,2)  ! sets it to 0-2

! below are parameters set in the UKCA panel for MODE

      IF (firstcall .AND. printstatus >= prstatus_normal) THEN
        WRITE(6,'(A20,I10)') 'i_mode_setup =      ', i_mode_setup
        WRITE(6,'(A20,L10)') 'L_UKCA_PRIMSU =     ', l_ukca_primsu
        WRITE(6,'(A20,F10.3)') 'mode_parfrac =      ', mode_parfrac
        WRITE(6,'(A20,L10)') 'L_UKCA_PRIMSS =     ', l_ukca_primss
        WRITE(6,'(A20,I10)') 'i_mode_ss_scheme =  ', i_mode_ss_scheme
        WRITE(6,'(A20,L10)') 'L_UKCA_PRIMBCOC =   ', l_ukca_primbcoc
        WRITE(6,'(A20,L10)') 'L_BCOC_ff =         ', l_bcoc_ff
        WRITE(6,'(A20,L10)') 'L_BCOC_bf =         ', l_bcoc_bf
        WRITE(6,'(A20,L10)') 'L_BCOC_bm =         ', l_bcoc_bm
        WRITE(6,'(A20,L10)') 'L_UKCA_PRIMDU =     ', l_ukca_primdu
        WRITE(6,'(A20,I10)') 'rl*r*ml     =       ', row_length*rows*   &
                                                     model_levels
      END IF

! set ISO2EMS according to i_mode_setup in UMUI
      IF (i_mode_setup == 1) iso2ems=1       ! SUSS_4mode
      IF (i_mode_setup == 2) iso2ems=3       ! SUSSBCOC_5mode
      IF (i_mode_setup == 3) iso2ems=3       ! SUSSBCOC_4mode
      IF (i_mode_setup == 4) iso2ems=3       ! SUSSBCOCSO_5mode
      IF (i_mode_setup == 5) iso2ems=3       ! SUSSBCOCSO_4mode
      IF (i_mode_setup == 6) iso2ems=0       ! DUonly_2mode
      IF (i_mode_setup == 7) iso2ems=0       ! DUonly_3mode
      IF (i_mode_setup == 8) iso2ems=3       ! SUSSBCOCDU_7mode
      IF (i_mode_setup == 9) iso2ems=3       ! SUSSBCOCDU_4mode
 
! set PRIMSU_ON according to L_UKCA_PRIMSU in UMUI
      IF (l_ukca_primsu) THEN
        primsu_on = 1
      ELSE
        primsu_on = 0
      END IF
 
! set PARFRAC according to mode_parfrac in UMUI
      parfrac=mode_parfrac/100.0
 
! set PRIMSS_ON according to L_UKCA_PRIMSS in UMUI
      IF (l_ukca_primss) THEN
        primss_on = 1
      ELSE
        primss_on = 0
      END IF
 
      IF (firstcall .AND. primss_on == 1) THEN
        IF (i_mode_ss_scheme == 1 .AND. printstatus >= prstatus_normal) &
          WRITE(6,'(A35)') 'Gong-Monahan SS scheme selected'
        IF (i_mode_ss_scheme >= 2) THEN
          cmessage='Bad value of I_MODE_SS_SCHEME'
          WRITE(6,'(A72,I10)') cmessage,i_mode_ss_scheme
          ierr = 1
          CALL EREPORT('UKCA_MODE_EMS_UM',ierr,cmessage)
        END IF
      END IF
 
! set PRIMBCOC_ON according to L_UKCA_PRIMBCOC in UMUI
      IF (l_ukca_primbcoc) THEN
        primbcoc_on = 1
      ELSE
        primbcoc_on = 0
      END IF

! set PRIMDU_ON according to L_UKCA_PRIMDU in UMUI
      IF (l_ukca_primdu) THEN
        primdu_on = 1
      ELSE
        primdu_on = 0
      END IF

      IF (L_UKCA_PRIMDU) THEN
        IF (printstatus >= prstatus_normal) THEN 
          WRITE(6,'(A30)') 'Dust emissions switched on'
        END IF
      ELSE
        IF (printstatus >= prstatus_normal)                             &
          WRITE(6,'(A30)') 'Dust emissions switched off'
      END IF
 

      IF (printstatus >= prstatus_diag) THEN

        WRITE(6,'(A18)') 'UKCA_MODE_EMS_UM:'

        WRITE(6,'(A18)') 'UKCA_MODE EMS_UM inputs: '
        WRITE(6,'(A18,I10)') 'i_mode_setup: ',i_mode_setup
        WRITE(6,'(A18,L10)') 'L_UKCA_PRIMSU: ',l_ukca_primsu
        WRITE(6,'(A18,F10.3)') 'mode_parfrac: ',mode_parfrac
        WRITE(6,'(A18,L10)') 'L_UKCA_PRIMSS=',l_ukca_primss
        WRITE(6,'(A18,I10)') 'i_mode_ss_scheme=',i_mode_ss_scheme
        WRITE(6,'(A18,L10)') 'L_UKCA_PRIMBCOC=',l_ukca_primbcoc
        WRITE(6,'(A18,L10)') 'L_BCOC_ff=',l_bcoc_ff
        WRITE(6,'(A18,L10)') 'L_BCOC_bf=',l_bcoc_bf
        WRITE(6,'(A18,L10)') 'L_BCOC_bm=',l_bcoc_bm
        WRITE(6,'(A18,L10)') 'L_UKCA_PRIMDU=',l_ukca_primdu

        WRITE(6,'(A15,I10)') 'imonth: ',imonth
        WRITE(6,'(A15,I10)') 'iday: ',iday
        WRITE(6,'(A15,I10)') 'row_length: ',row_length
        WRITE(6,'(A15,I10)') 'rows: ',rows
        WRITE(6,'(A15,I10)') 'model_levels: ',model_levels
        WRITE(6,'(A15,I10)') 'n_mode_tracers: ',n_mode_tracers

        WRITE(6,'(A15,3E12.4)') 'SO2_volc_3d: ',MINVAL(so2_volc_3d),    &
                                   MAXVAL(so2_volc_3d),                 &
                                   SUM(so2_volc_3d)/SIZE(so2_volc_3d)
        WRITE(6,'(A15,3E12.4)') 'BC_biom_3D,: ',MINVAL(BC_biom_3D),     &
                                   MAXVAL(BC_biom_3D),                  &
                                   SUM(BC_biom_3D)/SIZE(BC_biom_3D)
        WRITE(6,'(A15,3E12.4)') 'OC_biom_3D,: ',MINVAL(OC_biom_3D),     &
                                   MAXVAL(OC_biom_3D),                  &
                                   SUM(OC_biom_3D)/SIZE(OC_biom_3D)
        WRITE(6,'(A15)') 'all_emissions:'
        DO i=1,n_chem_emissions
          WRITE(6,'(3E12.4)')  MINVAL(all_emissions(:,:,i)),            &
                   MAXVAL(all_emissions(:,:,i)),                        &
                   SUM(all_emissions(:,:,i))/SIZE(all_emissions(:,:,i))
        END DO

      END IF    ! IF (printst >=....)

      IF (firstcall) THEN

! .. Initialise indices
        lso2emlo = imdi
        lso2emhi = imdi
        lbcff    = imdi
        lbcbf    = imdi
        locff    = imdi
        locbf    = imdi

! .. Set indices for emissions array
        DO itra=1,n_chem_emissions
          IF (em_index(itra) == 58)  lso2emlo = itra
          IF (em_index(itra) == 126) lso2emhi = itra
          IF (em_index(itra) == 310) lbcff    = itra
          IF (em_index(itra) == 311) lbcbf    = itra
          IF (em_index(itra) == 312) locff    = itra
          IF (em_index(itra) == 313) locbf    = itra
        END DO

! .. Check that these are all set now
        errcode = 0
        IF (lso2emlo == imdi) errcode = errcode + 1
        IF (lso2emhi == imdi) errcode = errcode + 2
        IF (lbcff == imdi) errcode = errcode + 4
        IF (lbcbf == imdi) errcode = errcode + 8
        IF (locff == imdi) errcode = errcode + 16
        IF (locbf == imdi) errcode = errcode + 32

        IF (errcode > 0) THEN
          cmessage = ' Indices for emission array not identified'
          WRITE(6,'(A72,I10)') cmessage,errcode
          CALL EREPORT('UKCA_MODE_EMS_UM',errcode,cmessage)
        END IF

        IF(firstcall) firstcall=.FALSE.
      END IF   ! firstcall

! Derived quantities
! ==================
      aird(:,:,:) = p_theta_levels(:,:,:)/(t_theta_levels(:,:,:)*       &
                                   zboltz*1E6)           ! No. density (/cm^3)

! Surface Emissions
! =================
      emanso2(:,:,:) = 0.0
      emvolconso2(:,:,:) = 0.0
      emvolexpso2(:,:,:) = 0.0
      embiomso2(:,:,:) = 0.0
      IF (primsu_on == 1) THEN      ! do primary S emissions
! .. Now required in units of kgSO2/m2/s
        emanso2(:,:,2) = all_emissions(:,:,lso2emlo)

! .. Add high-level primary SO4 emissions to type 4
        emanso2(:,:,4) = all_emissions(:,:,lso2emhi)

! .. Add volcanic SO4 emissions
        emvolconso2(:,:,:) = SO2_volc_3D(:,:,:)

      END IF ! if primsu_on = 1

      emc(:,:,:)=0.0
      emcbm(:,:,:,:)=0.0
      IF (primbcoc_on == 1) THEN
! Set carbonaceous aerosol emissions:
! ..  emissions(:,lBCbf/lBCff/lOCbf/lOCff/BC_biom_3D/OC_biom_3D) are in kgC/m2/s

        IF(L_bcoc_bf) THEN
         IF(printstatus >= prstatus_oper)                               &
           WRITE(6,'(A40)') 'Setting biofuel BC/OC emissions'
         emc(:,:,1) = all_emissions(:,:,lbcbf)           ! 2D BC biofuel emiss.
         emc(:,:,3) = all_emissions(:,:,locbf)           ! 2D OC biofuel emiss.
        ELSE
         IF(printstatus >= prstatus_oper)                               &
           WRITE(6,'(A40)') 'Not set biofuel BC/OC emissions'
        END IF
        IF(L_bcoc_ff) THEN
         IF(printstatus >= prstatus_oper)                               &
           WRITE(6,'(A40)') 'Setting fossil-fuel BC/OC emissions'
         emc(:,:,2) = all_emissions(:,:,lbcff)       ! 2D BC fossil fuel emiss.
         emc(:,:,4) = all_emissions(:,:,locff)       ! 2D OC fossil fuel emiss.
        ELSE
         IF(printstatus >= prstatus_oper)                               &
           WRITE(6,'(A40)') 'Not set fossil-fuel BC/OC emissions'
        END IF
        IF(L_bcoc_bm) THEN
         IF(printstatus >= prstatus_oper)                               &
           WRITE(6,'(A40)') 'Setting biomass burning BC/OC ems'
         emcbm(:,:,:,1) = BC_biom_3D(:,:,:)
         emcbm(:,:,:,2) = OC_biom_3D(:,:,:)/emfactor

        ELSE
         IF(printstatus >= prstatus_oper)                               &
           WRITE(6,'(A40)') 'Not set biomass burning BC/OC ems'
        END IF

      END IF  ! primbcoc_on = 1

      emcdu(:,:,:) = 0.0
      IF (primdu_on == 1) THEN

        emcdu(:,:,:) = dust_flux(:,:,:)

      END IF ! PRIMDU_ON=1


!  Call routine to return number and mass emission fluxes 
!  ------------------------------------------------------

      CALL ukca_mode_ems(row_length, rows, model_levels, ndiv,          &
          parfrac, mass, land_fraction, area, rough_length,             &
          seaice_frac, aird, u_scalar_10m,                              &
          emanso2, embiomso2, emvolconso2, emvolexpso2, so2_high_level, &
          emc, emcbm, emcdu,                                            &
          iso2ems, i_mode_ss_scheme,                                    &
          primsu_on, primbcoc_on, primss_on, primdu_on,verbose,         &
          aer_mas_primsu, aer_mas_primbc,                               &
          aer_mas_primoc, aer_mas_primss, aer_mas_primdu,               &
          aer_num_primsu, aer_num_primcar,                              &
          aer_num_primss, aer_num_primdu)

       itra=0

       em_field_mode(:,:,:,:) = 0.0

       DO imode = 1,nmodes
         IF (mode(imode)) THEN
           itra = ii_nd(imode)

! .. First set em_field_mode for total number flux from all primary emissions
            em_field_mode(:,:,:,itra) =                                 &
                                     aer_num_primsu(:,:,:,imode) +      &
                                     aer_num_primcar(:,:,:,imode) +     &
                                     aer_num_primss(:,:,:,imode) +      &
                                     aer_num_primdu(:,:,:,imode)

! .. Then set em_field_mode for component mass flux emissions
!     convert mass flux to be mass of component rather than mass of air
!     (units are mass mixing ratios not volume mixing ratios)
         DO icp = 1,ncp
           IF (component(imode,icp)) THEN
             itra = ii_md(imode,icp)

! .. Set emission flux for primary sulphate emissions
             IF (icp == cp_su) em_field_mode(:,:,:,itra) =              &
                        aer_mas_primsu(:,:,:,imode)*mm(cp_su)/mm_da

! .. Set emission flux for primary BC emissions to insoluble mode
             IF (icp == cp_bc) em_field_mode(:,:,:,itra) =              &
                        aer_mas_primbc(:,:,:,imode)*mm(cp_bc)/mm_da

! .. Set emission flux for primary OM emissions
             IF (icp == cp_oc) em_field_mode(:,:,:,itra) =              &
                        aer_mas_primoc(:,:,:,imode)*mm(cp_oc)/mm_da

! .. Set emission flux for primary sea-spray emissions
             IF (icp == cp_cl) em_field_mode(:,:,:,itra) =              &
                        aer_mas_primss(:,:,:,imode)*mm(cp_cl)/mm_da

! .. Set emission flux for primary dust emissions
             IF (icp == cp_du) em_field_mode(:,:,:,itra) =              &
                        aer_mas_primdu(:,:,:,imode)*mm(cp_du)/mm_da

          END IF   ! if component
         END DO   ! loop over components
        END IF   ! if MODE
       END DO   ! loop over modes

! Diagnostics for emitted SO4 mass (mol/gridbox/s)
! ------------------------------------------------

      section = mode_diag_sect
      DO item1 = 201,203
        IF (sf(item1,section)) THEN
          imode = item1 - 199     ! Aitken, accumulation and coarse (SOL)
! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork38(si(item1,section,im_index)),     &
               aer_mas_primsu(:,:,:,imode)*area(:,:,:)/mm(cp_su),       &
               row_length,rows,model_levels,0,0,0,0, at_extremity,      &
               stlist(1,stindex(1,item1,section,im_index)),len_stlist,  &
               stash_levels,num_stash_levels+1,                         &
               atmos_im,section,item1,ierr,cmessage)

          IF (ierr >  0)                                               &
            CALL EREPORT('UKCA_MODE_EMS_UM',section*1000+item1,cmessage)
        END IF
      END DO

! Diagnostics for emitted sea-salt mass (mol/gridbox/s)
! -----------------------------------------------------

      section = mode_diag_sect
      DO item1 = 204,205
        IF (sf(item1,section)) THEN
          imode = item1 - 201      ! accumulation and coarse (SOL)
! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork38(si(item1,section,im_index)),     &
               aer_mas_primss(:,:,:,imode)*area(:,:,:)/mm(cp_cl),       &
               row_length,rows,model_levels,0,0,0,0, at_extremity,      &
               stlist(1,stindex(1,item1,section,im_index)),len_stlist,  &
               stash_levels,num_stash_levels+1,                         &
               atmos_im,section,item1,ierr,cmessage)

          IF (ierr >  0)                                               &
            CALL EREPORT('UKCA_MODE_EMS_UM',section*1000+item1,cmessage)
        END IF
      END DO

! Diagnostics for emitted BC mass (mol/gridbox/s)
! -----------------------------------------------

      section = mode_diag_sect
      DO item1 = 206,207
        IF (sf(item1,section)) THEN
          IF (item1 == 206) imode = 2                ! Aitken SOL
          IF (item1 == 207) imode = 5                ! Aitken INSOL
! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork38(si(item1,section,im_index)),     &
             aer_mas_primbc(:,:,:,imode)*area(:,:,:)/mm(cp_bc),         &
             row_length,rows,model_levels,0,0,0,0, at_extremity,        &
             stlist(1,stindex(1,item1,section,im_index)),len_stlist,    &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,section,item1,ierr,cmessage)

          IF (ierr >  0)                                               &
            CALL EREPORT('UKCA_MODE_EMS_UM',section*1000+item1,cmessage)
        END IF
      END DO

! Diagnostics for emitted OC mass (mol/gridbox/s)
! -----------------------------------------------

      section = mode_diag_sect
      DO item1 = 208,209
        IF (sf(item1,section)) THEN
          IF (item1 == 208) imode = 2                ! Aitken SOL
          IF (item1 == 209) imode = 5                ! Aitken INSOL
! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork38(si(item1,section,im_index)),     &
             aer_mas_primoc(:,:,:,imode)*area(:,:,:)/mm(cp_oc),         &
             row_length,rows,model_levels,0,0,0,0, at_extremity,        &
             stlist(1,stindex(1,item1,section,im_index)),len_stlist,    &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,section,item1,ierr,cmessage)

          IF (ierr >  0)                                               &
           CALL EREPORT('UKCA_MODE_EMS_UM',section*1000+item1,cmessage)
        END IF
      END DO

! Diagnostics for emitted dust mass
! ---------------------------------

      section = mode_diag_sect
      DO item1 = 210,212
        IF (sf(item1,section)) THEN
          imode = item1 - 205       ! Aitken, accummulation, and coarse INSOL
! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork38(si(item1,section,im_index)),     &
               aer_mas_primdu(:,:,:,imode)*area(:,:,:)/mm(cp_du),       &
               row_length,rows,model_levels,0,0,0,0, at_extremity,      &
               stlist(1,stindex(1,item1,section,im_index)),len_stlist,  &
               stash_levels,num_stash_levels+1,                         &
               atmos_im,section,item1,ierr,cmessage)

          IF (ierr >  0)                                               &
            CALL EREPORT('UKCA_MODE_EMS_UM',section*1000+item1,cmessage)
        END IF
      END DO


      IF (lhook) CALL dr_hook('UKCA_MODE_EMS_UM',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE UKCA_MODE_EMS_UM
END MODULE UKCA_MODE_EMS_UM_MOD
