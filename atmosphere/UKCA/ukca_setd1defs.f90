! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
! Description:
!  Subroutine to define fields required from D1 by UKCA.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!  Set D1 section and item codes depending on the selected
!  chemistry.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 vn8.2 programming standards.
!
! ----------------------------------------------------------------------
!
      SUBROUTINE UKCA_SETD1DEFS(row_length,rows,n_rows,model_levels,   &
                 bl_levels,tr_levels,wet_levels,land_pts,sm_levels,    &
                 ntiles,tr_ukca)

      USE ukca_d1_defs
      USE ukca_option_mod, ONLY: L_ukca_rado3, L_ukca_radch4,          &
                                 L_ukca_mode, L_ukca_chem, L_ukca_dust,&
                                 L_ukca_qch4inter, i_ukca_photol,       &
                                 L_ukca_arg_act, L_ukca_radn2o,        &
                                 L_ukca_ageair, L_ukca_achem,          &
                                 L_ukca_aerchem, L_ukca_tropisop,      &
                                 L_ukca_trop, L_ukca_std_trop,         &
                                 L_ukca_stratcfc, L_ukca_trophet,      &
                                 L_ukca_strattrop, L_ukca_raq,         &
                                 L_ukca_strat, jpctr, l_ukca 
      USE ukca_photo_scheme_mod, ONLY: i_ukca_fastjx
      USE asad_mod,        ONLY: advt, nadvt, speci
      USE nstypes,         ONLY: ntype, npft
      USE parkind1,        ONLY: jprb, jpim
      USE yomhook,         ONLY: lhook, dr_hook
      USE ereport_mod,     ONLY: ereport
      USE PrintStatus_mod
      USE Control_Max_Sizes
      USE dec_spec
      USE spec_sw_lw
      USE um_input_control_mod,  ONLY: l_so2, l_dms, l_nh3, l_sulpc_so2, &
                                       l_so2_surfem, l_so2_hilem,        &
                                       l_so2_natem, l_dms_em, l_nh3_em
      USE rad_input_mod, ONLY: l_sw_radiate_prog

! version_mod items required by cstash.h and model.h
      USE version_mod,     ONLY: nproftp, nprofdp, nprofup,            &
                                 ndiagpm, ntimep, NTimSerP,            &
                                 nlevp, npslevp, npslistp,             &
                                 outfile_s, outfile_e, nsectp

      USE switches,        ONLY: l_ctile
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
! MODEL Defines model-dependent quantities used by data addressing and
! STASH
!
! submodel_mod must be used before this one
! VERSION_MOD module is required for nsectp, outfile_s and outfile_e
!
      INTEGER, PARAMETER :: AASSETS    = 9
      INTEGER, PARAMETER :: MEAD_TYPES = 4
      INTEGER, PARAMETER :: A_MAX_TRVARS=150 !Max.no.of tracers allowed
      INTEGER, PARAMETER :: A_MAX_UKCAVARS=150 ! Max.no.of UKCA allowed
      INTEGER, PARAMETER :: MAX_AOBS=100

      REAL :: H_A_EWSPACE
      REAL :: H_A_NSSPACE
      REAL :: H_A_FIRSTLAT
      REAL :: H_A_FIRSTLONG
      REAL :: H_A_POLELAT
      REAL :: H_A_POLELONG

      INTEGER :: H_A_GROUP
      INTEGER :: H_OROG_ROUGH
      INTEGER :: A_ASSMGRPS
      INTEGER :: NUM_PVPR

      LOGICAL :: A_RECON
      LOGICAL :: H_OROG_GRAD
      LOGICAL :: ATMODS
      LOGICAL :: CMODS
      LOGICAL :: LMESO

      LOGICAL :: TRACER_A (0:A_MAX_TRVARS)
      LOGICAL :: TR_UKCA_A (0:A_MAX_UKCAVARS)
      LOGICAL :: AASSET   (AASSETS)
      INTEGER :: AASPF    (AASSETS)
      INTEGER :: AASPL    (AASSETS)
      INTEGER :: AOBINC   (MAX_AOBS)
      INTEGER :: AOBGRP   (MAX_AOBS)
      INTEGER :: RUN_TARGET_END( 6)

      COMMON/MODELA/ H_A_EWSPACE,H_A_NSSPACE,H_A_FIRSTLAT,H_A_FIRSTLONG,&
     &  H_A_POLELAT,H_A_POLELONG,A_ASSMGRPS,NUM_PVPR ,A_RECON,H_A_GROUP,&
     &  H_OROG_GRAD,ATMODS,CMODS,LMESO,TRACER_A,TR_UKCA_A,              &
     &  AASSET,AASPF,AASPL

!Total data length for primary fields for each submodel data partition
      INTEGER      LPRIM(N_SUBMODEL_PARTITION_MAX)
! Global (ie. dump on disk) version of LPRIM
      INTEGER      global_LPRIM(N_SUBMODEL_PARTITION_MAX)
!Total data length for primary fields for each internal model
      INTEGER      LPrimIM(N_INTERNAL_MODEL_MAX)
!Total data length for diagnostic flds for each submodel data partition
! Global (ie. dump on disk) version of LPrimIM
      INTEGER      global_LPrimIM(N_INTERNAL_MODEL_MAX)
      INTEGER      LDUMP(N_SUBMODEL_PARTITION_MAX)
! Global (ie. dump on disk) version of LDUMP
      INTEGER      global_LDUMP(N_SUBMODEL_PARTITION_MAX)
!Total data length for diagnostic flds for each internal model
      INTEGER      LDumpIM(N_INTERNAL_MODEL_MAX)
! Global (ie. dump on disk) version of LDumpIM
      INTEGER      global_LDumpIM(N_INTERNAL_MODEL_MAX)
!Total data length for secondary flds for each submodel data partition
      INTEGER      LSECD(N_SUBMODEL_PARTITION_MAX)
!Total data length for secondary flds for each internal model
      INTEGER      LSecdIM(N_INTERNAL_MODEL_MAX)
!Total workspace length for each submodel data partition
      INTEGER      LWORK(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each submodel data partition
      INTEGER      NHeadSub(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each internal model
      INTEGER      NHEAD(N_INTERNAL_MODEL_MAX)
!Total length of extra space for each submod. data part.
      INTEGER      LEXTRA(N_SUBMODEL_PARTITION_MAX)
!Data length for dual-time level ocean fields
      INTEGER      LPRIM_O2
      INTEGER      ITEM_MAX_REQ
      INTEGER      ITEM_MAX_ALL

      INTEGER      NRECS_S
      INTEGER      NTIMES_S
      INTEGER      NSERBLK_S
      INTEGER      NSERREC_S
      INTEGER      NLEVL_S
      INTEGER      NMAXLEV_S
      INTEGER      NPSLISTS_S
      INTEGER      NMAXPSL_S
      INTEGER      NHEAD_FILE(OUTFILE_S:OUTFILE_E)
      LOGICAL      LSTUSER

      COMMON/STRET/                                                     &
     &  LPRIM,LDUMP,LSECD,LWORK,NHEAD,LEXTRA,LPRIM_O2,LPrimIM,LDumpIM,  &
     &  LSecdIM,NHeadSub,ITEM_MAX_REQ,ITEM_MAX_ALL,NSERBLK_S,NSERREC_S, &
     &  NLEVL_S,NMAXLEV_S,NPSLISTS_S,NMAXPSL_S,LSTUSER,NRECS_S,NTIMES_S,&
     &  NHEAD_FILE,                                                     &
     &  global_LPRIM,global_LPrimIM,global_LDUMP,global_LDumpIM
      CHARACTER(LEN=1)  H_ATMOS
      CHARACTER(LEN=1)  H_FLOOR
      CHARACTER(LEN=1)  H_STRAT
      CHARACTER(LEN=1)  H_GLOBAL(N_INTERNAL_MODEL_MAX         )
      INTEGER      H_VERS  (N_INTERNAL_MODEL_MAX,0:NSECTP)

      COMMON/CHOICE/ H_ATMOS,H_GLOBAL,H_FLOOR,H_STRAT

      COMMON/HVERS/ H_VERS

! These are set in SETMODL:
      INTEGER MEAN_NUMBER(N_INTERNAL_MODEL_MAX)
      COMMON/MODLMEAN/ MEAN_NUMBER


! Variables read in by namelist and used in SETMODL
      INTEGER      OCAAA 
      REAL         EWSPACEA,NSSPACEA
      REAL         FRSTLATA,FRSTLONA

      LOGICAL      ZonAvOzone
      LOGICAL      ZonAvTppsOzone
      REAL         LATS
      REAL         LONS
      INTEGER      LWBND
      INTEGER      OCALB
      REAL         POLELATA
      REAL         POLELONA
      INTEGER      SWBND
      INTEGER      TCA(A_MAX_TRVARS)
      INTEGER      TCA_LBC(A_MAX_TRVARS)  ! =1 if tracer in lbc file 
      INTEGER      TC_UKCA(A_MAX_UKCAVARS)
      INTEGER      TC_LBC_UKCA(A_MAX_UKCAVARS) ! =1 if tr in lbc file 
      INTEGER      StLevGWdrag
      INTEGER      BotVDiffLev
      INTEGER      TopVDiffLev


      COMMON/STSHCOMM/                                                  &
     &  RUN_TARGET_END,                                                 &
     &  OCAAA,EWSPACEA,POLELATA,FRSTLATA,LATS,                          &
     &  NSSPACEA,POLELONA,FRSTLONA,LONS,                                &
     &  SWBND,LWBND,                                                    &
     &  ZonAvOzone ,ZonAvTppsOzone, AOBINC,  AOBGRP,                    &
     &  StLevGWdrag, BotVDiffLev,TopVDiffLev,                           &
     &  OCALB,TCA,TCA_LBC,TC_UKCA,TC_LBC_UKCA


      CHARACTER(LEN=1) :: LFLOOR
      CHARACTER(LEN=1) :: OROGR
      CHARACTER(LEN=1) :: SWMCR
      CHARACTER(LEN=1) :: MESO

      COMMON/STSHCHAR/                                                  &
     &     LFLOOR,                                                      &
     &  OROGR,   SWMCR, MESO

      NAMELIST/STSHCOMP/                                                &
        RUN_TARGET_END,                                                 &
        OCAAA       ,EWSPACEA    ,POLELATA ,FRSTLATA  ,LATS   ,         &
                     NSSPACEA    ,POLELONA ,FRSTLONA  ,LONS   ,         &
        SWBND       ,LWBND                            ,OROGR  ,         &
        ZonAvOzone  ,SWMCR       ,MESO     ,                            &
        OCALB       ,LFLOOR      ,AOBINC   ,TCA,                        &
        TCA_LBC     ,TC_UKCA     ,TC_LBC_UKCA   ,AOBGRP          
  
! MODEL end

      INTEGER, INTENT(IN) :: row_length     ! length of row
      INTEGER, INTENT(IN) :: rows           ! number of rows
      INTEGER, INTENT(IN) :: n_rows         ! number of rows on v grid
      INTEGER, INTENT(IN) :: model_levels   ! number of levels
      INTEGER, INTENT(IN) :: bl_levels      ! number of levels in BL
      INTEGER, INTENT(IN) :: tr_levels      ! number of tracer levels
      INTEGER, INTENT(IN) :: wet_levels     ! number of wet levels
      INTEGER, INTENT(IN) :: land_pts       ! no of land points
      INTEGER, INTENT(IN) :: sm_levels      ! no of soil moisture levels
      INTEGER, INTENT(IN) :: ntiles         ! no of land tile types
      INTEGER, INTENT(IN) :: tr_ukca        ! no of activated tracers

      INTEGER :: I,J,idiag                  ! counters
      INTEGER :: errcode                    ! Variable passed to ereport

      CHARACTER (LEN=70) :: cmessage        ! Error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_SETD1DEFS',zhook_in,zhook_handle)

!     Set max fluxdiags to zero initially
      nmax_strat_fluxdiags = 0
      nmax_mode_diags      = 0
      n_nonchem_tracers    = 0
      n_chem_diags         = 13   ! constant for all N-R schemes

! as lbc_spec and lbc_mmr are passed to emission_ctl they must be allocated
     ALLOCATE(lbc_mmr(n_boundary_vals))
     ALLOCATE(lbc_spec(n_boundary_vals))

     lbc_mmr(:) = rmdi
!  Species with potential lower boundary conditions (same for all flavours)
     lbc_spec =                                                        &
          (/'N2O       ','CF2Cl2    ',                                 &
            'CFCl3     ','MeBr      ',                                 &
            'H2        ','CH4       ',                                 &
            'COS       '/)

      IF (L_UKCA_trop) THEN

! Standard tropospheric chemistry for B-E solver
! ==============================================
        n_chem_emissions = 8
        n_3d_emissions = 1       ! aircraft NOX
        n_aero_tracers = 0
        ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
        n_chem_tracers =  26
        n_chem_diags   =  22
        nr_therm       = 102        ! thermal reactions
        nr_phot        = 27         ! photolytic ---"---
        nmax_strat_fluxdiags = n_chem_tracers
        em_chem_spec =                                                 &
        (/'NO        ','CH4       ','CO        ','HCHO      ',         &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',         &
          'NO_aircrft'/)
      ELSEIF (L_ukca_std_trop) THEN

! Std tropospheric chemistry for N-R solver
! =========================================
        n_chem_emissions = 8
        n_3d_emissions = 1       ! aircraft NOX
        n_aero_tracers = 0
        ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
        n_chem_tracers =  39        ! O1D and O3P in SS
        nr_therm       = 102        ! thermal reactions
        nr_phot        = 27         ! photolytic ---"---
        nmax_strat_fluxdiags = n_chem_tracers
        em_chem_spec =                                                 &
        (/'NO        ','CH4       ','CO        ','HCHO      ',         &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',         &
          'NO_aircrft'/)
      ELSEIF (L_ukca_tropisop .AND. L_ukca_achem) THEN

! Std tropospheric chemistry + MIM with aerosol scheme (N-R)
! ==========================================================
        n_chem_emissions = 19    ! 2D emission fields
        n_3d_emissions = 4       ! SO2_nat, BC & OC biomass, aircraft NOX
        n_aero_tracers = 9       ! DMS, SO2... aerosol precursor species
        n_chem_tracers = 51
        ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
        nr_therm       = 113
        nr_phot        = 37
        nmax_strat_fluxdiags = n_chem_tracers
! Table refers to emissions, more species are emitted using surrogates
        em_chem_spec =                                                  &
        (/'NO        ','CH4       ','CO        ','HCHO      ',          &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',          &
          'C5H8      ','BC_fossil ','BC_biofuel','OC_fossil ',          &
          'OC_biofuel','Monoterp  ','NVOC      ','SO2_low   ',          &
          'SO2_high  ','DMS       ','NH3       ','SO2_nat   ',          &
          'BC_biomass','OC_biomass','NO_aircrft'/)
      ELSEIF (L_ukca_tropisop .AND. .NOT. L_ukca_achem) THEN

! Std tropospheric chemistry with MIM isoprene scheme (N-R)
! =========================================================
        n_chem_emissions = 9
        n_3d_emissions = 1       ! aircraft NOX
        n_aero_tracers = 0
        ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
        n_chem_tracers = 49
        nr_therm       = 132
        nr_phot        = 35
        nmax_strat_fluxdiags = n_chem_tracers
        em_chem_spec =                                                 &
        (/'NO        ','CH4       ','CO        ','HCHO      ',         &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',         &
          'C5H8      ','NO_aircrft'/)
      ELSEIF (L_ukca_aerchem) THEN

! Std trop chem with SO2, DMS, NH3, and monoterpene (BE)
! ======================================================
        n_chem_emissions = 18       ! Surface/ high-level emissions
        n_3d_emissions = 4          ! SO2_nat, aircraft NOX, OC & BC Biomass
        n_chem_tracers = 26         ! advected chemical tracers
        n_chem_diags   = 22         ! Number of non-advected prognostics
        n_aero_tracers =  7         ! advected aerochem ---"---
        nr_therm       = 137        ! thermal reactions
        nr_phot        = 27         ! photolytic ---"---
        nmax_strat_fluxdiags = n_chem_tracers
        ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
        em_chem_spec =                                                  &
        (/'NO        ','CH4       ','CO        ','HCHO      ',          &
          'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',          &
          'C5H8      ','BC_fossil ','BC_biofuel','OC_fossil ',          &
          'Monoterp  ','MeOH      ','SO2_low   ','SO2_high  ',          &
          'DMS       ','NH3       ','SO2_nat   ','BC_biomass',          &
          'OC_biomass','NO_aircrft'/)
      ELSEIF (L_UKCA_RAQ) THEN

! Regional air quality chemistry (RAQ), based on STOCHEM
! ========================================================
        n_chem_emissions  = 16
        n_3d_emissions    = 1       ! aircraft NOx
        n_chem_tracers    = 40      ! advected chemical tracers
        n_chem_diags      = 22      ! diagnosed species and other prognostics
        n_aero_tracers    = 0       ! advected aerochem tracers
        nr_therm          = 192     ! thermal reactions
        nr_phot           = 23      ! photolytic reacs 
        nmax_strat_fluxdiags = n_chem_tracers
        ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
        em_chem_spec =                                                  &
          (/'NO        ','CH4       ','CO        ','HCHO      ',        &
            'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ',        &
            'C5H8      ','C4H10     ','C2H4      ','C3H6      ',        &
            'TOLUENE   ','oXYLENE   ','CH3OH     ','H2        ',        &
            'NO_aircrft' /)     
      ELSEIF (L_UKCA_strat .OR. L_ukca_strattrop .OR. L_UKCA_stratcfc) &
        THEN

! Stratospheric chemistry
! =======================
          n_nonchem_tracers  = 2            ! Age of air and passive O3
     
          IF (L_ukca_strat) THEN
             IF (.NOT. L_ukca_achem) THEN ! NOT using aerosol chemistry
                ! emissions:  
                n_chem_emissions = 4
                n_3d_emissions = 1       ! aircraft NOX
                n_aero_tracers =  0
                ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
                em_chem_spec =                                          &
                     (/'NO        ','CH4       ',                       &
                       'CO        ','HCHO      ',                       &
                       'NO_aircrft'/)
                ! tracers and reactions:
                n_chem_tracers  = 37   ! CCMVal !!No H2OS, but does have H2O
                nr_therm       = 135        
                nr_phot        = 34         
             ELSE ! USING AEROSOL CHEMISTRY
                ! emissions:  
                n_chem_emissions = 7       ! em_chem_spec below
                n_3d_emissions   = 2       ! volc SO2 & aircraft NOX
                ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
                em_chem_spec =                                          &
                     (/'NO        ','CH4       ',                       &
                       'CO        ','HCHO      ',                       &
                       'SO2_low   ','SO2_high  ',                       &
                       'DMS       ','SO2_nat   ',                       &
                       'NO_aircrft'/)
                ! tracers and reactions:
                n_chem_tracers  = 45 
                n_aero_tracers  = 8
                nr_therm        = 149
                nr_phot         = 38
             END IF
          ELSE IF (L_ukca_strattrop .AND. .NOT. L_ukca_achem) THEN
             n_chem_emissions = 9
             n_3d_emissions = 1       ! aircraft NOX
             n_aero_tracers =  0
             ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
             em_chem_spec =                                             &
                 (/'NO        ','CH4       ','CO        ','HCHO      ', &
                   'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ', &
                   'C5H8      ','NO_aircrft'/)
             n_chem_tracers = 71         ! No chem tracers
             nr_therm       = 220        ! thermal reactions
             nr_phot        = 55         ! photolytic (ATA)

          ELSE IF (L_ukca_strattrop .AND. L_ukca_achem) THEN
             n_chem_emissions = 21      ! em_chem_spec below
             n_3d_emissions   = 2       ! volc SO2 & aircraft NOX
             ALLOCATE(em_chem_spec(n_chem_emissions+n_3d_emissions))
             em_chem_spec =                                             &
                 (/'NO        ','CH4       ','CO        ','HCHO      ', &
                   'C2H6      ','C3H8      ','Me2CO     ','MeCHO     ', &
                   'C5H8      ','BC_fossil ','BC_biofuel','OC_fossil ', &
                   'OC_biofuel','Monoterp  ','NVOC      ','SO2_low   ', &
                   'SO2_high  ','DMS       ','NH3       ','SO2_nat   ', &
                   'BC_biomass','OC_biomass','NO_aircrft'/)
             n_aero_tracers = 12
             n_chem_tracers = 71         ! No chem tracers
             IF (L_ukca_trophet) THEN
               nr_therm     = 241        ! thermal reactions
             ELSE
               nr_therm     = 239        ! thermal reactions
             ENDIF
             nr_phot        = 59         ! photolytic (ATA)
    
          ELSE IF (L_ukca_stratcfc) THEN
            n_chem_tracers = 43
          END IF
          n_chem_diags   = 13         ! Always the same for all N-R schemes
          nr_therm       = 102        !
          nr_phot        = 27         !
          nmax_strat_fluxdiags = n_chem_tracers
      ELSE

! No chemistry
! ============
        n_chem_emissions = 0
        n_chem_tracers   = 0
        n_chem_diags     = 0
      ENDIF

      IF (L_ukca_ageair .AND. .NOT. L_ukca_chem) THEN
! Age of air only
        n_nonchem_tracers = 1
      ENDIF

      IF (L_UKCA_dust) THEN
        n_dust_emissions = 6
        n_dust_tracers   = 6
        ALLOCATE(em_dust_spec(n_dust_emissions))
        em_dust_spec(1:n_dust_emissions)=                              &
        (/'Dust_div_1','Dust_div_2','Dust_div_3','Dust_div_4',         &
          'Dust_div_5','Dust_div_6'/)
      ELSE
        n_dust_emissions = 0
        n_dust_tracers   = 0
      ENDIF

      IF (L_UKCA_MODE) THEN
        IF (.NOT. (L_ukca_aerchem .OR. L_ukca_achem)) THEN
          cmessage=' L_ukca_aerchem or L_ukca_tropisop is required'//   &
                   ' to run UKCA_MODE'
          errcode=1
          CALL EREPORT('UKCA_SETD1DEFS',errcode,cmessage)
        ENDIF
        n_MODE_emissions = 0          ! See aerosol chemistry section
        n_MODE_tracers   = 31         ! For SU, SS, BC, OC, DU, SO
        nmax_mode_diags  = 300
      ELSE
        n_MODE_emissions = 0
        n_MODE_tracers   = 0
      ENDIF

      n_use_tracers   = n_chem_tracers + n_mode_tracers +               &
                        n_aero_tracers + n_dust_tracers +               &
                        n_nonchem_tracers
      n_use_emissions = n_chem_emissions + n_mode_emissions +           &
                        n_3d_emissions   + n_dust_emissions

      n_in_progs     =  30     ! max no of progs reqd other than tracers/ems
      n_in_diags0    =   4     ! max no of diags (sect 0) reqd
      n_in_diags1    =   2     ! max no of diags (sect 1) reqd
      n_in_diags2    =   1     ! max no of diags (sect 1) reqd
      n_in_diags3    =  17     ! max no of diags (sect 3) reqd
      n_in_diags4    =   5     ! max no of diags (sect 4) reqd
      n_in_diags5    =   4     ! max no of diags (sect 5) reqd
      n_in_diags8    =   1     ! max no of diags (sect 8) reqd
      n_in_diags15   =   1     ! max no of diags (sect 15) reqd
      n_in_diags30   =   1     ! max no of diags (sect 30) reqd
      n_in_diags34   = n_chem_diags
!                              ! max no UKCA diags (sect 34) reqd
      n_in_diags38   = nmax_mode_diags + nmax_strat_fluxdiags
!                              ! max no UKCA diags (sect 38) reqd


      n_emiss_first  = 301     ! first stash no for UKCA emissions
      n_emiss_last   = 320     !  last stash no for UKCA emissions

!     Create species names array

      ALLOCATE(nm_spec(n_all_tracers))
      IF (L_UKCA_RAQ) THEN
          !This list of tracers is valid for the RAQ chemistry.
          !If MODE aerosols are used with it but their positions  
          !change in the array then the list needs to be updated.
          nm_spec(1:n_all_tracers) = (/                                     &
         'O3        ','NO        ','NO3       ','NO2       ','N2O5      ',  &
         'HO2NO2    ','HONO2     ','H2O2      ','CH4       ','CO        ',  &
         'HCHO      ','MeOOH     ','HONO      ','C2H6      ','ETOOH     ',  &
         'MeCHO     ','PAN       ','C3H8      ','N-PrOOH   ','I-PrOOH   ',  &
         'EtCHO     ','Me2CO     ','MeCOCH2OOH','PPAN      ','MeONO2    ',  &
         'O3S       ','C5H8      ','ISOOH     ','ISON      ','MACR      ',  &
         'MACROOH   ','MPAN      ','HACET     ','MGLY      ','NALD      ',  &
         'HCOOH     ','MeCO3H    ','MeCO2H    ','MVK       ','MVKOOH    ',  &
         'Cl        ','ClO       ','Cl2O2     ','OClO      ','Br        ',  &
         'BrO       ','BrCl      ','BrONO2    ','N2O       ','HCl       ',  &
         'HOCl      ','HBr       ','HOBr      ','ClONO2    ','CFCl3     ',  &
         'CF2Cl2    ','MeBr      ','N         ','O(3P)     ','ORGNIT    ',  &
         'MeCl      ','CF2ClBr   ','CCl4      ','CF2ClCFCl2','CHF2Cl    ',  &
         'MeCCl3    ','CF3Br     ','H2OS      ','CH3OH     ','H2        ',  &
         'SO2       ','H2SO4     ','DMS       ','MSA       ','DMSO      ',  &
         'NH3       ','CS2       ','COS       ','H2S       ','H         ',  &
         'OH        ','HO2       ','MeOO      ','EtOO      ','MeCO3     ',  &
         'n-PrOO    ','i-PrOO    ','EtCO3     ','MeCOCH2OO ','RNC2H4    ',  &
         'RNC3H6    ','C2H4      ','C3H6      ','C4H10     ','C4H9OOH   ',  &
         'MEK       ','TOLUENE   ','MEMALD    ','GLYOXAL   ','oXYLENE   ',  &
         'ND_Nuc_SOL','Nuc_SOL_SU','ND_Ait_SOL','Ait_SOL_SU','Ait_SOL_BC',  &
         'Ait_SOL_OC','ND_Acc_SOL','Acc_SOL_SU','Acc_SOL_BC','Acc_SOL_OC',  &
         'Acc_SOL_SS','Acc_SOL_DU','ND_Cor_SOL','Cor_SOL_SU','Cor_SOL_BC',  &
         'Cor_SOL_OC','Cor_SOL_SS','Cor_SOL_DU','ND_Ait_INS','Ait_INS_BC',  &
         'Ait_INS_OC','ND_Acc_INS','Acc_INS_DU','ND_Cor_INS','Cor_INS_Du',  &
         'Nuc_SOL_OC','Ait_SOL_SS','Nuc_SOL_OZ','Ait_SOL_OZ','Acc_SOL_OZ',  &
         'Cor_SOL_OZ','Nuc_SOL_NH','Ait_SOL_NH','Acc_SOL_NH','Cor_SOL_NH',  &
         'Nuc_SOL_NT','Ait_SOL_NT','Acc_SOL_NT','Cor_SOL_NT','XXX       ',  &
         'Dust_Div_1','Dust_Div_2','Dust_Div_3','Dust_Div_4','Dust_Div_5',  &
         'Dust_Div_6','Rn-222    ','Pb-210    ','XXX       ','XXX       '   &
          /)
      ELSE
! Tracers 98,99 & 100 are for lumped Nitrogen, Br and Cl for stratospheric chemistry,
!  but can only be renamed in STASHmaster file not in advt or nm_spec.
          nm_spec(1:n_all_tracers) = (/                                     &
         'O3        ','NO        ','NO3       ','NO2       ','N2O5      ',  &
         'HO2NO2    ','HONO2     ','H2O2      ','CH4       ','CO        ',  & !10
         'HCHO      ','MeOOH     ','HONO      ','C2H6      ','EtOOH     ',  &
         'MeCHO     ','PAN       ','C3H8      ','n-PrOOH   ','i-PrOOH   ',  & !20
         'EtCHO     ','Me2CO     ','MeCOCH2OOH','PPAN      ','MeONO2    ',  &
         'O3_S      ','C5H8      ','ISOOH     ','ISON      ','MACR      ',  & !30
         'MACROOH   ','MPAN      ','HACET     ','MGLY      ','NALD      ',  &
         'HCOOH     ','MeCO3H    ','MeCO2H    ','H2O       ','ISO2      ',  & !40
         'Cl        ','ClO       ','Cl2O2     ','OClO      ','Br        ',  &
         'BrO       ','BrCl      ','BrONO2    ','N2O       ','HCl       ',  & !50
         'HOCl      ','HBr       ','HOBr      ','ClONO2    ','CFCl3     ',  &
         'CF2Cl2    ','MeBr      ','N         ','O(3P)     ','MACRO2    ',  & !60
         'MeCl      ','CF2ClBr   ','CCl4      ','CF2ClCFCl2','CHF2Cl    ',  &
         'MeCCl3    ','CF3Br     ','H2OS      ','CH2Br2    ','H2        ',  & !70
         'DMS       ','SO2       ','H2SO4     ','MSA       ','DMSO      ',  &
         'NH3       ','CS2       ','COS       ','H2S       ','H         ',  & !80
         'OH        ','HO2       ','MeOO      ','EtOO      ','MeCO3     ',  &
         'n-PrOO    ','i-PrOO    ','EtCO3     ','MeCOCH2OO ','MeOH      ',  & !90
         'Monoterp  ','Sec_Org   ','SESQUITERP','SO3       ','AROM      ',  &
         'O(3P)_S   ','O(1D)_S   ','NO2       ','BrO       ','HCl       ',  & !100
         'ND_Nuc_SOL','Nuc_SOL_SU','ND_Ait_SOL','Ait_SOL_SU','Ait_SOL_BC',  &
         'Ait_SOL_OC','ND_Acc_SOL','Acc_SOL_SU','Acc_SOL_BC','Acc_SOL_OC',  & !110
         'Acc_SOL_SS','Acc_SOL_DU','ND_Cor_SOL','Cor_SOL_SU','Cor_SOL_BC',  &
         'Cor_SOL_OC','Cor_SOL_SS','Cor_SOL_DU','ND_Ait_INS','Ait_INS_BC',  & !120
         'Ait_INS_OC','ND_Acc_INS','Acc_INS_DU','ND_Cor_INS','Cor_INS_Du',  &
         'Nuc_SOL_OC','Ait_SOL_SS','Nuc_SOL_OZ','Ait_SOL_OZ','Acc_SOL_OZ',  & !130
         'Cor_SOL_OZ','Nuc_SOL_NH','Ait_SOL_NH','Acc_SOL_NH','Cor_SOL_NH',  &
         'Nuc_SOL_NT','Ait_SOL_NT','Acc_SOL_NT','Cor_SOL_NT','XXX       ',  & !140
         'Anth_Prec ','Bio_Prec  ','Anth_Cond ','Bio_Cond  ','XXX       ',  &
         'XXX       ','XXX       ','XXX       ','PASSIVE O3','AGE OF AIR'   & !150
          /)
      END IF
! Mode components: Su: sulphate, BC: black carbon, OC: organic carbon
!                  SS: sea-salt, Du: dust,         OZ: organic carbon 2
!                  NH: ammonium, NT: nitrate,      ND: number density

! If radiative feedback is specified, check that the tracer array addresses 
! correspond to the named tracer
      IF (L_ukca_rado3) THEN
        IF (nm_spec(i_ukca_grg_o3) /= 'O3        ' .AND.                &
            ANY(advt(:) /= 'O3        ')) THEN
          errcode = i_ukca_grg_o3
          cmessage = ' Tracer address for O3 radiation feedback '//     &
                     ' does not correspond with O3'
          write(6,*) cmessage,i_ukca_grg_o3,nm_spec(i_ukca_grg_o3),     &
                     advt(i_ukca_grg_o3)
          CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
        END IF
      END IF
      IF (L_ukca_radch4) THEN
        IF (nm_spec(i_ukca_grg_ch4) /= 'CH4       ' .AND.               &
            ANY(advt(:) /= 'CH4       ')) THEN
          errcode = i_ukca_grg_ch4
          cmessage = ' Tracer address for CH4 radiation feedback '//    &
                     ' does not correspond with CH4'
          write(6,*) cmessage,i_ukca_grg_ch4,nm_spec(i_ukca_grg_ch4),   &
                     advt(i_ukca_grg_ch4)
          CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
        END IF
      END IF
      IF (L_ukca_radn2o) THEN
        IF (nm_spec(i_ukca_grg_n2o) /= 'N2O       ' .AND.               &
            ANY(advt(:) /= 'N2O       ')) THEN
          errcode = i_ukca_grg_n2o
          cmessage = ' Tracer address for N2O radiation feedback '//    &
                     ' does not correspond with CH4'
          write(6,*) cmessage,i_ukca_grg_ch4,nm_spec(i_ukca_grg_ch4),   &
                     advt(i_ukca_grg_ch4)
          CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
        END IF
      END IF

!     Specify the section and item codes of prognostics and diagnostics
!     required from D1.  Set array dimensions, but ignore halos as these
!     are set in the call to UKCA_SET_ARRAY_BOUNDS from the halosize
!     array.  Set %prognostic and %required t/f.
!     Other components of UkcaD1Codes are read in from D1 address array
!     in UKCA_MAIN1.

      Nukca_D1items = n_use_tracers  + n_use_emissions +                &
                      n_in_progs     + n_in_diags0     +                &
                      n_in_diags1    + n_in_diags2     +                &
                      n_in_diags3    +                                  &
                      n_in_diags4    + n_in_diags5     +                &
                      n_in_diags8    + n_in_diags15    +                &
                      n_in_diags30   + n_in_diags34    +                &
                      n_in_diags38

      ALLOCATE(UkcaD1Codes(Nukca_D1items))

      UkcaD1Codes(:)%section=IMDI
      UkcaD1Codes(:)%item=IMDI
      UkcaD1Codes(:)%n_levels=IMDI
      UkcaD1Codes(:)%address=IMDI
      UkcaD1Codes(:)%length=IMDI
      UkcaD1Codes(:)%halo_type=IMDI
      UkcaD1Codes(:)%grid_type=IMDI
      UkcaD1Codes(:)%field_type=IMDI
      UkcaD1Codes(:)%len_dim1=IMDI
      UkcaD1Codes(:)%len_dim2=IMDI
      UkcaD1Codes(:)%len_dim3=IMDI
      UkcaD1Codes(:)%prognostic=.TRUE.
      UkcaD1Codes(:)%required=.FALSE.

! If using Newton-Raphson solver, re-order tracers so that they 
!  exist in the same order as in ASAD

! Only retain names for active tracers and set required logical

      IF (printstatus >= Prstatus_oper)                                 &
        WRITE(6,*) ' UKCA: The following tracers were selected:'

      IF (.NOT.((L_UKCA_trop .OR. L_UKCA_raq) .OR. L_ukca_aerchem)) THEN
         DO i=1,jpctr
            SELECT CASE (advt(i))
               ! Link UM H2O tracer into chemistry scheme if selected in ASAD.
            CASE ('H2O       ')
               UkcaD1Codes(i)%section    = 0
               UkcaD1Codes(i)%item       = 10
               UkcaD1Codes(i)%len_dim1   = row_length
               UkcaD1Codes(i)%len_dim2   = rows
               UkcaD1Codes(i)%len_dim3   = wet_levels
               UkcaD1Codes(i)%prognostic = .TRUE.
               UkcaD1Codes(i)%required   = .TRUE. 
               UkcaD1Codes(J)%name       = advt(i)
               IF (printstatus >= prstatus_oper) WRITE(6,*) advt(i)
           CASE DEFAULT    
               ! Link tracers into UKCA in order of selection in ASAD.
       loop_j: DO j=1,100
                  IF (nm_spec(j) == advt(i)) THEN
!                    ! Check whether tracer is switched on in UMUI, but exclude these 3
!                    !  cases as there are duplicates in nm_spec to handle lumped cases.
                     IF ((tc_ukca(j) /= 1) .AND.                        &
                        (nm_spec(j) == 'NO2       ' .OR.                &
                         nm_spec(j) == 'BrO       ' .OR.                &
                         nm_spec(j) == 'HCl       ')) CYCLE loop_j 
                     IF (tc_ukca(j) /= 1) THEN
                        cmessage =                                      &
                 'Tracer '//advt(i)//' not selected in UMUI.'
                        errcode = 1
                        CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
                     ELSE
                        ! Check whether tracer is found only once here.
                        IF (UkcaD1Codes(i)%section /= IMDI) THEN
                           cmessage = &
                 'Tracer '//advt(i)//' multiply defined in SETD1DEFS.'
                           errcode = 1
                           CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
                        END IF
                        UkcaD1Codes(i)%section    = UKCA_sect
                        UkcaD1Codes(i)%item       = j
                        UkcaD1Codes(i)%len_dim1   = row_length
                        UkcaD1Codes(i)%len_dim2   = rows
                        UkcaD1Codes(i)%len_dim3   = tr_levels
                        UkcaD1Codes(i)%prognostic = .TRUE.
                        UkcaD1Codes(i)%required   = .TRUE.
                        UkcaD1Codes(J)%name       = advt(i)
                        IF (printstatus >= prstatus_oper) WRITE(6,*) advt(i)
                     END IF
                  END IF
               END DO loop_j

! Check whether tracer was found here.
               IF (UkcaD1Codes(i)%section==IMDI) THEN
                  cmessage = 'Tracer '//advt(i)//' unknown.'
                  errcode = 1
                  CALL ereport('UKCA_SETD1DEFS',errcode,cmessage)
               END IF
            END SELECT
         END DO

         ! Link remaining non-gasphase chemistry tracers
         J=jpctr
         DO I=101,n_all_tracers
            IF (TC_UKCA(I) /= 1) THEN
               nm_spec(I)='XXX       '
            ELSE
               J = J+1
               UkcaD1Codes(J)%section    = UKCA_sect
               UkcaD1Codes(J)%item       = I
               UkcaD1Codes(J)%len_dim1   = row_length
               UkcaD1Codes(J)%len_dim2   = rows
               UkcaD1Codes(J)%len_dim3   = tr_levels
               UkcaD1Codes(J)%prognostic = .TRUE.
               UkcaD1Codes(J)%required   = .TRUE.
               UkcaD1Codes(J)%name       = nm_spec(i)
               IF (printstatus >= prstatus_oper) WRITE(6,*) nm_spec(I)
            ENDIF
         ENDDO

      ELSE

         J=0
         DO I=1,n_all_tracers
            IF (TC_UKCA(I) /= 1) THEN
               nm_spec(I)='XXX       '
            ELSE
               J = J+1
               UkcaD1Codes(J)%section    = UKCA_sect
               UkcaD1Codes(J)%item       = I
               UkcaD1Codes(J)%len_dim1   = row_length
               UkcaD1Codes(J)%len_dim2   = rows
               UkcaD1Codes(J)%len_dim3   = tr_levels
               UkcaD1Codes(J)%prognostic = .TRUE.
               UkcaD1Codes(J)%required   = .TRUE.
               UkcaD1Codes(J)%name       = nm_spec(i)
               IF (printstatus >= prstatus_oper) WRITE(6,*) nm_spec(I)
            ENDIF
         ENDDO

         IF (J /= n_use_tracers) THEN
            cmessage='UKCA: wrong number of tracers'
            WRITE(6,*) 'Expected: ',n_use_tracers,' Encountered: ',J
            WRITE(6,*) 'tr_ukca: ',tr_ukca
            errcode = 1
            CALL EREPORT('UKCA_SETD1DEFS',errcode,cmessage)
         ENDIF

      END IF


      IF (J /= n_use_tracers) THEN
        cmessage='UKCA: wrong number of tracers'
        errcode = 1 
        WRITE(6,*) 'Expected: ',n_use_tracers,' Encountered: ',J
        WRITE(6,*) 'n_chem: ',n_chem_tracers,' n_mode: ',n_mode_tracers
        WRITE(6,*) 'n_aero: ',n_aero_tracers,' n_dust: ',n_dust_tracers
        WRITE(6,*) 'n_nonchem: ',n_nonchem_tracers
        WRITE(6,*) 'tr_ukca: ',tr_ukca
        CALL EREPORT('UKCA_SETD1DEFS',errcode,cmessage)
      ENDIF

! Prognostics from section 0, emissions required set above in em_chem_spec
! Surface emissions use stash codes from 301
! SO2 emissions from UM aerosol scheme
! 3D aerosol biomass emissions use stash 322-323, aircraft NOx= 340

      J = n_use_tracers
      IF (n_chem_emissions+n_3d_emissions+n_mode_emissions > 0) THEN
        DO i=1,n_chem_emissions + n_3d_emissions
          UkcaD1Codes(J+i)%section    = 0
          UkcaD1Codes(J+i)%item       = n_emiss_first+i-1  ! trop chemistry
          UkcaD1Codes(J+i)%len_dim1   = row_length         ! uses stash codes
          UkcaD1Codes(J+i)%len_dim2   = rows               ! 301-309 for
          UkcaD1Codes(J+i)%required   = .TRUE.             ! surface emissions
          UkcaD1Codes(J+i)%prognostic = .TRUE.             ! from Section 0
! Special cases, emissions already available in UM
          IF (em_chem_spec(i)(1:7) == 'SO2_low') THEN
            UkcaD1Codes(J+i)%item     = 58
            IF (.NOT. L_SO2_SURFEM .AND. (L_ukca_aerchem .OR.           &
                                          L_ukca_achem)) THEN
              cmessage='SO2 surface emissions from UM are not flagged'
              errcode=58

              CALL EREPORT('UKCA_SETD1DEFS',errcode,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:7) == 'SO2_nat') THEN
            UkcaD1Codes(J+i)%item     = 121
            UkcaD1Codes(J+i)%len_dim3 = tr_levels
            IF (.NOT. L_SO2_NATEM .AND. (L_ukca_aerchem .OR.            &
                      L_ukca_achem)) THEN
              cmessage='SO2 natural emissions from UM are not flagged'
              errcode=121

              CALL EREPORT('UKCA_SETD1DEFS',errcode,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:8) == 'SO2_high') THEN
            UkcaD1Codes(J+i)%item     = 126
            IF (.NOT. L_SO2_HILEM .AND. (L_ukca_aerchem .OR.            &
                                         L_ukca_achem)) THEN
              cmessage='SO2 high-level emissions are not flagged'
              errcode = UkcaD1Codes(J+i)%item
              CALL EREPORT('UKCA_SETD1DEFS',errcode,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:3) == 'NH3') THEN
            UkcaD1Codes(J+i)%item     = 127
            IF (.NOT. L_NH3_EM .AND. (L_ukca_aerchem .OR.               &
                                      L_ukca_achem)) THEN
              cmessage='NH3 surface emissions from UM are not flagged'
              errcode = UkcaD1Codes(J+i)%item
              CALL EREPORT('UKCA_SETD1DEFS',errcode,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i) == 'BC_fossil ') THEN
            UkcaD1Codes(J+i)%item = 310
          ELSEIF (em_chem_spec(i) == 'BC_biofuel') THEN
            UkcaD1Codes(J+i)%item = 311
          ELSEIF (em_chem_spec(i) == 'OC_fossil ') THEN
            UkcaD1Codes(J+i)%item = 312
          ELSEIF (em_chem_spec(i) == 'OC_biofuel') THEN
            UkcaD1Codes(J+i)%item = 313
          ELSEIF (em_chem_spec(i) == 'Monoterp  ') THEN
            UkcaD1Codes(J+i)%item = 314
          ELSEIF (em_chem_spec(i) == 'NVOC      ') THEN
            UkcaD1Codes(J+i)%item = 315
          ELSEIF (em_chem_spec(i) == 'BC_biomass') THEN
            UkcaD1Codes(J+i)%item = 322
            UkcaD1Codes(J+i)%len_dim3 = tr_levels
          ELSEIF (em_chem_spec(i) == 'OC_biomass') THEN
            UkcaD1Codes(J+i)%item = 323
            UkcaD1Codes(J+i)%len_dim3 = tr_levels
          ELSEIF (em_chem_spec(i) == 'SO2_biomas') THEN
            UkcaD1Codes(J+i)%item = 324
            UkcaD1Codes(J+i)%len_dim3 = tr_levels
          ELSEIF (em_chem_spec(i)(1:3) == 'DMS') THEN
            UkcaD1Codes(J+i)%section  = 17
            UkcaD1Codes(J+i)%item     = 205
            UkcaD1Codes(J+i)%prognostic = .FALSE.
            IF (.NOT. L_DMS_EM .AND. (L_ukca_aerchem .OR.               &
                      L_ukca_achem)) THEN
              cmessage='DMS surface emissions from UM are not flagged'
              errcode = UkcaD1Codes(J+i)%section*1000 +                   &
                       UkcaD1Codes(J+i)%item
              CALL EREPORT('UKCA_SETD1DEFS',errcode,cmessage)
            ENDIF
          ELSEIF (em_chem_spec(i)(1:7) == 'NO_airc') THEN
            UkcaD1Codes(J+i)%item     = 340
            UkcaD1Codes(J+i)%len_dim3 = tr_levels
          ENDIF
        ENDDO
      ENDIF

      J = J + n_chem_emissions + n_3d_emissions
      IF (L_UKCA_dust) THEN
        DO i=1,n_dust_emissions
          UkcaD1Codes(J+i)%section    = 0
          UkcaD1Codes(J+i)%item       = 311+i-1          ! Dust experiment
          UkcaD1Codes(J+i)%len_dim1   = row_length       ! uses stash code
          UkcaD1Codes(J+i)%len_dim2   = rows             ! 311-316 in Section 0
          UkcaD1Codes(J+i)%required   = .TRUE.           ! surface emissions
          UkcaD1Codes(J+i)%prognostic = .TRUE.           ! of dust
        ENDDO
      ENDIF

! Prognostic fields
      J = J + n_dust_emissions

      UkcaD1Codes(J+1:J+n_in_progs)%section    = 0
      UkcaD1Codes(J+1:J+n_in_progs)%prognostic = .TRUE.
      UkcaD1Codes(J+1:J+n_in_progs)%required   = .TRUE.
      UkcaD1Codes(J+1)%item=4              ! Potential Temperature
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+1)%len_dim3=model_levels
      UkcaD1Codes(J+2)%item=9              ! Soil Moisture
      UkcaD1Codes(J+2)%len_dim1=land_pts
      UkcaD1Codes(J+2)%len_dim2=sm_levels
      UkcaD1Codes(J+3)%item=10             ! Q
      UkcaD1Codes(J+3)%len_dim1=row_length
      UkcaD1Codes(J+3)%len_dim2=rows
      UkcaD1Codes(J+3)%len_dim3=wet_levels
      UkcaD1Codes(J+4)%item=12             ! QCF
      UkcaD1Codes(J+4)%len_dim1=row_length
      UkcaD1Codes(J+4)%len_dim2=rows
      UkcaD1Codes(J+4)%len_dim3=wet_levels
      UkcaD1Codes(J+5)%item=16             ! Conv cloud liquid water path
      UkcaD1Codes(J+5)%len_dim1=row_length
      UkcaD1Codes(J+5)%len_dim2=rows
      UkcaD1Codes(J+6)%item=24             ! Surface temperature
      UkcaD1Codes(J+6)%len_dim1=row_length
      UkcaD1Codes(J+6)%len_dim2=rows
      UkcaD1Codes(J+7)%item=25             ! Boundary layer height
      UkcaD1Codes(J+7)%len_dim1=row_length
      UkcaD1Codes(J+7)%len_dim2=rows
      UkcaD1Codes(J+8)%item=26             ! Roughness length
      UkcaD1Codes(J+8)%len_dim1=row_length
      UkcaD1Codes(J+8)%len_dim2=rows
      UkcaD1Codes(J+9)%item=30            ! Land sea mask
      UkcaD1Codes(J+9)%len_dim1=row_length
      UkcaD1Codes(J+9)%len_dim2=rows
      UkcaD1Codes(J+10)%item=31            ! Sea ice fraction
      UkcaD1Codes(J+10)%len_dim1=row_length
      UkcaD1Codes(J+10)%len_dim2=rows
      UkcaD1Codes(J+11)%item=60            ! Climatological ozone
      IF (ZONAVOZONE) THEN 
        UkcaD1Codes(J+11)%len_dim1=1 
      ELSE 
        UkcaD1Codes(J+11)%len_dim1=row_length 
      ENDIF
      UkcaD1Codes(J+11)%len_dim2=rows
      UkcaD1Codes(J+11)%len_dim3=model_levels
      UkcaD1Codes(J+12)%item=103           ! SO4 Aitken Mode
      UkcaD1Codes(J+12)%len_dim1=row_length
      UkcaD1Codes(J+12)%len_dim2=rows
      UkcaD1Codes(J+12)%len_dim3=tr_levels
      IF (.NOT. L_SULPC_SO2) UkcaD1Codes(J+12)%required=.FALSE.
      UkcaD1Codes(J+13)%item=104           ! SO4 accumulation Mode
      UkcaD1Codes(J+13)%len_dim1=row_length
      UkcaD1Codes(J+13)%len_dim2=rows
      UkcaD1Codes(J+13)%len_dim3=tr_levels
      IF (.NOT. L_SULPC_SO2) UkcaD1Codes(J+13)%required=.FALSE.
      UkcaD1Codes(J+14)%item=150           ! W component of wind
      UkcaD1Codes(J+14)%len_dim1=row_length
      UkcaD1Codes(J+14)%len_dim2=rows
      UkcaD1Codes(J+14)%len_dim3=model_levels+1
      IF (.NOT. L_ukca_arg_act) UkcaD1Codes(J+14)%required=.FALSE.
      UkcaD1Codes(J+15)%item=211           ! Conv cloud amount
      UkcaD1Codes(J+15)%len_dim1=row_length
      UkcaD1Codes(J+15)%len_dim2=rows
      UkcaD1Codes(J+15)%len_dim3=wet_levels
      UkcaD1Codes(J+16)%item=216           ! Fraction of surface types
      UkcaD1Codes(J+16)%len_dim1=land_pts
      UkcaD1Codes(J+16)%len_dim2=ntype
      UkcaD1Codes(J+17)%item=217           ! LAI of PFTs
      UkcaD1Codes(J+17)%len_dim1=land_pts
      UkcaD1Codes(J+17)%len_dim2=npft
      UkcaD1Codes(J+18)%item=218           ! Canopy heights of PFTs
      UkcaD1Codes(J+18)%len_dim1=land_pts
      UkcaD1Codes(J+18)%len_dim2=npft
      UkcaD1Codes(J+19)%item=229           ! Canopy water content on tiles
      UkcaD1Codes(J+19)%len_dim1=land_pts
      UkcaD1Codes(J+19)%len_dim2=ntiles
      UkcaD1Codes(J+20)%item=233           ! Surface temperature on tiles
      UkcaD1Codes(J+20)%len_dim1=land_pts
      UkcaD1Codes(J+20)%len_dim2=ntiles
      UkcaD1Codes(J+21)%item=234           ! Surface roughness lengths on t
      UkcaD1Codes(J+21)%len_dim1=land_pts
      UkcaD1Codes(J+21)%len_dim2=ntiles
      UkcaD1Codes(J+22)%item=240           ! Snow depth on tiles
      UkcaD1Codes(J+22)%len_dim1=land_pts
      UkcaD1Codes(J+22)%len_dim2=ntiles
      UkcaD1Codes(J+23)%item=253           ! Rho_r2
      UkcaD1Codes(J+23)%len_dim1=row_length
      UkcaD1Codes(J+23)%len_dim2=rows
      UkcaD1Codes(J+23)%len_dim3=model_levels
      UkcaD1Codes(J+24)%item=254           ! QCL
      UkcaD1Codes(J+24)%len_dim1=row_length
      UkcaD1Codes(J+24)%len_dim2=rows
      UkcaD1Codes(J+24)%len_dim3=wet_levels
      UkcaD1Codes(J+25)%item=255           ! Exner pressure on rho levels
      UkcaD1Codes(J+25)%len_dim1=row_length
      UkcaD1Codes(J+25)%len_dim2=rows
      UkcaD1Codes(J+25)%len_dim3=model_levels+1
      UkcaD1Codes(J+26)%item=265           ! Area cloud fraction in each la
      UkcaD1Codes(J+26)%len_dim1=row_length
      UkcaD1Codes(J+26)%len_dim2=rows
      UkcaD1Codes(J+26)%len_dim3=wet_levels
      UkcaD1Codes(J+27)%item=266           ! Bulk Cloud fraction
      UkcaD1Codes(J+27)%len_dim1=row_length
      UkcaD1Codes(J+27)%len_dim2=rows
      UkcaD1Codes(J+27)%len_dim3=wet_levels
      UkcaD1Codes(J+28)%item=267           ! Cloud Liquid fraction
      UkcaD1Codes(J+28)%len_dim1=row_length
      UkcaD1Codes(J+28)%len_dim2=rows
      UkcaD1Codes(J+28)%len_dim3=wet_levels
      UkcaD1Codes(J+29)%item=505           ! Land fraction
      UkcaD1Codes(J+29)%len_dim1=land_pts
      IF (.NOT. L_CTILE) UkcaD1Codes(J+29)%required = .FALSE.
      UkcaD1Codes(J+30)%item=510           ! Land albedo
      UkcaD1Codes(J+30)%len_dim1=row_length
      UkcaD1Codes(J+30)%len_dim2=rows
      IF (.NOT. L_SW_Radiate_prog .OR. .NOT. L_CTILE) THEN
! Required only on radiation TS, not available if coastal tiling off
        UkcaD1Codes(J+30)%required = .FALSE.
      ENDIF

! Diagnostics from section zero
      J = J + n_in_progs

      UkcaD1Codes(J+1:J+n_in_diags0)%section    = 0
      UkcaD1Codes(J+1:J+n_in_diags0)%prognostic = .FALSE.
      UkcaD1Codes(J+1:J+n_in_diags0)%required   = .TRUE.
      UkcaD1Codes(J+1)%item=406           ! Exner Press on theta levels
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+1)%len_dim3=model_levels
      UkcaD1Codes(J+2)%item=407           ! P on Rho Levels
      UkcaD1Codes(J+2)%len_dim1=row_length
      UkcaD1Codes(J+2)%len_dim2=rows
      UkcaD1Codes(J+2)%len_dim3=model_levels
      UkcaD1Codes(J+3)%item=408           ! P on Theta Levels
      UkcaD1Codes(J+3)%len_dim1=row_length
      UkcaD1Codes(J+3)%len_dim2=rows
      UkcaD1Codes(J+3)%len_dim3=model_levels
      UkcaD1Codes(J+4)%item=409           ! Surface Pressure
      UkcaD1Codes(J+4)%len_dim1=row_length
      UkcaD1Codes(J+4)%len_dim2=rows

! Diagnostic items from section 1 (SW radiation)
      J = J + n_in_diags0

      UkcaD1Codes(J+1:J+n_in_diags1)%section    = 1
      UkcaD1Codes(J+1:J+n_in_diags1)%prognostic = .FALSE.   ! always needed
      UkcaD1Codes(J+1:J+n_in_diags1)%required   = .TRUE.
      UkcaD1Codes(J+1)%item=201           ! Net downward surface SW flux
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+2)%item=235           ! Total downward surface SW flux
      UkcaD1Codes(J+2)%len_dim1=row_length
      UkcaD1Codes(J+2)%len_dim2=rows

! Diagnostic items from section 2 (LW radiation)
      J = J + n_in_diags1
      UkcaD1Codes(J+1:J+n_in_diags2)%section    = 2
      UkcaD1Codes(J+1:J+n_in_diags2)%prognostic = .FALSE.
      UkcaD1Codes(J+1)%item=284                     ! Sulphate optical depth
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+1)%len_dim3=sw_spectrum(1)%n_band  ! sw wave band pseudo levels
      write(6,*) 'n_band: ', sw_spectrum(1)%n_band
      IF ((i_ukca_photol == i_ukca_fastjx) .AND. L_SULPC_SO2 .AND.    &
        .NOT. L_ukca_mode) THEN
          UkcaD1Codes(J+1)%required=.TRUE.
      ELSE
        UkcaD1Codes(J+1)%required=.FALSE.
      END IF

! Diagnostic variables in section 3 (Boundary Layer)
      J = J + n_in_diags2

      UkcaD1Codes(J+1:J+n_in_diags3)%section    = 3
      UkcaD1Codes(J+1:J+n_in_diags3)%prognostic = .FALSE.
      UkcaD1Codes(J+1:J+n_in_diags3)%required   = .TRUE.
      UkcaD1Codes(J+1)%item=217           ! Surface heat flux
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+2)%item=462            ! Stomatal conductance
      UkcaD1Codes(J+2)%len_dim1=row_length
      UkcaD1Codes(J+2)%len_dim2=rows
      Ukcad1codes(J+2)%len_dim3=npft
      UkcaD1Codes(J+3)%item=465            ! Surface friction velocity
      UkcaD1Codes(J+3)%len_dim1=row_length
      UkcaD1Codes(J+3)%len_dim2=rows
      UkcaD1Codes(J+4)%item=60            ! rhokh_mix
      UkcaD1Codes(J+4)%len_dim1=row_length
      UkcaD1Codes(J+4)%len_dim2=rows
      Ukcad1codes(J+4)%len_dim3=bl_levels
      UkcaD1Codes(J+5)%item=64            ! dtrdz_charney_grid
      UkcaD1Codes(J+5)%len_dim1=row_length
      UkcaD1Codes(J+5)%len_dim2=rows
      Ukcad1codes(J+5)%len_dim3=bl_levels
      UkcaD1Codes(J+6)%item=65            ! kent
      UkcaD1Codes(J+6)%len_dim1=row_length
      UkcaD1Codes(J+6)%len_dim2=rows
      UkcaD1Codes(J+7)%item=66            ! we_lim
      UkcaD1Codes(J+7)%len_dim1=row_length
      UkcaD1Codes(J+7)%len_dim2=rows
      Ukcad1codes(J+7)%len_dim3=npft
      UkcaD1Codes(J+8)%item=67            ! t_frac
      UkcaD1Codes(J+8)%len_dim1=row_length
      UkcaD1Codes(J+8)%len_dim2=rows
      Ukcad1codes(J+8)%len_dim3=npft
      UkcaD1Codes(J+9)%item=68            ! zrzi
      UkcaD1Codes(J+9)%len_dim1=row_length
      UkcaD1Codes(J+9)%len_dim2=rows
      Ukcad1codes(J+9)%len_dim3=npft
      UkcaD1Codes(J+10)%item=69            ! kent_dsc
      UkcaD1Codes(J+10)%len_dim1=row_length
      UkcaD1Codes(J+10)%len_dim2=rows
      UkcaD1Codes(J+11)%item=70            ! we_lim_dsc
      UkcaD1Codes(J+11)%len_dim1=row_length
      UkcaD1Codes(J+11)%len_dim2=rows
      Ukcad1codes(J+11)%len_dim3=npft
      UkcaD1Codes(J+12)%item=71            ! t_frac_dsc
      UkcaD1Codes(J+12)%len_dim1=row_length
      UkcaD1Codes(J+12)%len_dim2=rows
      Ukcad1codes(J+12)%len_dim3=npft
      UkcaD1Codes(J+13)%item=72            ! zrzi_dsc
      UkcaD1Codes(J+13)%len_dim1=row_length
      UkcaD1Codes(J+13)%len_dim2=rows
      Ukcad1codes(J+13)%len_dim3=npft
      UkcaD1Codes(J+14)%item=73            ! zhsc
      UkcaD1Codes(J+14)%len_dim1=row_length
      UkcaD1Codes(J+14)%len_dim2=rows
      UkcaD1Codes(J+15)%item=25            ! ml_depth
      UkcaD1Codes(J+15)%len_dim1=row_length
      UkcaD1Codes(J+15)%len_dim2=rows
      UkcaD1Codes(J+16)%item=230           ! 10 m wind speed on C grid
      UkcaD1Codes(J+16)%len_dim1=row_length
      UkcaD1Codes(J+16)%len_dim2=rows
      IF (.NOT. L_UKCA_MODE) UkcaD1Codes(J+16)%required=.FALSE.
      UkcaD1Codes(J+17)%item=473           ! Turbulent Kinetic Energy
      UkcaD1Codes(J+17)%len_dim1=row_length
      UkcaD1Codes(J+17)%len_dim2=rows
      UkcaD1Codes(J+17)%len_dim3=bl_levels
      IF (.NOT. L_UKCA_ARG_ACT) UkcaD1Codes(J+18)%required=.FALSE.

! Diagnostic variables in section 4 (LS Precipitation)
! We currently do not use 4.206 and 4.227 so set required to .FALSE.
      J = J + n_in_diags3

      UkcaD1Codes(J+1:J+n_in_diags4)%section    = 4
      UkcaD1Codes(J+1:J+n_in_diags4)%prognostic = .FALSE.
      UkcaD1Codes(J+1:J+n_in_diags4)%required   = .TRUE.
      UkcaD1Codes(J+1)%item=205           ! Cloud Liquid Water after LS
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+1)%len_dim3=wet_levels
      UkcaD1Codes(J+2)%item=206           ! Cloud Ice Content after LS
      UkcaD1Codes(J+2)%len_dim1=row_length
      UkcaD1Codes(J+2)%len_dim2=rows
      UkcaD1Codes(J+2)%len_dim3=wet_levels
      UkcaD1Codes(J+2)%required=.FALSE.
      UkcaD1Codes(J+3)%item=222           ! Rainfall rate out of model levs
      UkcaD1Codes(J+3)%len_dim1=row_length
      UkcaD1Codes(J+3)%len_dim2=rows
      UkcaD1Codes(J+3)%len_dim3=wet_levels
      UkcaD1Codes(J+4)%item=223           ! Snowfall rate out of model levs
      UkcaD1Codes(J+4)%len_dim1=row_length
      UkcaD1Codes(J+4)%len_dim2=rows
      UkcaD1Codes(J+4)%len_dim3=wet_levels
      UkcaD1Codes(J+5)%item=227           ! Rain fraction (3C ppn only)
      UkcaD1Codes(J+5)%len_dim1=row_length
      UkcaD1Codes(J+5)%len_dim2=rows
      UkcaD1Codes(J+5)%len_dim3=wet_levels
      UkcaD1Codes(J+5)%required=.FALSE.

! Diagnostic variables in section 5 (Convection)
      J = J + n_in_diags4

      UkcaD1Codes(J+1:J+n_in_diags5)%section    = 5
      UkcaD1Codes(J+1:J+n_in_diags5)%prognostic = .FALSE.
      UkcaD1Codes(J+1:J+n_in_diags5)%required   = .TRUE.
      UkcaD1Codes(J+1)%item=227           ! 3D Convective rainfall rate
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows
      UkcaD1Codes(J+1)%len_dim3=model_levels
      UkcaD1Codes(J+2)%item=228           ! 3D Convective snowfall rate
      UkcaD1Codes(J+2)%len_dim1=row_length
      UkcaD1Codes(J+2)%len_dim2=rows
      UkcaD1Codes(J+2)%len_dim3=model_levels
! REMOVED FROM PROGNOSTICS (00014 & 00015) AND NOW TAKE DIAGNOSTICS 05218 & 05219
      UkcaD1Codes(J+3)%item=218             ! Conv cloud base level
      UkcaD1Codes(J+3)%len_dim1=row_length
      UkcaD1Codes(J+3)%len_dim2=rows
      UkcaD1Codes(J+4)%item=219             ! Conv cloud top level
      UkcaD1Codes(J+4)%len_dim1=row_length
      UkcaD1Codes(J+4)%len_dim2=rows


! Diagnostic variables in section 8 (Hydrology)
      J = J + n_in_diags5

      UkcaD1Codes(J+1)%section    = 8
      UkcaD1Codes(J+1)%prognostic = .FALSE.
      UkcaD1Codes(J+1)%required   = L_ukca_qch4inter
      UkcaD1Codes(J+1)%item       = 242        ! CH4 wetland flux
      UkcaD1Codes(J+1)%len_dim1   = row_length
      UkcaD1Codes(J+1)%len_dim2   = rows

! Diagnostic variables in section 15 (Processed Climate Diagnostics)
      J = J + n_in_diags8

      UkcaD1Codes(J+1)%section    = 15
      UkcaD1Codes(J+1)%prognostic = .FALSE.
      UkcaD1Codes(J+1)%required   = .TRUE.
      UkcaD1Codes(J+1)%item       = 218           ! PV on theta levels
      UkcaD1Codes(J+1)%len_dim1   = row_length
      UkcaD1Codes(J+1)%len_dim2   = rows
      UkcaD1Codes(J+1)%len_dim3   = model_levels

! Diagnostic variables in section 30 (Processed dynamics diagnostics)
      J = J + n_in_diags15

! Take tropopause height here. Needed only for volcanic SO2 emissions into
! the stratosphere. (Always required as in call to emission_ctl)
      UkcaD1Codes(J+1)%section = 30
      UkcaD1Codes(J+1)%prognostic = .FALSE.
      UkcaD1Codes(J+1)%required   = .TRUE.
      UkcaD1Codes(J+1)%item=453           ! Tropopause height
      UkcaD1Codes(J+1)%len_dim1=row_length
      UkcaD1Codes(J+1)%len_dim2=rows


! Diagnostic variables (UKCA), i.e. non-transported chemical species only
! These are held in start and dump files and may be used to initialise the solver
      J = J + n_in_diags30
      idiag_first = J+1
      idiag_last  = J+n_chem_diags

      DO I=1,n_chem_diags
        UkcaD1Codes(J+I)%section    = UKCA_sect
        UkcaD1Codes(J+I)%item       = item1_chem_diags+I-1
        UkcaD1Codes(J+I)%len_dim1   = row_length
        UkcaD1Codes(J+I)%len_dim2   = rows
        UkcaD1Codes(J+I)%len_dim3   = model_levels
        UkcaD1Codes(J+I)%required   = .TRUE.
        UkcaD1Codes(J+I)%prognostic = .TRUE.
      ENDDO

      IF (L_ukca_ageair .AND. .NOT. L_ukca_chem) THEN
        UkcaD1Codes(n_use_tracers+1:idiag_last)%required = .FALSE.
      ENDIF

! Diagnostic variables for stratospheric fluxes in section 38
      J = J + n_chem_diags
      istrat_first = J+1
      istrat_last  = J+nmax_strat_fluxdiags

      n_strat_fluxdiags = 0
      DO I=1,nmax_strat_fluxdiags
        DO idiag=1,ndiag     ! search for stash requests
          IF (modl_b(idiag) == SUBMODEL_FOR_SM(atmos_im) .AND.         &
              isec_b(idiag) == mode_diag_sect .AND.                    &
              item_b(idiag) == item1_stratflux+I-1) THEN
            n_strat_fluxdiags = n_strat_fluxdiags + 1
            UkcaD1Codes(J+I)%section    = UKCA_sect
            UkcaD1Codes(J+I)%item       = item_b(idiag)
            UkcaD1Codes(J+I)%len_dim1   = row_length
            UkcaD1Codes(J+I)%len_dim2   = rows
            UkcaD1Codes(J+I)%len_dim3   = model_levels
            UkcaD1Codes(J+I)%required   = .FALSE.
            UkcaD1Codes(J+I)%prognostic = .FALSE.
            IF (.NOT.(L_ukca_stratflux)) L_ukca_stratflux = .TRUE.
            EXIT
          ENDIF
        ENDDO
      ENDDO

! Diagnostic variables for UKCA-MODE in section 34 or 38
      J = J + nmax_strat_fluxdiags
      imode_first = J+1
      imode_last  = J+nmax_mode_diags

      IF (L_ukca_mode) THEN
        n_mode_diags = 0
        DO I=1,nmax_mode_diags
          DO idiag=1,ndiag     ! search for stash requests
            IF (modl_b(idiag) == SUBMODEL_FOR_SM(atmos_im) .AND.        &
              isec_b(idiag) == MODE_diag_sect .AND.                     &
              item_b(idiag) == item1_mode_diags+I-1) THEN
              n_mode_diags = n_mode_diags + 1
              UkcaD1Codes(J+I)%section    = MODE_diag_sect
              UkcaD1Codes(J+I)%item       = item_b(idiag)
              UkcaD1Codes(J+I)%len_dim1   = row_length
              UkcaD1Codes(J+I)%len_dim2   = rows
              UkcaD1Codes(J+I)%len_dim3   = model_levels
              UkcaD1Codes(J+I)%required   = .FALSE.
              UkcaD1Codes(J+I)%prognostic = .FALSE.
              IF (.NOT.(L_ukca_mode_diags)) L_ukca_mode_diags = .TRUE.
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      IF ( printstatus >= prstatus_diag ) THEN         
        WRITE(6,*) ' '
        WRITE(6,*) 'n_use_tracers: ', n_use_tracers
        WRITE(6,*) 'n_use_emissions: ', n_use_emissions
        WRITE(6,*) 'Total no of items required = ', Nukca_D1items
        WRITE(6,*) ' '
        WRITE(6,*) ' UKCA: UkcaD1Codes required items:'
        WRITE(6,*) 'section,item,prognostic,required,dim1,dim2,dim3'
        DO I=1,Nukca_D1items
          IF (UkcaD1Codes(I)%required)                                  &
          WRITE(6,*) I, UkcaD1Codes(I)%section,                         &
                        UkcaD1Codes(I)%item,                            &
                        UkcaD1Codes(I)%prognostic,                      &
                        UkcaD1Codes(I)%required,                        &
                        UkcaD1Codes(I)%len_dim1,                        &
                        UkcaD1Codes(I)%len_dim2,                        &
                        UkcaD1Codes(I)%len_dim3
        ENDDO
      ENDIF

      IF (lhook) CALL dr_hook('UKCA_SETD1DEFS',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE UKCA_SETD1DEFS
