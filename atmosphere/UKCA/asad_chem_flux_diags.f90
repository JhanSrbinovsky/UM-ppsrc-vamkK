! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to deal with the chemical flux diagnostics from ASAD
!
!  Part of the UKCA model. UKCA is a community model supported
!  by The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      MODULE ASAD_CHEM_FLUX_DIAGS

      USE UKCA_D1_DEFS,     ONLY : code, asad_diag_sect, item1_asad_diags
      USE ASAD_FLUX_DAT,    ONLY : asad_flux_defn,                      &
                                   asad_chemical_fluxes
      USE ASAD_MOD
      USE UKCA_TROPOPAUSE,  ONLY: L_troposphere
      USE UKCA_CSPECIES
      USE parkind1,         ONLY: jprb, jpim
      USE yomhook,          ONLY: lhook, dr_hook
      USE ereport_mod,      ONLY: ereport
      USE PrintStatus_mod
      USE UM_ParVars
      USE um_input_control_mod,  ONLY: h_sect
! version_mod items required by cstash.h and model.h
      USE version_mod,      ONLY: nproftp, nprofdp, nprofup,            &
                                  ndiagpm, ntimep, NTimSerP,            &
                                  nlevp, npslevp, npslistp,             &
                                  outfile_s, outfile_e, nsectp
      USE Submodel_Mod

      IMPLICIT NONE

      PRIVATE

! public types
      PUBLIC :: chemdiag
! public variables
      PUBLIC :: asad_chemdiags
      PUBLIC :: chemD1codes
      PUBLIC :: nmax_chemdiags
      PUBLIC :: n_chemdiags
      PUBLIC :: L_asad_use_chem_diags
      PUBLIC :: L_asad_use_air_ems
      PUBLIC :: L_asad_use_light_ems
      PUBLIC :: L_asad_use_light_diags_tot
      PUBLIC :: L_asad_use_light_diags_c2g
      PUBLIC :: L_asad_use_light_diags_c2c
      PUBLIC :: L_asad_use_light_diags_N
      PUBLIC :: L_asad_use_volc_ems
      PUBLIC :: L_asad_use_sulp_ems
      PUBLIC :: L_asad_use_surf_ems
      PUBLIC :: L_asad_use_flux_rxns
      PUBLIC :: L_asad_use_rxn_rates
      PUBLIC :: L_asad_use_wetdep
      PUBLIC :: L_asad_use_drydep
      PUBLIC :: L_asad_use_tendency      
      PUBLIC :: L_asad_use_mass_diagnostic
      PUBLIC :: L_asad_use_LiN_diagnostic
      PUBLIC :: L_asad_use_psc_diagnostic
      PUBLIC :: L_asad_use_trop_mask
      PUBLIC :: L_asad_use_STE
      PUBLIC :: L_asad_use_output_tracer
      PUBLIC :: aircraft_emissions
      PUBLIC :: lightning_emissions
      PUBLIC :: total_lightning_flashes
      PUBLIC :: lightning_flashes_cloud2ground
      PUBLIC :: lightning_flashes_cloud2cloud
      PUBLIC :: total_N_2D
      PUBLIC :: volcanic_emissions
      PUBLIC :: so2_emchar
      PUBLIC :: calculate_tendency
      PUBLIC :: calculate_STE
! public routines
      PUBLIC :: asad_chemical_diagnostics
      PUBLIC :: asad_setstash_chemdiag
      PUBLIC :: asad_init_chemdiag   
      PUBLIC :: asad_emissions_diagnostics
      PUBLIC :: asad_3D_emissions_diagnostics
      PUBLIC :: asad_tendency_STE
      PUBLIC :: asad_tropospheric_mask
      PUBLIC :: asad_mass_diagnostic
      PUBLIC :: asad_flux_put_stash
      PUBLIC :: asad_allocate_chemdiag
      PUBLIC :: asad_psc_diagnostic
      PUBLIC :: asad_output_tracer
      PUBLIC :: asad_lightning_diagnostics
      PUBLIC :: asad_LiN_diagnostic
      


      TYPE chemdiag
        INTEGER :: location          ! location of reaction in
                                     ! y/rk/prk/em_field/dpd/dpw array
        INTEGER :: stash_number      ! Stash number
        INTEGER :: find_rxn_loc      ! Reaction location
        INTEGER :: num_reactants     ! No of reactants
        INTEGER :: num_products      ! No of products
        LOGICAL :: tropospheric_mask ! are we only outputting the troposphere?
        REAL    :: molecular_mass    ! needed for emissions, tendency, STE
        REAL    :: c_vmr_to_mmr      ! conversion factor from vmr to mmr
        INTEGER :: num_levs          ! No of levels
        LOGICAL :: can_deallocate    ! can the throughput array be deallocated ?
        LOGICAL :: output_on_chem_tstep ! is the data in the throughput array
                                        ! only available on chemical timesteps?
        CHARACTER(LEN=3) :: diag_type ! Type: RXN,DDP,WDP,EMS,NET,STE,TPM
        CHARACTER(LEN=1) :: rxn_type  ! RXN only: B,J,H,T 
                                      ! (DDP==D,WDP==W,EMS==E),X
        CHARACTER(LEN=10) :: species  ! For DDP, WDP, EMS, NET only
        ! fix number of reactants at 2, but allow number of products to change
        CHARACTER(LEN=10), DIMENSION(1:2) :: reactants ! for use with RXN only
        ! for use with RXN only
        CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE :: products
        REAL, DIMENSION(:,:,:), ALLOCATABLE :: throughput 
      END TYPE chemdiag

      TYPE stashsumdiag
        INTEGER :: stash_value               ! section*1000 + item
        INTEGER :: stash_section             ! section
        INTEGER :: stash_item                ! item
        INTEGER :: chemd1codes_location      ! D1 location
        INTEGER :: number_of_fields          ! No of fields
        INTEGER :: len_dim1                  ! Array dimension
        INTEGER :: len_dim2                  ! Array dimension
        INTEGER :: len_dim3                  ! Array dimension
        LOGICAL :: output_on_chem_tstep      ! False for o/p at every timestep
        INTEGER, DIMENSION(:), ALLOCATABLE :: chemdiags_location
      END TYPE stashsumdiag

      INTEGER :: n_stashsumdiag              ! No of requests

      CHARACTER(LEN=10) :: blank = '          '  ! Represents undefined string

! Array to define diagnostic properties
      TYPE(chemdiag), DIMENSION(:), ALLOCATABLE, SAVE :: asad_chemdiags

! Array to define D1 lengths, address, etc.
      TYPE(code), DIMENSION(:), ALLOCATABLE, SAVE :: chemD1codes

      TYPE(stashsumdiag), DIMENSION(:), ALLOCATABLE, SAVE ::            &
                                                   stash_handling

      INTEGER, PARAMETER :: nmax_chemdiags = 205      ! For CJ
      INTEGER, SAVE :: n_chemdiags

      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: temp_chemdiag
      REAL, ALLOCATABLE, DIMENSION(:,:)   :: temp_chemdiag_2d

      INTEGER, SAVE :: N_chemdiags_first

      CHARACTER(LEN=3), SAVE :: end_label='EOF'

! This is set to .FALSE. initially, and changed in setstash
      LOGICAL, SAVE :: L_asad_use_chem_diags=.FALSE.

! Initially set to these - if not required then will be turned off.
      LOGICAL, SAVE :: L_asad_use_air_ems=.FALSE.   ! For aircraft emission
      LOGICAL, SAVE :: L_asad_use_light_ems=.FALSE. ! For lightning emission
      LOGICAL, SAVE :: L_asad_use_surf_ems=.FALSE.  ! For surface emission
      LOGICAL, SAVE :: L_asad_use_volc_ems=.FALSE.  ! For volcanic emission
      LOGICAL, SAVE :: L_asad_use_sulp_ems=.FALSE.  ! For 3D SO2 emission
      LOGICAL, SAVE :: L_asad_use_flux_rxns=.FALSE. ! For reaction flux
      LOGICAL, save :: L_asad_use_rxn_rates=.FALSE.
      LOGICAL, SAVE :: L_asad_use_wetdep=.FALSE.    ! For wet deposition
      LOGICAL, SAVE :: L_asad_use_drydep=.FALSE.    ! For dry deposition
      LOGICAL, SAVE :: L_asad_use_tendency=.FALSE.  ! For tendency
      LOGICAL, SAVE :: L_asad_use_STE=.FALSE.       ! For strat-trop exchange
      
      ! For air mass diagnostic
      LOGICAL, SAVE :: L_asad_use_mass_diagnostic=.FALSE.
      
      ! For psc diagnostic
      LOGICAL, SAVE :: L_asad_use_psc_diagnostic=.FALSE.
      
      ! To output tropospheric mask
      LOGICAL, SAVE :: L_asad_use_trop_mask_output=.FALSE.
      
      ! To output tracer
      LOGICAL, SAVE :: L_asad_use_output_tracer=.FALSE.
      
      ! For o/p in troposphere only
      LOGICAL, SAVE :: L_asad_use_trop_mask=.FALSE.
      
      LOGICAL, save :: L_asad_use_LiN_diagnostic=.FALSE.     
      LOGICAL, save :: L_asad_use_light_diags_tot=.FALSE.
      LOGICAL, save :: L_asad_use_light_diags_c2g=.FALSE.
      LOGICAL, save :: L_asad_use_light_diags_c2c=.FALSE.
      LOGICAL, save :: L_asad_use_light_diags_N=.FALSE.


      CHARACTER(LEN=3), PARAMETER :: cdrxn='RXN'  ! chars for reaction
      CHARACTER(LEN=3), PARAMETER :: cddep='DEP'  ! chars for deposition
      CHARACTER(LEN=3), PARAMETER :: cdems='EMS'  ! chars for emission
      CHARACTER(LEN=3), PARAMETER :: cdnet='NET'  ! chars for tendency
      CHARACTER(LEN=3), PARAMETER :: cdste='STE'  ! chars for strattrop exchange
      CHARACTER(LEN=3), PARAMETER :: cdmas='MAS'  ! chars for air mass
      CHARACTER(LEN=3), PARAMETER :: cdpsc='PSC'  ! chars for polar strat cloud
      CHARACTER(LEN=3), PARAMETER :: cdtpm='TPM'  ! chars for tropospheric mask
      CHARACTER(LEN=3), PARAMETER :: cdout='OUT'  ! chars for tracer output
      CHARACTER(LEN=1), PARAMETER :: cdbimol='B'  ! char for bimolecular rxn
      CHARACTER(LEN=1), PARAMETER :: cdtermol='T' ! char for termolecular rxn
      CHARACTER(LEN=1), PARAMETER :: cdphot='J'   ! char for photolytic rxn
      CHARACTER(LEN=1), PARAMETER :: cdhetero='H' ! char for heterogenous rxn
      CHARACTER(LEN=1), PARAMETER :: cdair='A'    ! char for aircraft emission
      CHARACTER(LEN=1), PARAMETER :: cdsurf='S'   ! char for surface emission
      CHARACTER(LEN=1), PARAMETER :: cdsulp='T'   ! char for 3D SO2 emission
      CHARACTER(LEN=1), PARAMETER :: cdvolc='V'   ! char for volcanic emission
      CHARACTER(LEN=1), PARAMETER :: cdlight='L'  ! char for lightning emission
      CHARACTER(LEN=1), PARAMETER :: cdwet='W'    ! char for wet deposition
      CHARACTER(LEN=1), PARAMETER :: cddry='D'    ! char for dry deposition
      CHARACTER(LEN=1), PARAMETER :: cdpsc_typ1='1' ! char for PSC type 1
      CHARACTER(LEN=1), PARAMETER :: cdpsc_typ2='2' ! char for PSC type 2
      CHARACTER(LEN=1), PARAMETER :: so2_emchar='S' !    "         "
      CHARACTER(LEN=3), PARAMETER :: calculate_tendency=cdnet
      CHARACTER(LEN=3), PARAMETER :: calculate_STE=cdste
      CHARACTER(LEN=3), PARAMETER :: cdrte='RTE'
      CHARACTER(LEN=3), PARAMETER :: cdlin='LIN'
      CHARACTER(LEN=3), PARAMETER :: cdlgt='LGT'
      CHARACTER(LEN=1), PARAMETER :: cdlgttot='T'
      CHARACTER(LEN=1), PARAMETER :: cdlgtc2g='G'
      CHARACTER(LEN=1), PARAMETER :: cdlgtc2c='C'
      CHARACTER(LEN=1), PARAMETER :: cdlgtN='N'      
      CHARACTER(LEN=1), PARAMETER :: cdignore='X'
      CHARACTER(LEN=1), PARAMETER :: aircraft_emissions='A'
      CHARACTER(LEN=1), PARAMETER :: lightning_emissions='L'
      CHARACTER(LEN=1), PARAMETER :: volcanic_emissions='V'
      CHARACTER(LEN=1), PARAMETER :: total_lightning_flashes=cdlgttot
      CHARACTER(LEN=1), PARAMETER :: lightning_flashes_cloud2ground=cdlgtc2g
      CHARACTER(LEN=1), PARAMETER :: lightning_flashes_cloud2cloud=cdlgtc2c
      CHARACTER(LEN=1), PARAMETER :: total_N_2D=cdlgtN
      


      INTEGER, PARAMETER :: number_of_reactants=2       ! No of reactants
      INTEGER :: icode                                  ! error code

! Interface section
      
      INTERFACE ASAD_SETSTASH_CHEMDIAG
        MODULE PROCEDURE ASAD_SETSTASH_CHEMDIAG
      END INTERFACE ASAD_SETSTASH_CHEMDIAG
      INTERFACE ASAD_INIT_CHEMDIAG
        MODULE PROCEDURE ASAD_INIT_CHEMDIAG
      END INTERFACE ASAD_INIT_CHEMDIAG
      INTERFACE ASAD_ALLOCATE_CHEMDIAG
        MODULE PROCEDURE ASAD_ALLOCATE_CHEMDIAG
      END INTERFACE ASAD_ALLOCATE_CHEMDIAG
      INTERFACE ASAD_CHEMICAL_DIAGNOSTICS
        MODULE PROCEDURE ASAD_CHEMICAL_DIAGNOSTICS
      END INTERFACE ASAD_CHEMICAL_DIAGNOSTICS
      INTERFACE ASAD_EMISSIONS_DIAGNOSTICS
        MODULE PROCEDURE ASAD_EMISSIONS_DIAGNOSTICS
      END INTERFACE ASAD_EMISSIONS_DIAGNOSTICS
      INTERFACE ASAD_3D_EMISSIONS_DIAGNOSTICS
        MODULE PROCEDURE ASAD_3D_EMISSIONS_DIAGNOSTICS
      END INTERFACE ASAD_3D_EMISSIONS_DIAGNOSTICS
      INTERFACE ASAD_TENDENCY_STE
        MODULE PROCEDURE ASAD_TENDENCY_STE
      END INTERFACE ASAD_TENDENCY_STE
      INTERFACE ASAD_MASS_DIAGNOSTIC
        MODULE PROCEDURE ASAD_MASS_DIAGNOSTIC
      END INTERFACE ASAD_MASS_DIAGNOSTIC
      INTERFACE ASAD_PSC_DIAGNOSTIC
        MODULE PROCEDURE ASAD_PSC_DIAGNOSTIC
      END INTERFACE ASAD_PSC_DIAGNOSTIC
      INTERFACE ASAD_OUTPUT_TRACER
        MODULE PROCEDURE ASAD_OUTPUT_TRACER
      END INTERFACE ASAD_OUTPUT_TRACER
      INTERFACE ASAD_FLUX_PUT_STASH
        MODULE PROCEDURE ASAD_FLUX_PUT_STASH
      END INTERFACE ASAD_FLUX_PUT_STASH


      CONTAINS

! #####################################################################

      SUBROUTINE ASAD_SETSTASH_CHEMDIAG(row_length,rows,model_levels)
      IMPLICIT NONE

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

      INTEGER, INTENT(IN) :: row_length     ! length of row
      INTEGER, INTENT(IN) :: rows           ! number of rows
      INTEGER, INTENT(IN) :: model_levels   ! number of levels

      INTEGER :: i,j,ierr,idiag   ! counters
      
      CHARACTER (LEN=70) :: cmessage          ! Error message

      LOGICAL, SAVE :: firsttime=.TRUE.     ! T for first call
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('ASAD_SETSTASH_CHEMDIAG',zhook_in,        &
                              zhook_handle)
      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED asad_setstash_chemdiag'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,IERR)
      END IF

      IF (firsttime) THEN

         ALLOCATE(chemD1codes(nmax_chemdiags))

! These are initialised after being read in from D1
         chemD1codes(:)%section    = IMDI
         chemD1codes(:)%item       = IMDI
         chemD1codes(:)%n_levels   = IMDI
         chemD1codes(:)%address    = IMDI
         chemD1codes(:)%length     = IMDI
         chemD1codes(:)%halo_type  = IMDI
         chemD1codes(:)%grid_type  = IMDI
         chemD1codes(:)%field_type = IMDI
! These are standard
         chemD1codes(:)%len_dim1   = row_length
         chemD1codes(:)%len_dim2   = rows
! Dry dep is single level - change when assigned, this will be re-set in 
! init_chemdiag
         chemD1codes(:)%len_dim3   = model_levels
         chemD1codes(:)%prognostic = .TRUE.
! Required set after D1 reading
         chemD1codes(:)%required   = .FALSE.
         
         ! this should now set up the model to all the stash requests that are 
         ! needed
         n_stashsumdiag = 0
         
         ! should have already been set above, but re-set to .false. here
         L_asad_use_chem_diags = .FALSE.
         
         DO I=1,nmax_chemdiags
            diag_loop: DO idiag=1,ndiag     ! search for stash requests
               IF ((modl_b(idiag) == submodel_for_sm(atmos_im))     &
                    .AND. ((isec_b(idiag) == asad_diag_sect)        & 
                    .AND. (item_b(idiag) == item1_asad_diags+I-1))) &
                    THEN
                  Chemd1codes(I)%section    = asad_diag_sect
                  Chemd1codes(I)%item       = item_b(idiag)
                  Chemd1codes(I)%len_dim1   = row_length
                  Chemd1codes(I)%len_dim2   = rows
                  ! check on size of model levels - currently only 
                  ! contiguous theta levels or surface allowed
                  IF (iopl_d(idom_b(idiag)) == 2) THEN  ! theta levels
                     Chemd1codes(I)%len_dim3 = model_levels
                  ELSE IF (iopl_d(idom_b(idiag)) == 5) THEN ! single level
                     Chemd1codes(I)%len_dim3 = 1
                  ELSE
                     WRITE(cmessage,'(A,5(i6,1X),A)') 'ERROR: ',        &
                          idiag,modl_b(idiag),                          &
                          isec_b(idiag),item_b(idiag),                  &
                          iopl_d(idom_b(idiag)),                        &
                          ' incorrect model levels'
                     ierr = isec_b(idiag)*1000 + item_b(idiag)
                     CALL EREPORT('ASAD_SETSTASH_CHEMDIAG',ierr,        &
                                  cmessage) 
                  END IF 
                  ! required is set to true here, but in the context of a
                  ! UKCA framework, it would be set to .false.
                  Chemd1codes(I)%required   = .true.
                  IF (.NOT. L_asad_use_chem_diags) &
                       L_asad_use_chem_diags = .TRUE.
                  Chemd1codes(I)%prognostic = .false.
                  n_stashsumdiag              = n_stashsumdiag + 1
                  EXIT diag_loop
               END IF
            END DO diag_loop
           IF (PrintStatus >= Prstatus_normal .AND.                     &
                    Chemd1codes(I)%section /= IMDI) THEN
               WRITE(6,*) 'ASAD_FLUXES: ',I,Chemd1codes(I)%section,&
                    Chemd1codes(I)%item,Chemd1codes(I)%len_dim1,&
                    Chemd1codes(I)%len_dim2,Chemd1codes(I)%len_dim3,&
                    Chemd1codes(I)%required             
            END IF
         END DO       ! I=1,nmax_chemdiags
         
         firsttime = .FALSE.
      END IF ! firstime
      
      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT asad_setstash_chemdiag'
         ! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF
      
    IF (lhook) CALL dr_hook('ASAD_SETSTASH_CHEMDIAG',zhook_out,         &
                            zhook_handle)
    RETURN
    END SUBROUTINE ASAD_SETSTASH_CHEMDIAG
   
! #####################################################################

      SUBROUTINE ASAD_ZERO_CHEMDIAG(diagnostic)
      IMPLICIT NONE


      TYPE(chemdiag), INTENT(INOUT) :: diagnostic

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('asad_zero_chemdiag',zhook_in,            &
                              zhook_handle)
      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED asad_zero_chemdiag'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF
      
      !     "zero" all fields
      diagnostic%location             = 0
      diagnostic%stash_number         = 0
      diagnostic%find_rxn_loc         = 0
      diagnostic%num_reactants        = 0
      diagnostic%num_products         = 0
      diagnostic%num_levs             = 0
      diagnostic%can_deallocate       = .FALSE.
      diagnostic%output_on_chem_tstep = .FALSE.
      diagnostic%tropospheric_mask    = .FALSE.
      diagnostic%molecular_mass       = 0.0
      diagnostic%c_vmr_to_mmr         = 0.0
      diagnostic%species              = blank
      diagnostic%diag_type            = blank
      diagnostic%rxn_type             = blank
      diagnostic%reactants(:)         = blank
      IF (ALLOCATED(diagnostic%products)) &
           diagnostic%products(:)          = blank
      IF (ALLOCATED(diagnostic%throughput)) &
           diagnostic%throughput(:,:,:)    = 0.0
      

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*)'mype=',mype,'LEFT asad_zero_chemdiag'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      IF (lhook) CALL dr_hook('ASAD_ZERO_CHEMDIAG',zhook_out,           &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_ZERO_CHEMDIAG

! #####################################################################

      SUBROUTINE ASAD_DEALLOC_CHEMDIAG(diagnostic)
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

      TYPE(chemdiag), INTENT(INOUT) :: diagnostic

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_DEALLOC_CHEMDIAG',zhook_in,         &
                              zhook_handle)

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED asad_dealloc_chemdiag'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      
      ! "deallocate" all fields
      ! set to integer missing data
      diagnostic%location             = IMDI
      diagnostic%stash_number         = IMDI
      diagnostic%find_rxn_loc         = IMDI
      diagnostic%num_reactants        = IMDI
      diagnostic%num_products         = IMDI
      diagnostic%num_levs             = IMDI
      diagnostic%can_deallocate       = .FALSE.
      diagnostic%output_on_chem_tstep = .FALSE.
      diagnostic%tropospheric_mask    = .FALSE.
      diagnostic%molecular_mass       = RMDI
      diagnostic%c_vmr_to_mmr         = RMDI
      diagnostic%species              = blank
      diagnostic%diag_type            = blank
      diagnostic%rxn_type             = blank
      diagnostic%reactants(:)         = blank
      IF (ALLOCATED(diagnostic%products)) &
           DEALLOCATE(diagnostic%products)
      IF (ALLOCATED(diagnostic%throughput)) &
           DEALLOCATE(diagnostic%throughput)
      

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT asad_dealloc_chemdiag'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      IF (lhook) CALL dr_hook('ASAD_DEALLOC_CHEMDIAG',zhook_out,        &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_DEALLOC_CHEMDIAG

! #####################################################################
      SUBROUTINE ASAD_INIT_CHEMDIAG(ierr)

      USE Control_Max_Sizes
      USE Submodel_Mod
      IMPLICIT NONE

! TYPSIZE start
!   Description:
!     This file contains sizes needed for dynamic allocation of
!   main data arrays within the model. Sizes read in from the user
!   interface via NAMELISTs are passed by /COMMON/. Other control
!   sizes that are fundamental in the definition of data structures
!   are assigned by PARAMETER statements.
!
!   Declarations for the NLSIZES namelist are also held in the module
!   nlsizes_namelist_mod. That module is currently only used by the
!   reconfiguration, while the UM uses this include file.
!
! All sizes
! Not dependent on sub-model
! DATA IN NAMLST#x MEMBER OF THE JOB LIBRARY
! ATMOS START
! Main sizes of fields for each submodel
! Grid-related sizes for ATMOSPHERE submodel.
INTEGER :: ROW_LENGTH           ! No of points per local row
INTEGER :: global_ROW_LENGTH    ! Points per global row
INTEGER :: ROWS                 ! No of local (theta) rows
INTEGER :: global_ROWS          ! No of global (theta) rows
INTEGER :: MODEL_LEVELS         ! No of model levels
INTEGER :: LAND_FIELD           ! No of land points in field
INTEGER :: NTILES               ! No of land surface tiles
INTEGER :: NICE                 ! No. of sea ice thickness categories
INTEGER :: NICE_USE             ! No. of sea ice categories used fully
                                !  in surface exchange and radiation
                                !  (If nice>1 & nice_use=1, categories only 
                                !  partially used in surface exchange)

! Physics-related sizes for ATMOSPHERE submodel
INTEGER :: WET_LEVELS          ! No of moist-levels
INTEGER :: CLOUD_LEVELS        ! No of cloud-levels
INTEGER :: ST_LEVELS           ! No of soil temperature levels
INTEGER :: SM_LEVELS           ! No of soil moisture levels
INTEGER :: BL_LEVELS           ! No of boundary-layer-levels
INTEGER :: OZONE_LEVELS        ! No of ozone-levels
INTEGER :: TPPS_OZONE_LEVELS   ! No of tropopause-ozone-levels
INTEGER :: RIVER_ROWS          ! No of rows for river routing
INTEGER :: RIVER_ROW_LENGTH    ! Row length for river routing
! Dynamics-related sizes for ATMOSPHERE submodel

INTEGER :: TR_LEVELS            ! No of tracer-levels
INTEGER :: TR_VARS              ! No of passive tracers
INTEGER :: TR_LBC_VARS          ! No of tracers in lbcs 
INTEGER :: TR_UKCA              ! No of UKCA tracers
INTEGER :: TR_LBC_UKCA          ! No of UKCA tracer lbcs 

! For Small executables

! Grid related sizes for data structure
! Data structure sizes for ATMOSPHERE submodel
INTEGER :: A_PROG_LOOKUP     ! No of prognostic fields
INTEGER :: A_PROG_LEN        ! Total length of prog fields
INTEGER :: A_LEN_INTHD       ! Length of INTEGER header
INTEGER :: A_LEN_REALHD      ! Length of REAL header
INTEGER :: A_LEN2_LEVDEPC    ! No of LEVEL-dependent arrays
INTEGER :: A_LEN2_ROWDEPC    ! No of ROW-dependent arrays
INTEGER :: A_LEN2_COLDEPC    ! No of COLUMN-dependent arrays
INTEGER :: A_LEN2_FLDDEPC    ! No of FIELD arrays
INTEGER :: A_LEN_EXTCNST     ! No of EXTRA scalar constants
INTEGER :: A_LEN_CFI1        ! Length of compressed fld index 1
INTEGER :: A_LEN_CFI2        ! Length of compressed fld index 2
INTEGER :: A_LEN_CFI3        ! Length of compressed fld index 3
! atmos end

! Data structure sizes for ATMOSPHERE ANCILLARY file control
! routines
INTEGER :: NANCIL_LOOKUPSA  ! Max no of fields to be read

! Data structure sizes for ATMOSPHERE INTERFACE file control
! routines
INTEGER :: N_INTF_A          ! No of atmosphere interface areas
INTEGER :: MAX_INTF_MODEL_LEVELS ! Max no of model levs in all areas
INTEGER :: MAX_LBCROW_LENGTH ! Max no of lbc row length in all areas
INTEGER :: MAX_LBCROWS ! Max no of lbc rows in all areas

!  Data structure sizes for ATMOSPHERE BOUNDARY file control
! routines

! Sizes applicable to all configurations (DUMPS/FIELDSFILE)

INTEGER :: PP_LEN_INTHD   ! Length of PP file integer header
INTEGER :: PP_LEN_REALHD  ! Length of PP file real    header


      ! Grid related sizes for COUPLING between ATMOS and OCEAN
      ! submodels [For MPP, sizes are 'global' values over all
      ! PEs.Also needed for river routing]
      INTEGER:: AOCPL_IMT                ! Ocean rowlength
      INTEGER:: AOCPL_JMT                ! Ocean no. of rows
      INTEGER:: AOCPL_ROW_LENGTH         ! Atmos rowlength
      INTEGER:: AOCPL_P_ROWS             ! Atmos no. of p rows

      COMMON/SIZE_AOCPL/                                                &
        AOCPL_IMT, AOCPL_JMT, AOCPL_ROW_LENGTH, AOCPL_P_ROWS

! Other sizes passed from namelist into common blocks
! Any additions to this common block must be mirrored in nlsizes_namelist_mod.
COMMON/NLSIZES/                                                     &
    ROW_LENGTH,global_ROW_LENGTH,ROWS,global_ROWS,                  &
    LAND_FIELD,MODEL_LEVELS,WET_LEVELS,                             &
    NTILES, NICE, NICE_USE,                                         &
    CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,           &
    OZONE_LEVELS,TPPS_OZONE_LEVELS,TR_VARS,TR_LBC_VARS,             &
    TR_UKCA,TR_LBC_UKCA,RIVER_ROWS,RIVER_ROW_LENGTH,                &
    A_PROG_LOOKUP,A_PROG_LEN,                                       &
    A_LEN_INTHD,A_LEN_REALHD,                                       &
    A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,                   &
    A_LEN2_FLDDEPC,A_LEN_EXTCNST,                                   &
    A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,                               &    
    NANCIL_LOOKUPSA,                                                &    
    N_INTF_A, MAX_INTF_MODEL_LEVELS, MAX_LBCROW_LENGTH,             &
    MAX_LBCROWS, PP_LEN_INTHD,PP_LEN_REALHD

!-----------------------------------------------------------------
! data in STASHC#x member of the job library

! Data structure sizes for ATMOSPHERE submodel (config dependent)
INTEGER :: A_LEN2_LOOKUP   ! Total no of fields (incl diags)
INTEGER :: A_LEN_DATA      ! Total no of words of data
INTEGER :: A_LEN_D1        ! Total no of words in atmos D1

! Size of main data array for this configuration

INTEGER :: LEN_TOT             ! Length of D1 array
INTEGER :: N_OBJ_D1_MAX         ! No of objects in D1 array

COMMON/STSIZES/                                                     &
    A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                              &
    LEN_TOT,N_OBJ_D1_MAX
! global (ie. dump version) of *_LEN_DATA
INTEGER :: global_A_LEN_DATA

COMMON /MPP_STSIZES_extra/ global_A_LEN_DATA
! Sizes of Stash Auxillary Arrays and associated index arrays
! Initialised in UMINDEX and UMINDEX_A/O/W
INTEGER :: LEN_A_IXSTS
INTEGER :: LEN_A_SPSTS

COMMON /DSIZE_STS/                                                  &
    LEN_A_IXSTS, LEN_A_SPSTS
!     The number of land points is computed for each PE
!     before the addressing section. All prognostics on land
!     points in the D1 space are now dimensioned by the local
!     no of land points rather than the global no of land points.

      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field
      ! ----------------------------------------------------------------
      ! extra variables not passed through user interface

      ! fundamental data sizes :
      ! Fundamental parameter  sizes of data structure
      ! Sizes applicable to all configurations (HISTORY FILE)

      ! Length of history file in dump
      INTEGER, PARAMETER :: LEN_DUMPHIST = 0

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      ! Length of dump fixed header
      INTEGER, PARAMETER :: LEN_FIXHD = 256

      ! Size of a single LOOKUP header
      INTEGER, PARAMETER :: LEN1_LOOKUP  = 64
      INTEGER, PARAMETER :: MPP_LEN1_LOOKUP= 2

      ! Size of compressed LBC LOOKUP (only used internally and
      ! contains just the items which change between each set of LBCs
      INTEGER, PARAMETER :: LEN1_LBC_COMP_LOOKUP = 8

      ! Sizes applicable to all configurations (STASH)
      ! Moved to typstsz.h

      INTEGER:: INTF_LEN2_LEVDEPC !1st dim of interface out lev dep cons
      INTEGER:: INTF_LEN2_ROWDEPC !2nd dim of interface out Row dep cons
      INTEGER:: INTF_LEN2_COLDEPC !2nd dim of interface out Col dep cons
      
      COMMON /DSIZE/                                                    &
        INTF_LEN2_LEVDEPC,INTF_LEN2_ROWDEPC,INTF_LEN2_COLDEPC
      ! sub-model atmosphere   :
      ! Data structure sizes derived from grid size
      INTEGER:: A_LEN1_LEVDEPC ! IN: 1st dim of level  dep const
      INTEGER:: A_LEN1_ROWDEPC ! IN: 1st dim of row    dep const
      INTEGER:: A_LEN1_COLDEPC ! IN: 1st dim of column dep const
      INTEGER:: A_LEN1_FLDDEPC ! IN: 1st dim of field  dep const

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER :: INTF_LOOKUPSA        ! No of interface lookups.
      COMMON /DSIZE_A/                                                  &
        A_LEN1_LEVDEPC,A_LEN1_FLDDEPC,A_LEN1_ROWDEPC,A_LEN1_COLDEPC,    &
        INTF_LOOKUPSA

      ! sub-model atmosphere   : derived sizes
      ! derived from model grid/levels. Arakawa B-grid

                                  ! Size of fields on THETA grid:
      INTEGER :: THETA_FIELD_SIZE     ! IN: with no halos
      INTEGER :: THETA_OFF_SIZE       ! IN: with simple halos
      INTEGER :: THETA_HALO_SIZE      ! IN: with extended halos

                                  ! Size of fields on U grid:
      INTEGER :: U_FIELD_SIZE         ! IN: with no halos
      INTEGER :: U_OFF_SIZE           ! IN: with simple halos
      INTEGER :: U_HALO_SIZE          ! IN: with extended halos

                                  ! Size of fields on V grid
      INTEGER :: V_FIELD_SIZE         ! IN: with no halos
      INTEGER :: V_OFF_SIZE           ! IN: with simple halos
      INTEGER :: V_HALO_SIZE          ! IN: with extended halos

      INTEGER :: N_ROWS               ! IN: No of V-rows
      INTEGER :: N_CCA_LEV            ! IN: No of CCA Levels
      COMMON/DRSIZE_A/                                                  &
        N_ROWS,N_CCA_LEV,THETA_FIELD_SIZE,U_FIELD_SIZE,V_FIELD_SIZE,    &
        THETA_OFF_SIZE,U_OFF_SIZE,V_OFF_SIZE,                           &
        THETA_HALO_SIZE,U_HALO_SIZE,V_HALO_SIZE
      ! boundary updating      : derived values
      ! Variables describing the Atmosphere Lateral Boundary Conditions
      ! Local (per processor) information


! TYPSIZE end
! TYPD1 Common block containing the ALT_N_SUBMODEL_PARTITION variables
! CALTSUBM
! TYPD1 needs access to N_SUBMODEL_PARTITION/_MAX in CSUBMODL. However,
! they are not always called in the same decks and in the right order.
! Therefore, copy the values to another file and include it from TYPD1

      INTEGER ALT_N_SUBMODEL_PARTITION

      INTEGER, PARAMETER :: ALT_N_SUBMODEL_PARTITION_MAX=1

      COMMON/CALTSUBM/ALT_N_SUBMODEL_PARTITION
! CALTSUBM end
! This file needs TYPSIZE included first

      REAL    ::  D1(LEN_TOT)       ! IN/OUT: Main data array
      LOGICAL :: LD1(LEN_TOT)       ! IN/OUT: Main data array (logical)
      INTEGER :: ID1(LEN_TOT)       ! I/OUT: Main data array (integer)

! D1_ADDR start
      ! Information for accessing D1 addressing array
      ! Number of items of info needed for each object and maximum
      ! number of objects in D1 -

      ! Number of items of information in D1 addressing array
      INTEGER,PARAMETER:: D1_LIST_LEN=17

! Names of items in D1 addressing array. Update D1_LIST_LEN above if
! items added

      ! Prognostic, Diagnostic, Secondary or other
      INTEGER,PARAMETER:: d1_object_type    = 1 ! Internal model id
      INTEGER,PARAMETER:: d1_imodl          = 2  ! Internal model id
      INTEGER,PARAMETER:: d1_section        = 3  ! Section
      INTEGER,PARAMETER:: d1_item           = 4  ! Item
      INTEGER,PARAMETER:: d1_address        = 5  ! Address in D1
      INTEGER,PARAMETER:: d1_length         = 6  ! Record length
      INTEGER,PARAMETER:: d1_grid_type      = 7  ! Grid type
      INTEGER,PARAMETER:: d1_no_levels      = 8  ! Number of levels

      ! Stash list number for diags. -1 for progs
      INTEGER,PARAMETER:: d1_stlist_no      = 9

      ! Pointer to dump header lookup table
      INTEGER,PARAMETER:: d1_lookup_ptr     = 10

      INTEGER,PARAMETER:: d1_north_code     = 11 ! Northern row
      INTEGER,PARAMETER:: d1_south_code     = 12 ! Southern row
      INTEGER,PARAMETER:: d1_east_code      = 13 ! Eastern row
      INTEGER,PARAMETER:: d1_west_code      = 14 ! Western row
      INTEGER,PARAMETER:: d1_gridpoint_code = 15 ! gridpoint info
      INTEGER,PARAMETER:: d1_proc_no_code   = 16 ! Processing Code
      INTEGER,PARAMETER:: d1_halo_type      = 17 ! Halo width type

      ! Types of items for d1_type

      INTEGER,PARAMETER:: prognostic = 0
      INTEGER,PARAMETER:: diagnostic = 1
      INTEGER,PARAMETER:: secondary  = 2
      INTEGER,PARAMETER:: other      = 3

! D1_ADDR end
      ! D1 addressing array and number of objects in each submodel
      INTEGER :: D1_ADDR(D1_LIST_LEN,N_OBJ_D1_MAX,                      &
     &  ALT_N_SUBMODEL_PARTITION)

      INTEGER :: NO_OBJ_D1(ALT_N_SUBMODEL_PARTITION_MAX)

      COMMON/common_D1_ADDRESS/ NO_OBJ_D1
! TYPD1 end
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

      INTEGER, INTENT(INOUT) :: ierr       ! Error code

      INTEGER :: i,j,k,cderr

      INTEGER    :: len               ! local dimension
      INTEGER    :: levs              ! number of levels
      INTEGER    :: section           ! stash section
      INTEGER    :: item              ! stash item
      INTEGER    :: addr                   ! stash address
      INTEGER    :: field_typ         ! Field type
      INTEGER    :: halo_typ          ! halo type
      INTEGER    :: grid_typ          ! grid type
      INTEGER    :: tag               ! stash tag
      INTEGER    :: m_atm_modl        ! sub model

      CHARACTER(LEN=10), DIMENSION(1:6) :: chems      ! Unused ?

      ! file handling
      CHARACTER(LEN=255) :: cdtmp,cdend
      CHARACTER(LEN=1) :: rxntyp
      CHARACTER(LEN=3) :: dgtyp
      INTEGER        :: stashval
      INTEGER        :: numrec
      INTEGER        :: cdunit
      INTEGER        :: rxnloc
      ! counters
      INTEGER :: icount,jdiag,istash

      INTEGER, EXTERNAL  :: GET_FLD_TYPE 

      CHARACTER(LEN=80) :: cmessage
      CHARACTER(LEN=72) :: outformat
      INTEGER       :: errcode

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('ASAD_INIT_CHEMDIAG',zhook_in,            &
                              zhook_handle)

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED asad_init_chemdiag'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      
      ! initially set to these - if not required then will be turned off.
      ! set up above, but also re-iterate here.
      L_asad_use_air_ems         = .FALSE.
      L_asad_use_light_ems       = .FALSE.
      L_asad_use_surf_ems        = .FALSE.
      L_asad_use_volc_ems        = .FALSE.
      L_asad_use_sulp_ems        = .FALSE.
      L_asad_use_flux_rxns       = .FALSE.
      L_asad_use_wetdep          = .FALSE.
      L_asad_use_drydep          = .FALSE.
      L_asad_use_tendency        = .FALSE.
      L_asad_use_mass_diagnostic = .FALSE.     
      L_asad_use_trop_mask       = .FALSE.


      ALLOCATE(stash_handling(1:n_stashsumdiag))
      stash_handling(:)%stash_value   = 0
      stash_handling(:)%stash_section = 0
      stash_handling(:)%stash_item    = 0
      stash_handling(:)%chemd1codes_location = 0
      stash_handling(:)%len_dim1 = 0
      stash_handling(:)%len_dim2 = 0
      stash_handling(:)%len_dim3 = 0
      stash_handling(:)%number_of_fields = 0
! Some fields are only output on chemical timesteps
      ! .FALSE. means all timesteps
      stash_handling(:)%output_on_chem_tstep = .FALSE. 
         
! Allocate diags array since we know how large it will
! need to be from STASH
! get stash item numbers - used to check/ch
! check to make sure they are required
      n_chemdiags = 0
      icount=0
      DO i=1,nmax_chemdiags
        IF (chemD1codes(i)%required) THEN
          icount=icount+1
          stash_handling(icount)%stash_section = &
                    chemd1codes(i)%section
          stash_handling(icount)%stash_item = &
                    chemd1codes(i)%item
          stash_handling(icount)%stash_value = &
                    ((stash_handling(icount)%stash_section)*1000) + &
                    stash_handling(icount)%stash_item
          stash_handling(icount)%chemd1codes_location = i
          stash_handling(icount)%len_dim1 = chemd1codes(i)%len_dim1
          stash_handling(icount)%len_dim2 = chemd1codes(i)%len_dim2
          stash_handling(icount)%len_dim3 = &
                    MAX(chemd1codes(i)%len_dim3,&
                   stash_handling(icount)%len_dim3)
          ! calculate how many diagnostics from flux_dat we are using
          DO j=1,SIZE(asad_chemical_fluxes)
            IF (stash_handling(icount)%stash_value ==             &
                       asad_chemical_fluxes(j)%stash_number) THEN
              n_chemdiags = n_chemdiags + 1
              ! NOTE: do not exit after first check, since we may have 
              !       more than one diagnostic associated with each  
              !       stash number.
            END IF
          END DO
          ! debug output
          IF (PrintStatus >= Prstatus_diag) THEN
                  WRITE(6,*) 'ASAD_STASH: ',icount,                   &
                    stash_handling(icount)%stash_section,               &
                    stash_handling(icount)%stash_item,                  &
                    stash_handling(icount)%len_dim1,                    &
                    stash_handling(icount)%len_dim2,                    &
                    stash_handling(icount)%len_dim3
          END IF
        END IF      ! chemD1codes(i)%required 
      END DO        ! i=1,nmax_chemdiags
         
      IF (L_asad_use_chem_diags) THEN

        IF (PrintStatus >= Prstatus_diag) WRITE(6,*)                  &
                 'mype= ',mype,' getting diagnostics'

            ! allocate diagnostics array
            ALLOCATE(asad_chemdiags(1:n_chemdiags))
            DO i=1,n_chemdiags
               CALL asad_zero_chemdiag(asad_chemdiags(i))
            END DO 

! Now we do a rather complicated do-loop to set up the asad_chemdiags array
            jdiag = 0
            icount=0
            DO istash=1,nmax_chemdiags
               ! check to see if stash is requested
               IF (chemD1codes(istash)%required) THEN
                  icount=icount+1
                  DO j=1,SIZE(asad_chemical_fluxes)
                     ! check if stash code has been requested, and if 
                     ! we have a corresponding flux specified
                     IF (stash_handling(icount)%stash_value ==          &
                          asad_chemical_fluxes(j)%stash_number) THEN
                        ! increment jdiag. We may have more than one 
                        ! diagnostic associated with each stash number
                        jdiag = jdiag + 1              
                        ! STASH number
                        asad_chemdiags(jdiag)%stash_number=&
                             asad_chemical_fluxes(j)%stash_number
                        ! diagnostic type
                        asad_chemdiags(jdiag)%diag_type = &
                             asad_cd_uppercase(&
                             asad_chemical_fluxes(j)%diag_type)
                        ! reaction (or otherwise) type
                        asad_chemdiags(jdiag)%rxn_type = &
                             asad_cd_uppercase(&
                             asad_chemical_fluxes(j)%rxn_type)
                        ! do we use a tropospheric mask only output the 
                        ! troposphere?
                        IF (asad_chemical_fluxes(&
                             j)%tropospheric_mask) THEN
                           asad_chemdiags(jdiag)%tropospheric_mask &
                                = .TRUE.
                           ! set logical for global value
                           IF (.NOT. L_asad_use_trop_mask) &
                                L_asad_use_trop_mask = .TRUE.
                        ELSE
                           asad_chemdiags(jdiag)%tropospheric_mask &
                                = .FALSE.
                        END IF
                        ! if we have two reactions with the exactly the same 
                        ! reactants and products, but different rates, 
                        ! find_reaction will only find the first one, this 
                        ! flag over-rides that
                        IF (asad_chemical_fluxes(j)%rxn_location &
                             <= 0) THEN
                           asad_chemdiags(jdiag)%find_rxn_loc = 1
                        ELSE
                           asad_chemdiags(jdiag)%find_rxn_loc = &
                                asad_chemical_fluxes(j)%rxn_location
                        END IF
                        ! fill in number of reactants and products
                        IF (asad_chemical_fluxes(j)%num_species  &
                            == 1) THEN 
                           ! are not dealing with a reaction
                           asad_chemdiags(jdiag)%num_reactants = 0
                           asad_chemdiags(jdiag)%num_products  = 0
                        ELSE IF (asad_chemical_fluxes(j)%num_species &
                             >= 3) THEN
                           ! we have a reaction (probably, or user error!)
                           asad_chemdiags(jdiag)%num_reactants = &
                                number_of_reactants ! = 2 currently
                           asad_chemdiags(jdiag)%num_products  = &
                                asad_chemical_fluxes(j)%num_species - &
                                asad_chemdiags(jdiag)%num_reactants
                        ELSE ! incorrect numbers specified!
                           WRITE(cmessage,*) &
                                'WRONG NUMBER OF SPECIES: ', &
                                asad_chemical_fluxes(j)%num_species, &
                                ' FOUND'
                           errcode=asad_chemdiags(jdiag)%stash_number
                           WRITE(6,*) cmessage, ' Error code: ', &
                                errcode,' PE: ',mype
                           CALL EREPORT('ASAD_CHEMICAL_DIAGNOSTICS', &
                                errcode,cmessage)
                        END IF
                        ! initialise reactants and products
                        asad_chemdiags(jdiag)%reactants(:) = blank
                        IF (ALLOCATED(asad_chemdiags(&
                             jdiag)%products)) &
                             DEALLOCATE(asad_chemdiags(&
                             jdiag)%products)
                        ALLOCATE(asad_chemdiags(jdiag)%products(1:&
                             asad_chemdiags(jdiag)%num_products))
                        asad_chemdiags(jdiag)%products(:) = blank
                        ! fill in from list - rxn specific
                        SELECT CASE (asad_chemdiags(jdiag)%diag_type)
                        CASE (cdrxn,cdrte)
                           ! copy over reacants
                           DO i=1,number_of_reactants
                              asad_chemdiags(jdiag)%reactants(i) = &
                                   trim(adjustl(&
                                   asad_chemical_fluxes(&
                                   j)%reactants(i)))
                           END DO
                           ! copy over products
                           DO i=1,asad_chemdiags(jdiag)%num_products
                              asad_chemdiags(jdiag)%products(i) = &
                                   trim(adjustl(&
                                   asad_chemical_fluxes(&
                                   j)%products(i)))
                           END DO
                        CASE (cddep,cdems,cdnet,cdste,cdout)
                           ! only one species in this case
                           asad_chemdiags(jdiag)%species = &
                                trim(adjustl(&
                                asad_chemical_fluxes(j)%reactants(1)))
                        CASE (cdpsc,cdmas,cdtpm,cdlgt,cdlin)
                           ! no need to do anything - no species
                        CASE DEFAULT
                           cmessage='DIAGNOSTIC TYPE '// &
                                asad_chemdiags(jdiag)%diag_type &
                                //' NOT FOUND'
                           errcode=asad_chemdiags(jdiag)%stash_number
                           WRITE(6,*) cmessage, ' Error code: ', &
                                errcode,' PE: ',mype
                        CALL EREPORT('ASAD_INIT_CHEMDIAG',              &
                                errcode,cmessage)
                        END SELECT

! Check on emissions, deposition etc, and set allocate/output logicals
                        SELECT CASE (asad_chemdiags(jdiag)%diag_type)
                        CASE (cdems)
                           ! emissions only on chemical timesteps, and can
                           ! deallocate at other times
                           asad_chemdiags(jdiag)%can_deallocate &
                                = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .TRUE.
                           SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
                           CASE (cdair) ! aircraft
                              L_asad_use_air_ems=.TRUE.
                           CASE (cdlight) ! lightning
                              L_asad_use_light_ems=.TRUE.
                           CASE (cdvolc) ! volcanic
                              L_asad_use_volc_ems=.TRUE.
                            CASE (cdsulp) ! volcanic
                              L_asad_use_sulp_ems = .TRUE.
                           CASE (cdsurf) ! surface
                              L_asad_use_surf_ems=.TRUE.
                           CASE DEFAULT
                              cmessage='INCORRECT EMISSION TYPE: '//&
                                   asad_chemdiags(jdiag)%diag_type//&
                                  ':'//asad_chemdiags(jdiag)%rxn_type
                              errcode=&
                                   asad_chemdiags(jdiag)%stash_number
                              WRITE(6,*) cmessage, ' Error code: ',   &
                                   errcode,' PE: ',mype
                              CALL EREPORT('ASAD_INIT_CHEMDIAG',        &
                                   errcode,cmessage)                           
                           END SELECT ! asad_chemdiags(jdiag)%rxn_type
                        CASE (cdrxn)
                           ! reactions only on chemical timesteps, and can
                           ! deallocate at other times
                          asad_chemdiags(jdiag)%can_deallocate = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .TRUE.
                           SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
                           CASE (cdbimol,cdtermol,cdphot,cdhetero)
                              L_asad_use_flux_rxns=.TRUE.
                           CASE DEFAULT
                              cmessage='INCORRECT REACTION TYPE: '// &
                                   asad_chemdiags(jdiag)%diag_type//&
                                   ':'// &
                                   asad_chemdiags(jdiag)%rxn_type
                              errcode=asad_chemdiags(&
                                   jdiag)%stash_number
                              WRITE(6,*) cmessage, ' Error code: ', &
                                   errcode,' PE: ',mype
                              CALL EREPORT('ASAD_CHEMICAL_DIAGNOSTICS',&
                                   errcode,cmessage)                           
                           END SELECT ! asad_chemdiags(jdiag)%rxn_type
                        CASE (cdrte)
                           ! rates only on chemical timesteps, and can
                           ! deallocate at other times
                           asad_chemdiags(jdiag)%can_deallocate &
                                = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .TRUE.
                           SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
                           CASE (cdbimol,cdtermol,cdphot,cdhetero)
                              L_asad_use_rxn_rates=.TRUE.
                           CASE DEFAULT
                              cmessage='INCORRECT RATE TYPE: '// &
                                   asad_chemdiags(jdiag)%diag_type//&
                                   ':'// &
                                   asad_chemdiags(jdiag)%rxn_type
                              errcode=asad_chemdiags(&
                                   jdiag)%stash_number
                              WRITE(6,*) cmessage, ' Error code: ', &
                                   errcode,' PE: ',mype
                              CALL EREPORT('ASAD_INIT_CHEMDIAG',        &
                                   errcode,cmessage)                           
                           END SELECT ! asad_chemdiags(jdiag)%rxn_type
                        CASE (cddep)
                           ! deposition only on chemical timesteps, and can
                           ! deallocate at other times
                          asad_chemdiags(jdiag)%can_deallocate = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .TRUE.
                           SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
                           CASE (cdwet)
                              L_asad_use_wetdep=.TRUE.
                           CASE (cddry)
                              L_asad_use_drydep=.TRUE.
                           CASE DEFAULT
                              cmessage='INCORRECT DEPOSTITION TYPE: '// &
                                   asad_chemdiags(jdiag)%diag_type//&
                                   ':'// &
                                   asad_chemdiags(jdiag)%rxn_type
                              errcode=asad_chemdiags(&
                                   jdiag)%stash_number
                              WRITE(6,*) cmessage, ' Error code: ',&
                                   errcode,' PE: ',mype
                              CALL EREPORT('ASAD_INIT_CHEMDIAG',        &
                                   errcode,cmessage)                           
                           END SELECT ! asad_chemdiags(jdiag)%rxn_type
                        CASE (cdnet)
                           ! tendency only on chemical timesteps, and can
                           ! deallocate at other times
                          asad_chemdiags(jdiag)%can_deallocate = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .TRUE.
                           ! tendency
                           L_asad_use_tendency = .TRUE.
                        CASE (cdste)
                           ! STE on all timesteps, and cannot be
                           ! deallocated at any time
                          asad_chemdiags(jdiag)%can_deallocate = .FALSE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .FALSE.
                           ! tendency
                           L_asad_use_STE = .TRUE.
                        CASE (cdmas)
                           ! air mass on all timesteps, and can
                           ! deallocate at other times
                          asad_chemdiags(jdiag)%can_deallocate = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .FALSE.
                           L_asad_use_mass_diagnostic=.TRUE.
                        CASE (cdpsc)
                           ! air mass on all timesteps, and can
                           ! deallocate at other times
                          asad_chemdiags(jdiag)%can_deallocate = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .TRUE.
                           L_asad_use_psc_diagnostic=.TRUE.
                        CASE (cdtpm)
                           ! Tropospheric mask on all timesteps, and can
                           ! deallocate at other times
                           asad_chemdiags(jdiag)%can_deallocate = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .FALSE.
                           ! set logical value to T for calculation of 
                           !tropospheric mask
                                L_asad_use_trop_mask = .TRUE.
                           ! set logical to true for outputting trop mask
                           L_asad_use_trop_mask_output=.TRUE.
                        CASE (cdout)
                           ! tracer on all timesteps, and cannot be
                           ! deallocated at any time
                           asad_chemdiags(jdiag)%can_deallocate &
                                = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .FALSE.
                           ! tracer
                           L_asad_use_output_tracer = .TRUE.
                        CASE (cdlgt)
                           ! lightning diagnostics on all timesteps, and cannot 
                           ! be deallocated at any time
                           asad_chemdiags(jdiag)%can_deallocate &
                                = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .FALSE.
                           ! check which diagnostics
                           SELECT CASE (asad_chemdiags(jdiag)%rxn_type)
                           CASE (cdlgttot)
                              L_asad_use_light_diags_tot = .TRUE.
                           CASE (cdlgtc2g)
                              L_asad_use_light_diags_c2g = .TRUE.
                           CASE (cdlgtc2c)
                              L_asad_use_light_diags_c2c = .TRUE.
                           CASE (cdlgtN)
                              L_asad_use_light_diags_N = .TRUE.
                           CASE default
                              cmessage=&
                                   'INCORRECT LIGHTNING DIAG TYPE: '//&
                                   asad_chemdiags(jdiag)%diag_type//&
                                   ':'//&
                                   asad_chemdiags(jdiag)%rxn_type
                              errcode=&
                                   asad_chemdiags(jdiag)%stash_number
                              WRITE(6,*) cmessage, ' Error code: ',&
                                   errcode,' PE: ',mype
                              CALL EREPORT(&
                                   'ASAD_CHEMICAL_DIAGNOSTICS',&
                                   errcode,cmessage)                           
                           END SELECT
                        CASE (cdlin)
                           ! Lightning N in kg(N)/km^2/s 
                           ! on all timesteps, and can
                           ! deallocate at other times
                           asad_chemdiags(jdiag)%can_deallocate &
                                = .TRUE.
                           asad_chemdiags(jdiag)%output_on_chem_tstep &
                                = .FALSE.
                           L_asad_use_LiN_diagnostic=.TRUE.
                        END SELECT ! asad_chemdiags(jdiag)%diag_type
                     END IF ! stash_handling(icount)%stash_value ==
                            ! asad_chemical_fluxes(j)%stash_number
                  END DO !j,SIZE(asad_chemical_diagnostics)
               END IF ! chemd1codes(i)%required
            END DO ! i,nmax_chemdiags

! Check on summing options - simplifies STASH output
            DO k=1,n_stashsumdiag
               DO i=1,n_chemdiags
                  ! calcs number of fields with the same stash code
                  IF (stash_handling(k)%stash_value ==  &
                       asad_chemdiags(i)%stash_number) THEN
                     stash_handling(k)%number_of_fields = &
                          stash_handling(k)%number_of_fields + 1
                  END IF
               END DO
               ! creates an array of this size
               ALLOCATE(stash_handling(k)%chemdiags_location(1: &
                    stash_handling(k)%number_of_fields))
               stash_handling(k)%chemdiags_location(:) = 0
               ! allocate %throughput array for each flux
               icount=0
               DO i=1,n_chemdiags
                  IF (stash_handling(k)%stash_value ==   &
                       asad_chemdiags(i)%stash_number) THEN
                     icount=icount+1
                     stash_handling(k)%chemdiags_location(icount) = i
                     ! get number of levels (needed in putstash)
                     asad_chemdiags(i)%num_levs = &
                          stash_handling(k)%len_dim3
                     ! set-up if outputting on chemical timesteps only 
                     ! - .FALSE. (i.e. all timesteps) initially
                     IF (asad_chemdiags(i)%output_on_chem_tstep) THEN
                        stash_handling(k)%output_on_chem_tstep = .TRUE.
                        ! i.e. if true for 1 diagnostic summed through 
                        ! STASH, then it's true for all
                     END IF                          
                  END IF
               END DO           ! i,n_chemdiags
            END DO              ! k,n_stashsumdiag
            ! error check
            DO k=1,n_stashsumdiag
           IF (stash_handling(k)%number_of_fields == 0) THEN
             cmessage=' STASH CODE NOT FOUND IN'//                      &
                      ' stash_handling ARRAY'
                  errcode = -1*stash_handling(k)%stash_value
                  WRITE(6,*) cmessage, ' Error code: ', &
                       errcode,' PE: ',mype,' K: ',k
             CALL EREPORT('ASAD_INIT_CHEMDIAG',errcode,cmessage)
               END IF
            END DO
            
         END IF                 ! L_asad_use_chem_diags
         


         ! print out diagnostics requested
         IF ((PrintStatus >= 2) .AND. (mype == 0)) THEN
            WRITE(6,'(A,A)') 'THE FOLLOWING UKCA FLUX DIAGNOSTICS ',&
                 'HAVE BEEN REQUESTED'
            DO i=1,n_chemdiags
               IF (asad_chemdiags(i)%num_products > 0) THEN
                  ! fancy formatting to cope with variable number of prods
                  IF (asad_chemdiags(i)%num_products < 10) THEN
                     WRITE(outformat,FMT='(A38,I1,A7)') &
                       '(I5,1X,A3,1X,A1,1X,A,1X,2(A,1X),A2,1X,',&
                       asad_chemdiags(i)%num_products,&
                       '(A,1X))'
                  ELSE
                     WRITE(outformat,FMT='(A38,I2,A7)') &
                          '(I5,1X,A3,1X,A1,1X,A,1X,2(A,1X),A2,1X,',&
                          asad_chemdiags(i)%num_products,&
                          '(A,1X))'                     
                  END IF
                  WRITE(6,FMT=trim(adjustl(outformat))) &
                       asad_chemdiags(i)%stash_number,&
                       asad_chemdiags(i)%diag_type,&
                       asad_chemdiags(i)%rxn_type,&
                       TRIM(ADJUSTL(asad_chemdiags(i)%species)),&
                       (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))), &
                       j=1,2),'->',&
                       (TRIM(ADJUSTL(asad_chemdiags(i)%products(j))), &
                       j=1,asad_chemdiags(i)%num_products)
               ELSE
                  ! not a reaction (or rather, has no products)
                  WRITE(6,FMT='(I5,1X,A3,1X,A1,1X,A,1X,2(A,1X))') &
                       asad_chemdiags(i)%stash_number,&
                       asad_chemdiags(i)%diag_type,&
                       asad_chemdiags(i)%rxn_type,&
                       TRIM(ADJUSTL(asad_chemdiags(i)%species)),&
                       (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))), &
                       j=1,2)
               END IF
            END DO
         END IF

         ! print out status of logicals
       IF (PrintStatus >= Prstatus_diag) THEN
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
                 'L_asad_use_air_ems         = ',L_asad_use_air_ems
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
                 'L_asad_use_light_ems       = ',L_asad_use_light_ems      
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
                 'L_asad_use_surf_ems        = ',L_asad_use_surf_ems       
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
                 'L_asad_use_volc_ems        = ',L_asad_use_volc_ems       
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
                 'L_asad_use_sulp_ems        = ',L_asad_use_sulp_ems
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',            &
                 'L_asad_use_flux_rxns       = ',L_asad_use_flux_rxns      
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
                 'L_asad_use_wetdep          = ',L_asad_use_wetdep         
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
            'L_asad_use_drydep          = ',L_asad_use_drydep         
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
                 'L_asad_use_tendency        = ',L_asad_use_tendency       
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
                 'L_asad_use_mass_diagnostic = ',&
                 L_asad_use_mass_diagnostic
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
                 'L_asad_use_trop_mask       = ',L_asad_use_trop_mask      
            WRITE(6,*) 'mype=',mype,' asad_init_chemdiag: ',&
                 'L_asad_use_chem_diags      = ',L_asad_use_chem_diags          
       END IF ! PrintStatus >= Prstatus_diag

      ierr = 0
      
      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT asad_init_chemdiag'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      IF (lhook) CALL dr_hook('ASAD_INIT_CHEMDIAG',zhook_out,           &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_INIT_CHEMDIAG

! #####################################################################
      FUNCTION ASAD_CD_UPPERCASE(string)
! To convert any string into uppercase letters

        IMPLICIT NONE     
      

      CHARACTER(len=*), intent(in) :: string
      CHARACTER(len=30)            :: asad_cd_uppercase
      INTEGER :: ia,iz,ishift
      INTEGER :: i,ic,lstring 

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_CD_UPPERCASE',zhook_in,zhook_handle)
      
      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED asad_cd_uppercase'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      ! required for shift
      ia     = ichar('a')
      iz     = ichar('z') 
      ishift = ichar('A')-ia
      ! will need to strip leading dashes
      lstring = len(string)
      
      asad_cd_uppercase = string       
      
! uppercase the string a letter at a time using the shift values calculated above
      DO i=1,30
         IC = ICHAR(asad_cd_uppercase(i:i))
         IF ((ic >= ia) .AND. (ic <= iz)) asad_cd_uppercase(i:i) =      &
              char(ishift+ic) 
      END DO
      
      asad_cd_uppercase = TRIM(adjustl(asad_cd_uppercase))
      
      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT asad_cd_uppercase'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      IF (lhook) CALL dr_hook('ASAD_CD_UPPERCASE',zhook_out,            &
                              zhook_handle)
      RETURN
      END FUNCTION ASAD_CD_UPPERCASE
      
! #####################################################################
      SUBROUTINE ASAD_ALLOCATE_CHEMDIAG(row_length,rows)
      IMPLICIT NONE
        
        
      INTEGER, INTENT(IN) :: row_length
      INTEGER, INTENT(IN) :: rows
        
      INTEGER :: i


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_ALLOCATE_CHEMDIAG',zhook_in,        &
                               zhook_handle)

      IF (PrintStatus >= Prstatus_diag) THEN
        WRITE(6,*) 'mype=',mype,'ENTERED asad_allocate_chemdiag'
! DEPENDS ON: UM_FORT_FLUSH
        CALL UM_FORT_FLUSH(6,ICODE)
      END IF
        
! Allocate flux array and initialise
      DO i=1,n_chemdiags
        IF (.NOT. ALLOCATED(asad_chemdiags(i)%throughput)) THEN
          ALLOCATE(asad_chemdiags(i)%throughput( &
                   1:row_length, &
                   1:rows, &
                   1:asad_chemdiags(i)%num_levs))
          asad_chemdiags(i)%throughput(:,:,:) = 0.0
        END IF
      END DO

      IF (PrintStatus >= Prstatus_diag) THEN
        WRITE(6,*)'mype=',mype,'LEFT asad_allocate_chemdiag'
! DEPENDS ON: UM_FORT_FLUSH
        CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      IF (lhook) CALL dr_hook('ASAD_ALLOCATE_CHEMDIAG',zhook_out,       &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_ALLOCATE_CHEMDIAG
      
! #####################################################################
      SUBROUTINE ASAD_CHEMICAL_DIAGNOSTICS(row_length, rows,            &
                  model_levels, klevel, secs_per_step, volume, ierr)

      USE UKCA_CONSTANTS,        ONLY: avogadro
      USE ukca_option_mod, ONLY: jpspec, jpbk, jptk, jphk, jppj
      USE Control_Max_Sizes
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: row_length       ! length of row
      INTEGER, INTENT(IN) :: rows             ! number of rows
      INTEGER, INTENT(IN) :: model_levels     ! the level that we are at
      INTEGER, INTENT(IN) :: klevel           ! the level that we are at
      REAL, INTENT(IN)    :: secs_per_step    ! timestep
      REAL, INTENT(IN)    :: volume(row_length, rows, model_levels)   ! cell volume
      INTEGER, INTENT(INOUT) :: ierr          ! error code

      REAL, PARAMETER :: convfac = 1.0D6/avogadro   ! conversion factor to convert
                                                    ! volume into cm^3 and result in mol

      LOGICAL, SAVE :: firstcall=.TRUE.

      INTEGER :: i,j,k,l,m,idep                     ! counters

      INTEGER, ALLOCATABLE, DIMENSION(:) :: findrxn_tmp   ! array to hold locations

      CHARACTER(LEN=72) :: cmessage
      CHARACTER(LEN=72) :: outformat
      INTEGER       :: errcode

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_CHEMICAL_DIAGNOSTICS',zhook_in,     &
                              zhook_handle)

      IF (PrintStatus >= Prstatus_diag) THEN
        WRITE(6,*) 'mype=',mype,'ENTERED asad_chemical_diagnostics,', &
                     ' Level: ',klevel
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      ierr = -1

      IF (firstcall) THEN

! Set up locations of reactions using (an altered version of) asad_findreaction
         DO i=1,n_chemdiags
            SELECT CASE (asad_chemdiags(i)%diag_type)
            CASE (cdrxn,cdrte)
               SELECT CASE (asad_chemdiags(i)%rxn_type)
               CASE (cdbimol) ! bimolecular reacions/rate
                  ALLOCATE(findrxn_tmp(1:(jpbk+1)))                 
                  findrxn_tmp = cd_findreaction(               &
                       asad_chemdiags(i)%num_products,         &
                       asad_chemdiags(i)%reactants,            &
                       asad_chemdiags(i)%products,             &
                       spb, nbrkx, (jpbk+1), jpspb )
                  asad_chemdiags(i)%location =                      &
                       findrxn_tmp(asad_chemdiags(i)%find_rxn_loc)
                  DEALLOCATE(findrxn_tmp)
                  IF (asad_chemdiags(i)%location <= 0) THEN
                     IF (asad_chemdiags(i)%num_products < 10)    &
                          THEN
                        WRITE(outformat,FMT='(A28,I1,A7)')          &
                             '(A13,1X,A1,1X,2(A,1X),A2,1X,',        &
                             asad_chemdiags(i)%num_products,        &
                             '(A,1X))'
                     ELSE
                        WRITE(outformat,FMT='(A28,I2,A7)')          &
                             '(A13,1X,A1,1X,2(A,1X),A2,1X,',        &
                             asad_chemdiags(i)%num_products,        &
                             '(A,1X))'
                     END IF
                     WRITE(6,*) outformat
                     WRITE(cmessage,fmt=trim(adjustl(outformat)))       &
                       'RXN NOT FOUND', cdbimol,                        &
                       (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))),  &
                        j=1,2),'->',&
                       (TRIM(ADJUSTL(asad_chemdiags(i)%products(j))),   &
                        j=1,asad_chemdiags(i)%num_products)
                     errcode=asad_chemdiags(i)%stash_number
                     WRITE(6,*) cmessage, ' Error code: ',errcode,    &
                          ' PE: ',mype
                     CALL EREPORT('ASAD_CHEMICAL_DIAGNOSTICS', &
                          errcode,cmessage)
                  END IF
               CASE (cdhetero) ! Heterogenous reactions/rates
                  allocate(findrxn_tmp(1:(jphk+1)))
                  findrxn_tmp = cd_findreaction( &
                       asad_chemdiags(i)%num_products, &
                       asad_chemdiags(i)%reactants, &
                       asad_chemdiags(i)%products, &
                       sph, nhrkx, (jphk+1), jpsph )
                  asad_chemdiags(i)%location = &
                       findrxn_tmp(asad_chemdiags(i)%find_rxn_loc)
                  deallocate(findrxn_tmp)
                  IF (asad_chemdiags(i)%location <= 0) THEN
                     IF (asad_chemdiags(i)%num_products < 10) THEN
                        WRITE(outformat,FMT='(A28,I1,A7)')            &
                             '(A13,1X,A1,1X,2(A,1X),A2,1X,',          &
                             asad_chemdiags(i)%num_products,          &
                             '(A,1X))'
                     ELSE
                        WRITE(outformat,FMT='(A28,I2,A7)')            &
                             '(A13,1X,A1,1X,2(A,1X),A2,1X,',          &
                             asad_chemdiags(i)%num_products,          &
                             '(A,1X))'
                     END IF
                     WRITE(cmessage,fmt=trim(adjustl(outformat)))       &
                        'RXN NOT FOUND', cdhetero,                      &
                        (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))), &
                        j=1,2),'->',&
                        (TRIM(ADJUSTL(asad_chemdiags(i)%products(j))),  &
                        j=1,asad_chemdiags(i)%num_products)
                     errcode=asad_chemdiags(i)%stash_number
                     CALL EREPORT('ASAD_CHEMICAL_DIAGNOSTICS',          &
                          errcode,cmessage)
                  END IF
               CASE (cdtermol) ! termolecular reactions/rates
                  ALLOCATE(findrxn_tmp(1:(jptk+1)))
                  findrxn_tmp = cd_findreaction(                        &
                       asad_chemdiags(i)%num_products,                  &
                       asad_chemdiags(i)%reactants,                     &
                       asad_chemdiags(i)%products,                      &
                       spt, ntrkx, (jptk+1), jpspt )
                  asad_chemdiags(i)%location =                          &
                       findrxn_tmp(asad_chemdiags(i)%find_rxn_loc)
                  DEALLOCATE(findrxn_tmp)
                  IF (asad_chemdiags(i)%location <= 0) THEN
                     IF (asad_chemdiags(i)%num_products < 10)           &
                          THEN
                        WRITE(outformat,FMT='(A28,I1,A7)')              &
                             '(A13,1X,A1,1X,2(A,1X),A2,1X,',            &
                             asad_chemdiags(i)%num_products,            &
                             '(A,1X))'
                     ELSE
                        WRITE(outformat,FMT='(A28,I2,A7)')              &
                             '(A13,1X,A1,1X,2(A,1X),A2,1X,',            &
                             asad_chemdiags(i)%num_products,            &
                             '(A,1X))'
                     END IF
                     WRITE(cmessage,fmt=trim(adjustl(outformat)))        &
                        'RXN NOT FOUND', cdtermol,                       &
                        (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))),  &
                        j=1,2),'->',                                     &
                        (TRIM(ADJUSTL(asad_chemdiags(i)%products(j))),   &
                        j=1,asad_chemdiags(i)%num_products)
                     errcode=asad_chemdiags(i)%stash_number
                     WRITE(6,*) cmessage, ' Error code: ',errcode,     &
                          ' PE: ',mype
                     CALL EREPORT('ASAD_CHEMICAL_DIAGNOSTICS',           &
                          errcode,cmessage)
                  END IF
               CASE (cdphot) ! photolysis reactions/rates
                  ALLOCATE(findrxn_tmp(1:(jppj+1)))
                  findrxn_tmp = cd_findreaction(                         &
                       asad_chemdiags(i)%num_products,                   &
                       asad_chemdiags(i)%reactants,                      &
                       asad_chemdiags(i)%products,                       &
                       spj, nprkx, (jppj+1), jpspj )
                  asad_chemdiags(i)%location =                           &
                       findrxn_tmp(asad_chemdiags(i)%find_rxn_loc)
                  DEALLOCATE(findrxn_tmp)
                  IF (asad_chemdiags(i)%location <= 0) THEN
                     IF (asad_chemdiags(i)%num_products < 10)            &
                          THEN
                        WRITE(outformat,FMT='(A28,I1,A7)')               &
                             '(A13,1X,A1,1X,2(A,1X),A2,1X,',             &
                             asad_chemdiags(i)%num_products,             &
                             '(A,1X))'
                     ELSE
                        WRITE(outformat,FMT='(A28,I2,A7)')               &
                             '(A13,1X,A1,1X,2(A,1X),A2,1X,',             &
                             asad_chemdiags(i)%num_products,             &
                             '(A,1X))'
                     END IF
                     WRITE(cmessage,fmt=trim(adjustl(outformat)))        &
                        'RXN NOT FOUND', cdphot,                         &
                        (TRIM(ADJUSTL(asad_chemdiags(i)%reactants(j))),  &
                        j=1,2),'->',                                     &
                        (TRIM(ADJUSTL(asad_chemdiags(i)%products(j))),   &
                        j=1,asad_chemdiags(i)%num_products)
                     errcode=asad_chemdiags(i)%stash_number
                     WRITE(6,*) cmessage, ' Error code: ',errcode,     &
                          ' PE: ',mype
                     CALL EREPORT('ASAD_CHEMICAL_DIAGNOSTICS',           &
                          errcode,cmessage)
                  END IF
               CASE DEFAULT
                  cmessage='REACTION TYPE '//                            &
                       asad_chemdiags(i)%rxn_type//' NOT FOUND'
                  errcode=asad_chemdiags(i)%stash_number
                  WRITE(6,*) cmessage, ' Error code: ',errcode,        &
                       ' PE: ',mype
                  CALL EREPORT('ASAD_CHEMICAL_DIAGNOSTICS',              &
                       errcode,cmessage)
               END SELECT       ! asad_chemdiags(i)%rxn_type
            CASE (cddep)
               SELECT CASE (asad_chemdiags(i)%rxn_type)
               CASE (cddry) ! DRY DEP
               ! %species contains species to be dry deposited
               ! %location contains location in requisite array
                  idep=1
                  DO j=1,jpspec
                    IF (ldepd(j) .AND. (asad_chemdiags(i)%species ==    &
                         speci(j))) THEN
                        ! it is being deposited logically and in diagnostic sense
                        asad_chemdiags(i)%location = j 
                        idep=idep+1
                     ELSE IF (.NOT. ldepd(j) .AND. &
                        (asad_chemdiags(i)%species == speci(j))) THEN
                      ! diagnostics of deposition requested, but species is not
                      ! labelled as being dry deposited in the chch_defs file
                      cmessage='SPECIES '//asad_chemdiags(i)%species    &
                             //' NOT DRY DEPOSITED'
                        errcode=asad_chemdiags(i)%stash_number
                        WRITE(6,*) cmessage, ' Error code: ',         &
                             errcode,' PE: ',mype
                        CALL EREPORT('ASAD_CHEMICAL_DIAGNOSTICS',       &
                             errcode,cmessage)
                     END IF
                  END DO               
               CASE (cdwet) ! WET DEP
                  DO j=1,jpspec
                   IF (ldepw(j) .AND. (asad_chemdiags(i)%species ==     &
                        speci(j))) THEN
                        asad_chemdiags(i)%location = j 
                        idep=idep+1
                     ELSE IF (.NOT. ldepw(j) .AND.                      &
                        (asad_chemdiags(i)%species == speci(j))) THEN
                        ! diagnostics of deposition requested, but species is not
                     ! labelled as being wet deposited in the chch_defs file
                     cmessage='SPECIES '//asad_chemdiags(i)%species     &
                             //' NOT WET DEPOSITED'
                        errcode=asad_chemdiags(i)%stash_number
                        WRITE(6,*) cmessage, ' Error code: ',         &
                             errcode,' PE: ',mype
                        CALL EREPORT('ASAD_CHEMICAL_DIAGNOSTICS',       &
                             errcode,cmessage)
                     END IF
                  END DO               
                  ! %species contains species to be wet deposited
                  ! %location contains location in requisite array
             END SELECT    ! asad_chemdiags(i)%rxn_type
           END SELECT    ! asad_chemdiags(i)%diag_type
            
               
            ! assumes 2 reactants
           IF (PrintStatus >= Prstatus_diag)  THEN
               SELECT CASE (asad_chemdiags(i)%diag_type)
               CASE (cdrxn,cdrte)
                  WRITE(6,*) 'mype=',mype,' CD:FIND: ',              &
                       asad_chemdiags(i)%stash_number,' ',             &
                       asad_chemdiags(i)%rxn_type,' : ',               &
                       asad_chemdiags(i)%reactants(1),                 &
                       ' + ',asad_chemdiags(i)%reactants(2),           &
                       ' rxn = ',                                      &
                       asad_chemdiags(i)%location
               CASE (cddep)
                  WRITE(6,*) 'mype=',mype,' CD:FIND: ',              &
                       asad_chemdiags(i)%stash_number,' ',             &
                       asad_chemdiags(i)%rxn_type,' : ',               &
                       asad_chemdiags(i)%species,                      &
                       ' deposition: ',                                &
                       asad_chemdiags(i)%location
               END SELECT
           END IF           ! Printstatus
           firstcall=.FALSE.
         END DO ! i
         
      END IF    ! firstcall

! Go through and pick up fluxes from ASAD arrays
! prk is in units of molecules.cm^-3.s^-1
      DO i=1,n_chemdiags
        IF (PrintStatus >= Prstatus_diag) THEN
          WRITE(6,*) 'Calculating Tput for i: ',i,                    &
                       asad_chemdiags(i)%stash_number,                  &
                       asad_chemdiags(i)%rxn_type,                      &
                       asad_chemdiags(i)%location
          CALL UM_FORT_FLUSH(6,ICODE)
        END IF

        SELECT CASE (asad_chemdiags(i)%diag_type)
         CASE (cdrxn)
            asad_chemdiags(i)%throughput(:,:,klevel) =                  &
              RESHAPE(prk(:,asad_chemdiags(i)%location),                &
                      (/row_length,rows/))*volume(:,:,klevel)*convfac
            IF (asad_chemdiags(i)%tropospheric_mask) THEN
              WHERE (.NOT. L_troposphere(:,:,klevel))
                asad_chemdiags(i)%throughput(:,:,klevel) = 0.0
              ENDWHERE
            END IF
         CASE (cddep)
            SELECT CASE (asad_chemdiags(i)%rxn_type)
            CASE (cddry) ! DRY DEP
               IF (klevel <=                                            &
                 SIZE(asad_chemdiags(i)%throughput(:,:,:),DIM=3)) THEN
              ! if 2D then only take lowest level, otherwise will be 3D
                  asad_chemdiags(i)%throughput(:,:,klevel)=             &
                       RESHAPE(dpd(:,asad_chemdiags(i)%location),       &
                              (/row_length,rows/))*                     &
                       RESHAPE(y(:,asad_chemdiags(i)%location),         &
                              (/row_length,rows/))*                     &
                              volume(:,:,klevel)*convfac
                  IF (asad_chemdiags(i)%tropospheric_mask) THEN
                     WHERE (.NOT. L_troposphere(:,:,klevel))
                        asad_chemdiags(i)%throughput(:,:,klevel) = 0.0
                     ENDWHERE
                  END IF
              END IF
! Not needed (?)
            CASE (cdwet) ! WET DEP
               asad_chemdiags(i)%throughput(:,:,klevel)=             &
                    RESHAPE(dpw(:,asad_chemdiags(i)%location),       &
                           (/row_length,rows/))*                     &
                    RESHAPE(y(:,asad_chemdiags(i)%location),         &
                           (/row_length,rows/))*                     &
                           volume(:,:,klevel)*convfac
! Not needed (?)
               IF (asad_chemdiags(i)%tropospheric_mask) THEN
                  WHERE (.NOT. L_troposphere(:,:,klevel))
                    asad_chemdiags(i)%throughput(:,:,klevel) = 0.0
                  ENDWHERE
                END IF
            END SELECT
         END SELECT
      END DO                ! i=1,n_chemdiags
         

      ierr = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT asad_chemical_diagnostics'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      IF (lhook) CALL dr_hook('ASAD_CHEMICAL_DIAGNOSTICS',zhook_out,    &
                              zhook_handle)
      RETURN

      CONTAINS

! #####################################################################

      FUNCTION CD_FINDREACTION( numprods, reactants, products,          &
                               reactions, r_index, nr, ns )

!     cd_findreaction   - does a search for a specified reaction in the given
!                      ratefile. Note that this function is clever enough to
!                      look for all the various combinations in which the
!                      species might appear.

!     Glenn Carver, Centre for Atmospheric Science, University of Cambridge
!
!     ASAD: cd_findreaction            Version: cd_findreaction.f 1.2 12/21/01

!     Purpose
!     -------
!     Given a reaction, this function looks for a matching entry in the ratefile.
!     The reactants and products given to this routine are not order dependent.
!     This function will look for the reactants in any order, ie. react1 / react2
!     or react2 / react1. And likewise for the products. If a reaction only has
!     two products, then specify the third product as a blank string.
!
!     Interface
!     ---------
!     react1 & react2   - character strings holding the two reactants of the
!                         reaction to look for.
!     prod1, prod2 & prod3 - character strings holding the 3 products of the
!                         reaction to look for. If there are less than 3 products
!                         then give a blank string.
!     reactions         - character array holding the array to search through.
!                         Dimensioned as (Nreacts,Entries). Nreacts (nr)is the no. 
!                         of reactions. Entries (ns) will be the total no. of 
!                         reactants and products (e.g. 2 reactants + 3 products = 5 )
!     index             - this is an integer array dimensioned as (nreacts). ASAD
!                         reorders the reactions from the external ratefiles for
!                         computational efficiency so the order they are stored
!                         is not the same order they appear in the ratefile. Pass
!                         the appropriate indexing array in order that the correct
!                         position of this reaction in the ascii ratefile can be
!                         returned. e.g. use one of the arrays, nbrkx, ntrkx etc
!
!     NOTE!!  Because of the way the string matching is done in this routine, you
!     **MUST** pass the strings with trailing spaces. e.g. if the string length
!     allowed for species names is 6, then the argument list must look like:
!
!     i = cd_findreaction( 'OH    ', 'CO2   ', 'HOCO2 ', '      ', '      ', .... )
!
!     If this function finds a matching reaction, it will return the index
!     to that reaction in the array reactions. If it doesn't find a match it
!     returns 0.
!
!----------------------------------------------------------------------
!
      IMPLICIT NONE
      
      ! number of products in product array
      INTEGER, INTENT(in)           :: numprods
      INTEGER, INTENT(IN)           :: nr
      INTEGER, INTENT(IN)           :: ns
      CHARACTER(LEN=10), INTENT(IN) :: reactants(number_of_reactants)
      CHARACTER(LEN=10), INTENT(IN) :: products(numprods)
      ! array containing reactions
      CHARACTER(len=10), intent(in) :: reactions(nr,ns)
      INTEGER, INTENT(IN)           :: r_index(nr)

      !Local
      CHARACTER(LEN=10), ALLOCATABLE :: prods(:)     ! found products
      CHARACTER(LEN=10), ALLOCATABLE :: inprods(:)   ! input products
      LOGICAL, ALLOCATABLE           :: found(:)     ! T if prods==inprods
      LOGICAL                        :: final        ! T if all prods match inprods

      INTEGER :: j, jp, jp2, jr
      INTEGER :: nreacts, nprods, iprods, ipargs, jprod

      INTEGER, DIMENSION(1:nr) :: cd_findreaction 
      INTEGER :: find_count, i

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('CD_FINDREACTION',zhook_in,zhook_handle)

      cd_findreaction(:) = -1 ! set to -1 to see if haven't got it at end
      find_count = 0

!     -----------------------------------------------------------------
!          Search for reactants first.
!          ------ --- --------- ------
!
      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED cd_findreaction'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      nreacts = size(reactions,DIM=1)
      nprods  = size(reactions,DIM=2) - number_of_reactants ! always 2 reactants
      
      IF (nprods < numprods) THEN ! if number of products is not equal to number of 
                                  ! number of products in reactions array then exit
         cmessage='wrong number of products in find reaction'
         errcode=nprods
         WRITE(6,*) cmessage, ' Error code: ',errcode,' PE: ',mype
         CALL EREPORT('cd_findreaction',errcode,cmessage)
      END IF

      ALLOCATE(prods(nprods), inprods(nprods), found(nprods) )

      !  Get the non-blank products that we're looking for.
      prods(:) = blank
      inprods(:) = blank
      ipargs = 0
      inprods = blank
      DO jprod=1,numprods ! cycle over number of products
         IF ( products(jprod) /= blank ) THEN
            ipargs = ipargs + 1
            inprods(ipargs) = products(jprod)
         END IF
      END DO

      j = 1
      find_loop: DO
        jr = r_index(j)            ! loop over reactions

        IF (PrintStatus >= Prstatus_diag) THEN 
         WRITE(6,*) 'mype=',mype,' cd_findreaction ',            &
              jr,j,(i,':',reactions(j,i),i=1,ns)
         WRITE(6,*) 'mype=',mype,' cd_findreaction ',            &
              jr,j,(i,':',reactants(i),i=1,number_of_reactants)
         WRITE(6,*) 'mype=',mype,' cd_findreaction ',            &
              jr,j,(i,':',products(i),i=1,numprods)
         ! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      ! ASSUME 2 REACTANTS HERE!!
      IF ((reactions(j,1)==reactants(1) .AND.                      &
           reactions(j,2)==reactants(2)) .OR.                      &
           (reactions(j,1)==reactants(2) .AND.                     & 
           reactions(j,2)==reactants(1))) THEN

!          Found a reaction with matching reactants, now check all the
!          possible cases that the products could match. Need to allow for
!          rate files which have varying no. of products.

                prods = blank
                found = .false.
                iprods = 0
                DO jp = 1, nprods   ! only copy the non-blank products
                  IF ( reactions(j,2+jp) /= blank ) THEN
                    iprods = iprods + 1
                    prods(iprods) = reactions(j,2+jp)
                  END IF
                END DO

! If no. of nonblank products doesn't match then try next reaction.
                IF ( iprods /= ipargs ) THEN
                    j = j + 1
                    IF ( j > nreacts ) THEN
                      IF (lhook) CALL dr_hook('CD_FINDREACTION',        &
                                          zhook_out,zhook_handle)
                      RETURN
                    END IF
                    CYCLE
                END IF

! Otherwise check the names of the products

                DO jp = 1, iprods
                  DO jp2 = 1, iprods
                    IF (.NOT. found(jp) .AND. inprods(jp)==prods(jp2))  &
                      THEN
                      found(jp) = .true.
                      prods(jp2) = blank
                    END IF
                  END DO
                END DO

! Check to see if we have found all the products. If we have then exit

                final = .TRUE.
                DO jp = 1, iprods
                  final = final .AND. found(jp)
                END DO

                IF ( final ) THEN
                   find_count = find_count+1
                   cd_findreaction(find_count) = jr
           ! don't return from here since need to see if there is 
           ! another reaction with the same reactants & products
                END IF
              END IF

! next reaction

        j = j + 1
        IF ( j > nreacts ) THEN
          EXIT find_loop
        END IF
      END DO find_loop

      IF (ALLOCATED(prods)) DEALLOCATE( prods )
      IF (ALLOCATED(inprods)) DEALLOCATE( inprods )
      IF (ALLOCATED(found)) DEALLOCATE(found)

      IF (PrintStatus >= Prstatus_diag) THEN
           WRITE(6,*) 'mype=',mype,'LEFT cd_findreaction'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      IF (lhook) CALL dr_hook('CD_FINDREACTION',zhook_out,zhook_handle)
      RETURN
      END FUNCTION CD_FINDREACTION

! #####################################################################

      INTEGER FUNCTION CD_FINDSPECIESLOC(react)
      IMPLICIT NONE

      CHARACTER(LEN=10), INTENT(IN) :: react

      INTEGER :: i

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('CD_FINDSPECIESLOC',zhook_in,zhook_handle)

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED cd_findspeciesloc'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      cd_findspeciesloc = 0

      DO i=1,jpspec
        IF (react == speci(i)) THEN
          cd_findspeciesloc = i
          IF (lhook) CALL dr_hook('CD_FINDSPECIESLOC',zhook_out,        &
                                  zhook_handle)
          RETURN
        END IF
      END DO

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT cd_findspeciesloc'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      IF (lhook) CALL dr_hook('CD_FINDSPECIESLOC',zhook_out,            &
                              zhook_handle)
      RETURN
      END FUNCTION cd_findspeciesloc
      END SUBROUTINE ASAD_CHEMICAL_DIAGNOSTICS

! #####################################################################
      SUBROUTINE ASAD_EMISSIONS_DIAGNOSTICS(row_length, rows,           &
                                       model_levels, n_chem_tracers,    &
                                       em_field, surf_area,             &
                                       timestep, n_emissions,           &
                                       n_boundary_vals, em_spec,        &
                                       lbc_spec, molmass,               &
                                       lbc_molmass, ierr)

      USE UKCA_CONSTANTS,    ONLY: c_no, c_no2
      USE ukca_option_mod, ONLY: jpctr
      USE Control_Max_Sizes
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: row_length                        ! array dimensions
      INTEGER, INTENT(IN) :: rows
      INTEGER, INTENT(IN) :: model_levels
      INTEGER, INTENT(IN) :: n_chem_tracers                    ! number of tracers
      REAL, INTENT(IN)    :: em_field(row_length,rows,n_chem_tracers)
      REAL, INTENT(IN)    :: surf_area(row_length,rows)        ! cell area
      REAL, INTENT(IN)    :: timestep
      INTEGER, INTENT(IN) :: n_emissions                       ! no of emitted species
      INTEGER, INTENT(IN) :: n_boundary_vals
      CHARACTER(LEN=10), INTENT(IN) :: em_spec(n_emissions)
      CHARACTER(LEN=10), INTENT(IN) :: lbc_spec(n_boundary_vals)
      REAL, DIMENSION(n_emissions), INTENT(IN) :: molmass
      REAL, DIMENSION(n_boundary_vals), INTENT(IN) :: lbc_molmass

      INTEGER, INTENT(OUT) :: ierr                               ! error code

      INTEGER :: i,j,l,ntracer,mspecies                          ! loop variables

      LOGICAL, DIMENSION(jpctr) :: L_cd_emitted                  ! T for emission diags 
      
      CHARACTER(LEN=72) :: cmessage
      INTEGER       :: errcode

      LOGICAL, SAVE :: firstcall=.TRUE.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_EMISSIONS_DIAGNOSTICS',zhook_in,    &
                              zhook_handle)

      ierr = -1

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED ', &
              'asad_emissions_diagnostics'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF



! Set up logicals and molmass for use in later calculations
      IF (firstcall) THEN

         L_cd_emitted(:) =.FALSE.
         DO i=1,n_chemdiags
          IF (asad_chemdiags(i)%diag_type == cdems) THEN
               DO ntracer=1,jpctr
                  DO mspecies=1,n_emissions
                     ! work out which molmass to take!
                IF ((advt(ntracer) == em_spec(mspecies)) .AND.          &
                   (asad_chemdiags(i)%species==em_spec(mspecies))) THEN 
                  asad_chemdiags(i)%molecular_mass = molmass(mspecies)
                  asad_chemdiags(i)%c_vmr_to_mmr = c_species(ntracer)
                        L_cd_emitted(ntracer) = .TRUE.
                  IF (asad_chemdiags(i)%species == 'NO        ')        &
                             asad_chemdiags(i)%molecular_mass = &
                             asad_chemdiags(i)%molecular_mass*c_no/c_no2
                     END IF
                  END DO        ! mspecies,n_emissions
                  DO mspecies=1,n_boundary_vals
                     ! work out which molmass to take!, can overwrite with
                     ! lbc_molmass, since emmissions_ctl does this anyway
                IF ((advt(ntracer) == lbc_spec(mspecies)) .AND.         &
                    (asad_chemdiags(i)%species == lbc_spec(mspecies)))  &
                           THEN     
                  asad_chemdiags(i)%molecular_mass=lbc_molmass(mspecies)
                  asad_chemdiags(i)%c_vmr_to_mmr = c_species(ntracer)
                        L_cd_emitted(ntracer) = .TRUE.
                        IF (asad_chemdiags(i)%species == 'NO        ') &
                             asad_chemdiags(i)%molecular_mass = &
                             asad_chemdiags(i)%molecular_mass*c_no/c_no2
                     END IF
                  END DO        ! mspecies,n_boundary_vals

                  IF (L_cd_emitted(ntracer)) THEN
                IF (advt(ntracer) == asad_chemdiags(i)%species)         &
                          asad_chemdiags(i)%location = ntracer
                  ELSE IF ((.NOT. L_cd_emitted(ntracer)) .AND.          &
                     (advt(ntracer) == asad_chemdiags(i)%species)) THEN
                    cmessage='SPECIES '// asad_chemdiags(i)%species     &
                          //' NOT EMITTED'
                    errcode=asad_chemdiags(i)%stash_number
                    WRITE(6,*) cmessage, ' Error code: ',             &
                          errcode,' PE: ',mype
                    CALL EREPORT('ASAD_EMISSIONS_DIGNOSTICS',           &
                          errcode,cmessage)    
                  END IF 
               END DO           ! ntracer,jpctr
          END IF              ! %diag_type == cdems
         END DO                 ! i,n_chemdiags
         firstcall=.FALSE.
      END IF                    ! first

      DO i=1,n_chemdiags
        IF (asad_chemdiags(i)%diag_type == cdems) THEN
          ! emissions are 2D and converted to (kg/m2/s)
          asad_chemdiags(i)%throughput(:,:,1) =                         &
              em_field(:,:,asad_chemdiags(i)%location)*surf_area(:,:)*  &
              (1000.0/asad_chemdiags(i)%molecular_mass)
        END IF                 ! %diag_type
      END DO                    ! i,n_chemdiags
      
            
      ierr = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT ', &
              'asad_emissions_diagnostics'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      IF (lhook) CALL dr_hook('ASAD_EMISSIONS_DIAGNOSTICS',zhook_out,   &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_EMISSIONS_DIAGNOSTICS

! #####################################################################
      SUBROUTINE ASAD_3D_EMISSIONS_DIAGNOSTICS( row_length, rows,       &
           model_levels, i_spec_emiss, em_field_3d_in, surf_area,       &
           total_number_density, volume, mass, em_molmass, timestep,    &
           emission_type, ierr)

      USE UKCA_CONSTANTS,      ONLY: avogadro
      IMPLICIT NONE


      INTEGER, INTENT(IN) :: row_length                      ! Array dimensions
      INTEGER, INTENT(IN) :: rows
      INTEGER, INTENT(IN) :: model_levels
      INTEGER, INTENT(IN) :: i_spec_emiss
      REAL, INTENT(IN) :: em_field_3d_in(1:row_length,1:rows,           &
                                         1:model_levels)
      REAL, INTENT(IN) :: surf_area(1:row_length,1:rows)
      REAL, INTENT(IN) :: total_number_density(1:row_length,1:rows,     &
                                               1:model_levels)
      REAL, INTENT(IN) :: volume(1:row_length,1:rows,1:model_levels)
      REAL, INTENT(IN) :: mass(1:row_length,1:rows,1:model_levels)
      REAL, INTENT(IN) :: em_molmass
      REAL, INTENT(IN) :: timestep
      CHARACTER(LEN=1), INTENT(IN) :: emission_type
      INTEGER, INTENT(OUT) :: ierr

      REAL, DIMENSION(1:row_length,1:rows,1:model_levels) :: em_field_3d


      integer :: i,j,l,klevel

      CHARACTER(LEN=72) :: cmessage
      integer       :: errcode

      LOGICAL, SAVE :: firstcall=.TRUE.
      LOGICAL       :: spec_emitted

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_3D_EMISSIONS_DIAGNOSTICS',zhook_in, &
                              zhook_handle)

      ierr = -1

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED ', &
              'asad_3D_emissions_diagnostics'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF



      IF (firstcall) THEN
         spec_emitted = .FALSE.
         DO i=1,n_chemdiags
          IF (asad_chemdiags(i)%diag_type == cdems .AND.                &
              asad_chemdiags(i)%rxn_type == emission_type .AND.         &
              asad_chemdiags(i)%species == advt(i_spec_emiss)) THEN
                     asad_chemdiags(i)%location = i_spec_emiss
                     asad_chemdiags(i)%molecular_mass = em_molmass
            asad_chemdiags(i)%c_vmr_to_mmr   = c_species(i_spec_emiss)
                     spec_emitted = .TRUE.
          END IF         ! diag_type == cdems &etc
         END DO
         IF (.NOT. spec_emitted) THEN 
            cmessage='SPECIES '//asad_chemdiags(i)%species              &
                 //' NOT EMITTED IN 3D'
            errcode=asad_chemdiags(i)%stash_number
            WRITE(6,*) cmessage, ' Error code: ',errcode,             &
                 ' PE: ',mype
            CALL EREPORT('ASAD_3D_EMISSIONS_DIGNOSTICS',                &
                 errcode,cmessage)    
         END IF            
         firstcall = .FALSE.
      END IF                    ! firstcall

      ! need to do some conversion because of the way that things are done in 
      ! UKCA for 3D emissions
      em_field_3d(:,:,:) = -1.0d0 ! i.e. if values are -ve, have a problem
      ! sanity check, just in case
      !em_field_3d(:,:,:) = em_field_3d_in(:,:,:)
      DO i=1,n_chemdiags
         SELECT CASE (asad_chemdiags(i)%diag_type)
         CASE (cdems)
            IF ((asad_chemdiags(i)%rxn_type == aircraft_emissions)     &
                 .AND. (emission_type == aircraft_emissions)) THEN
               ! these are in kg(NO2)/m2/s 
               ! - need to multiply by the surface area and / air mass 
              DO klevel=1,model_levels
                DO j=1,rows
                  DO l=1,row_length
                    em_field_3d(l,j,klevel) =                           &
                     (em_field_3d_in(l,j,klevel)/mass(l,j,klevel))*     &
                      surf_area(l,j)
                  END DO ! l
                END DO ! j
              END DO ! klevel
            ELSE IF ((asad_chemdiags(i)%rxn_type ==                     &
                 lightning_emissions) .AND. (emission_type ==           &
                 lightning_emissions)) THEN
               ! do nothing - in kg(NO)/kg(air)/gridcell/s
               em_field_3d(:,:,:) = em_field_3d_in(:,:,:)
            ELSE IF ((asad_chemdiags(i)%rxn_type ==                     &
                 volcanic_emissions) .AND. (emission_type ==            &
                 volcanic_emissions)) THEN
               ! volcanic emissions are kg(SO2)/kg(air)/gridcell/timestep
               ! - so divide by timestep to get /s
               em_field_3d(:,:,:) = em_field_3d_in(:,:,:)/timestep
            ELSE
               ! do nothing - assume kg(Species)/kg(air)/gridcell/s
               em_field_3d(:,:,:) = em_field_3d_in(:,:,:)
            END IF
         END SELECT
      END DO


! Volume calculation method, converts to mol per second
      DO i=1,n_chemdiags
        IF (asad_chemdiags(i)%diag_type == cdems) THEN
          IF ((advt(i_spec_emiss) == asad_chemdiags(i)%species) .AND.   &
                   (emission_type == asad_chemdiags(i)%rxn_type)) THEN
            asad_chemdiags(i)%throughput(:,:,:)= (em_field_3d(          &
                           1:row_length,1:rows,1:model_levels)          &
                           *total_number_density(1:row_length,1:rows,   &
                           1:model_levels)*volume(1:row_length,1:rows,  &
                           1:model_levels)                              &
                           /(avogadro*asad_chemdiags(i)%c_vmr_to_mmr))  &
                           /timestep
            IF (asad_chemdiags(i)%tropospheric_mask) THEN
              WHERE (.NOT. L_troposphere)
                asad_chemdiags(i)%throughput(:,:,:) = 0.0
              END WHERE
            END IF
            END IF              ! advt==%species
        END IF             ! %diag_type == cdems
      END DO                    ! i,n_chemdiags
      
      ierr = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT ', &
              'asad_3D_emissions_diagnostics'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

 
      IF (lhook) CALL dr_hook('asad_3D_emissions_diagnostics',          &
                              zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE asad_3D_emissions_diagnostics
      
! #####################################################################
      SUBROUTINE ASAD_TROPOSPHERIC_MASK(rows,row_length,model_levels,   &
                                        ierr)

      USE UKCA_TROPOPAUSE , ONLY : L_troposphere
      IMPLICIT NONE



      INTEGER, INTENT(IN)  :: row_length     ! length of row
      INTEGER, INTENT(IN)  :: rows           ! number of rows
      INTEGER, INTENT(IN)  :: model_levels   ! number of levels
      INTEGER, INTENT(OUT) :: ierr

      INTEGER :: i,j,k,l

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_TROPOSPHERIC_MASK',zhook_in,        &
                               zhook_handle)
      
      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED ', &
              'asad_tropospheric_mask'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      ierr = 1

      IF (L_asad_use_trop_mask_output) THEN
         DO l=1,n_chemdiags
            SELECT CASE (asad_chemdiags(l)%diag_type)
            CASE (cdtpm)
               asad_chemdiags(l)%throughput(:,:,:) = 1.e0
               WHERE (.NOT. L_Troposphere(:,:,:))
                  asad_chemdiags(l)%throughput(:,:,:) = 0.e0
               ENDWHERE
            END SELECT ! %diag_type
         END DO ! l
      END IF
      
      ierr = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT ', &
              'asad_tropospheric_mask'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      IF (lhook) CALL dr_hook('ASAD_TROPOSPHERIC_MASK',zhook_out,       &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_TROPOSPHERIC_MASK

! #####################################################################
      SUBROUTINE ASAD_TENDENCY_STE(row_length, rows, model_levels,      &
           n_chem_tracers,all_tracers, total_number_density, volume,    &
           mass,timestep, which_diagnostic, L_store_value, ierr)

      USE UKCA_CONSTANTS,   ONLY: avogadro, m_air
      USE UKCA_CSPECIES , ONLY : c_species   
      IMPLICIT NONE


      INTEGER, INTENT(IN)  :: row_length     ! length of row
      INTEGER, INTENT(IN)  :: rows           ! number of rows
      INTEGER, INTENT(IN)  :: model_levels   ! number of levels
      INTEGER, INTENT(IN)  :: n_chem_tracers ! number of chemical tracers

! tracer array
      REAL, INTENT(IN)    :: all_tracers(1-Offx:row_length+Offx,        &
                         1-Offy:rows+Offy,model_levels,n_chem_tracers)
! number density array
      REAL, INTENT(IN)    :: total_number_density(row_length, rows,     &
                                                  model_levels)
! cell volume array
      REAL, INTENT(IN)    :: volume(row_length,rows,model_levels)
! cell mass array
      REAL, INTENT(IN) :: mass(row_length,rows,model_levels)
      REAL, INTENT(IN) :: timestep

      CHARACTER(LEN=3), INTENT(IN) :: which_diagnostic

      LOGICAL, INTENT(IN) :: L_store_value    ! T to store value, F to make difference

      INTEGER, INTENT(OUT) :: ierr            ! error code

      INTEGER :: i,j,k,l,m                    ! counters

      LOGICAL, SAVE :: firstcall_tendency=.TRUE.
      LOGICAL, SAVE :: firstcall_STE=.TRUE.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_TENDENCY_STE',zhook_in,zhook_handle)

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED ','asad_tendency_STE'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      ierr = 1

      ! set up location and molecular mass from c_species array
      IF (firstcall_tendency .OR. firstcall_STE) THEN
         
         DO m=1,n_chem_tracers
            chemdiags_loop_first: DO l=1,n_chemdiags
            IF (which_diagnostic /= asad_chemdiags(l)%diag_type)        &
                    CYCLE chemdiags_loop_first
               SELECT CASE (asad_chemdiags(l)%diag_type)
               CASE (cdnet,cdste)
                IF (advt(m) == asad_chemdiags(l)%species) THEN
                     asad_chemdiags(l)%location = m
                  asad_chemdiags(l)%molecular_mass = c_species(m)*m_air
                  asad_chemdiags(l)%c_vmr_to_mmr = c_species(m)
                  IF (PrintStatus >= Prstatus_diag) WRITE(6,*)     &
                          'mype= ',mype,' TENDENCY/STE: ', &
                          asad_chemdiags(l)%diag_type,' : ', &
                          asad_chemdiags(l)%species,' : ', &
                          asad_chemdiags(l)%location,' : ', &
                          asad_chemdiags(l)%molecular_mass
                  END IF
               END SELECT
            END DO chemdiags_loop_first
         END DO
         
         SELECT CASE (asad_chemdiags(l)%diag_type)
         CASE (cdnet)
            firstcall_tendency=.FALSE.
         CASE (cdste)
            firstcall_STE=.FALSE.
         END SELECT
      END IF         ! firstcall_tendency .OR. firstcall_STE

      ! calculate tendency/STE
      chemdiags_loop: DO l=1,n_chemdiags
        IF (which_diagnostic /= asad_chemdiags(l)%diag_type)            &
              CYCLE chemdiags_loop 
         SELECT CASE (asad_chemdiags(l)%diag_type)
         CASE(cdnet,cdste)
            IF (L_store_value) THEN         ! store value in %throughput
! Volume calculation method, converts to mol
              asad_chemdiags(l)%throughput(:,:,:) =                     &
                    all_tracers(1:row_length,1:rows,1:model_levels,     &
                    asad_chemdiags(l)%location)       &
                    *total_number_density(1:row_length,&
                    1:rows,1:model_levels)*&
                    volume(1:row_length,1:rows,1:model_levels) &
                           /(avogadro*asad_chemdiags(l)%c_vmr_to_mmr)
            ELSE
! Take difference from previously stored value
              asad_chemdiags(l)%throughput(:,:,:) =                     &
                  ((all_tracers(1:row_length,1:rows,1:model_levels,     &
                           asad_chemdiags(l)%location)       &
                           *total_number_density(1:row_length,1:rows,   &
                           1:model_levels)*&
                           volume(1:row_length,1:rows,1:model_levels)   &
                           /(avogadro *asad_chemdiags(l)%c_vmr_to_mmr)) &
                           -asad_chemdiags(l)%throughput(:,:,:))        &
                           /timestep
! Mask out stratospheric values if requested
              IF (asad_chemdiags(l)%tropospheric_mask) THEN
                WHERE (.NOT. L_troposphere)
                  asad_chemdiags(l)%throughput(:,:,:) = 0.0
                ENDWHERE
              END IF
            END IF              ! L_store_value
         END SELECT             ! asad_chemdiags(l)%diag_type
      END DO chemdiags_loop     ! l=1,n_chemdiags
      
      ierr = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT ','asad_tendency_STE'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      IF (lhook) CALL dr_hook('ASAD_TENDENCY_STE',zhook_out,            &
                               zhook_handle)
      RETURN
      END SUBROUTINE ASAD_TENDENCY_STE
      
! #####################################################################
      SUBROUTINE ASAD_MASS_DIAGNOSTIC(row_length, rows, model_levels,   &
                                      mass, ierr)
      IMPLICIT NONE


      INTEGER, INTENT(IN)  :: row_length                           ! array dimension
      INTEGER, INTENT(IN)  :: rows                                 ! array dimension
      INTEGER, INTENT(IN)  :: model_levels                         ! array dimension
      REAL,    INTENT(IN)  :: mass(row_length,rows,model_levels)   ! mass of cell
      INTEGER, INTENT(OUT) :: ierr                                 ! error code

      INTEGER :: i                      ! counter

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_MASS_DIAGNOSTIC',zhook_in,          &
                               zhook_handle)

      ierr = -1


      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED ',                          &
              'asad_mass_diagnostic'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF



      DO i=1,n_chemdiags
        IF (asad_chemdiags(i)%diag_type == cdmas) THEN
          asad_chemdiags(i)%throughput(:,:,:) = mass(:,:,:)
          IF (asad_chemdiags(i)%tropospheric_mask) THEN
            WHERE (.NOT. L_Troposphere)
              asad_chemdiags(i)%throughput(:,:,:) = 0.0
            ENDWHERE
          END IF
        END IF
      END DO                    ! i,n_chemdiags
      
      ierr = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT ','asad_mass_diagnostic'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      IF (lhook) CALL dr_hook('ASAD_MASS_DIAGNOSTIC',zhook_out,         &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_MASS_DIAGNOSTIC

! #####################################################################
      SUBROUTINE ASAD_LIN_DIAGNOSTIC(row_length, rows, model_levels,    &
                                   N_field, ierr)
      IMPLICIT NONE


      INTEGER, INTENT(IN) :: row_length          ! Field dimension
      INTEGER, INTENT(IN) :: rows                !  "      "
      INTEGER, INTENT(IN) :: model_levels        ! No of levels
      REAL, INTENT(IN) :: N_field(1:row_length,1:rows,1:model_levels)
      INTEGER, INTENT(OUT) :: ierr

      integer :: i

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_LIN_DIAGNOSTIC',zhook_in,           &
                              zhook_handle)

      ierr = -1

      IF (PrintStatus >= Prstatus_diag) THEN
        WRITE(6,*) 'mype=',mype,'ENTERED ','asad_LiN_diagnostic'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      DO i=1,n_chemdiags
         SELECT CASE (asad_chemdiags(i)%diag_type)
         CASE (cdlin)
            asad_chemdiags(i)%throughput(:,:,:) = N_field(:,:,:)
            IF (asad_chemdiags(i)%tropospheric_mask) THEN
               WHERE (.NOT. L_Troposphere(:,:,:))
                  asad_chemdiags(i)%throughput(:,:,:) = 0.0
               ENDWHERE
            END IF
         END SELECT             ! %diag_type
      END DO                    ! i,n_chemdiags
      
      ierr = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT ','asad_LiN_diagnostic'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      IF (lhook) CALL dr_hook('ASAD_LIN_DIAGNOSTIC',zhook_out,          &
                               zhook_handle)
      RETURN
      END SUBROUTINE ASAD_LIN_DIAGNOSTIC

! #####################################################################
      SUBROUTINE ASAD_PSC_DIAGNOSTIC(row_length, rows, klevel, ierr)

      IMPLICIT NONE


      INTEGER, INTENT(IN)  :: row_length             ! array dimension
      INTEGER, INTENT(IN)  :: rows                   ! array dimension
      INTEGER, INTENT(IN)  :: klevel                 ! level number
      INTEGER, INTENT(OUT) :: ierr                   ! error code

      INTEGER :: i,j,l,m                             ! counters

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_PSC_DIAGNOSTIC',zhook_in,           &
                              zhook_handle)
      ierr = -1


      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED ','asad_psc_diagnostic'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      DO i=1,n_chemdiags
         SELECT CASE (asad_chemdiags(i)%diag_type)
         CASE (cdpsc)
            SELECT CASE (asad_chemdiags(i)%rxn_type)
            CASE (cdpsc_typ1)
              asad_chemdiags(i)%throughput(:,:,klevel) =                &
                          RESHAPE(fpsc1(:),(/row_length,rows/))
              IF (asad_chemdiags(i)%tropospheric_mask) THEN
                WHERE (.NOT. L_Troposphere(:,:,klevel))
                  asad_chemdiags(i)%throughput(:,:,klevel) = 0.0
                ENDWHERE
              END IF
            CASE (cdpsc_typ2)
              asad_chemdiags(i)%throughput(:,:,klevel) =                &
                          RESHAPE(fpsc2(:),(/row_length,rows/))
              IF (asad_chemdiags(i)%tropospheric_mask) THEN
                WHERE (.NOT. L_Troposphere(:,:,klevel))
                  asad_chemdiags(i)%throughput(:,:,klevel) = 0.0
                ENDWHERE
              END IF
            END SELECT
         END SELECT
      END DO
         
      ierr = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT ','asad_psc_diagnostic'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      IF (lhook) CALL dr_hook('ASAD_PSC_DIAGNOSTIC',zhook_out,          &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_PSC_DIAGNOSTIC

! #####################################################################
      SUBROUTINE ASAD_OUTPUT_TRACER(row_length, rows, model_levels,     &
                                    n_chem_tracers, all_tracers, ierr)

      USE UKCA_CSPECIES,  ONLY : c_species   
      USE UKCA_CONSTANTS, ONLY : avogadro, m_air
      IMPLICIT NONE


      INTEGER, INTENT(IN)  :: row_length     ! length of row
      INTEGER, INTENT(IN)  :: rows           ! number of rows
      INTEGER, INTENT(IN)  :: model_levels   ! number of levels
      INTEGER, INTENT(IN)  :: n_chem_tracers ! number of chemical tracers

! Tracer array
      REAL, INTENT(IN)     :: all_tracers(1-Offx:row_length+Offx,       &
                    1-Offy:rows+Offy,1:model_levels,1:n_chem_tracers)

      INTEGER, INTENT(OUT) :: ierr           ! error code
      INTEGER :: i,j,k,l,m                   ! counters

      LOGICAL, SAVE :: firstcall=.TRUE.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_OUTPUT_TRACER',zhook_in,            &
                              zhook_handle)

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED ','asad_output_tracer'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      ierr = 1

      ! set up location and molecular mass from c_species array
      IF (firstcall) THEN
         
         DO m=1,n_chem_tracers
            DO l=1,n_chemdiags
               SELECT CASE (asad_chemdiags(l)%diag_type)
               CASE (cdout)
                IF (advt(m) == asad_chemdiags(l)%species) THEN
                     asad_chemdiags(l)%location = m
                  asad_chemdiags(l)%molecular_mass = c_species(m)*m_air
                  asad_chemdiags(l)%c_vmr_to_mmr = c_species(m)
                  IF (PrintStatus >= Prstatus_diag) WRITE(6,*)        &
                          'mype= `',mype,' OUTPUT TRACER: ',            &
                          asad_chemdiags(l)%diag_type,' : ',            &
                          asad_chemdiags(l)%species,' : ',              &
                          asad_chemdiags(l)%location,' : ',             &
                          asad_chemdiags(l)%molecular_mass,' : ',       &
                          asad_chemdiags(l)%c_vmr_to_mmr
                  END IF
               END SELECT
            END DO
         END DO

         firstcall=.FALSE.         
      END IF

      ! take copy of tracer, and apply tropospheric mask if requested
      DO l=1,n_chemdiags
         SELECT CASE (asad_chemdiags(l)%diag_type)
         CASE(cdout)
            asad_chemdiags(l)%throughput(:,:,:) =                       &
                 all_tracers(:,:,:,asad_chemdiags(l)%location)
            IF (asad_chemdiags(l)%tropospheric_mask) THEN
               WHERE (.NOT. L_Troposphere(:,:,:))
                  asad_chemdiags(l)%throughput(:,:,:) = 0.0
               ENDWHERE
            END IF
         END SELECT             ! asad_chemdiags(l)%diag_type
      END DO      ! l=1,n_chemdiags
      
      ierr = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT ','asad_output_tracer'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      IF (lhook) CALL dr_hook('ASAD_OUTPUT_TRACER',zhook_out,           &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_OUTPUT_TRACER

! #####################################################################
      SUBROUTINE ASAD_LIGHTNING_DIAGNOSTICS( row_length, rows,          &
           which_diagnostic, diagnostic_field, ierr)

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: row_length                   ! Field dimension
      INTEGER, INTENT(IN) :: rows                         ! Field dimension
      CHARACTER(LEN=1), INTENT(IN) :: which_diagnostic
      REAL, INTENT(IN) :: diagnostic_field(1:row_length,1:rows)
      INTEGER, INTENT(OUT) :: ierr

      INTEGER :: l

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_LIGHTNING_DIAGNOSTICS',zhook_in,    &
                              zhook_handle)

      ierr = -1

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED ',                          &
              'asad_lightning_diagnostics: ', which_diagnostic
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      DO l=1,n_chemdiags
         SELECT CASE (asad_chemdiags(l)%diag_type)
         CASE (cdlgt)
            IF (asad_chemdiags(l)%rxn_type == which_diagnostic) THEN  ! only copy into diagnostic requested
               asad_chemdiags(l)%throughput(:,:,1) = &
                          diagnostic_field(:,:)

! never used
               IF (asad_chemdiags(l)%tropospheric_mask) THEN
                  WHERE (.NOT. L_Troposphere(:,:,1))
                     asad_chemdiags(l)%throughput(:,:,1) = 0.0
                  ENDWHERE
               END IF
            END IF
         END SELECT
      END DO
         
      ierr = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT ', &
              'asad_lightning_diagnostics'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      IF (lhook) CALL dr_hook('ASAD_LIGHTNING_DIAGNOSTICS',zhook_out,   &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_LIGHTNING_DIAGNOSTICS


! #####################################################################
      SUBROUTINE ASAD_FLUX_PUT_STASH(                                   &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
         stashwork_size,STASHwork,                                      &
         row_length,rows,model_levels,do_chemistry)

      USE Submodel_Mod
      IMPLICIT NONE

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

      INTEGER, INTENT(IN) :: stashwork_size
      REAL, DIMENSION(1:stashwork_size), INTENT(INOUT) :: STASHwork
      INTEGER, INTENT(IN) :: row_length       
      INTEGER, INTENT(IN) :: rows 
      INTEGER, INTENT(IN) :: model_levels  
      LOGICAL, INTENT(IN) :: do_chemistry

! Flux arrays
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: upload_array_3D
      REAL, DIMENSION(:,:),   ALLOCATABLE :: upload_array_2D

      INTEGER       :: item                             ! stash item
      INTEGER       :: section                          ! stash section
      INTEGER       :: im_index                         ! atmos model index
      INTEGER       :: i,j,k,l                          ! Counters
      INTEGER       :: icode                            ! Error code
      CHARACTER(LEN=72) :: cmessage                         ! Error return message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ASAD_FLUX_PUT_STASH',zhook_in,           &
                              zhook_handle)

      icode = 0

      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'ENTERED ','asad_flux_put_stash'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF

      im_index = internal_model_index(atmos_im) 

      ALLOCATE(upload_array_2D(1:row_length,1:rows))
      ALLOCATE(upload_array_3D(1:row_length,1:rows,1:model_levels))
      
      ! copy items into STASHwork array here, which is then uploaded with
      ! the call to STASH at the end of UKCA_MAIN.
      DO l=1,n_stashsumdiag 
         item = stash_handling(l)%stash_item 
         section = stash_handling(l)%stash_section
         upload_array_2D(:,:)   = 0.0
         upload_array_3D(:,:,:) = 0.0
         IF (PrintStatus >= Prstatus_diag) THEN 
            WRITE(6,*) 'PUTSTASH:',mype,' section = ',section
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' item = ',item
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' row_length = ',row_length
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' rows = ',rows
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' model_levels = ',&
                 model_levels
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' at_extremity = ',&
                 at_extremity
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' len_stlist = ',len_stlist
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' stash_levels = ',&
                 stash_levels
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' num_stash_levels+1 = ',&
                 num_stash_levels+1
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' atmos_im = ',atmos_im
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' im_index = ',im_index
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' stindex = ',&
                 stindex(1,item,section,im_index)
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            WRITE(6,*) 'PUTSTASH:',mype,' stlist = ',&
                 stlist(1,stindex(1,item,section,im_index))
            ! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
         END IF
! Now check if can output diagnostic, since some of these
!  will only be allowed on chemical timesteps
         IF (stash_handling(l)%len_dim3 == 1) THEN                   ! 2D field
           IF (PrintStatus >= Prstatus_diag) WRITE(6,*)               &
                    'PUTSTASH: CALLING COPYDIAG ',section,item
             DO k=1,stash_handling(l)%number_of_fields
                  ! will only be summing surface fields in this case
               upload_array_2D(:,:) = upload_array_2D(:,:) +            &
                       asad_chemdiags(&
                       stash_handling(l)%chemdiags_location(k)          &
                       )%throughput(:,:,1)
             END DO
             IF (sf(item,section)) THEN 
             IF (PrintStatus >= Prstatus_diag) WRITE(6,*)             &
               'PUTSTASH: CALLING COPYDIAG ',section,item
! DEPENDS ON: copydiag 
                  Call copydiag (stashwork(si(item,section,im_index)),  & 
                       upload_array_2D(:,:),                            & 
                       row_length,rows,0,0,0,0, at_extremity,           & 
                       atmos_im,section,item,icode,cmessage) 
               END IF
            ELSE ! len_dim3 /= 1
              IF (PrintStatus >= Prstatus_diag) WRITE(6,*)            &
                    'PUTSTASH: CALLING COPYDIAG_3D ',section,item
              DO k=1,stash_handling(l)%number_of_fields
                  ! may be summing both 3D and surface fields into the same 3D array
                DO j=1,asad_chemdiags(                                  &
                  stash_handling(l)%chemdiags_location(k))%num_levs
                  upload_array_3D(:,:,j) = upload_array_3D(:,:,j) +     &
                          asad_chemdiags(                               &
                          stash_handling(l)%chemdiags_location(k)       &
                          )%throughput(:,:,j)
                  END DO
               END DO
               IF (sf(item,section)) THEN 
               IF (PrintStatus >= Prstatus_diag) WRITE(6,*)           &
               'PUTSTASH: CALLING COPYDIAG_3D ',section,item
! DEPENDS ON: copydiag_3d 
                  Call copydiag_3d (stashwork(si(                       &
                       item,section,im_index)),                         &
                       upload_array_3D(:,:,:),                          &
                       row_length,rows,model_levels,0,0,0,0,            &
                       at_extremity,                                    &
                       stlist(1,stindex(1,item,section,im_index)),      &
                       len_stlist,                                      & 
                       stash_levels,num_stash_levels+1,                 & 
                       atmos_im,section,item,icode,cmessage) 
               END IF
            END IF ! len_dim3
         IF (icode >  0) THEN 
           WRITE(cmessage,'(A,A,5(I6,1X))') 'ASAD_FLUX_PUT_STASH: ',    &
                 'ERROR WITH copydiag ',                                &
                 item,section,im_index,atmos_im,icode
           WRITE(6,*) cmessage 
! DEPENDS ON: UM_FORT_FLUSH
            CALL UM_FORT_FLUSH(6,ICODE)
            icode=(1000*section) + item
            CALL EREPORT('ASAD_FLUX_PUT_STASH',icode,cmessage) 
         END IF
      END DO       ! l,n_stashsumdiag

! Deallocate arrays, if are able to
      DO i=1,n_chemdiags
         IF (asad_chemdiags(i)%can_deallocate) &
              DEALLOCATE(asad_chemdiags(i)%throughput)
      END DO


      IF (PrintStatus >= Prstatus_diag) THEN
         WRITE(6,*) 'mype=',mype,'LEFT ','asad_flux_put_stash'
! DEPENDS ON: UM_FORT_FLUSH
         CALL UM_FORT_FLUSH(6,ICODE)
      END IF


      IF (lhook) CALL dr_hook('ASAD_FLUX_PUT_STASH',zhook_out,          &
                              zhook_handle)
      RETURN
      END SUBROUTINE ASAD_FLUX_PUT_STASH

  END MODULE ASAD_CHEM_FLUX_DIAGS
