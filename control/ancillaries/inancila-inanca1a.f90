! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!   Subroutine INANCILA
!
!   Purpose : Takes as input,the code defining the frequency of update
!             of ancillary fields as set by the user interface.
!             Converts them into a list of numbers of timesteps after
!             which each field must be updated, and calculates the
!             frequency with which this list must be interrogated.
!             Where the update interval is in months or years,
!             the check will be carried out each day. The physical
!             files required are also determined by input code,
!             and the headers and lookup tables are read into
!             the arguments FIXHD,INTHD,LOOKUP which are in
!             COMMON/ANCILHDA/ of calling routine INANCCTL.
!             Indexes for each possible ancillary field are set up in
!             COMMON/IXANCILA/
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Ancillaries

SUBROUTINE inancila(len_fixhd,len_inthd,len_realhd,                            &
                                                     !Intent (In)
                    len1_levdepc,len2_levdepc,                                 &
                    fixhd,inthd,realhd,lookup,                                 &
                    a_fixhd,a_realhd,a_levdepc,                                &
                    ndatasets,nlookups,ftnancil,                               &
                    lookup_start,len1_lookup,row_length,                       &
                    p_rows,u_rows,p_levels,                                    &
                    tr_levels,st_levels,sm_levels,ntiles,                      &
                    ozone_levels,tpps_ozone_levels,title,                      &
                    si_atmos,silen,                                            &
                    ancillary_steps,steps_per_hr,                              &
                    icode,cmessage,lcal360)         ! Intent (Out)


USE nstypes
USE check_iostat_mod
USE ancilcta_namelist_mod   ! All
USE clmchfcg_scenario_mod, ONLY: nsulpat

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IO, ONLY :       &
    file_open,       &
    setpos,          &
    ioLocalReplicate
USE model_file
USE ereport_mod, ONLY : ereport
USE PrintStatus_mod
USE UM_ParVars
USE lookup_addresses
USE science_fixes_mod, ONLY : l_error_ancil_struct
! version_mod items required by model.h
USE version_mod, ONLY: nsectp, outfile_e, outfile_s
USE Submodel_Mod

IMPLICIT NONE

LOGICAL lcal360  ! Logical switch for 360-day calendar
LOGICAL :: Found ! Used in DO loops, replacing GO TO

INTEGER                                                                        &
        len_fixhd,                                                             &
                         ! Length of header blocks in ancillary
!                              ! data sets
        len_inthd,                                                             &
                         !
        len_realhd,                                                            &
                         !
        len1_levdepc,                                                          &
                         ! Dimension of LEVDEPC in model
      len2_levdepc                                                             &
     ,ancillary_steps,                                                         &
        steps_per_hr


INTEGER                                                                        &
        ndatasets,                                                             &
                         ! No of physical files
        nlookups,                                                              &
                         ! No of lookups required(set by User I.)
        iounit,                                                                &
                         ! I/O unit number
        ftnancil(ndatasets),                                                   &
                         ! Fortran nos of physical files
        lookup_start(ndatasets),                                               &
                         ! Start of each individual lookup
                         !        in overall LOOKUP array
        len1_lookup,                                                           &
                         ! Length of PP header
        row_length,                                                            &
                         ! Atmosphere model dimensions
        p_rows,                                                                &
                         ! No. of rows for pressure-type variables
        u_rows,                                                                &
                         ! No. of rows for wind-type variables
        p_levels,                                                              &
                         ! No. of pressure levels
        tr_levels,                                                             &
                         ! No. of tracer levels
        file_levels,                                                           &
                         ! Number of levels of data in files
!                              ! contining multi-level data.
        st_levels,                                                             &
                         ! No. of soil temperature levels
        sm_levels,                                                             &
                         ! No. of soil moisture levels
        ntiles,                                                                &
                         ! No. of surface tiles.
        ozone_levels,                                                          &
                           ! No. of ozone levels
        tpps_ozone_levels                                                      &
                           ! No of ozone levs in TppsOzon dataset

!      For atmos only runs SI_SLAB is a copy of SI_ATMOS
!      SI_SLAB is only used in SLAB runs.

       ,silen                                                                  &
                          ! Length for SI_ATMOS array
       ,si_atmos(silen) 
                          ! ) STASHin addresses of atmos
CHARACTER(LEN=80) title(ndatasets) ! Titles of each dataset

LOGICAL :: l_vert_mismatch    ! T : Vertical levels mismatch

INTEGER                                                                        &
        fixhd(len_fixhd,ndatasets),                                            &
                                   ! Overall Fixed header array
        a_fixhd(len_fixhd),                                                    &
                            ! Fixed header for Dump
        inthd(len_inthd,ndatasets),                                            &
                                   ! Overall Integer header array
        lookup(len1_lookup,nlookups),                                          &
                                     !Overall Lookup array
        icode,                                                                 &
                        ! Return code =0 Normal Exit  >0 Error
        errorstatus                
          

REAL                                                                           &
      realhd(len_realhd,ndatasets),                                            &
                                   !
      a_realhd(len_realhd),                                                    &
                           !
      a_levdepc(len1_levdepc,len2_levdepc),                                    &
      levdepc( (p_levels+1)*4 ) ! Space to hold level dependent
                                ! constants from Anc File

REAL :: coldepc(row_length+1)
REAL :: rowdepc(p_rows+1)

CHARACTER(LEN=80)                                                              &
        cmessage         ! Out error message if I>0
CHARACTER (LEN=*), PARAMETER :: routinename = 'INANCILA'

! Comdecks:----------------------------------------------------------
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
!*L--------------------COMDECK  CANCILA ---------------------------
!
! Purpose : Contains index blocks for control of update of
!           ancillary fields.
!
! -------------------------------------------------------------------
!
! CANCMAXA Store maximum total no of atmosphere/slab ancillary fields.
! -------------------------------------------------------------------
! Type Declarations

      INTEGER NANCIL_FIELDS  ! No of Atmosphere & Slab Ancillary fields
PARAMETER (NANCIL_FIELDS = 196)
! CANCMAXA end
! Type Declarations

      INTEGER                                                           &
     &  FILEANCIL,                                                      &
                         ! File number associated with ancillary fields
     &  NLOOKUP,                                                        &
                         ! Position of ancillary field in lookup tables.
     &  LOOKUP_STEP,                                                    &
                         ! Interval between PP Headers refering to
!                        ! to the same ancillary fields at diferent time
     &  LEVELS,                                                         &
                         ! Number of levels of data in each ancillary
!                        ! field.
     &  STASHANCIL,                                                     &
                         ! Stash codes for ancillary files
     &  D1_ANCILADD      ! Address of ancillary field in main data block


      COMMON/IXANCILA/ FILEANCIL(NANCIL_FIELDS),                        &
     &           NLOOKUP(NANCIL_FIELDS),                                &
     &           LOOKUP_STEP(NANCIL_FIELDS),                            &
     &           LEVELS(NANCIL_FIELDS),                                 &
     &           STASHANCIL(NANCIL_FIELDS),                             &
     &           D1_ANCILADD(NANCIL_FIELDS)

!*L---------- Control data calculated from NAMELIST-------------------
      LOGICAL                                                           &
     &         UPDATE

      INTEGER  FIELDCODE,                                               &
     &         STEPS
!*----------------------------------------------------------------------
      COMMON/CTANCILA/                                                  &
     &         FIELDCODE(2,NANCIL_FIELDS),                              &
     &         STEPS(NANCIL_FIELDS),UPDATE(NANCIL_FIELDS)
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
!*L --------------------- Comdeck: CENVIR   ----------------------------
!
!   Purpose: COMDECK defining Character enviroment variables used
!            by portable IO to open and close files
!
!
! Type declarations
!
      CHARACTER(LEN=8) FT_ENVIRON(199)  ! Array holding enviroment variables
!                                  for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
!
!
!Common Blocks for character and integer arrays
!
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
!


! Comdecks for ancillary files/fields.
! CANCFLDA List of Ancillary Fields - Atmosphere Stash Codes,
! Model Codes and Logical file Numbers
!
! -------------------------------------------------------------------
      INTEGER :: ITEM_CODES_ANCIL(NANCIL_FIELDS)  ! Stash Codes
      INTEGER :: MODEL_CODES_ANCIL(NANCIL_FIELDS) ! Model Codes
      INTEGER :: ANCIL_FILE_NO(NANCIL_FIELDS)     ! Logical file numbers

! -----------------------------------------
! Note 127:154 not used in UM ; set to zero
! -----------------------------------------

      DATA ITEM_CODES_ANCIL(1:100)/                                     &
     &  30,  33,  34,  35,  36,  37,  60,   0,  23,  20,                &
     &  40,  41, 190,  43,  44,   0,  46,  47,   0,   0,                &
     &   0,   0,   0,   0,   0,  26,  31,  24,  32,  28,                &
     &  29,  93, 274, 275,  48,   9,   0,   0,  58,  59,                &
     &  88,  87,  85,  57,  90,  17,  18, 301, 302, 303,                &
     & 304, 305, 306, 307, 308, 309, 310, 311, 312, 313,                &
     & 314, 315, 316, 317, 318, 319, 320, 127, 128, 129,                &
     &   0, 121, 122, 123, 124, 125, 126, 251, 207,   0,                &
     &   0, 160, 216, 217, 218, 213, 219, 220, 223, 321,                &
     & 322, 323, 324, 325, 326, 327, 328, 329, 330, 331/

      DATA ITEM_CODES_ANCIL(101:NANCIL_FIELDS)/                         &
     & 332, 333, 334, 335, 336, 337, 338, 339, 340, 341,                &
     & 505, 418, 419, 420, 421, 422, 423, 424, 425, 426,                &
     & 130, 131, 132, 153, 151, 152,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   5,   6, 351, 352, 353, 354,                &
     & 355, 356, 357, 358, 359, 360, 361, 362, 363, 364,                &
     & 365, 366, 367, 368, 369, 370, 371, 480, 481, 482,                &
     & 483, 484, 485, 486, 487, 134, 135,   7,   0,   0,                &
     &   0,   0,   0, 243, 244, 245/

      DATA MODEL_CODES_ANCIL(1:100) /                                   &
     &   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,                &
     &   1,   1,   1,   0,   1,   0,   1,   1,   0,   0,                &
     &   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   0,   0,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/

      DATA MODEL_CODES_ANCIL(101:NANCIL_FIELDS) /                       &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   0,   1,   1,   1,   1,   1,   1,   1,   0,   0,                &
     &   0,   0,   0,   1,   1,   1/
 
      DATA ANCIL_FILE_NO(1:100) /                                       &
     &   9,  10,  10,  10,  10,  10,   1,   0,   2,   3,                &
     &   4,   4,  11,   4,   4,   0,   4,   4,   0,   0,                &
     &   0,   0,   0,   0,   0,   5,   7,   6,   7,   8,                &
     &   8,   9,  98,  99,   4,   2,   0,   0,  12,  12,                &
     &  13,  13,  13,  14,  14,  10,  10,  15,  15,  15,                &
     &  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,                &
     &  15,  15,  15,  15,  15,  15,  15,  12,  23,  23,                &
     &   0,  17,  18,  18,  18,  18,  12,  24,   4,   0,                &
     &   0,  19,  20,  21,  21,  21,  22,   4,   4,  16,                &
     &  16,  16,  16,  16,  16,  16,  16,  16,  16,  16/

      DATA ANCIL_FILE_NO(101:NANCIL_FIELDS) /                           &
     &  16,  16,  16,  16,  16,  16,  16,  16,  16,  25,                &
     &  26,  27,  27,  27,  27,  27,  27,  27,  27,  27,                &
     &  28,  28,  29,  30,  31,  31,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,  10,  10,  38,  39,  39,  39,                &
     &  40,  40,  41,  41,  42,  42,  42,  43,  43,  43,                &
     &  43,  43,  43,  44,  44,  44,  45,  46,  46,  46,                &
     &  46,  46,  46,  46,  46,  47,  47,  10,   0,   0,                &
     &   0,   0,   0,   5,   5,   5/
! CANCFLDA end

!  Namelist input
!     UPANCA Namelist
INTEGER                                                                        &
   anc_ref_no                                                                  &
                    ! Ancil Ref. No : See comdeck CANCFLDA
  ,period                                                                      &
                    ! Period of Updating Interval (Y/M/D/H)
  ,interval         ! Updating Interval

NAMELIST /upanca/ anc_ref_no,period,interval

! Local Variables

INTEGER                                                                        &
        i,                                                                     &
                         !
        item,                                                                  &
                         !
        j,                                                                     &
                         !
        j1,                                                                    &
                         !
        k,                                                                     &
                         !
        len_io,                                                                &
                         !
        lookups,                                                               &
                         !
        nftin,                                                                 &
                         ! Current FTN number for ancillary data
        start_block                                                            &
                         !
       ,stash_code                                                             &
                         ! Stash item code
       ,nrec_a                                                                 &
                         ! No of atmos records
       ,stash_addr                                                             &
                         ! Stash address
       ,dummy                                                                  &
                         !
       ,n_anc_upd        ! No of ancillaries to be updated
DATA dummy /1/

CHARACTER(LEN=8) cperiod      ! PERIOD in characters.
LOGICAL                                                                        &
        lfile            !

!     SoilDepths : position of soil depths in level dependent constants
INTEGER, PARAMETER :: soildepths = 4

REAL p1,p2
LOGICAL lner

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

lner(p1,p2) = ((ABS(p1-p2))  >   (1.e-6*ABS(p1+p2)))

!  Internal Structure

IF (lhook) CALL dr_hook('INANCILA',zhook_in,zhook_handle)
icode=0
cmessage=' '
iounit=0

!   1.  Initialisation for atmosphere model

DO i=1,nancil_fields
  fileancil(i) =ancil_file_no(i)
  stashancil(i)=item_codes_ancil(i)
END DO
fieldcode(:,:) = 0


! Set default values

l_sstanom=.FALSE.
lamipii=.FALSE.

!   Read in control information from namelist

REWIND 5
READ (UNIT=5, NML=ancilcta, IOSTAT=ErrorStatus) 
CALL check_iostat(errorstatus, "namelist ancilcta")

!     Initialise FIELDCODE from Namelist UPANCA
n_anc_upd = 0
DO i=1,nancil_fields
  READ (5,upanca,err=101,END=101)
  fieldcode(1,anc_ref_no) = period
  fieldcode(2,anc_ref_no) = interval
  n_anc_upd = n_anc_upd+1
END DO

 101  CONTINUE
WRITE (6,*) ' '
WRITE (6,*) n_anc_upd,' Atmos & Slab Ancillaries to be updated.'
DO i=1,nancil_fields
  IF (fieldcode(1,i) >  0) THEN
  IF (fieldcode(1,i) == 1) cperiod=' Years'
  IF (fieldcode(1,i) == 2) cperiod=' Months'
  IF (fieldcode(1,i) == 3) cperiod=' Days'
  IF (fieldcode(1,i) == 4) cperiod=' Hours'
  WRITE (6,*) 'Anc Ref No ',i,' Stash code ',item_codes_ancil(i),              &
  ' Interval ',fieldcode(2,i),cperiod
  END IF
END DO
WRITE (6,*) ' '

! Check that ancillary field has valid address (>1) before proceding
!  to try and update it.  If not, switch off updating via FIELDCODE.
DO i=1,nancil_fields
IF (stashancil(i)  >   0) THEN
  stash_addr = si_atmos(stashancil(i))
ELSE
  stash_addr=0
END IF
  IF (stash_addr  <=  1) THEN
    IF (fieldcode(1,i) >  0) THEN
           WRITE(6,*)' INANCILA: update requested for item ',i,                &
     &     ' STASHcode ',stashancil(i),' but prognostic address not set'
      WRITE(6,*)' FIELDCODE values reset to zeroes'
      fieldcode(1,i) = 0
      fieldcode(2,i) = 0
    END IF
  END IF
END DO

!   1.1 Set number of steps after which each ancillary field is updated
!       Zero is used for fields not to be updated

DO i=1,nancil_fields
  steps(i)=0
  IF (fieldcode(1,i) == 4)THEN
    steps(i)=fieldcode(2,i)*steps_per_hr
  END IF
  IF (fieldcode(1,i) == 3) THEN
    steps(i)=fieldcode(2,i)*24*steps_per_hr
  END IF

  IF (lcal360) THEN
    IF (fieldcode(1,i) == 2) THEN
      steps(i)=fieldcode(2,i)*30*24*steps_per_hr
    END IF
    IF (fieldcode(1,i) == 1) THEN
      steps(i)=fieldcode(2,i)*360*24*steps_per_hr
    END IF
  ELSE
! Gregorian calender:
! If update interval is months or years, test each day. Further testing
! done in REPLANCA.

  IF (fieldcode(1,i) == 1.OR.fieldcode(1,i) == 2)THEN
    steps(i)=24*steps_per_hr
  END IF
END IF

END DO

!   1.2 Set master number of steps ANCILLARY_STEPS at which
!       individual switches are tested.

!   Find first active field

Found = .FALSE.
DO i=1,nancil_fields
  IF (steps(i) >  0) THEN
    ancillary_steps=steps(i)
    Found=.TRUE.
    item = i
    EXIT
  END IF
END DO

! No above fields found

IF (.NOT.Found) THEN
  ancillary_steps=0
  GO TO 9999
END IF

!       Set ANCILLARY_STEPS to lowest common denominater of
!       frequencies for active fields

DO i=item+1,nancil_fields
  IF (steps(i) <  ancillary_steps                                              &
      .AND. steps(i) >  0) THEN
    IF (MOD(ancillary_steps,steps(i)) == 0) THEN
      ancillary_steps=steps(i)
    ELSE
      j1=steps(i)-1
      DO j=j1,1,-1
        IF ((MOD(ancillary_steps,j) == 0).AND.                                 &
            (MOD(steps(i),j)        == 0)) THEN
          ancillary_steps = j
          EXIT
        END IF
      END DO
    END IF
  END IF
END DO

!  1.2.4 Sea surface temperature must be updated when sea ice is update

IF (steps(27) >  0.AND.steps(28) <= 0) THEN
   steps(28)=1
END IF


!  1.3 Set number of headers for each ancillary field

DO i=1,nancil_fields
  levels(i)=1
!   Multilayer hydrology
  IF(i == 36)levels(i)=sm_levels
!   Multilayer aerosols
  IF(i >= 41.AND.i <= 43) levels(i)=tr_levels
!   Multilayer murk concentration and source
  IF(i >= 44.AND.i <= 45) levels(i)=p_levels
!   Multilayer user ancillaries
  IF(i >= 90.AND.i <= 109) levels(i)=p_levels
!   Multi-level ancillaries for sulphur cycle
  IF (i == 72) levels(i) = p_levels
  IF (i == 73) levels(i) = p_levels
  IF (i == 74) levels(i) = p_levels
  IF (i == 75) levels(i) = p_levels
  IF (i == 76) levels(i) = p_levels
  IF (i == 82) levels(i) = nsulpat
  IF (i == 83) levels(i) = ntype
  IF (i == 84) levels(i) = npft
  IF (i == 85) levels(i) = npft
!   Multi-level ancillaries aerosol climatology
  IF (i >= 157.AND.i <= 177) levels(i)=p_levels
  IF (i >= 178.AND.i <= 185) levels(i)=p_levels
  IF (i == 193) levels(i) = ntiles
END DO

levels(7)=ozone_levels
! consider do a check for l_use_tpps_ozone and if set then
! set levels(110)=tpps_ozone_levels
levels(10)=st_levels


!  1.4 Read headers

lookups=0

DO i=1,ndatasets

!  Initialise LOOKUP_START (=0 implies file I not required)
  lookup_start(i)=0

!  Check whether each physical file is needed

  lfile=.FALSE.
  DO j=1,nancil_fields
    IF (fileancil(j) == i.AND.steps(j) >  0) THEN
      lfile=.TRUE.
    END IF
  END DO

  IF(lfile) THEN

    WRITE(6,*) ' '
    WRITE(6,*) ' Ancillary data file ',I,', unit no ',FTNANCIL(I),             &
   &           ', ',TITLE(I)

!   Read headers for physical files required

    nftin=ftnancil(i)

!   1.4.1 Buffer in fixed length header record



    CALL model_file_open(nftin,ft_environ(nftin),                           &
        len_ft_envir(nftin),0,0,icode,ioLocality=ioLocalReplicate)
    IF(icode /= 0)THEN
      cmessage='INANCLA: Error opening file'
      WRITE(6,*) 'INANCILA: Error opening file on unit ',nftin,                &
               ' accessed from env.var.: ',ft_environ(nftin)
      IF (lhook) CALL dr_hook('INANCILA',zhook_out,zhook_handle)
      RETURN
    END IF

    CALL setpos(nftin,0,icode)

!   Read in fixed header to get array dimensions
! DEPENDS ON: read_flh
    CALL read_flh(nftin,fixhd(1,i),len_fixhd,icode,cmessage)
    IF (icode >  0) THEN
      WRITE (6,*) ' Error in reading fixed header for file ',i
      GO TO 9999   !  Return
    END IF

!       Check for negative dimensions
    IF (fixhd(101,i) <= 0) fixhd(101,i)=1
    IF (fixhd(106,i) <= 0) fixhd(106,i)=1
    IF (fixhd(111,i) <= 0) fixhd(111,i)=1
    IF (fixhd(112,i) <= 0) fixhd(112,i)=1
    IF (fixhd(151,i) <= 0) fixhd(151,i)=1
    IF (fixhd(152,i) <= 0) fixhd(152,i)=1
    IF (fixhd(161,i) <= 0) fixhd(161,i)=1

! Set start position of boundary fields for file
  lookup_start(i)=lookups+1

    IF (lookups+fixhd(152,i) >  nlookups) THEN
      WRITE (6,*) 'No room in LOOKUP table for Ancillary File ',i
      cmessage='INANCILA: Insufficient space for LOOKUP headers'
      icode=14
      GO TO 9999   !  Return
    END IF


    CALL setpos(nftin,0,icode)
    IF (icode >  0) THEN
      WRITE (6,*) ' ERROR in SETPOS called from INANCA1A'
      WRITE (6,*) ' SETPOS attempted with Unit No ',nftin
      cmessage = 'INANCA1A : ERROR in SETPOS'
      GO TO 9999    !   Return
    END IF

! DEPENDS ON: readhead
    CALL readhead(nftin,                                                       &
                fixhd(1,i),len_fixhd,                                          &
                inthd(1,i),fixhd(101,i),                                       &
                realhd(1,i),fixhd(106,i),                                      &
                levdepc,fixhd(111,i),fixhd(112,i),                             &
                rowdepc,fixhd(116,i),fixhd(117,i),                             &
                coldepc,fixhd(121,i),fixhd(122,i),                             &
                dummy,dummy,dummy,                                             &
                dummy,dummy,                                                   &
                dummy,dummy,                                                   &
                dummy,dummy,                                                   &
                dummy,dummy,                                                   &
                dummy,dummy,                                                   &
                lookup(1,lookups+1),fixhd(151,i),fixhd(152,i),                 &
                fixhd(161,i),                                                  &
                start_block,icode,cmessage)

    IF (icode >  0) THEN
       WRITE(6,*) 'ERROR in READHEAD for Ancillary File ',i
       WRITE(6,*) 'Unit Number ',nftin
       GO TO 9999   !   Return
    END IF

!   Check calendar indicator is correct if this is set in the ancillary
    IF (fixhd(8,i) /=  imdi) THEN
      IF ((     lcal360 .AND. fixhd(8,i) /= 2) .OR.                            &
           (.NOT.lcal360 .AND. fixhd(8,i) /= 1) ) THEN
        icode=100+i
        cmessage='INANCILA : Wrong calendar set in Ancillary File'
        WRITE (6,'(a)') ' ******** Error in INANCILA ********'
        WRITE (6,'(a,i4)') ' Wrong calendar setting in Ancillary File ', i
        IF (lcal360) THEN
          WRITE (6,'(a)') ' Model run is set up for 360 day calendar.'
          WRITE (6,'(a)') ' Ancillary File is for 365 day calendar.'
        ELSE
          WRITE (6,'(a)') ' Model run is set up for 365 day calendar.'
          WRITE (6,'(a)') ' Ancillary File is for 360 day calendar.'
        END IF
        WRITE (6,'(a)') ' Rerun with correct ancillary file.'
        GO TO 9999   !  Return
      END IF
    ELSE
      WRITE (6,'(a,i4)') 'Unspecified calendar type in ancillary file ', i
    END IF
      

    file_levels=1

    IF(i == 1) THEN
      file_levels=ozone_levels
    ELSE IF(i == 2) THEN
      file_levels=sm_levels
! This is the maximum value that might be present on the ancillary
! file if it includes soil moisture in layers; otherwise only single
! level data is present and PR_FIXHD will not check value since
! FIXHD(110) will be zero
    ELSE IF(i == 3) THEN
        file_levels=st_levels
    ELSE IF(i == 13) THEN   ! for multilevel aerosols
        file_levels=tr_levels
    ELSE IF(i == 14.OR.i == 16) THEN   ! for murk and user ancil.
        file_levels=p_levels
    ELSE IF(i == 17.OR.i == 18) THEN
!           multi-level sulphur cycle ancillary files.
        file_levels=p_levels
    ELSE IF (i == 25) THEN
! tropopause-based ozone file with tpps_ozone_levels
       file_levels=tpps_ozone_levels
    END IF


!  1.4.2 Buffer in integer constants

    IF(fixhd(100,i) >  0) THEN

! Check for error in file pointers

! Check validity of integer data and print out information
! All files except ozone should contain full fields

      IF(inthd(6,i) /= row_length) THEN
! Ozone may contain zonal mean data
! also applies to tropopause-based ozone -- no 25.
        IF(.NOT.((i == 1) .OR. (i  ==  25))                                    &
         .OR. inthd(6,i) /= 1) THEN
          icode=4
          cmessage='INANCILA:integer header error'
          WRITE(6,*) ' INTHD(6) : ',inthd(6,i),' ?'
          IF (lhook) CALL dr_hook('INANCILA',zhook_out,zhook_handle)
          RETURN
        END IF
      END IF

      IF(inthd(7,i) /= p_rows.AND.(i == 9.AND.inthd                            &
        (7,i) /= u_rows)) THEN
        icode=5
        cmessage='INANCILA:integer header error'
        WRITE(6,*) ' INTHD(7) : ',inthd(7,i),' ?'
        IF (lhook) CALL dr_hook('INANCILA',zhook_out,zhook_handle)
        RETURN
      END IF

      IF (i == 1 .OR. i  ==  25) THEN  ! Ozone or tpps-ozone file
        WRITE (6,*) ' '
        IF (inthd(6,i) == 1)THEN
          WRITE (6,*) ' OZONE file contains zonal mean data for ',             &
          inthd(6,i),' points x ',inthd(7,i),' rows'
        ELSE IF (inthd(6,i) == row_length)THEN
          WRITE (6,*) ' OZONE file contains full fields for ',                 &
          inthd(6,i),' points x ',inthd(7,i),' rows'
        END IF
! Check that correct ozone file has been provided.
        IF ((zonavozone .AND. i  ==  1) .OR.                                   &
          (zonavtppsozone .AND. i  ==  25)) THEN
!! Where is ZonAvOzone (and ZonAvTppsOzone) defined.
          IF (inthd(6,i) /= 1) THEN
            WRITE (6,*) ' Zonal Ozone Data is expected',                       &
            ' for 1 point x ',p_rows,' rows'
            icode = 51
            cmessage = 'INANCA1A : Wrong Ozone data provided.'
            GO TO 9999   !  Return
          END IF
        ELSE
          IF (inthd(6,i) /= row_length) THEN
            WRITE (6,*) ' Ozone Data is expected for ',                        &
            row_length,' points x ',p_rows,' rows.'
            icode = 52
            cmessage = 'INANCA1A : Wrong Ozone data provided.'
            GO TO 9999   !  Return
          END IF
        END IF
      END IF

    END IF

!  1.4.3 Buffer in real constants
    
    IF(fixhd(105,i) >  0) THEN

! Check validity of real header and print out information

! Only perform this check if ancillary and model is on C grid (with P at poles)
     IF (fixhd(9,i) == 3 .AND. a_fixhd(9) == 3 ) THEN
       DO j=1,6
         IF(realhd(j,i) >  (a_realhd(j)+0.1).OR.                               &
           realhd(j,i) <  (a_realhd(j)-0.1))THEN
         IF(i /= 1.OR.(j /= 1.AND.j /= 4))THEN
           WRITE(6,*)(realhd(k,i),k=1,6),(a_realhd(k),k=1,6)
           icode=8
           cmessage='INANCILA: REAL header Error.'
           IF (lhook) CALL dr_hook('INANCILA',zhook_out,zhook_handle)
           RETURN
         END IF
         END IF
       END DO
     END IF

   END IF

!  1.4.4 Buffer in level dependent constants if required
!        Not retained in model after initial check

   IF(fixhd(110,i) >  0) THEN

!  Only files 1 (Ozone), and 3 (Soil temperature)should contain multi
!  level data. File 2 (Soil moisture,snow depth,fractional snow time
!  and soil moisture in layers) may possibly also have multi level data.
!  FILES 13,14,16 (aerosols, murkiness, user ancil.) may also have
!   multi level data.
!! Files 25 (TppsOzon) may also contain multi-level data
!! File 46 has multilevel data for cariolle ozone scheme

     IF (i == 1   .OR.                                                         &
                          !  Ozone File
         i == 14  .OR.                                                         &
                          !  Murkiness File
         i == 16  .OR.                                                         &
                          !  User Ancillary File
         i == 46) THEN    !  cariolle ozone File

! Check that ancillary file is set up for correct vertical levels

      IF (fixhd(111,i)-1 /= p_levels) THEN
        icode=110
        WRITE (cmessage,*) ' Ancillary File set up for wrong',                 &
        ' no of model levels. Anc ',fixhd(111,i)-1,                            &
        ' Model ',p_levels

        CALL ereport ( routinename, icode, cmessage )
      END IF

      l_vert_mismatch = .FALSE.

! Check eta_theta and eta_rho

      DO j=1,p_levels+1
        IF (lner( levdepc(j), a_levdepc(j,1) )) THEN
          l_vert_mismatch = .TRUE.
          EXIT
        END IF
      END DO

      DO j=1,p_levels
        IF (lner( levdepc(fixhd(111,i)+j), a_levdepc(j,2) )) THEN
          l_vert_mismatch = .TRUE.
          EXIT
        END IF
      END DO

! Abort if there is a mis-match

      IF (l_vert_mismatch) THEN
        WRITE (6,*) 'Mismatch in vertical levels between model ',              &
                    'and Ancillary File.'
        WRITE (6,*) 'Anc File : ',title(i)
        WRITE (6,*) 'Eta_Theta - Model'
        WRITE (6,'(5F10.7)') (a_levdepc(k,1),k=1,p_levels+1)
        WRITE (6,*) 'Eta_Theta - Anc File'
        WRITE (6,'(5F10.7)') (levdepc(k),k=1,p_levels+1)
        WRITE (6,*) 'Eta_Rho   - Model'
        WRITE (6,'(5F10.7)') (a_levdepc(k,2),k=1,p_levels)
        WRITE (6,*) 'Eta_Rho   - Anc File'
        WRITE (6,'(5F10.7)') (levdepc(p_levels+1+k),k=1,p_levels)
             icode=11
        WRITE (cmessage,*) 'Mismatch in LEVDEPC ',                             &
        'between model and Ancillary File.'

        CALL ereport ( routinename, icode, cmessage )
      END IF

     ELSE IF (i  ==  25) THEN !! tropopause-based ozone
       !! no checks to run....

     ELSE IF (i == 2) THEN  !  Soil Moisture File

       IF (printstatus >= prstatus_diag .AND. mype == 0 )THEN
         WRITE (6,*)
         WRITE (6,*) 'SoilDepths = ',soildepths
         WRITE (6,*) 'SM_Levels  = ',sm_levels
         DO j=1,sm_levels
           WRITE (6,*) 'model ',a_levdepc(j,soildepths),                       &
                       ' anc ',levdepc(fixhd(111,i)*3+j)
         END DO
       END IF

! Check Soil moisture levels

       DO j=1,sm_levels
         IF (lner(levdepc(fixhd(111,i)*3+j),                                   &
                  a_levdepc(j,soildepths))) THEN
           icode=12
           cmessage='INANCILA: error in LEVDEPC.'
           IF (lhook) CALL dr_hook('INANCILA',zhook_out,zhook_handle)
           RETURN
         END IF
       END DO

     ELSE IF (i == 3) THEN  !  Deep Soil Temperature File

       IF (printstatus >= prstatus_diag .AND. mype == 0) THEN
         WRITE (6,*)
         WRITE (6,*) 'SoilDepths = ',soildepths
         WRITE (6,*) 'st_levels  = ',st_levels
         DO j=1,st_levels
           WRITE (6,*) 'model ',a_levdepc(j,soildepths),                       &
                       ' anc ',levdepc(fixhd(111,i)*3+j)
         END DO
       END IF

! Check Deep Soil levels

       DO j=1,st_levels
         IF (lner(levdepc(fixhd(111,i)*3+j),                                   &
                  a_levdepc(j,soildepths))) THEN
           icode=13
           cmessage='INANCILA: error in LEVDEPC.'
           IF (lhook) CALL dr_hook('INANCILA',zhook_out,zhook_handle)
           RETURN
         END IF
       END DO

!  If aerosol file, check against model levels

     ELSE IF (i == 13) THEN

       DO j=1,tr_levels
         DO j1=1,4
           IF(lner(levdepc(j+(j1-1)*fixhd(111,i)),a_levdepc                    &
                   (j,j1))) THEN
             WRITE(6,*)'Error in level dependent constants:Level=',j
             WRITE(6,*)'Position=',j1
             WRITE(6,*)'Value in model =',a_levdepc(j,j1)
             WRITE(6,*)'Value in ancillary data =',levdepc(j+                  &
                             (j1-1)*fixhd(111,i))
             icode=16
             cmessage='INANCILA: error in LEVDEPC.'
             IF (lhook) CALL dr_hook('INANCILA',zhook_out,zhook_handle)
             RETURN
           END IF
         END DO
       END DO

     END IF  !  If I

   END IF  !  If Fixhd(110,I) > 0

!  1.4.5 Buffer in lookup table
! Set start position of boundary fields for file

   IF(fixhd(150,i) >  0) THEN

     nrec_a = 0
     DO j = 1,fixhd(152,i)
       IF (lookup(model_code,lookups+j)  ==  0 .OR.                            &
           lookup(model_code,lookups+j)  ==  imdi) THEN
         stash_code = lookup(item_code,lookups+j)
         lookup(model_code,lookups+j) = atmos_im
         nrec_a = nrec_a+1
       END IF
     END DO
     IF (nrec_a >  0) THEN
       WRITE (6,*) ' '
       WRITE (6,*) ' INANCA1A : submodel_id in ',nrec_a,                       &
       ' records set to atmos_im in ancillary file ',i
     END IF

   END IF

   lookups=lookups+fixhd(152,i)

 ELSE

!   If file not required, zero fixed length header
   DO j=1,len_fixhd
fixhd(j,i)=0
   END DO

   lookup_start(i)=lookups+1
 END IF

END DO

!  1.5 Set positions in main data blocks


DO i=1,nancil_fields
  IF (stashancil(i)  >   0) THEN
    d1_anciladd(i)=si_atmos(stashancil(i))
  ELSE
    d1_anciladd(i)=0
  END IF
END DO

!  1.51 If a request is made to update a field, ensure that space for
!      that field has been allocted in D1.

DO i=1,nancil_fields
  IF((fieldcode(1,i) >  0).AND.(d1_anciladd(i) <= 1)) THEN
    WRITE(6,*)                                                                 &
    ' An address in D1 has not been set for ancillary field number ',i
    icode=30
    WRITE(cmessage,*) 'INANCILA: updating for ancillary field is requested',   &
                      ' but no space has been allocated in D1'
    IF (lhook) CALL dr_hook('INANCILA',zhook_out,zhook_handle)
    RETURN
  END IF
END DO

!  1.6 Set positions of data

DO i=1,nancil_fields
nlookup(i) =0
lookup_step(i)=0

! If LOOKUP_START=0 for file FILEANCIL(I), no fields required.
  IF   (stashancil(i)  >   0) THEN
  IF (lookup_start(fileancil(i)) >  0) THEN

  DO j=lookup_start(fileancil(i)),lookups

    IF (lookup(item_code,j) == stashancil(i)) THEN
      nlookup(i)=j-lookup_start(fileancil(i))+1
      lookup_step(i)=0
      EXIT
    END IF

  END DO

! Find second occurence of data to set LOOKUP_STEP

  IF(j <  lookups) THEN

    DO j1=j+levels(i),lookups
      IF (lookup(item_code,j1) == stashancil(i)) THEN
        lookup_step(i)=j1-nlookup(i)-lookup_start(fileancil(i))+1
        EXIT
      END IF
    END DO

  END IF

! Check Ancillary files are consistent since we assume:
! for month in month_list:
!   for STASH in STASH_list:
!     for level in level_psuedolevel_list:
! First check if its periodic since we can only check these at
! this point and whether it has multiple fields to step through.
  IF (fixhd(10,fileancil(i)) == 2 .AND. lookup_step(i) > 0) THEN
    ! Loop over the fields which should only differ by time.
    DO j = lookup_start(fileancil(i)) + nlookup(i)-1,               & 
           lookup_start(fileancil(i)) + fixhd(152,fileancil(i)) - 1, &
           lookup_step(i)
      ! If the STASH is different we have a problem.
      IF (lookup(item_code,j) /= stashancil(i) ) THEN
        IF (l_error_ancil_struct) THEN
          icode = 17
        ELSE
          icode = -17
        END IF
        WRITE (cmessage,'(A,I4)')                                  &
          'Incorrect structure for ancillary file ', fileancil(i)
        CALL ereport ( routinename, icode, cmessage )
        EXIT
      END IF
    END DO
  END IF


  END IF
  END IF

END DO

!  SET LEVELS=2 FOR ICE FRACTION AND SNOW DEPTH, TO INDICATE PRESCENCE
!  fractional time fields

levels(27)=2
levels(9)=2


9999 CONTINUE
IF (lhook) CALL dr_hook('INANCILA',zhook_out,zhook_handle)
RETURN
END SUBROUTINE inancila

