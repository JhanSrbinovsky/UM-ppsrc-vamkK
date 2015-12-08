! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Data module defining STSHCOMP namelist

MODULE Rcf_Model_Mod

! Description:
!   Data module to define the STSHCOMP namelist and related addressing
!   variables.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.

USE Submodel_Mod, ONLY :      &
    N_Submodel_partition_max, &
    N_Internal_Model_Max

USE version_mod, ONLY :       &
    NsectP,                   &
    NitemP,                   &
    outfile_E,                &
    outfile_S

IMPLICIT NONE

INTEGER, PARAMETER  ::   A_Max_UKCAVars=150 ! Max.no.of tracers allowed
INTEGER, PARAMETER  ::   A_Max_TrVars = 150 ! Max.no.of tracers allowed

INTEGER, PARAMETER  ::   AASSETS    = 9
INTEGER, PARAMETER  ::   MEAD_TYPES = 4

REAL, SAVE          ::   H_A_EWSPACE ,H_A_NSSPACE
REAL, SAVE          ::   H_A_FIRSTLAT,H_A_FIRSTLONG
REAL, SAVE          ::   H_A_POLELAT ,H_A_POLELONG

INTEGER, SAVE       ::   H_A_GROUP
INTEGER, SAVE       ::   H_OROG_ROUGH
INTEGER, SAVE       ::   A_ASSMGRPS
INTEGER, SAVE       ::   NUM_PVPR

LOGICAL, SAVE       ::   A_RECON
LOGICAL, SAVE       ::   H_OROG_GRAD
LOGICAL, SAVE       ::   ATMODS
LOGICAL, SAVE       ::   CMODS
LOGICAL, SAVE       ::   LMESO

LOGICAL, SAVE       ::   TRACER_A    (0:A_Max_TrVars)
LOGICAL, SAVE       ::   TR_UKCA_A    (0:A_Max_UKCAVars)
LOGICAL, SAVE       ::   AASSET   (AASSETS)
INTEGER, PARAMETER  ::   MAX_AOBS  = 100
INTEGER, SAVE       ::   AOBINC   (MAX_AOBS)
INTEGER, SAVE       ::   AOBGRP   (MAX_AOBS)
INTEGER, SAVE       ::   AASPF    (AASSETS)
INTEGER, SAVE       ::   AASPL    (AASSETS)
INTEGER, SAVE       ::   RUN_TARGET_END( 6)

!Total data length for primary fields for each submodel data partition
INTEGER, SAVE       ::   LPRIM(N_SUBMODEL_PARTITION_MAX)
!Total data length for primary fields for each internal model
INTEGER, SAVE       ::   LPrimIM(N_INTERNAL_MODEL_MAX)
!Total data length for diagnostic flds for each submodel data partition
INTEGER, SAVE       ::   LDUMP(N_SUBMODEL_PARTITION_MAX)
!Total data length for diagnostic flds for each internal model
INTEGER, SAVE       ::   LDumpIM(N_INTERNAL_MODEL_MAX)
!Total data length for secondary flds for each submodel data partition
INTEGER, SAVE       ::   LSECD(N_SUBMODEL_PARTITION_MAX)
!Total data length for secondary flds for each internal model
INTEGER, SAVE       ::   LSecdIM(N_INTERNAL_MODEL_MAX)
!Total workspace length for each submodel data partition
INTEGER, SAVE       ::   LWORK(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each submodel data partition
INTEGER, SAVE       ::   NHeadSub(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each internal model
INTEGER, SAVE       ::   NHEAD(N_INTERNAL_MODEL_MAX)
!Total length of extra space for each submod. data part.
INTEGER, SAVE       ::   LEXTRA(N_SUBMODEL_PARTITION_MAX)
!Data length for dual-time level ocean fields
INTEGER, SAVE       ::   LPRIM_O2
INTEGER, SAVE       ::   ITEM_MAX_REQ
INTEGER, SAVE       ::   ITEM_MAX_ALL

INTEGER, SAVE       ::   NRECS_S
INTEGER, SAVE       ::   NTIMES_S
INTEGER, SAVE       ::   NSERBLK_S
INTEGER, SAVE       ::   NSERREC_S
INTEGER, SAVE       ::   NLEVL_S
INTEGER, SAVE       ::   NMAXLEV_S
INTEGER, SAVE       ::   NPSLISTS_S
INTEGER, SAVE       ::   NMAXPSL_S
INTEGER, SAVE       ::   NHEAD_FILE(OUTFILE_S:OUTFILE_E)
LOGICAL, SAVE       ::   LSTUSER
INTEGER, SAVE       ::   ReconItems(NITEMP)

CHARACTER (LEN = 1), SAVE  ::  H_ATMOS
CHARACTER (LEN = 1), SAVE  ::  H_FLOOR
CHARACTER (LEN = 1), SAVE  ::  H_STRAT
CHARACTER (LEN = 1), SAVE  ::  H_GLOBAL(N_INTERNAL_MODEL_MAX)
INTEGER, SAVE              ::  H_VERS  &
                               (N_INTERNAL_MODEL_MAX,0:NSECTP)

! These are set in SETMODL:
INTEGER, SAVE       ::   MEAN_NUMBER(N_INTERNAL_MODEL_MAX)

REAL, SAVE          ::   H_W_EWSPACE ,H_W_NSSPACE
REAL, SAVE          ::   H_W_FIRSTLAT,H_W_FIRSTLONG

! Variables read in by namelist and used in SETMODL
INTEGER, SAVE       ::   OCAAA
INTEGER, SAVE       ::   NWTRAIN
INTEGER, SAVE       ::   LWBND
INTEGER, SAVE       ::   OCALB
INTEGER, SAVE       ::   SWBND
INTEGER, SAVE       ::   TCA(A_Max_TrVars)
INTEGER, SAVE       ::   TC_UKCA(A_Max_UKCAVars)
INTEGER, SAVE       ::   TCA_LBC(A_Max_TrVars)
INTEGER, SAVE       ::   TC_LBC_UKCA(A_Max_UKCAVars)
INTEGER, SAVE       ::   StLevGWdrag
INTEGER, SAVE       ::   BotVDiffLev
INTEGER, SAVE       ::   TopVDiffLev

REAL, SAVE          ::   EWSPACEA, NSSPACEA
REAL, SAVE          ::   FRSTLATA, FRSTLONA
REAL, SAVE          ::   POLELATA, POLELONA
REAL, SAVE          ::   LATS, LONS

LOGICAL, SAVE       ::   ZonAvOzone
LOGICAL, SAVE       ::   ZonAvTppsOzone

CHARACTER (LEN = 1), SAVE  ::  lfloor
CHARACTER (LEN = 1), SAVE  ::  OROGR
CHARACTER (LEN = 1), SAVE  ::  SWMCR
CHARACTER (LEN = 1), SAVE  ::  MESO

! Copies of namelist variables with more user-friendly names
REAL, SAVE                 :: delta_lon !  = ewspacea
REAL, SAVE                 :: delta_lat !  = nsspacea

NAMELIST/STSHCOMP/                                               &
  RUN_TARGET_END,                                                &
  OCAAA,     EWSPACEA,    POLELATA, FRSTLATA,  LATS,             &
             NSSPACEA,    POLELONA, FRSTLONA,  LONS,             &
  SWBND,     LWBND,                            OROGR,            &
  ZONAVOZONE,             SWMCR,    MESO,                        &
  OCALB,     LFLOOR,      AOBINC,                                &
  TCA,       TCA_LBC,     TC_UKCA,  TC_LBC_UKCA, AOBGRP      
                                         


END MODULE Rcf_Model_Mod
