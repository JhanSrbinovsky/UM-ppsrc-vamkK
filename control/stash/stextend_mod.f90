! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stash

MODULE stextend_mod

USE version_mod, ONLY: NPROFDP, NSECTP, NITEMP,  &
                       NRECDP, NTIMEP, NPROFTP, NELEMP,    &
                       NLEVP_S, NPSLISTP

IMPLICIT NONE

! Description:
! Contains variables and arrays involved in STASH processing in the UM.
! NOTE: CSUBMODEL and VERSION must be included before this file

!   Output levels lists
!     List type (real/int)
      CHARACTER(LEN=1) LLISTTY(NPROFDP)
!     Real levels
      REAL        RLEVLST_S(NLEVP_S, NPROFDP)
!     Integer (i.e. model) levels
      INTEGER      LEVLST_S(NLEVP_S, NPROFDP)

!   STASH list array (extra row only for internal indexing in
!                   processing routines)
      INTEGER LIST_S  (NELEMP+1, NRECDP)
!   Output times tables
      INTEGER ITIM_S  (NTIMEP, 2*NPROFTP+2)

! CSUBMAX start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Declarations:
!
!
!  1. Maximum internal model/submodel array sizes for this version.
!
      ! Max no. of internal models
      INTEGER,PARAMETER:: N_INTERNAL_MODEL_MAX=2

      ! Max no. of submodel dump partitions
      INTEGER,PARAMETER:: N_SUBMODEL_PARTITION_MAX=1

      ! Max value of internal model id
      INTEGER,PARAMETER:: INTERNAL_ID_MAX=N_INTERNAL_MODEL_MAX

      ! Max value of submodel dump id
      INTEGER,PARAMETER:: SUBMODEL_ID_MAX=N_SUBMODEL_PARTITION_MAX

! CSUBMAX end

!   STASH lengths and addresses
      INTEGER IN_S    (2,N_INTERNAL_Model_MAX,0:NSECTP,NITEMP)
!   STASH list index
      INTEGER INDX_S  (2,N_INTERNAL_Model_MAX,0:NSECTP,NITEMP)

!   Start addresses for pp headers
      INTEGER PPIND_S (N_INTERNAL_Model_MAX,NITEMP)
!   Time series block information
!     No. of records in a block
      INTEGER NRECS_TS(NPROFDP)
!     Start position of block
      INTEGER NPOS_TS (NPROFDP)
!   lengths of pseudo-levels lists
      INTEGER LENPLST (NPSLISTP)

!     Set up preliminary array for addressing D1:
!     Number of items of info needed for each object and likely maximum
!     number of objects in D1 - this can be increased if necessary

      INTEGER, PARAMETER :: D1_Items_Prel = 5
      INTEGER, PARAMETER :: Max_D1_Len    = 1500

      ! Names of items

      INTEGER, PARAMETER :: d1_type       = 1 ! Prognostic, diagnostic
                                              ! or other
      INTEGER, PARAMETER :: d1_im         = 2 ! Internal model id
      INTEGER, PARAMETER :: d1_extra_info = 3 ! Progs and other :-
                                              ! PPXREF item no
                                              ! Diags :-
                                              ! Stash list item no
      INTEGER, PARAMETER :: d1_levs       = 4 ! No of levels
      INTEGER, PARAMETER :: d1_sect       = 5 ! Section No

      ! Types of items for d1_type

      INTEGER, PARAMETER :: Prog     = 0
      INTEGER, PARAMETER :: Diag     = 1
      INTEGER, PARAMETER :: Seco     = 2
      INTEGER, PARAMETER :: Extra_d1 = 3

      ! Stores number of objects in D1
      INTEGER :: N_Obj_D1(N_SUBMODEL_Partition_MAX)

!     Preliminary array for addressing D1. Holds the minimum amount of
!     info required for order of objects of D1; this amount of info is
!     enough to obtain any other required info from stashlist or ppxref

      INTEGER :: D1_PAddr(D1_Items_Prel,Max_D1_Len,N_SUBMODEL_Partition_MAX)

END MODULE stextend_mod
