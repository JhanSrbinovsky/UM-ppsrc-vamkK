! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE ADD_ENG_CORR
!  
!    PURPOSE : PART OF ENERGY CORRECTION SUITE OF ROUTINES
!              - TO ADD IN TEMPERATURE CORRECTION TO
!                GLOBAL TEMPERATURE FIELD SO TO
!                CONSERVE TOTAL ENERGY GLOBALLY
!  
!    NOT SUITABLE FOR SINGLE COLUMN MODEL USE
!  
!    PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!  
!    DOCUMENTATION :
!  
!----------------------------------------------------------------------
!
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Energy Correction
      MODULE add_eng_corr_mod
      IMPLICIT NONE

      CONTAINS

      SUBROUTINE add_eng_corr (energy_corr,t,                           &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                               STASHwork14)

      USE yomhook,               ONLY: lhook, dr_hook
      USE parkind1,              ONLY: jprb, jpim
      USE ereport_mod,           ONLY: ereport
      USE proc_info_mod,         ONLY: at_extremity
      USE Submodel_Mod
      USE atm_fields_bounds_mod, ONLY: tdims
      USE timestep_mod,          ONLY: tstep => timestep

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

!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      REAL :: ENERGY_CORR         ! IN ENERGY CORRECTION
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE IN AND OUT
!----------------------------------------------------------------------
!
! INOUT sum of temperature increments
      REAL :: T(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                            1:tdims%k_end)
!  diagnostics  out
      REAL :: STASHwork14(*)   ! STASH workspace
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
      REAL ::                                                           &
        work(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
                         1:tdims%k_end)

      INTEGER :: i,j,k,                                                 &
                                   ! LOOP COUNTERs
        icode,                                                          &
                         ! return code
        im_index         ! model index

      CHARACTER(LEN=80) ::                                                   &
        cmessage        ! return message
      CHARACTER(LEN=*), PARAMETER :: routinename ='Add_Eng_corr'


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!----------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('ADD_ENG_CORR',zhook_in,zhook_handle)
      icode=0      ! initialise return code to zero

!----------------------------------------------------------------------
! CORRECT TEMPERATURE FOR ERROR IN ENERGY BUDGET OF THE
! PREVIOUS DAY
!----------------------------------------------------------------------
!
      DO k=1,tdims%k_end
        DO j=tdims%j_start,tdims%j_end
          DO i=tdims%i_start,tdims%i_end
            t(i,j,k) = t(i,j,k) + energy_corr*tstep
          END DO
        END DO
      END DO

!-----------------------------------------------------------------------
! Stash diagnostics
!-----------------------------------------------------------------------
! 14 181  T increment on model levels

      im_index = internal_model_index(atmos_im)

      IF (sf(181,14)) THEN
        DO k=1,tdims%k_end
          DO j=tdims%j_start,tdims%j_end
            DO i=tdims%i_start,tdims%i_end

              work(I,J,K) = ENERGY_CORR*TSTEP

            END DO
          END DO
        END DO

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork14(si(181,14,im_index)),             &
             work,tdims%j_end-tdims%j_start+1,                          &
             tdims%j_end-tdims%j_start+1,                               &
             tdims%k_end,                                               &
             0,0,0,0, at_extremity,                                     &
             stlist(1,stindex(1,181,14,im_index)),len_stlist,           &
             stash_levels,num_stash_levels+1,atmos_im,14,181,           &
             icode,cmessage)

        IF (icode >  0) THEN
          cmessage=":error in copydiag_3d(item 181)"//cmessage
        END IF

      END IF

! Note old 14201 would be expensive to calculate here as it is
! now defined as cv*dt*(column integral of dry mass) and is therefore
! calculated and output from section 30.

      IF (icode /= 0) THEN
         CALL ereport(RoutineName,icode,cmessage)
      END IF
!----------------------------------------------------------------------
! ADD ENERGY CORRECTION INTO SUM OF DIABATIC FLUXES
!----------------------------------------------------------------------
! From UM 5.1 Moved elsewhere in code

      IF (lhook) CALL dr_hook('ADD_ENG_CORR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ADD_ENG_CORR

      END MODULE
