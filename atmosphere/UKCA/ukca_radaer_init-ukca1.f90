! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!+ Gather D1 information about UKCA-MODE inputs required by UKCA_RADAER
!
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
SUBROUTINE ukca_radaer_init(                                      &
                            ierr                                  &
 ,                          cmessage,                             &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                            ukca_radaer                           &
)

USE parkind1, ONLY: jpim, jprb 
USE yomhook,  ONLY: lhook, dr_hook
USE ukca_mode_setup
USE ukca_radaer_precalc
USE ukca_radaer_struct_mod
USE ukca_option_mod, ONLY: i_mode_setup
USE Submodel_Mod
USE Control_Max_Sizes

IMPLICIT NONE

!
! Arguments
!
!
! Error indicator (0 is OK, >0 error)
!
INTEGER ierr

!
! Error message if ierr is larger than 0
!
CHARACTER (len=256) :: cmessage

!
! Structure for UKCA/radiation interaction
!
TYPE (ukca_radaer_struct) :: ukca_radaer
!
! STASH flags and other related variables
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
!
! Local variables
!
! Loop variables
!
INTEGER i, j

!     
! In-loop copy of mode names
!
CHARACTER(len=7) :: this_name

!
! In-loop mode type
!
INTEGER this_type

!
! Local counters for number of modes and components
!
INTEGER n_loc_mode
INTEGER n_loc_cpnt

!
! STASH codes for the diagnostics given the modal dry and wet
! diameters, the component mass-mixing ratios, and the modal
! densities.
! Those arrays should be obtained through UKCA_MODE_SETUP for
! maximum flexibility.
!
! The following prognostics are expected in section 34.
!
INTEGER, DIMENSION(nmodes), PARAMETER :: stashc_nbr = (/ &
           ! nuc_sol ait_sol acc_sol cor_sol ait_ins acc_ins cor_ins
               101,     103,    107,    113,    119,    122,     124 &
           /)

INTEGER, DIMENSION(nmodes, ncp), PARAMETER :: stashc_mmr = &
  reshape( source = (/ &
           ! nuc_sol ait_sol acc_sol cor_sol ait_ins acc_ins cor_ins
    ! h2so4
                102,    104,    108,    114,     -1,     -1,     -1, &
    ! bcarbon       
                 -1,    105,    109,    115,    120,     -1,     -1, &
    ! ocarbon
                126,    106,    110,    116,    121,     -1,     -1, &
    ! nacl
                 -1,    127,    111,    117,     -1,     -1,     -1, &
    ! dust
                 -1,     -1,    112,    118,     -1,    123,    125, &
    ! sec_org
                128,    129,    130,    131,     -1,     -1,     -1  &
                    /), & 
           shape = (/ nmodes, ncp /) )
!
! The following diagnostics are expected in section 38.
!
INTEGER, DIMENSION(nmodes), PARAMETER :: stashc_dry = (/ &
           ! nuc_sol ait_sol acc_sol cor_sol ait_ins acc_ins cor_ins
                401,    402,    403,    404,    405,   406,      407 &
           /)
INTEGER, DIMENSION(nmodes), PARAMETER :: stashc_wet = (/ &
           ! nuc_sol ait_sol acc_sol cor_sol ait_ins acc_ins cor_ins
                408,    409,    410,    411,     -1,     -1,      -1 &
           /)
INTEGER, DIMENSION(nmodes), PARAMETER :: stashc_rho = (/ &
           ! nuc_sol ait_sol acc_sol cor_sol ait_ins acc_ins cor_ins
                430,    431,    432,    433,    434,    435,     436 &
           /)

!
! Diagnostics that will have to be taken in also include:
! 
!
! Component volumes (including water as a component)
!
INTEGER, DIMENSION(nmodes, ncp+1), PARAMETER :: stashc_cvl = &
  reshape( source = (/ &
           ! nuc_sol ait_sol acc_sol cor_sol ait_ins acc_ins cor_ins
    ! h2so4
                442,    446,    451,    458,     -1,     -1,     -1, &
    ! bcarbon       
                 -1,    447,    452,    459,    465,     -1,     -1, &
    ! ocarbon
                443,    448,    453,    460,    466,     -1,     -1, &
    ! nacl
                 -1,     -1,    454,    461,     -1,     -1,     -1, &
    ! dust
                 -1,     -1,    455,    462,     -1,    467,    468, &
    ! sec_org
                444,    449,    456,    463,     -1,     -1,     -1, &
    ! water
                445,    450,    457,    464,     -1,     -1,     -1  &
                    /), & 
           shape = (/ nmodes, ncp+1 /) )

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UKCA_RADAER_INIT', zhook_in, zhook_handle)

ierr = 0

!
! Check that fixed array dimensions are the same in ukca_setup and
! ukca_radaer.
!
IF (nmodes /= npd_ukca_mode) THEN
  ierr = 700
  cmessage = &
    'ukca_radaer_init: Maximum number of mode inconsistent.'
  IF (lhook) CALL dr_hook('UKCA_RADAER_INIT', zhook_out, zhook_handle)
  RETURN
END IF

IF (ncp /= npd_ukca_cpnt_in_mode) THEN
  ierr = 700
  cmessage = &
    'ukca_radaer_init: Maximum number of components inconsistent.'
  IF (lhook) CALL dr_hook('UKCA_RADAER_INIT', zhook_out, zhook_handle)
  RETURN
END IF

! Define UKCA-mode setup based on value of i_mode_setup from UMUI
! below are called from modules UKCA_SETUP_INDICES and UKCA_MODE_SETUP.
! Note that this piece of code must be consistent with that in 
! UKCA_MAIN().
IF (i_mode_setup == 1) THEN
  CALL ukca_mode_allcp_4mode
ELSE IF (i_mode_setup == 2) THEN
  CALL ukca_mode_sussbcoc_5mode
END IF
!
! Loop on all modes and components per mode, and only retain
! included modes and components.
!
n_loc_mode = 0
n_loc_cpnt = 0
DO i = 1, npd_ukca_mode

  IF (mode(i)) THEN
    
    this_name = mode_names(i)
    !
    ! Get the mode type. Since there is no direct information,
    ! it is obtained from the mode names.
    !
    SELECT CASE (this_name(1:3))
      
      CASE ('nuc')
        this_type = IP_UKCA_MODE_NUCLEATION
      
      CASE ('ait')
        this_type = IP_UKCA_MODE_AITKEN
      
      CASE ('acc')
        this_type = IP_UKCA_MODE_ACCUM
      
      CASE ('cor')
        this_type = IP_UKCA_MODE_COARSE
      
      CASE default
        ierr = 701
        cmessage = 'ukca_radaer_init: Unexpected mode name.'
        IF (lhook) CALL dr_hook('UKCA_RADAER_INIT', zhook_out, zhook_handle)
        RETURN
    
    END SELECT
    
    !
    ! Interaction of nucleation modes with radiation is
    ! neglected.
    !
    IF (this_type /= IP_UKCA_MODE_NUCLEATION) THEN
      
      n_loc_mode = n_loc_mode + 1
      
      ukca_radaer%i_mode_type(n_loc_mode) = this_type
      
      ukca_radaer%l_soluble(n_loc_mode) = modesol(i) == 1
      
      ukca_radaer%d0low(n_loc_mode) = ddplim0(i)
      
      ukca_radaer%d0up(n_loc_mode) = ddplim1(i)
      
      ukca_radaer%sigma(n_loc_mode) = sigmag(i)
      
      ukca_radaer%n_cpnt_in_mode(n_loc_mode) = 0
      
      ukca_radaer%stashcode_dry(n_loc_mode) = stashc_dry(i)
      ukca_radaer%stashcode_wet(n_loc_mode) = stashc_wet(i)
      ukca_radaer%stashcode_rho(n_loc_mode) = stashc_rho(i)
      ukca_radaer%stashcode_wtv(n_loc_mode) = stashc_cvl(i, IP_UKCA_WATER)
      ukca_radaer%stashcode_nbr(n_loc_mode) = stashc_nbr(i)
      
      !
      ! Loop on components within that mode.
      !
      DO j = 1, npd_ukca_cpnt_in_mode
      
        IF (component(i, j)) THEN
          
          n_loc_cpnt = n_loc_cpnt + 1
          
          !
          ! Update the number of components in that mode and
          ! retain the array index of the current component.
          !
          ukca_radaer%n_cpnt_in_mode(n_loc_mode) =                    &
            ukca_radaer%n_cpnt_in_mode(n_loc_mode) + 1
          ukca_radaer%i_cpnt_index(                                   &
            ukca_radaer%n_cpnt_in_mode(n_loc_mode), n_loc_mode) =     &
            n_loc_cpnt
          
          !
          ! Get the component type. Since there is no direct
          ! information, it is obtained from the component
          ! names.
          !
          SELECT CASE (component_names(j))
            
            CASE ('h2so4  ')
              ukca_radaer%i_cpnt_type(n_loc_cpnt) = IP_UKCA_SULPHATE
            
            CASE ('bcarbon')
              ukca_radaer%i_cpnt_type(n_loc_cpnt) = IP_UKCA_BLACKCARBON
            
            CASE ('ocarbon')
              ukca_radaer%i_cpnt_type(n_loc_cpnt) = IP_UKCA_ORGANICCARBON
            
            CASE ('nacl   ')
              ukca_radaer%i_cpnt_type(n_loc_cpnt) = IP_UKCA_SEASALT
            
            CASE ('dust   ')
              ukca_radaer%i_cpnt_type(n_loc_cpnt) = IP_UKCA_DUST
            
            CASE ('sec_org')
              ukca_radaer%i_cpnt_type(n_loc_cpnt) = IP_UKCA_SECONDORGANIC
            
            CASE default
              ierr = 702
              cmessage = &
                'ukca_radaer_init: Unexpected component name: ' // &
                component_names(j)
              IF (lhook) THEN
                CALL dr_hook('UKCA_RADAER_INIT', zhook_out, zhook_handle)
              END IF
              RETURN
          
          END SELECT
          
          
          ukca_radaer%density(n_loc_cpnt) = rhocomp(j)
          
          ukca_radaer%i_mode(n_loc_cpnt) = n_loc_mode
          
          ukca_radaer%stashcode_mmr(n_loc_cpnt) = stashc_mmr(i, j)
          
          ukca_radaer%stashcode_cvl(n_loc_cpnt) = stashc_cvl(i, j)
          
        END IF
      
      END DO ! j
    
    END IF
  
  END IF

END DO ! i

ukca_radaer%n_mode = n_loc_mode
ukca_radaer%n_cpnt = n_loc_cpnt

IF (ukca_radaer%n_mode == 0 .or. ukca_radaer%n_cpnt == 0) THEN
  ierr = 703
  cmessage = 'ukca_radaer_init: Setup includes no UKCA aerosols.'
  IF (lhook) CALL dr_hook('UKCA_RADAER_INIT', zhook_out, zhook_handle)
  RETURN
END IF

IF (lhook) CALL dr_hook('UKCA_RADAER_INIT', zhook_out, zhook_handle)
      
RETURN
END SUBROUTINE ukca_radaer_init
