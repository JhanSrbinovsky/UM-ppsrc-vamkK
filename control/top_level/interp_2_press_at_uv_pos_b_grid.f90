! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Interpolate vertically from model levels onto pressure levels and 
! horizontally from atmospheric u points on the c grid to uv positions 
! on the B grid.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

SUBROUTINE Interp_2_Press_At_UV_Pos_B_Grid(                       &
           im_index, item, section,                               &
           from, row_length, rows, n_rows, num_p_levels,          &
           model_levels, off_x,off_y, exner_theta_levels,         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                                             to)
  USE Submodel_Mod
  USE atmos_constants_mod, ONLY: kappa, p_zero

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  USE interpor_mod, ONLY: interp_order_linear
  USE uc_to_ub_mod, ONLY: uc_to_ub
  USE vert_interp2_mod, ONLY: vert_interp2
  IMPLICIT NONE

  INTEGER :: i,j,k
  INTEGER :: im_index
  INTEGER :: item, section
  INTEGER :: row_length, rows, n_rows
  INTEGER :: num_p_levels, model_levels
  INTEGER :: off_x, off_y

  REAL :: exner_theta_levels                                      &
     (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)

  REAL :: from(row_length,rows,model_levels)
  REAL :: exner_at_theta(row_length,rows,model_levels)
  REAL :: to(row_length, n_rows, num_p_levels) 
  REAL :: work(row_length, rows)

  INTEGER :: index

  REAL :: exner

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

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

  IF (lhook) CALL dr_hook('INTERP_2_PRESS_AT_UV_POS_B_GRID',zhook_in,zhook_handle)

!  Calculate exner at theta points. Store in exner_at_theta

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        exner_at_theta(i,j,k) = exner_theta_levels(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k

  index = -stlist(10,stindex(1,item,section,im_index))
  DO i = 1, num_p_levels
    exner = (((stash_levels(i+1,index)/1000.0) * 100.0)/p_zero)**kappa

    CALL vert_interp2 (from                                       &
      ,row_length, rows, model_levels                             &
      ,exner                                                      &
      ,0,0,0,0                                                    &
      ,exner_at_theta, interp_order_linear                        &
      ,work)

! Perform simple horizontal interpolation from 'C' to 'B' grid
    CALL  uC_to_uB(work,                                          &
      row_length,rows,n_rows,1,off_x,off_y,                       &
      to(1,1,i))
  END DO

  IF (lhook) CALL dr_hook('INTERP_2_PRESS_AT_UV_POS_B_GRID',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE Interp_2_Press_At_UV_Pos_B_Grid
