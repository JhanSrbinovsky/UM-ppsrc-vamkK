! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Interface:
      SUBROUTINE Diagnostics_solver(                                    &
     & row_length,rows,n_rows,model_levels,                             &
! wind field increments after difusion:
     & R_u,R_v,R_w,                                                     &
! wind field increments before diffusion (on stashflag):
     & u_incr_diagnostic,v_incr_diagnostic,                             &
     & w_incr_diagnostic,                                               &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork)

      USE atm_fields_bounds_mod

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Submodel_Mod
      IMPLICIT NONE
!
! Description:
!   Diagnostics_solver calculates diagnostics for the Helmholtz solver
!   for output by STASH routines for UM section 10.
!
! Method:
!   Simple copying of relevant fields into STASHwork array for output.
!   Vertical compression of levels is handled by copydiag_3D routine
!   using STASH control arrays passed via argsts.h.
!   Sequential processing of each diagnostic is performed, dependent
!   upon STASH flags being set.
!   List of diagnostics: (item,section)
!   (185,10) u increment          = delta(R_u) across solver
!   (186,10) v increment          = delta(R_v) across solver
!   (187,10) w increment          = delta(R_v) across solver
!
!
!   Note: no error trapping - cosmetic Errorstatus/cmessage needed for
!   compilation - no checks performed at lower levels.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: dynamics solver
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER, INTENT(IN) :: ExampleVariable  ! Description of variable
!
! Global variables (#include statements etc):

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

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER, INTENT(IN) :: row_length,rows                            &
                                              ! horizontal dimensions
     &,n_rows                                                           &
                        ! rows size for v-arrays
     &,model_levels     ! vertical levels

!   Array  arguments with intent(in):
      REAL, INTENT(IN) :: R_u    (udims_s%i_start:udims_s%i_end,        &
                                  udims_s%j_start:udims_s%j_end,        &
                                  udims_s%k_start:udims_s%k_end)
      REAL, INTENT(IN) :: R_v    (vdims_s%i_start:vdims_s%i_end,        &
                                  vdims_s%j_start:vdims_s%j_end,        &
                                  vdims_s%k_start:vdims_s%k_end) 
      REAL, INTENT(IN) :: R_w    (wdims%i_start:wdims%i_end,            &
                                  wdims%j_start:wdims%j_end,            &
                                  wdims%k_start:wdims%k_end)

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
      REAL, INTENT(INOUT) :: u_incr_diagnostic(udims%i_start:udims%i_end, &
                                               udims%j_start:udims%j_end, &
                                               udims%k_start:udims%k_end)
      REAL, INTENT(INOUT) :: v_incr_diagnostic(vdims%i_start:vdims%i_end, &
                                               vdims%j_start:vdims%j_end, &
                                               vdims%k_start:vdims%k_end)
      REAL, INTENT(INOUT) :: w_incr_diagnostic(wdims%i_start:wdims%i_end, &
                                               wdims%j_start:wdims%j_end, &
                                               wdims%k_start:wdims%k_end)

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL, INTENT(OUT) ::                                              &
     & STASHwork(*)   ! Output array holding diagnostic fields

! Local parameters:
      CHARACTER (LEN=*), PARAMETER :: RoutineName='Diagnostics_solver'

      INTEGER, PARAMETER :: sect =10  ! STASH section for diagnostics

! Local scalars:
      INTEGER  ::                                                       &
     & i,j,k                                                            &
                   !  loop indices
     &,im_index                                                         &
                   !  internal model index for STASH arrays
     &,item        !  STASH item of diagnostic
      INTEGER :: Errorstatus = 0  ! initial value for error code

      CHARACTER (LEN=80) :: CMessage !  Error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Local dynamic arrays: NONE at present
!      REAL  ::

!- End of header

!
! 1. Initialisation
!
      IF (lhook) CALL dr_hook('DIAGNOSTICS_SOLVER',zhook_in,zhook_handle)
      im_index    = internal_model_index(atmos_im)
      Cmessage    = ''
!
! 2. Extract diagnostic fields dependent on STASHflags sf
!
! u wind increment
      item = 185           ! u increment
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

        DO k= udims%k_start, udims%k_end
          DO j= udims%j_start, udims%j_end
            DO i= udims%i_start, udims%i_end
              u_incr_diagnostic(i,j,k) = R_u(i,j,k) -                   &
     &                                        u_incr_diagnostic(i,j,k)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        u_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      ENDIF ! sf(item,sect)

! v wind increment
      item = 186           ! v increment
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

        DO k= vdims%k_start, vdims%k_end
          DO j= vdims%j_start, vdims%j_end
            DO i= vdims%i_start, vdims%i_end
              v_incr_diagnostic(i,j,k) = R_v(i,j,k) -                   &
     &                                        v_incr_diagnostic(i,j,k)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        v_incr_diagnostic,                                        &
     &        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      ENDIF ! sf(item,sect)

! w  increment
      item = 187           ! w increment
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

        DO k= wdims%k_start, wdims%k_end
          DO j= wdims%j_start, wdims%j_end
            DO i= wdims%i_start, wdims%i_end
              w_incr_diagnostic(i,j,k) = R_w(i,j,k) -                   &
     &                             w_incr_diagnostic(i,j,k)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        w_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      ENDIF ! sf(item,sect)

!
! 3. Error handling
!
      IF (Errorstatus /= 0) THEN

        CALL Ereport(RoutineName,Errorstatus,Cmessage)
      ENDIF

      IF (lhook) CALL dr_hook('DIAGNOSTICS_SOLVER',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Diagnostics_solver
