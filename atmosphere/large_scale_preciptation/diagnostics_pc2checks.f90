! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: large scale precip
MODULE diagnostics_pc2checks_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE diagnostics_pc2checks(                                 &
                             row_length, rows, model_levels             &
      ,                      wet_model_levels                           &
      ,                      me, timestep                               &
      ,                      at_extremity                               &
      ,                      t_inc, q_inc, qcl_inc, qcf_inc             &
      ,                      cf_inc, cfl_inc, cff_inc                   &
      ,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
       stashwork                                                        &
        )

! Purpose:
!          Calculates diagnostics and outputs them.

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      USE Submodel_Mod
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      LOGICAL                                                           &
        at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Parameters

      INTEGER                                                           &
        row_length                                                      &
                            ! number of points on a row
      , rows                                                            &
                            ! number of rows in a theta field
      , model_levels                                                    &
                            ! number of model levels
      , wet_model_levels    ! number of model levels where moisture

      INTEGER                                                           &
        me                  ! IN. Processor number

      REAL                                                              &
        timestep

! Primary Arrays used in all models
      REAL                                                              &
        t_inc(row_length, rows, model_levels)                           &
      , q_inc(row_length,rows, wet_model_levels)                        &
      , qcl_inc(row_length, rows, wet_model_levels)                     &
      , qcf_inc(row_length, rows, wet_model_levels)                     &
      , cf_inc(row_length,rows, wet_model_levels)                       &
      , cfl_inc(row_length, rows, wet_model_levels)                     &
      , cff_inc(row_length, rows, wet_model_levels)

      ! 3d work array for calculating separate +/- PC2 increments
      ! only to be allocated if it is needed.
      REAL, ALLOCATABLE ::  work3d(:,:,:)

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

! Diagnostic variables
      REAL                                                              &
       stashwork(*)        ! STASH workspace for section 4 (PC2 checks)

! Local variables
      INTEGER                                                           &
       i, j, k                                                          &
      ,    icode           ! Return code  =0 Normal exit  >1 Error

      INTEGER sect,item    ! STASH section, item no.s
      PARAMETER (sect = 4) !  for microphysics - large scale rain

      CHARACTER(LEN=80)  cmessage

      CHARACTER(LEN=*) RoutineName
      PARAMETER ( RoutineName='diagnostics_pc2checks')

      INTEGER                                                           &
        im_index        ! internal model index

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('DIAGNOSTICS_PC2CHECKS',zhook_in,zhook_handle)
      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! Allocate work3d array if calculating +/- increments for cfl,cff,qcl,qcf
      IF (sf(130,sect) .OR. sf(131,sect) .OR.                           &
          sf(132,sect) .OR. sf(133,sect) .OR.                           &
          sf(136,sect) .OR. sf(137,sect) .OR.                           &
          sf(138,sect) .OR. sf(139,sect) ) THEN
        ALLOCATE ( work3d(row_length,rows,wet_model_levels) )
      END IF

! Copy diagnostic information to STASHwork for STASH processing

! increment diagnostics= modified - previous

      item = 141  ! temperature increment
      IF (icode.le.0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
         CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
              t_inc,                                                    &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         IF (icode.gt.0) THEN
            cmessage=": error in copydiag_3d(item 141)"//cmessage
         END IF

      END IF  !  sf(item,sect)

      item = 142  ! vapour increment
      IF (icode.le.0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
         CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
              q_inc,                                                    &
              row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         IF (icode.gt.0) THEN
            cmessage=": error in copydiag_3d(item 142)"//cmessage
         END IF

      END IF  !  sf(item,sect)

      item = 143  ! liquid increment net
      IF (icode.le.0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
         CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
              qcl_inc,                                                  &
              row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         IF (icode.gt.0) THEN
            cmessage=": error in copydiag_3d(item 143)"//cmessage
         END IF

      END IF  !  sf(item,sect)

      item = 130  ! liquid increment: positive
      IF (icode.le.0 .AND. sf(item,sect)) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MAX(0.0,qcl_inc(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work3d,                                                    &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode.gt.0) THEN
          cmessage=": error in copydiag_3d(item 130)"//cmessage
        END IF

      END IF  !  sf(item,sect)

      item = 131  ! liquid increment: negative
      IF (icode.le.0 .AND. sf(item,sect)) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MIN(0.0,qcl_inc(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work3d,                                                    &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode.gt.0) THEN
          cmessage=": error in copydiag_3d(item 131)"//cmessage
        END IF

      END IF  !  sf(item,sect)

      item = 144  ! ice increment: net
      IF (icode.le.0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
         CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
              qcf_inc,                                                  &
              row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         IF (icode.gt.0) THEN
            cmessage=": error in copydiag_3d(item 144)"//cmessage
         END IF

      END IF  !  sf(item,sect)

      item = 132  ! ice increment: positive
      IF (icode.le.0 .AND. sf(item,sect)) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MAX(0.0,qcf_inc(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work3d,                                                    &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode.gt.0) THEN
          cmessage=": error in copydiag_3d(item 132)"//cmessage
        END IF

      END IF  !  sf(item,sect)

      item = 133  ! ice increment: negative
      IF (icode.le.0 .AND. sf(item,sect)) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MIN(0.0,qcf_inc(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work3d,                                                    &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode.gt.0) THEN
          cmessage=": error in copydiag_3d(item 133)"//cmessage
        END IF

      END IF  !  sf(item,sect)

      item = 152  ! total cloud fraction increment
      IF (icode.le.0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
         CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
              cf_inc,                                                   &
              row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         IF (icode.gt.0) THEN
            cmessage=": error in copydiag_3d(item 152)"//cmessage
         END IF

      END IF  !  sf(item,sect)

      item = 153  ! liquid cloud fraction increment: net
      IF (icode.le.0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
         CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
              cfl_inc,                                                  &
              row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         IF (icode.gt.0) THEN
            cmessage=": error in copydiag_3d(item 153)"//cmessage
         END IF

      END IF  !  sf(item,sect)

      item = 136  ! liquid cloud fraction increment: positive 
      IF (icode.le.0 .AND. sf(item,sect)) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MAX(0.0,cfl_inc(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work3d,                                                    &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode.gt.0) THEN
          cmessage=": error in copydiag_3d(item 136)"//cmessage
        END IF

      END IF  !  sf(item,sect)

      item = 137  ! liquid cloud fraction increment: positive 
      IF (icode.le.0 .AND. sf(item,sect)) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MIN(0.0,cfl_inc(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work3d,                                                    &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode.gt.0) THEN
          cmessage=": error in copydiag_3d(item 137)"//cmessage
        END IF

      END IF  !  sf(item,sect)

      item = 154  ! ice cloud fraction increment: net
      IF (icode.le.0 .AND. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
         CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
              cff_inc,                                                  &
              row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         IF (icode.gt.0) THEN
            cmessage=": error in copydiag_3d(item 154)"//cmessage
         END IF

      END IF  !  sf(item,sect)

      item = 138  ! ice cloud fraction increment: positive 
      IF (icode.le.0 .AND. sf(item,sect)) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MAX(0.0,cff_inc(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work3d,                                                    &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode.gt.0) THEN
          cmessage=": error in copydiag_3d(item 138)"//cmessage
        END IF

      END IF  !  sf(item,sect)

      item = 139  ! ice cloud fraction increment: negative 
      IF (icode.le.0 .AND. sf(item,sect)) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              work3d(i,j,k) = MIN(0.0,cff_inc(i,j,k))
            END DO
          END DO
        END DO
! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work3d,                                                    &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode.gt.0) THEN
          cmessage=": error in copydiag_3d(item 139)"//cmessage
        END IF

      END IF  !  sf(item,sect)

! Allocate work3 array if calculating +/- increments for cfl,cff,qcl,qcf.
      IF (sf(130,sect) .OR. sf(131,sect) .OR.                           &
          sf(132,sect) .OR. sf(133,sect) .OR.                           &
          sf(136,sect) .OR. sf(137,sect) .OR.                           &
          sf(138,sect) .OR. sf(139,sect) ) THEN
        DEALLOCATE ( work3d )
      END IF

 9999 CONTINUE
      IF(icode.ne.0) THEN
        CALL ereport(RoutineName,icode,cmessage)
      END IF

      IF (lhook) CALL dr_hook('DIAGNOSTICS_PC2CHECKS',zhook_out,zhook_handle)
      RETURN
    END SUBROUTINE diagnostics_pc2checks
END MODULE diagnostics_pc2checks_mod
