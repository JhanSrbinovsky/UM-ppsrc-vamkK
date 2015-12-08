! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate diagnostic quantities from the initial atmosphere dump
!

      MODULE diagnostics_adv_mod

      IMPLICIT NONE

      CONTAINS

! Subroutine Interface:
      SUBROUTINE Diagnostics_adv(                                       &
     & row_length,rows,n_rows,model_levels,wet_levels,                  &
! primary wind fields:
     & u,v,theta,q,qcl,qcf,qrain,qgraup,qcf2,cf,cfl,cff,                &
! wind field increments after advection:
     & R_u,R_v,R_W,                                                     &
! wind field increments before advection (on stashflag):
     & u_incr_diagnostic,v_incr_diagnostic,                             &
     & T_incr_diagnostic,q_incr_diagnostic,                             &
     & qcl_incr_diagnostic,qcf_incr_diagnostic,                         &
       qrain_incr_diagnostic,qgraup_incr_diagnostic,                    &
       qcf2_incr_diagnostic,                                            &
     & cf_incr_diagnostic,cfl_incr_diagnostic,cff_incr_diagnostic,      &
     & theta_star,q_star,qcl_star,qcf_star,                             &
       qrain_star,qgraup_star,qcf2_star,                                & 
     & cf_star,cfl_star,cff_star,                                       &
     & exner_theta_levels,                                              &
! w departure point information
     & depart_lambda,depart_phi,depart_r,r_theta_levels,                &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork)

Use ac_diagnostics_mod, Only :                                          &
    qcl_adv

USE atm_fields_bounds_mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE UM_ParVars
USE Submodel_Mod
IMPLICIT NONE
!
! Description:
!   Diagnostics_adv extracts diagnostics of (N+1) time-level estimates
!   of primary fields after advection has been called, to be processed
!   by STASH routines for UM section 12 (advection).
! Method:
!   Simple copying of relevant fields into STASHwork array for output.
!   Vertical compression of levels is handled by copydiag_3D routine
!   using STASH control arrays passed via argsts.h.
!   Sequential processing of each diagnostic is performed, dependent
!   upon STASH flags being set.
!   List of diagnostics: (item,section)
!   (  2,12) u wind               = u + R_u
!   (  3,12) v wind               = v + R_v
!   (  4,12) temperature          = theta_star/exner_theta_levels
!   ( 10,12) specific humidity    = q_star
!   (254,12) qcl                  = qcl_star
!   ( 12,12) qcf                  = qcf_star
!   (185,12) u increment          = delta(R_u) across advection
!   (186,12) v increment          = delta(R_v) across advection
!   (187,12) w increment          = delta(w  ) across advection
!   (181,12) T increment          = delta(theta)/exner
!   (182,12) q   increment        = delta(q)   across advection
!   (183,12) qcl increment        = delta(qcl) across advection
!   (184,12) qcf increment        = delta(qcf) across advection
!   (189,12) qrain increment      = delta(qrain) across advection
!   (190,12) qgraup increment     = delta(qgraup) across advection
!   (191,12) qcf2 increment       = delta(qcf2) across advection
!   (188,12) cf  increment        = delta(cf)   across advection
!   (199,12) cfl increment        = delta(cfl) across advection
!   (190,12) cff increment        = delta(cff) across advection
!   (204,12) departure point (w)  = depart_lambda
!   (203,12) departure point (w)  = depart_phi
!   (205,12) model height diff    = depart_r - r_theta_levels
!
!   _star fields are estimates of N+1 time-level quantities; R_u/v are
!   wind increments from physics1.
!
!   Note: no error trapping - cosmetic Errorstatus/cmessage needed for
!   compilation - no checks performed at lower levels.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: dynamics advection
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
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
      INTEGER                                                           &
     & row_length,rows                                                  &
                        ! horizontal dimensions
     &,n_rows                                                           &
                        ! rows for last (N) row of pes
     &,model_levels                                                     &
                        ! vertical levels
     &,wet_levels       ! vertical levels with moisture

!   Array  arguments with intent(in):
      REAL :: u      (udims_s%i_start:udims_s%i_end,        &
                      udims_s%j_start:udims_s%j_end,        &
                      udims_s%k_start:udims_s%k_end)
      REAL :: v      (vdims_s%i_start:vdims_s%i_end,        &
                      vdims_s%j_start:vdims_s%j_end,        &
                      vdims_s%k_start:vdims_s%k_end)

      REAL :: theta  (tdims_s%i_start:tdims_s%i_end,        &
                      tdims_s%j_start:tdims_s%j_end,        &
                      tdims_s%k_start:tdims_s%k_end)

      REAL :: q      (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)
      REAL :: qcl    (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)       
      REAL :: qcf    (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)

      REAL :: qrain  (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)
      REAL :: qgraup (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)       
      REAL :: qcf2   (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)

      REAL :: cf     (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)
      REAL :: cfl    (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)       
      REAL :: cff    (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)

      REAL :: R_u    (udims_s%i_start:udims_s%i_end,        &
                      udims_s%j_start:udims_s%j_end,        &
                      udims_s%k_start:udims_s%k_end)
      REAL :: R_v    (vdims_s%i_start:vdims_s%i_end,        &
                      vdims_s%j_start:vdims_s%j_end,        &
                      vdims_s%k_start:vdims_s%k_end) 
      REAL :: R_w    (wdims%i_start:wdims%i_end,            &
                      wdims%j_start:wdims%j_end,            &
                      wdims%k_start:wdims%k_end)

      REAL :: u_incr_diagnostic    (udims%i_start:udims%i_end,         &
                                    udims%j_start:udims%j_end,         &
                                    udims%k_start:udims%k_end)
      REAL :: v_incr_diagnostic    (vdims%i_start:vdims%i_end,         &
                                    vdims%j_start:vdims%j_end,         &
                                    vdims%k_start:vdims%k_end)
      REAL :: T_incr_diagnostic    (tdims%i_start:tdims%i_end,         &
                                    tdims%j_start:tdims%j_end,         &
                                                1:tdims%k_end)       

      REAL :: q_incr_diagnostic    (qdims%i_start:qdims%i_end,         &
                                    qdims%j_start:qdims%j_end,         &
                                                1:qdims%k_end)
      REAL :: qcl_incr_diagnostic  (qdims%i_start:qdims%i_end,         &
                                    qdims%j_start:qdims%j_end,         &
                                                1:qdims%k_end)
      REAL :: qcf_incr_diagnostic  (qdims%i_start:qdims%i_end,         &
                                    qdims%j_start:qdims%j_end,         &
                                                1:qdims%k_end)

      REAL :: qrain_incr_diagnostic (qdims%i_start:qdims%i_end,        &
                                     qdims%j_start:qdims%j_end,        &
                                                 1:qdims%k_end)
      REAL :: qgraup_incr_diagnostic(qdims%i_start:qdims%i_end,        &
                                     qdims%j_start:qdims%j_end,        &
                                                 1:qdims%k_end)       
      REAL :: qcf2_incr_diagnostic  (qdims%i_start:qdims%i_end,        &
                                     qdims%j_start:qdims%j_end,        &
                                                 1:qdims%k_end)

      REAL :: cf_incr_diagnostic    (qdims%i_start:qdims%i_end,        &
                                     qdims%j_start:qdims%j_end,        &
                                                 1:qdims%k_end)
      REAL :: cfl_incr_diagnostic   (qdims%i_start:qdims%i_end ,       &
                                     qdims%j_start:qdims%j_end ,       &
                                                 1:qdims%k_end)
      REAL :: cff_incr_diagnostic   (qdims%i_start:qdims%i_end,        &
                                     qdims%j_start:qdims%j_end,        &
                                                 1:qdims%k_end)

      REAL :: theta_star  (tdims_s%i_start:tdims_s%i_end,        &
                           tdims_s%j_start:tdims_s%j_end,        &
                           tdims_s%k_start:tdims_s%k_end)

      REAL :: q_star      (qdims_s%i_start:qdims_s%i_end,        &
                           qdims_s%j_start:qdims_s%j_end,        &
                           qdims_s%k_start:qdims_s%k_end)
      REAL :: qcl_star    (qdims_s%i_start:qdims_s%i_end,        &
                           qdims_s%j_start:qdims_s%j_end,        &
                           qdims_s%k_start:qdims_s%k_end)       
      REAL :: qcf_star    (qdims_s%i_start:qdims_s%i_end,        &
                           qdims_s%j_start:qdims_s%j_end,        &
                           qdims_s%k_start:qdims_s%k_end)

      REAL :: qrain_star  (qdims_s%i_start:qdims_s%i_end,        &
                           qdims_s%j_start:qdims_s%j_end,        &
                           qdims_s%k_start:qdims_s%k_end)
      REAL :: qgraup_star (qdims_s%i_start:qdims_s%i_end,        &
                           qdims_s%j_start:qdims_s%j_end,        &
                           qdims_s%k_start:qdims_s%k_end)       
      REAL :: qcf2_star   (qdims_s%i_start:qdims_s%i_end,        &
                           qdims_s%j_start:qdims_s%j_end,        &
                           qdims_s%k_start:qdims_s%k_end) 
 
      REAL :: cf_star     (qdims_s%i_start:qdims_s%i_end,        &
                           qdims_s%j_start:qdims_s%j_end,        &
                           qdims_s%k_start:qdims_s%k_end)
      REAL :: cfl_star    (qdims_s%i_start:qdims_s%i_end,        &
                           qdims_s%j_start:qdims_s%j_end,        &
                           qdims_s%k_start:qdims_s%k_end)       
      REAL :: cff_star    (qdims_s%i_start:qdims_s%i_end,        &
                           qdims_s%j_start:qdims_s%j_end,        &
                           qdims_s%k_start:qdims_s%k_end) 

      REAL :: exner_theta_levels(tdims_s%i_start:tdims_s%i_end,        &
                                 tdims_s%j_start:tdims_s%j_end,        &
                                 tdims_s%k_start:tdims_s%k_end) 

      REAL &
     & depart_lambda(row_length, rows, model_levels)                    &
     &,depart_phi(row_length, rows, model_levels)                       &
     &,depart_r(row_length, rows, model_levels)

    REAL :: r_theta_levels      (tdims_l%i_start:tdims_l%i_end,        &
                                 tdims_l%j_start:tdims_l%j_end,        &
                                               0:tdims_l%k_end) 

    REAL, ALLOCATABLE :: qcl_incr_store(:,:,:),                         &
                         qcf_incr_store(:,:,:),                         &
                         cfl_incr_store(:,:,:),                         &
                         cff_incr_store(:,:,:)

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL                                                              &
     & STASHwork(*)   ! Output array holding diagnostic fields

! Local parameters:
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='Diagnostics_adv')

      INTEGER sect      ! STASH section for diagnostics
      PARAMETER ( sect = 12 )

! Local scalars:
      INTEGER                                                           &
     & i,j,k,ji                                                         &
                   !  loop indices
     &,im_index                                                         &
                   !  internal model index for STASH arrays
     &,item                                                             &
                   !  STASH item of diagnostic
     &,Errorstatus !  Error status

      CHARACTER(LEN=80) CMessage !  Error message

! Local dynamic arrays:
      REAL :: work_1(tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                                 1:tdims%k_end)
      REAL :: work_u(udims%i_start:udims%i_end,        &
                     udims%j_start:udims%j_end,        &
                     udims%k_start:udims%k_end)
      REAL :: work_v(vdims%i_start:vdims%i_end,        &
                     vdims%j_start:vdims%j_end,        &
                     vdims%k_start:vdims%k_end)
      REAL :: w_incr_diagnostic                        &
                    (wdims%i_start:wdims%i_end,        &
                     wdims%j_start:wdims%j_end,        &
                                 1:wdims%k_end)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of header

!
! 1. Initialisation
!
      IF (lhook) CALL dr_hook('DIAGNOSTICS_ADV',zhook_in,zhook_handle)
      im_index    = internal_model_index(atmos_im)
      Errorstatus = 0
      Cmessage    = ''

! 1.5. If calculating separate positive and negative increments for cloud
! fields, need to store the values into temporary arrays.

      IF ( (sf(170,sect) .AND. Errorstatus == 0) .OR.                   &
           (sf(171,sect) .AND. Errorstatus == 0)) THEN
        ALLOCATE ( qcl_incr_store(qdims%i_start:qdims%i_end,            &
                                  qdims%j_start:qdims%j_end,            &
                                              1:qdims%k_end) )
        DO k = 1, qdims%k_end 
          DO j = qdims%j_start,qdims%j_end 
           DO i = qdims%i_start, qdims%i_end 
              qcl_incr_store(i,j,k) = qcl_incr_diagnostic(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k
      END IF ! if requesting separate positive or negative incr for qcl

      IF ( (sf(172,sect) .AND. Errorstatus == 0) .OR.                   &
           (sf(173,sect) .AND. Errorstatus == 0)) THEN
        ALLOCATE ( qcf_incr_store(qdims%i_start:qdims%i_end,            &
                                  qdims%j_start:qdims%j_end,            &
                                              1:qdims%k_end) )
        DO k = 1, qdims%k_end 
          DO j = qdims%j_start,qdims%j_end 
           DO i = qdims%i_start, qdims%i_end 
              qcf_incr_store(i,j,k) = qcf_incr_diagnostic(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k
      END IF ! if requesting separate positive or negative incr for qcf

      IF ( (sf(176,sect) .AND. Errorstatus == 0) .OR.                   &
           (sf(177,sect) .AND. Errorstatus == 0)) THEN
        ALLOCATE ( cfl_incr_store(qdims%i_start:qdims%i_end,            &
                                  qdims%j_start:qdims%j_end,            &
                                              1:qdims%k_end) )
        DO k = 1, qdims%k_end 
          DO j = qdims%j_start,qdims%j_end 
           DO i = qdims%i_start, qdims%i_end 
              cfl_incr_store(i,j,k) = cfl_incr_diagnostic(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k
      END IF ! if requesting separate positive or negative incr for cfl

      IF ( (sf(178,sect) .AND. Errorstatus == 0) .OR.                   &
           (sf(179,sect) .AND. Errorstatus == 0)) THEN
        ALLOCATE ( cff_incr_store(qdims%i_start:qdims%i_end,            &
                                  qdims%j_start:qdims%j_end,            &
                                              1:qdims%k_end) )
        DO k = 1, qdims%k_end 
          DO j = qdims%j_start,qdims%j_end 
           DO i = qdims%i_start, qdims%i_end 
              cff_incr_store(i,j,k) = cff_incr_diagnostic(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k
      END IF ! if requesting separate positive or negative incr for cfl

!
! 2. Extract diagnostic fields dependent on STASHflags sf
!

! u wind estimate = u + physics1 increment
      item = 2             ! u
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

        DO k= udims%k_start, udims%k_end
          DO j= udims%j_start, udims%j_end
            DO i= udims%i_start, udims%i_end
              work_u(i,j,k) = u(i,j,k) + R_u(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),work_u,      &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! v wind estimate = v + physics1 increment
      item = 3             ! v
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

        DO k= vdims%k_start, vdims%k_end
          DO j= vdims%j_start, vdims%j_end
            DO i= vdims%i_start, vdims%i_end
              work_v(i,j,k) = v(i,j,k) + R_v(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),work_v,      &
     &        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)


! temperature estimate = theta / exner pressure
      item = 4             ! temperature
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

        DO k=             1, tdims%k_end
          DO j= tdims%j_start, tdims%j_end
            DO i= tdims%i_start, tdims%i_end
              work_1(i,j,k)=theta_star(i,j,k)*exner_theta_levels(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),work_1,      &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! specific humidity
      item = 10            ! q
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
              q_star(qdims_s%i_start:qdims_s%i_end,                     &
                     qdims_s%j_start:qdims_s%j_end,                     &
                                   1:qdims_s%k_end),                    &
     &        row_length,rows,wet_levels,0,0,1,1,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcl
      item = 254           ! qcl
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

      If (.not.allocated(qcl_adv)) Then
        allocate ( qcl_adv(row_length*rows,wet_levels) )
      End If

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
              qcl_star(qdims_s%i_start:qdims_s%i_end,                   &
                       qdims_s%j_start:qdims_s%j_end,                   &
                                     1:qdims_s%k_end),                  &
     &        row_length,rows,wet_levels,0,0,1,1,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

! Copy non-halo area of qcl_star into qcl_adv in module ac_diagnostics_mod
      DO k=             1, qdims%k_end
        DO j= qdims%j_start, qdims%j_end
          DO i= qdims%i_start, qdims%i_end
            ji = (j-1)*row_length+i
            qcl_adv(ji,k) = qcl_star(i,j,k)
          end do
        end do
      end do

      END IF ! sf(item,sect)

! qcf
      item = 12            ! qcf
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
              qcf_star(qdims_s%i_start:qdims_s%i_end,                   &
                       qdims_s%j_start:qdims_s%j_end,                   &
                                     1:qdims_s%k_end),                  &
     &        row_length,rows,wet_levels,0,0,1,1,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! u wind increment
      item = 185           ! u increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= udims%k_start, udims%k_end
          DO j= udims%j_start, udims%j_end
            DO i= udims%i_start, udims%i_end
              u_incr_diagnostic(i,j,k) = R_u(i,j,k) -                   &
     &                                        u_incr_diagnostic(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        u_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! v wind increment
      item = 186           ! v increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= vdims%k_start, vdims%k_end
          DO j= vdims%j_start, vdims%j_end
            DO i= vdims%i_start, vdims%i_end
              v_incr_diagnostic(i,j,k) = R_v(i,j,k) -                   &
     &                                        v_incr_diagnostic(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        v_incr_diagnostic,                                        &
     &        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! w wind increment
      item = 187           ! w increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=             1, wdims%k_end
          DO j= wdims%j_start, wdims%j_end
            DO i= wdims%i_start, wdims%i_end
              w_incr_diagnostic(i,j,k) = R_w(i,j,k)
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

      END IF ! sf(item,sect)

! T wind increment
! theta_star now holds theta+dtheta whereas t_incr holds dtheta before
! advection
      item = 181           ! T increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= 1, tdims%k_end
          DO j= tdims%j_start, tdims%j_end
!CDIR NOUNROLL
            DO i= tdims%i_start, tdims%i_end
              T_incr_diagnostic(i,j,k) = (theta_star(i,j,k) -           &
     &                 (theta(i,j,k) + T_incr_diagnostic(i,j,k)))       &
     &                                 *exner_theta_levels(i,j,k)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        T_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! q wind increment
      item = 182           ! q increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= 1, qdims%k_end
!CDIR NOUNROLL
          DO j= qdims%j_start, qdims%j_end
            DO i= qdims%i_start, qdims%i_end
              q_incr_diagnostic(i,j,k) = q_star(i,j,k) -                &
     &                     (q(i,j,k) + q_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        q_incr_diagnostic,                                        &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcl wind increment
      item = 183           ! qcl increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= 1, qdims%k_end
!CDIR NOUNROLL
          DO j= qdims%j_start, qdims%j_end
            DO i= qdims%i_start, qdims%i_end
              qcl_incr_diagnostic(i,j,k) = qcl_star(i,j,k) -            &
     &                     (qcl(i,j,k) + qcl_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        qcl_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcl wind increment: positive
      item = 170
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1, qdims%k_end
!CDIR NOUNROLL
          DO j=qdims%j_start, qdims%j_end
            DO i=qdims%i_start, qdims%i_end 
              qcl_incr_diagnostic(i,j,k) = max( 0.0, qcl_star(i,j,k) -  &
                (qcl(i,j,k) + qcl_incr_store(i,j,k)) )
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
              qcl_incr_diagnostic,                                      &
              row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcl wind increment: negative
      item = 171
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1, qdims%k_end
!CDIR NOUNROLL
          DO j=qdims%j_start, qdims%j_end
            DO i=qdims%i_start, qdims%i_end 
              qcl_incr_diagnostic(i,j,k) = min( 0.0, qcl_star(i,j,k) -  &
                (qcl(i,j,k) + qcl_incr_store(i,j,k)) )
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
              qcl_incr_diagnostic,                                      &
              row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcf increment
      item = 184           ! qcf increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= 1, qdims%k_end
!CDIR NOUNROLL
          DO j= qdims%j_start, qdims%j_end
            DO i= qdims%i_start, qdims%i_end
              qcf_incr_diagnostic(i,j,k) = qcf_star(i,j,k) -            &
     &                     (qcf(i,j,k) + qcf_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        qcf_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcf increment: positive
      item = 172
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1, qdims%k_end
!CDIR NOUNROLL
          DO j=qdims%j_start, qdims%j_end
            DO i=qdims%i_start, qdims%i_end
              qcf_incr_diagnostic(i,j,k) = max( 0.0, qcf_star(i,j,k) -   &
                (qcf(i,j,k) + qcf_incr_store(i,j,k)) )
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        qcf_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcf increment: negative
      item = 173
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1, qdims%k_end
!CDIR NOUNROLL
          DO j=qdims%j_start, qdims%j_end
            DO i=qdims%i_start, qdims%i_end
              qcf_incr_diagnostic(i,j,k) = min( 0.0, qcf_star(i,j,k) -   &
                (qcf(i,j,k) + qcf_incr_store(i,j,k)) )
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        qcf_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qrain increment
      item = 189           ! qrain increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= 1, qdims%k_end
!CDIR NOUNROLL
          DO j= qdims%j_start, qdims%j_end
            DO i= qdims%i_start, qdims%i_end
              qrain_incr_diagnostic(i,j,k) = qrain_star(i,j,k) -            &
     &                     (qrain(i,j,k) + qrain_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
              qrain_incr_diagnostic,                                    &
              row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qgraup increment
      item = 190           ! qrain increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= 1, qdims%k_end
!CDIR NOUNROLL
          DO j= qdims%j_start, qdims%j_end
            DO i= qdims%i_start, qdims%i_end
              qgraup_incr_diagnostic(i,j,k) = qgraup_star(i,j,k) -      &
                         (qgraup(i,j,k) + qgraup_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
              qgraup_incr_diagnostic,                                   &
              row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcf2 increment
      item = 191           ! qrain increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= 1, qdims%k_end
!CDIR NOUNROLL
          DO j= qdims%j_start, qdims%j_end
            DO i= qdims%i_start, qdims%i_end
              qcf2_incr_diagnostic(i,j,k) = qcf2_star(i,j,k) -      &
                         (qcf2(i,j,k) + qcf2_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
              qcf2_incr_diagnostic,                                     &
              row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              Errorstatus,cmessage)

      END IF ! sf(item,sect)

! cf increment
      item = 192           ! cf increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= 1, qdims%k_end
!CDIR NOUNROLL
          DO j= qdims%j_start, qdims%j_end
            DO i= qdims%i_start, qdims%i_end
              cf_incr_diagnostic(i,j,k) = cf_star(i,j,k) -              &
     &                     (cf(i,j,k) + cf_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        cf_incr_diagnostic,                                       &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! cfl increment
      item = 193           ! cfl increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= 1, qdims%k_end
!CDIR NOUNROLL
          DO j= qdims%j_start, qdims%j_end
            DO i= qdims%i_start, qdims%i_end
              cfl_incr_diagnostic(i,j,k) = cfl_star(i,j,k) -            &
     &                     (cfl(i,j,k) + cfl_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        cfl_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! cfl increment: positive
      item = 176
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1, qdims%k_end
!CDIR NOUNROLL
          DO j=qdims%j_start, qdims%j_end
            DO i=qdims%i_start, qdims%i_end
              cfl_incr_diagnostic(i,j,k) = max( 0.0, cfl_star(i,j,k) -  &
                (cfl(i,j,k) + cfl_incr_store(i,j,k)) )
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        cfl_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! cfl increment: positive
      item = 177
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1, qdims%k_end
!CDIR NOUNROLL
          DO j=qdims%j_start, qdims%j_end
            DO i=qdims%i_start, qdims%i_end
              cfl_incr_diagnostic(i,j,k) = min( 0.0, cfl_star(i,j,k) -  &
                (cfl(i,j,k) + cfl_incr_store(i,j,k)) )
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        cfl_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! cff increment
      item = 194           ! cff increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k= 1, qdims%k_end
!CDIR NOUNROLL
          DO j= qdims%j_start, qdims%j_end
            DO i= qdims%i_start, qdims%i_end
              cff_incr_diagnostic(i,j,k) = cff_star(i,j,k) -            &
     &                     (cff(i,j,k) + cff_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        cff_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! cff increment: positive
      item = 178
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1, qdims%k_end
!CDIR NOUNROLL
          DO j=qdims%j_start, qdims%j_end
            DO i=qdims%i_start, qdims%i_end
              cff_incr_diagnostic(i,j,k) = max( 0.0, cff_star(i,j,k) -   &
                (cff(i,j,k) + cff_incr_store(i,j,k)) )
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        cff_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! cff increment: negative
      item = 179
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1, qdims%k_end
!CDIR NOUNROLL
          DO j=qdims%j_start, qdims%j_end
            DO i=qdims%i_start, qdims%i_end
              cff_incr_diagnostic(i,j,k) = min( 0.0, cff_star(i,j,k) -   &
                (cff(i,j,k) + cff_incr_store(i,j,k)) )
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        cff_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! Departure point diagnostics for w
! (a) lambda
      item = 204
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &         depart_lambda,                                           &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)
! (b) phi
      item = 205
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &         depart_phi,                                              &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)
! (c) dr  difference from model height of departure point
      item = 203
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=  1, tdims%k_end
          DO j= tdims%j_start, tdims%j_end
            DO i= tdims%i_start, tdims%i_end
              work_1(i,j,k) = depart_r(i,j,k) -                         &
     &                                        r_theta_levels(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &         work_1,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! 2.5. Deallocate extra arrays used if separate positive and negative
!      increments were requested.

      IF ( (sf(178,sect) .AND. Errorstatus == 0) .OR.                   &
           (sf(179,sect) .AND. Errorstatus == 0)) THEN
        DEALLOCATE( cff_incr_store )
      END IF ! sf(item,sect)

      IF ( (sf(176,sect) .AND. Errorstatus == 0) .OR.                   &
           (sf(177,sect) .AND. Errorstatus == 0)) THEN
        DEALLOCATE( cfl_incr_store )
      END IF ! sf(item,sect)

      IF ( (sf(172,sect) .AND. Errorstatus == 0) .OR.                   &
           (sf(173,sect) .AND. Errorstatus == 0)) THEN
        DEALLOCATE( qcf_incr_store )
      END IF ! sf(item,sect)

      IF ( (sf(170,sect) .AND. Errorstatus == 0) .OR.                   &
           (sf(171,sect) .AND. Errorstatus == 0)) THEN
        DEALLOCATE( qcl_incr_store )
      END IF ! sf(item,sect)

! 3. Error handling
!
      IF (Errorstatus /= 0) THEN

        CALL Ereport(RoutineName,Errorstatus,Cmessage)
      END IF

      IF (lhook) CALL dr_hook('DIAGNOSTICS_ADV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE diagnostics_adv

      END MODULE diagnostics_adv_mod
