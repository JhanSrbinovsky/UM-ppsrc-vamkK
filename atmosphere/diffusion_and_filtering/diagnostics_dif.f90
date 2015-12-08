! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate quantities from diffusion and divergence damping (sect 13)
!
! Subroutine Interface:
      SUBROUTINE Diagnostics_dif(                                       &
     & row_length, rows, n_rows, model_levels, bl_levels,               &
! primary fields:
     & theta,q,                                                         &
! wind field increments after difusion:
     & R_u,R_v,                                                         &
! wind field increments before diffusion (on stashflag):
     & u_incr_diagnostic,v_incr_diagnostic,                             &
     & T_incr_diagnostic,q_incr_diagnostic,                             &
     & w_local_mask,                                                    &
! Current theta+dtheta and q+dq values
     & theta_star,q_star,                                               &
     & exner_theta_levels,                                              &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork)

      USE atm_fields_bounds_mod
      USE turb_diff_mod, ONLY: L_subfilter_horiz, L_subfilter_vert
      USE turb_diff_ctl_mod, ONLY: visc_m, visc_h, shear, rneutml_sq

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Submodel_Mod
      IMPLICIT NONE
!
! Description:
!   Diagnostics_fildif calculates diagnostics for divergence damping
!   and diffusion for output by STASH routines for UM section 13.
!
! Method:
!   Simple copying of relevant fields into STASHwork array for output.
!   Vertical compression of levels is handled by copydiag_3D routine
!   using STASH control arrays passed via argsts.h.
!   Sequential processing of each diagnostic is performed, dependent
!   upon STASH flags being set.
!   List of diagnostics: (item,section)
!   (181,13) T increment          = delta(theta)/exner
!   (182,13) q   increment        = delta(q)   across advection
!   (185,13) u increment          = delta(R_u) across advection
!   (186,13) v increment          = delta(R_v) across advection
!
!   _star fields are estimates of N+1 time-level quantities; R_u/v are
!   wind increments from all routines.
!
!   Note: no error trapping - cosmetic Errorstatus/cmessage needed for
!   compilation - no checks performed at lower levels.
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: diffusion and filtering
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
      INTEGER, INTENT(IN) :: row_length, rows                           &
                                              ! horizontal dimensions
     &,n_rows                                                           &
                        ! rows for last (N) row of pes
     &,model_levels                                                     &
                        ! vertical levels
     &,bl_levels        ! vertical levels with moisture

!   Array  arguments with intent(in):
      REAL, INTENT(IN) :: theta  (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN) :: q      (qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)
      REAL, INTENT(IN) :: R_u    (udims_s%i_start:udims_s%i_end,        &
                                  udims_s%j_start:udims_s%j_end,        &
                                  udims_s%k_start:udims_s%k_end)
      REAL, INTENT(IN) :: R_v    (vdims_s%i_start:vdims_s%i_end,        &
                                  vdims_s%j_start:vdims_s%j_end,        &
                                  vdims_s%k_start:vdims_s%k_end) 
      REAL, INTENT(IN) :: theta_star                                    &  
                                 (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(IN) :: q_star (qdims_s%i_start:qdims_s%i_end,        &
                                  qdims_s%j_start:qdims_s%j_end,        &
                                  qdims_s%k_start:qdims_s%k_end)
      REAL, INTENT(IN) :: exner_theta_levels                            &
                                 (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end)

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
      REAL, INTENT(INOUT) :: u_incr_diagnostic                         &   
                                   (udims%i_start:udims%i_end,         &
                                    udims%j_start:udims%j_end,         &
                                    udims%k_start:udims%k_end)
      REAL, INTENT(INOUT) :: v_incr_diagnostic                         & 
                                   (vdims%i_start:vdims%i_end,         &
                                    vdims%j_start:vdims%j_end,         &
                                    vdims%k_start:vdims%k_end)
      REAL, INTENT(INOUT) :: T_incr_diagnostic                         &  
                                   (tdims%i_start:tdims%i_end,         &
                                    tdims%j_start:tdims%j_end,         &
                                                1:tdims%k_end)       
      REAL, INTENT(INOUT) :: q_incr_diagnostic                         &
                                   (qdims%i_start:qdims%i_end,         &
                                    qdims%j_start:qdims%j_end,         &
                                                1:qdims%k_end)

      REAL, INTENT(INOUT) ::  w_local_mask(row_length, rows)

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL, INTENT(OUT) ::                                              &
     & STASHwork(*)   ! Output array holding diagnostic fields

! Local parameters:
      CHARACTER (LEN=*), PARAMETER :: RoutineName='Diagnostics_dif'

      INTEGER, PARAMETER :: sect =13  ! STASH section for diagnostics

! Local scalars:
      INTEGER  ::                                                       &
     & i,j,k                                                            &
                   !  loop indices
     &,im_index                                                         &
                   !  internal model index for STASH arrays
     &,item        !  STASH item of diagnostic
      INTEGER :: Errorstatus = 0  ! initial value for error code

      CHARACTER (LEN=80) :: CMessage !  Error message

      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     & work_visc,rneutml

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of header

!
! 1. Initialisation
!
      IF (lhook) CALL dr_hook('DIAGNOSTICS_DIF',zhook_in,zhook_handle)
      im_index    = internal_model_index(atmos_im)
      Cmessage    = ''
!
! 2. Extract diagnostic fields dependent on STASHflags sf
!

! T increment
      item = 181           ! T increment
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

        DO k= 1, tdims%k_end
          DO j= tdims%j_start, tdims%j_end
            DO i= tdims%i_start, tdims%i_end
              T_incr_diagnostic(i,j,k) = (theta_star(i,j,k) -           &
     &                        T_incr_diagnostic(i,j,k))                 &
     &                                  *exner_theta_levels(i,j,k)
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

      ENDIF ! sf(item,sect)

! q  increment
      item = 182           ! q increment
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

        DO k= 1, qdims%k_end
          DO j= qdims%j_start, qdims%j_end
            DO i= qdims%i_start, qdims%i_end
              q_incr_diagnostic(i,j,k) = q_star(i,j,k) -                &
     &                             q_incr_diagnostic(i,j,k)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        q_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      ENDIF ! sf(item,sect)

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
      ALLOCATE (work_visc(tdims%i_start:tdims%i_end,     &
                          tdims%j_start:tdims%j_end,     &
                                      1:tdims%k_end))

      If (.not. L_subfilter_vert .and. .not. L_subfilter_horiz) then
        sf(190,sect)=.false.
        sf(191,sect)=.false.
        sf(192,sect)=.false.
        sf(193,sect)=.false.
        sf(194,sect)=.false.
        sf(195,sect)=.false.
        sf(196,sect)=.false.
        sf(197,sect)=.false.
      Else If (.not. L_subfilter_vert) then
        sf(196,sect)=.false.
        sf(197,sect)=.false.  
      End if

      item = 190           ! momentum viscosity coeff
      If(sf(item,sect).AND.Errorstatus == 0) THEN

        Do k = 1, tdims%k_end
          Do j = tdims%j_start, tdims%j_end
            Do i = tdims%i_start, tdims%i_end
              work_visc(i,j,k) = visc_m(i,j,k)
            End Do
          End Do
        End Do

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_visc,                                                &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If

      item = 191           ! scalar viscosity coeff
      If(sf(item,sect).AND.Errorstatus == 0) THEN

        Do k = 1, tdims%k_end
          Do j = tdims%j_start, tdims%j_end
            Do i = tdims%i_start, tdims%i_end
              work_visc(i,j,k) = visc_h(i,j,k)
            End Do
          End Do
        End Do

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_visc,                                                &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If
      item = 192 ! shear
      If(sf(item,sect).AND.Errorstatus == 0) THEN

        Do k = 1, tdims%k_end
          Do j = tdims%j_start, tdims%j_end
            Do i = tdims%i_start, tdims%i_end
              work_visc(i,j,k) = shear(i,j,k)
            End Do
          End Do
        End Do

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_visc,                                                &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If

      item = 193 ! mixing length
      If(sf(item,sect).AND.Errorstatus == 0) THEN

        ALLOCATE (rneutml(row_length,rows,model_levels))
        DO k=1,model_levels-1
         DO j=1,rows
          DO i=1,row_length
            rneutml(i,j,k) = SQRT(rneutml_sq(i,j,k))
          END DO
         END DO
        END DO
        k=model_levels
        DO j=1,rows
          DO i=1,row_length
            rneutml(i,j,k) = 0.0
          END DO
        END DO
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        rneutml,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If

      DEALLOCATE (work_visc)

! Counter for occurances of local q diffusion
      item = 201           ! local q diffusion at a point
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

! DEPENDS ON: copydiag
        CALL copydiag(STASHwork(si(item,sect,im_index)),                &
     &        w_local_mask,                                             &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      ENDIF ! sf(item,sect)
!
! 3. Error handling
!
      IF(Errorstatus /= 0) THEN

        CALL Ereport(RoutineName,Errorstatus,Cmessage)
      ENDIF

      IF (lhook) CALL dr_hook('DIAGNOSTICS_DIF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Diagnostics_dif
