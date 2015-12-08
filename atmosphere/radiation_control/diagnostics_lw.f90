! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************
!
!+ Subroutine diagnostics_lw
!
! Purpose:
!   Calculates diagnostics and outputs them.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------

      Subroutine diagnostics_lw(                                        &
     &      row_length, rows, model_levels                              & 
     &,     wet_model_levels, ozone_levels                              &
     &,     cloud_levels                                                & 
     &,     n_rows, global_row_length, global_rows                      &
     &,     halo_i, halo_j, off_x, off_y, me                            &
     &,     n_proc, n_procx, n_procy                                    &
     &,     g_rows, g_row_length                                        &
     &,     at_extremity                                                &
     &,     timestep,i_off                                              & 
     &,     T_n, T_inc                                                  & 
     &,     q_n, qcl_n, cf_n, cfl_n                                     &
     &,     T_latest, q_latest, qcl_latest                              &
     &,     cf_latest, cfl_latest                                       &
     &,     surflw, OLR, lw_down                                        &
     &,     T_incr_diagnostic                                           &
     &,     n_channel                                                   &
     &,     ozone                                                       &
     &,     O3_trop_level                                               &
     &,     O3_trop_height                                              &
     &,     T_trop_level                                                &
     &,     T_trop_height                                               &
     &,     LW_incs, LWsea                                              &
     &,     n_aod_wavel                                                 &
     &,     j_lw                                                        &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &      STASHwork)

      Use lw_diag_mod, Only: LW_diag
      Use rad_input_mod, Only: l_rad_perturb
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE solinc_data, ONLY: sky, l_skyview
      USE ereport_mod, ONLY : ereport
      USE Submodel_Mod
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Logical ::  at_extremity(4)  
!        Indicates if this processor is at north,
!        south, east or west of the processor grid.

      Real :: timestep         
!         Atmosphere timestep

      Integer :: row_length       ! number of points on a row
      Integer :: rows             ! number of rows in a theta field
      Integer :: n_rows           ! number of rows in a v field
      Integer :: model_levels     ! number of model levels
      Integer :: cloud_levels     ! number of cloudy levels
      Integer :: wet_model_levels ! number of model levels where moisture
                                  ! variables are held
      Integer :: ozone_levels     ! number of levels where ozone is held
      Integer :: number_format    ! switch controlling number format diagnostics
                                  ! are written out in. See PP_WRITE for details.
      Integer :: model_domain     ! indicator as to model type, ie global, lam
      Integer :: n_channel        ! Number of satellite channels used

      Integer :: global_row_length   !IN. NUMBER OF points on a global row
      Integer :: global_rows         !IN. NUMBER OF global rows
      Integer :: me                  !IN. Processor number
      Integer :: halo_i              !IN. size of large halo in x direction
      Integer :: halo_j              !IN. size of large halo in y direction
      Integer :: off_x               !IN. size of small halo in x direction
      Integer :: off_y               !IN. size of small halo in y direction
      Integer :: n_proc
      Integer :: n_procx
      Integer :: n_procy
      Integer :: g_rows (0:n_proc-1)
      Integer :: g_row_length (0:n_proc-1)

      Integer, Intent(IN) :: i_off
!          Offset for diagnostics
      Integer, Intent(IN) ::  n_aod_wavel
!          Aerosol optical depth diagnostics
      Integer, Intent(IN) :: j_lw
!       Call to LW radiation

!
! Primary Arrays used in all models
!
      Real :: T_n(row_length, rows, model_levels)
      Real :: T_inc(row_length, rows, model_levels)
      Real :: q_n(row_length, rows, wet_model_levels)
      Real :: qcl_n(row_length, rows, wet_model_levels)
      Real :: cf_n(row_length, rows, wet_model_levels)
      Real :: cfl_n(row_length, rows, wet_model_levels)
      Real :: T_latest(row_length, rows, model_levels)
      Real :: q_latest(row_length, rows, wet_model_levels)
      Real :: qcl_latest(row_length, rows, wet_model_levels)
      Real :: cf_latest(row_length, rows, wet_model_levels)
      Real :: cfl_latest(row_length, rows, wet_model_levels)

      Real :: OLR (row_length, rows)
      Real :: lw_down(row_length, rows)
      Real :: surflw (row_length, rows)
      Real :: LWsea(row_length, rows)
      Real :: T_incr_diagnostic(row_length,rows,model_levels)
      Real :: ozone(row_length,rows,ozone_levels)
      Real :: O3_trop_level(row_length,rows)
      Real :: O3_trop_height(row_length,rows)
      Real :: T_trop_level(row_length,rows)
      Real :: T_trop_height(row_length,rows)
      Real :: LW_incs(row_length, rows, 0:model_levels)

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
      
      Real :: STASHwork(*)    ! STASH workspace

! Local array & variables
      
      Real :: work_3d(row_length, rows, model_levels)
      Real :: isccp_dummy_3d(row_length,rows,cloud_levels)

      Integer :: i, j, k
      Integer :: icode                ! Return code  =0 Normal exit  >1 Error
      Integer :: item                 ! STASH item
      Integer, Parameter :: sect=2    ! STASH LW Radiation section

      CHARACTER(LEN=80)  cmessage
      CHARACTER(LEN=*) RoutineName
      Parameter ( RoutineName='diagnostics_lw')

      Integer :: im_index        ! internal model index
      Integer :: ptr_stash       ! Pointer to position in STASH

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('DIAGNOSTICS_LW',zhook_in,zhook_handle)
      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1
! First treat diagnostics that are not contained in the LW_DIAG
! structure (and the ISCCP diagnostics). These can only be output once
! and are valid after the last call to radiation where j_lw==1
! ----------------------------------------------------------------------

  If (j_lw == 1) Then
!
! Pseudo temperature after lw radiation (diagnostic only)
!
      item =   4  ! T + LW T_increment
      If (icode <= 0 .and. sf_calc(item,sect)) Then

        Do k = 1,model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3d(i,j,k) = T_n(i,j,k) + T_inc(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       work_3d,                                                   &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 004)"//cmessage
        End if

      Endif  !  sf_calc(item,sect)

! Skyview factor
    IF (l_skyview) THEN
      item = 101
      IF (sf_calc(item,sect)) THEN
! DEPENDS ON: copydiag
        CALL copydiag (STASHwork(si(item,sect,im_index)), sky,          &
          row_length,rows,0,0,0,0, at_extremity,                        &
          atmos_im,sect,item,icode,cmessage)
        IF (icode > 0) THEN
          WRITE(cmessage,'(A,I3,A)') 'Error in copydiag(item ',item,')'
          GOTO 9999
        END IF
      END IF
    END IF

!
! Increment diagnostics= modified - previous
!
      item = 161  ! temperature increment
      If (icode <= 0 .and. sf_calc(item,sect)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &       T_incr_diagnostic,                                         &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 161)"//cmessage
         End if

      Endif  !  sf_calc(item,sect)

!
! Temperature increment including condensation
!
      item = 181  
      If (icode <= 0 .and. sf_calc(item,sect)) Then

      Do k = 1,model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = T_latest(i,j,k) - T_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &       work_3d,                                                   &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        & 
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 181)"//cmessage
         End if

      Endif  !  sf_calc(item,sect)

!
! Vapour increment
!
      item = 182  
      If (icode <= 0 .and. sf_calc(item,sect)) Then

      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = q_latest(i,j,k) - q_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &      work_3d,                                                    &
     &      row_length,rows,wet_model_levels,0,0,0,0, at_extremity,     &
     &      stlist(1,stindex(1,item,sect,im_index)),len_stlist,         &
     &      stash_levels,num_stash_levels+1,                            &
     &      atmos_im,sect,item,                                         &
     &      icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 182)"//cmessage
         End if

      Endif  !  sf_calc(item,sect)

!
! Liquid water content increment
!
      item = 183  
      If (icode <= 0 .and. sf_calc(item,sect)) Then

      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = qcl_latest(i,j,k) - qcl_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &       work_3d,                                                   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 183)"//cmessage
         End if

      Endif  !  sf_calc(item,sect)

!
! Total cloud fraction increment
!
      item = 192  
      If (icode <= 0 .and. sf_calc(item,sect)) Then

      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = cf_latest(i,j,k) - cf_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           & 
     &       work_3d,                                                   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 192)"//cmessage
         End if

      Endif  !  sf_calc(item,sect)

!
! Liquid cloud fraction increment
!
      item = 193  
      If (icode <= 0 .and. sf_calc(item,sect)) Then

      Do k = 1,wet_model_levels
        Do j = 1,rows
          Do i = 1,row_length
             work_3d(i,j,k) = cfl_latest(i,j,k) - cfl_n(i,j,k)
          Enddo ! i
        Enddo   ! j
      Enddo     ! k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &       work_3d,                                                   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        & 
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 193)"//cmessage
         End if

      Endif  !  sf_calc(item,sect)

!
! Liquid water content increment: positive
!
      item = 194 
      IF (icode <= 0 .AND. sf_calc(item,sect)) THEN

        DO k = 1,wet_model_levels
          DO j = 1,rows
            DO i = 1,row_length
              work_3d(i,j,k) = MAX(0.0,qcl_latest(i,j,k) - qcl_n(i,j,k))
            END DO ! i
          END DO ! j
        END DO ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work_3d,                                                   &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode >  0) THEN
          cmessage=": error in copydiag_3d(item 194)"//cmessage
        END IF 

      END IF !  sf_calc(item,sect)

!
! Liquid water content increment: negative
!
      item = 195 
      IF (icode <= 0 .AND. sf_calc(item,sect)) THEN

        DO k = 1,wet_model_levels
          DO j = 1,rows
            DO i = 1,row_length
              work_3d(i,j,k) = MIN(0.0,qcl_latest(i,j,k) - qcl_n(i,j,k))
            END DO ! i
          END DO ! j
        END DO ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work_3d,                                                   &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

         IF (icode >  0) THEN
            cmessage=": error in copydiag_3d(item 195)"//cmessage
         END IF

      END IF !  sf_calc(item,sect)

!
! Liquid cloud fraction increment: positive
!
      item = 198  
      IF (icode <= 0 .AND. sf_calc(item,sect)) THEN

        DO k = 1,wet_model_levels
          DO j = 1,rows
            DO i = 1,row_length
              work_3d(i,j,k) = MAX(0.0,cfl_latest(i,j,k) - cfl_n(i,j,k))
            END DO ! i
          END DO  ! j
        END DO   ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work_3d,                                                   &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        & 
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode >  0) THEN
            cmessage=": error in copydiag_3d(item 198)"//cmessage
        END IF

      END IF !  sf_calc(item,sect)

!
! Liquid cloud fraction increment: negative
!
      item = 199  
      IF (icode <= 0 .AND. sf_calc(item,sect)) THEN

        DO k = 1,wet_model_levels
          DO j = 1,rows
            DO i = 1,row_length
              work_3d(i,j,k) = MIN(0.0,cfl_latest(i,j,k) - cfl_n(i,j,k))
            END DO ! i
          END DO ! j
        END DO ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d (stashwork(si(item,sect,im_index)),            &
             work_3d,                                                   &
             row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,        & 
             stash_levels,num_stash_levels+1,                           &
             atmos_im,sect,item,                                        &
             icode,cmessage)

        IF (icode >  0) THEN
          cmessage=": error in copydiag_3d(item 199)"//cmessage
        END IF

      END IF !  sf_calc(item,sect)

!
! Surflw
!
      If (sf_calc(201,2)) Then
      
! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(201,2,im_index)),surflw,           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,201,                                            & 
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if

!
! OLR
!
      If (sf_calc(205,2)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(205,2,im_index)),OLR,              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,205,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 205)"
            goto 9999
         End if

      End if

!
! Surface Down Flux.
!
        If (sf_calc(207,2)) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(207,2,im_index)), lw_down,        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,207,                                         &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 207)"
            goto 9999
          End if

        End if

!
! LWsea : 'net down lw rad flux: open sea'
!
      If (sf_calc(203,2)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(203,2,im_index)),                  &
     &       LWsea,                                                     &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,203,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 203 = LWsea)"
            goto 9999
         End if

      End if

!
! Ozone & ozone Troposphere Diagnostics
!
      If (sf_calc(260,2)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(260,2,im_index)),               &
     &       ozone,                                                     &
     &       row_length,rows,ozone_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,260,2,im_index)),len_stlist,            & 
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,260,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 260)"//cmessage
            goto 9999
         End if

      End if

      If (sf_calc(280,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(280,2,im_index)),                  &
     &       O3_trop_level,                                             & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,280,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 280)"//cmessage
         End if
      End if

      If (sf_calc(281,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(281,2,im_index)),                  & 
     &       O3_trop_height,                                            &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,281,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 281)"//cmessage
         End if
      End if

      If (sf_calc(282,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(282,2,im_index)),                  &
     &       T_trop_level,                                              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,282,                                            &
     &       icode,cmessage)      

         If (icode  >   0) then
            cmessage=": error in copydiag(item 282)"//cmessage
         End if
      End if

      If (sf_calc(283,2)) then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(283,2,im_index)),                  &
     &       T_trop_height,                                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,283,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag(item 283)"//cmessage
         End if
      End if

!
! LW Heating = Lw radiation temp. incr. per timestep / timestep
!
      item = 232  ! LW heating
      If (icode <= 0 .and. sf_calc(item,sect)) Then

        Do k = 1,model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3d(i,j,k) = LW_incs(i,j,k)/timestep
            Enddo ! i
          Enddo   ! j
        Enddo     ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            & 
     &       work_3d,                                                   &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        & 
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        & 
     &       icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 232)"//cmessage
        End if

      Endif  !  sf_calc(item,sect)

!
! ISCCP Weights: Stash: sf_calc(269,2)
!
      If (sf_calc(269,2)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(269,2,im_index)),            &
     &       LW_diag(1)%isccp_weights,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,269,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 269)"
            goto 9999
         End if

      End if

! Copy ISCCP diagnostics by looping over 7 levels in call to copydiag.
! This is because copydiag_3d cannot handle ISCCP levels.

!
! ISCCP CF: Stash: sf_calc(270,2)
!
      If (sf_calc(270,2)) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(270,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                &
     &       LW_diag(1)%isccp_cf(1,1,k),                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,270,                                      &
     &       icode,cmessage)
        Enddo


         If (icode  >   0) then
            cmessage=": error in copydiag( item 270)"//cmessage
            goto 9999
         End if

      End if

!
! isccp_cf_tau_0_to_p3. Stash: sf_calc(271,2)
!
      If (sf_calc(271,2)) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(271,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                & 
     &       LW_diag(1)%isccp_cf_tau_0_to_p3(1,1,k),                 & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,271,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 271)"//cmessage
            goto 9999
         End if

      End if

!
! Isccp_cf_tau_p3_to_1p3. Stash: sf_calc(272,2)
!
      If (sf_calc(272,2)) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(272,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                & 
     &       LW_diag(1)%isccp_cf_tau_p3_to_1p3(1,1,k),               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,272,                                      & 
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 272)"//cmessage
            goto 9999
         End if

      End if

!
! Isccp_cf_tau_1p3_to_3p6. Stash: sf_calc(273,2)
!
      If (sf_calc(273,2)) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(273,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                & 
     &       LW_diag(1)%isccp_cf_tau_1p3_to_3p6(1,1,k),              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,273,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 273)"//cmessage
            goto 9999
         End if

      End if

!
! Isccp_cf_tau_3p6_to_9p4. Stash: sf_calc(274,2)
!
      If (sf_calc(274,2)) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(274,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                & 
     &       LW_diag(1)%isccp_cf_tau_3p6_to_9p4(1,1,k),              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,274,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 274)"//cmessage
            goto 9999
         End if

      End if

!
! Isccp_cf_tau_9p4_to_23. Stash: sf_calc(275,2)
!
      If (sf_calc(275,2)) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(275,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                & 
     &       LW_diag(1)%isccp_cf_tau_9p4_to_23(1,1,k),               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,275,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 275)"//cmessage
            goto 9999
         End if

      End if

!
!   isccp_cf_tau_23_to_60. Stash: sf_calc(276,2)
!
      If (sf_calc(276,2)) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(276,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                &
     &       LW_diag(1)%isccp_cf_tau_23_to_60(1,1,k),                &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,276,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 276)"//cmessage
            goto 9999
         End if

      End if

!
! Isccp_cf_tau_ge_60. Stash: sf_calc(277,2)
!
      If (sf_calc(277,2)) Then

        Do k = 1,7
        
! DEPENDS ON: copydiag  
         Call copydiag (STASHwork(si(277,2,im_index)              &
     &        +(row_length*rows*(k-1))),                                &
     &       LW_diag(1)%isccp_cf_tau_ge_60(1,1,k),                   &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,277,                                      &
     &       icode,cmessage)
        Enddo

         If (icode  >   0) then
            cmessage=": error in copydiag( item 277)"//cmessage
            goto 9999
         End if

      End if

!
! Mean Albedo Cloud. Stash: sf_calc(290,2)
!
      If (sf_calc(290,2)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(290,2,im_index)),            &
     &        LW_diag(1)%meanalbedocld,                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,2,290,                                     &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 290)"
            goto 9999
         End if

      End if

!
! Mean Tau Cloud. Stash: sf_calc(291,2)
!
      If (sf_calc(291,2)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(291,2,im_index)),            &
     &       LW_diag(1)%meantaucld,                                  &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,291,                                      & 
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 291)"
            goto 9999
         End if

      End if

!
! Mean Top: Stash: sf_calc(292,2)
!
      If (sf_calc(292,2)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(292,2,im_index)),            & 
     &       LW_diag(1)%meanptop,                                    &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,292,                                      & 
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 292)"
            goto 9999
         End if

      End if

!
! Total Cloud Area
!
      If (sf_calc(293,2)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(293,2,im_index)),            &
     &       LW_diag(1)%totalcldarea,                                &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,293,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 293)"
            goto 9999
         End if

      End if

  End if ! j_lw == 1


! ----------------------------------------------------------------------
! Section 2
! Now treat diagnostics contained in structure LW_DIAG. These may be 
! available for each of the radiation calls dependent on the scheme.
! ----------------------------------------------------------------------

  If (l_rad_perturb.and.(j_lw == 1)) Then

! ----------------------------------------------------------------------
! Section 2.1
! For the incremental time-stepping scheme many of the diagnostics are
! not calculated on the "cloud only" radiation calls. For these the 
! diagnostics from the last full radiation call are used.
! ----------------------------------------------------------------------

!   Clear-sky upward flux on levels. Stash: sf_calc(219,2)
    IF (sf_calc(219,2)) THEN
      CALL copydiag_3d(STASHwork(si(219,2,im_index)),                   &
           LW_diag(1)%flux_up_clear,                                    &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,219,2,im_index)),len_stlist,              &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,2,219,                                              &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 219)"
        GOTO 9999
      END IF
    END IF

!   Clear-sky downward flux on levels. Stash: sf_calc(220,2)
    IF (sf_calc(220,2)) THEN
      CALL copydiag_3d(STASHwork(si(220,2,im_index)),                   &
           LW_diag(1)%flux_down_clear,                                  &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,220,2,im_index)),len_stlist,              &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,2,220,                                              &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 220)"
        GOTO 9999
      END IF
    END IF

!
! Clear Olr. Stash: sf_calc(206,2)
!
        If (sf_calc(206,2)) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(206,2,im_index)),                 &
     &       LW_diag(1)%clear_olr,                                      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,206,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 206)"
            goto 9999
          End if

        End if

!
! Clear Sky Surface Down Flux. Stash: sf_calc(208,2) 
!
        If (sf_calc(208,2)) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(208,2,im_index)),                 &
     &       LW_diag(1)%surf_down_clr,                                  & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,208,                                         &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
          End if

        End if

!
! Clear-Sky heating Rates. Stash: sf_calc(233,2)
!
        If (sf_calc(233,2)) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(233,2,im_index)),              &
     &       LW_diag(1)%clear_hr,                                       & 
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,233,2,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,233,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 233)"//cmessage
            goto 9999
          End if

        End if

!
! Cloud Absorptivity. Stash: sf_calc(262,2)
!
        If (sf_calc(262,2)) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(262,2,im_index)),              &
     &       LW_diag(1)%cloud_absorptivity,                             &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,262,2,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,262,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 262)"//cmessage
            goto 9999
          End if

        End if

!
! Cloud Weight Absorptivity. Stash: sf_calc(263,2)
!
        If (sf_calc(263,2)) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(263,2,im_index)),              &
     &       LW_diag(1)%cloud_weight_absorptivity,                      &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,263,2,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,263,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 263)"//cmessage
            goto 9999
          End if

        End if

!
! LS. Cloud Absorptivity. Stash: sf_calc(264,2)
!
        If (sf_calc(264,2)) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(264,2,im_index)),              & 
     &       LW_diag(1)%ls_cloud_absorptivity,                          &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,264,2,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,264,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 264)"//cmessage
            goto 9999
          End if

        End if

!
! Ls. Cloud Weight Absorptivity. Stash: sf_calc(265,2)
!
        If (sf_calc(265,2)) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(265,2,im_index)),              & 
     &       LW_diag(1)%ls_cloud_weight_absorptivity,                   &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,265,2,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,265,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 265)"//cmessage
            goto 9999
          End if

        End if

!
! Convective Cloud Absorptivity. Stash: sf_calc(266,2)
!
        If (sf_calc(266,2)) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(266,2,im_index)),              &
     &       LW_diag(1)%cnv_cloud_absorptivity,                         &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,266,2,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,266,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 266)"//cmessage
            goto 9999
          End if

        End if

!
! Convective Cloud Weight Absorptivity. Stash: sf_calc(267,2)
!
        If (sf_calc(267,2)) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(267,2,im_index)),              & 
     &       LW_diag(1)%cnv_cloud_weight_absorptivity,                  &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,267,2,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,267,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 267)"//cmessage
            goto 9999
          End if

        End if

!   Aerosol optical depth diagnostics
!   (loop on wavelength)

!
! AOD Sulphate. Stash: sf_calc(284,2)
!
      If (sf_calc(284,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(284,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(1)%aod_sulphate(1,1,k),                            & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,284,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 284)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Dust. Stash: sf_calc(285,2)
!
      If (sf_calc(285,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(285,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(1)%aod_dust(1,1,k),                                & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,285,                                            & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 285)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Seasalt. Stash: sf_calc(286,2)
!
      If (sf_calc(286,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(286,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(1)%aod_seasalt(1,1,k),                             & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,286,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 286)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Soot. Stash: sf_calc(287,2)
!
      If (sf_calc(287,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(287,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(1)%aod_soot(1,1,k),                                &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,287,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 287)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Biomass. Stash: sf_calc(288,2)
!
      If (sf_calc(288,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(288,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(1)%aod_biomass(1,1,k),                             & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,288,                                            & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 288)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
      
!
! AOD Biogenic. Stash: sf_calc(289,2)
!
      If (sf_calc(289,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(289,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(1)%aod_biogenic(1,1,k),                            &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,289,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 289)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Fossil-fuel organic carbon. Stash: sf_calc(295,2)
!
      If (sf_calc(295,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(295,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(1)%aod_ocff(1,1,k),                                &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,295,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 295)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Delta aerosol. Stash: sf_calc(296,2)
!
      If (sf_calc(296,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(296,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(1)%aod_delta(1,1,k),                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,296,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 296)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Nitrate. Stash: sf_calc(297,2)
!
      If (sf_calc(297,2)) Then
        Do k = 1, n_aod_wavel

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(297,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(1)%aod_nitrate(1,1,k),                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,297,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 297)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! Total AOD in Radiaiton. Stash: sf_calc(298,2)
!
      If (sf_calc(298,2)) Then
        Do k = 1, n_aod_wavel

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(298,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(1)%aod_total_radn(1,1,k),                          &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,298,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 298)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
!
! Angstrom Exponent of the Total AOD in Radiation: Stash: sf_calc(299,2) 
!
      If (sf_calc(299,2)) Then
        Do k = 1, n_aod_wavel

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(299,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(1)%angst_total_radn(1,1,k),                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,299,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 299)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

! 
! UKCA AOD Aitken soluble. Stash: sf_calc(300,2) 
! 
      item = 300
      IF (sf_calc(item,2)) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(1)%aod_ukca_ait_sol(1,1,k),                        & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 300)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF  
! 
! UKCA AOD accum. soluble. Stash: sf_calc(301,2) 
! 
      item = 301
      IF (sf_calc(item,2)) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(1)%aod_ukca_acc_sol(1,1,k),                        & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 301)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF 
  
! 
! UKCA AOD coarse soluble. Stash: sf_calc(302,2) 
! 
      item = 302
      IF (sf_calc(item,2)) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(1)%aod_ukca_cor_sol(1,1,k),                        & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 302)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF 
  
! 
! UKCA AOD Aitken insoluble. Stash: sf_calc(303,2) 
! 
      item = 303
      IF (sf_calc(item,2)) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(1)%aod_ukca_ait_ins(1,1,k),                        & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 303)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF 
  
! 
! UKCA AOD accum. insoluble. Stash: sf_calc(304,2) 
! 
      item = 304
      IF (sf_calc(item,2)) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(1)%aod_ukca_acc_ins(1,1,k),                        & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 304)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF 
  
! 
! UKCA AOD coarse insoluble. Stash: sf_calc(305,2) 
! 
      item = 305
      IF (sf_calc(item,2)) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(1)%aod_ukca_cor_ins(1,1,k),                        & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 305)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF 
!
! Prognositc Sulphate AOD. Stash: sf_calc(421,2)
!
      If (sf_calc(421,2)) Then
        Do k = 1, n_aod_wavel

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(421,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(1)%aod_prog_sulphate(1,1,k),                       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,421,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 421)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! Prognostic Dust AOD. Stash: sf_calc(422,2)
!
      If (sf_calc(422,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(422,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(1)%aod_prog_dust(1,1,k),                           & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,422,                                            & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 422)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
!
! Diagnosed Seasalt AOD. Stash: sf_calc(423,2)
!
      If (sf_calc(423,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(423,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(1)%aod_prog_seasalt(1,1,k),                        & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,423,                                            & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 423)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
!
! Prognostic Soot AOD. Stash: sf_calc(424,2)
!
      If (sf_calc(424,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(424,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(1)%aod_prog_soot(1,1,k),                           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,424,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 424)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
!
! Prognostic Biomass AOD. Stash: sf_calc(425,2)
!
      If (sf_calc(425,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(425,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(1)%aod_prog_biomass(1,1,k),                        & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,425,                                            & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 425)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
!
! AOD Fossil-fuel organic carbon. Stash: sf_calc(426,2)
!
      If (sf_calc(426,2)) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(426,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(1)%aod_prog_ocff(1,1,k),                           &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,426,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 426)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
!
! Prognostic Nitrate AOD. Stash: sf_calc(427,2)
!
      If (sf_calc(427,2)) Then
        Do k = 1, n_aod_wavel

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(427,2,im_index)                   &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(1)%aod_prog_nitrate(1,1,k),                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,427,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 427)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

  Else ! l_rad_perturb.and.(j_lw==1)

! ----------------------------------------------------------------------
! Section 2.2
! The following diagnostics are available on "cloud only" radiation
! calls for the incremental time-stepping scheme, and all calls for
! other schemes.
! ----------------------------------------------------------------------

!   Upward flux on levels. Stash: sf_calc(217,2)
    IF (LW_diag(j_lw)%l_flux_up) THEN
      CALL copydiag_3d(STASHwork(si(217+i_off,2,im_index)),             &
           LW_diag(j_lw)%flux_up,                                       &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,217+i_off,2,im_index)),len_stlist,        &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,2,217+i_off,                                        &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 217)"
        GOTO 9999
      END IF
    END IF

!   Downward flux on levels. Stash: sf_calc(218,2)
    IF (LW_diag(j_lw)%l_flux_down) THEN
      CALL copydiag_3d(STASHwork(si(218+i_off,2,im_index)),             &
           LW_diag(j_lw)%flux_down,                                     &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,218+i_off,2,im_index)),len_stlist,        &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,2,218+i_off,                                        &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 218)"
        GOTO 9999
      END IF
    END IF

!
! Net LW Flux at Tropopause. Stash: sf_calc(237,2)
!
        If (LW_diag(j_lw)%L_net_flux_trop) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(237+i_off,2,im_index)),           &
     &       LW_diag(j_lw)%net_flux_trop,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,237+i_off,                                      &
     &       icode,cmessage) 

          If (icode  >   0) then
            cmessage="Error in copydiag( item 237)"
            goto 9999
          End if

        End if

!
! LW Down Flux at Tropopause. Stash: sf_calc(238,2)
!
        If (LW_diag(j_lw)%L_down_flux_trop) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(238+i_off,2,im_index)),           &
     &       LW_diag(j_lw)%down_flux_trop,                              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,238+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 238)"
            goto 9999
          End if
        End if

!
! Top of the Atmosphere Radiances. Stash: sf_calc(297+i,2)
!
        If (LW_diag(j_lw)%L_toa_radiance) Then

           Do i=1, n_channel

             ptr_stash = si(297+i_off,2,im_index)                      &
     &         + (i-1) * row_length * rows
     
! DEPENDS ON: copydiag     
             Call copydiag (STASHwork(ptr_stash),                      &
     &           LW_diag(j_lw)%toa_radiance(1, 1, i),                  &
     &           row_length,rows,0,0,0,0, at_extremity,                &
     &           atmos_im,2,297+i_off,                                 &
     &           icode,cmessage)

             If (icode  >   0) then
                cmessage="Error in copydiag( item 297+i_off)"
                goto 9999
             End if
           Enddo

        End if

!
! Total Cloud Cover. Stash: sf_calc(204,2)
!
      If (LW_diag(j_lw)%L_total_cloud_cover) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(204+i_off,2,im_index)),            &
     &       LW_diag(j_lw)%total_cloud_cover,                           & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,204+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 204)"
            goto 9999
         End if

      End if

!
! Total Cloud on Levels. Stash: sf_calc(261,2)
!
      If (LW_diag(j_lw)%L_total_cloud_on_levels) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(261+i_off,2,im_index)),         &
     &       LW_diag(j_lw)%total_cloud_on_levels,                       &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,261+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,261+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 261)"//cmessage
            goto 9999
         End if

      End if

!
!   Cloud water mixing ratios
!
      item = 308+i_off  ! LS cloud liquid water mixing ratio
      If (LW_diag(j_lw)%L_ls_qcl_rad) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%ls_qcl_rad,                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 308)"//cmessage
            goto 9999
         End if
      Endif

      item = 309+i_off  ! LS cloud ice water mixing ratio
      If (LW_diag(j_lw)%L_ls_qcf_rad) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%ls_qcf_rad,                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 309)"//cmessage
            goto 9999
         End if
      Endif

      item = 310+i_off  ! Convective cloud liquid water mixing ratio
      If (LW_diag(j_lw)%L_cc_qcl_rad) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%cc_qcl_rad,                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 310)"//cmessage
            goto 9999
         End if
      Endif

      item = 311+i_off  ! Convective cloud ice water mixing ratio
      If (LW_diag(j_lw)%L_cc_qcf_rad) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%cc_qcf_rad,                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 311)"//cmessage
            goto 9999
         End if
      Endif

!
!   Cloud amounts
!
      item = 312+i_off  ! LS cloud fraction of grdbox seen by radiation.
                        ! Liquid
      If (LW_diag(j_lw)%L_ls_cl_rad) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%ls_cl_rad,                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 312)"//cmessage
            goto 9999
         End if
      Endif

      item = 313+i_off  ! LS cloud fraction of grdbox seen by radiation.
                        ! Ice
      If (LW_diag(j_lw)%L_ls_cf_rad) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%ls_cf_rad,                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 313)"//cmessage
            goto 9999
         End if
      Endif

      item = 314+i_off  ! CONV cloud fraction of grdbox seen by radiation.
                        ! Liquid
      If (LW_diag(j_lw)%L_cc_cl_rad) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%cc_cl_rad,                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 314)"//cmessage
            goto 9999
         End if
      Endif

      item = 315+i_off  ! CONV cloud fraction of grdbox seen by radiation.
                        ! Ice
      If (LW_diag(j_lw)%L_cc_cf_rad) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              LW_diag(j_lw)%cc_cf_rad,                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
         If (icode.gt.0) Then
            cmessage=": error in copydiag_3d(item 315)"//cmessage
            goto 9999
         End if
      Endif

  End if  ! l_rad_perturb.and.(j_lw==1) - end of Else block



  If (.not.l_rad_perturb) Then

! ----------------------------------------------------------------------
! Section 2.3
! The following diagnostics have been output in section 2.1 for the
! incremental time-stepping scheme. They are now output for other
! schemes on any radiation call.
! ----------------------------------------------------------------------

!   Clear-sky upward flux on levels. Stash: sf_calc(219,2)
    IF (LW_diag(j_lw)%l_flux_up_clear) THEN
      CALL copydiag_3d(STASHwork(si(219+i_off,2,im_index)),             &
           LW_diag(j_lw)%flux_up_clear,                                 &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,219+i_off,2,im_index)),len_stlist,        &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,2,219+i_off,                                        &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 219)"
        GOTO 9999
      END IF
    END IF

!   Clear-sky downward flux on levels. Stash: sf_calc(220,2)
    IF (LW_diag(j_lw)%l_flux_down_clear) THEN
      CALL copydiag_3d(STASHwork(si(220+i_off,2,im_index)),             &
           LW_diag(j_lw)%flux_down_clear,                               &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,220+i_off,2,im_index)),len_stlist,        &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,2,220+i_off,                                        &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 220)"
        GOTO 9999
      END IF
    END IF

!
! Clear Olr. Stash: sf_calc(206,2)
!
        If (LW_diag(j_lw)%L_clear_olr) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(206+i_off,2,im_index)),           &
     &       LW_diag(j_lw)%clear_olr,                                   &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,206+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 206)"
            goto 9999
          End if

        End if

!
! Clear Sky Surface Down Flux. Stash: sf_calc(208,2) 
!
        If (icode <= 0 .and.LW_diag(j_lw)%L_surf_down_clr) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(208+i_off,2,im_index)),           &
     &       LW_diag(j_lw)%surf_down_clr,                               & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,sect,208+i_off,                                   &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
          End if

        End if

!
! Clear-Sky heating Rates. Stash: sf_calc(233,2)
!
        If (LW_diag(j_lw)%L_clear_hr) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(233+i_off,2,im_index)),        &
     &       LW_diag(j_lw)%clear_hr,                                    & 
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,233+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,233+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 233)"//cmessage
            goto 9999
          End if

        End if

!
! Cloud Absorptivity. Stash: sf_calc(262,2)
!
        If (LW_diag(j_lw)%L_cloud_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(262+i_off,2,im_index)),        &
     &       LW_diag(j_lw)%cloud_absorptivity,                          &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,262+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,262+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 262)"//cmessage
            goto 9999
          End if

        End if

!
! Cloud Weight Absorptivity. Stash: sf_calc(263,2)
!
        If (LW_diag(j_lw)%L_cloud_weight_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(263+i_off,2,im_index)),        &
     &       LW_diag(j_lw)%cloud_weight_absorptivity,                   &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,263+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,263+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 263)"//cmessage
            goto 9999
          End if

        End if

!
! LS. Cloud Absorptivity. Stash: sf_calc(264,2)
!
        If (LW_diag(j_lw)%L_ls_cloud_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(264+i_off,2,im_index)),        & 
     &       LW_diag(j_lw)%ls_cloud_absorptivity,                       &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,264+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,264+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 264)"//cmessage
            goto 9999
          End if

        End if

!
! Ls. Cloud Weight Absorptivity. Stash: sf_calc(265,2)
!
        If (LW_diag(j_lw)%L_ls_cloud_weight_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(265+i_off,2,im_index)),        & 
     &       LW_diag(j_lw)%ls_cloud_weight_absorptivity,                &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,265+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,265+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 265)"//cmessage
            goto 9999
          End if

        End if

!
! Convective Cloud Absorptivity. Stash: sf_calc(266,2)
!
        If (LW_diag(j_lw)%L_cnv_cloud_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(266+i_off,2,im_index)),        &
     &       LW_diag(j_lw)%cnv_cloud_absorptivity,                      &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,266+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,266+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 266)"//cmessage
            goto 9999
          End if

        End if

!
! Convective Cloud Weight Absorptivity. Stash: sf_calc(267,2)
!
        If (LW_diag(j_lw)%L_cnv_cloud_weight_absorptivity) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(267+i_off,2,im_index)),        & 
     &       LW_diag(j_lw)%cnv_cloud_weight_absorptivity,               &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,267+i_off,2,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,2,267+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag_3d( item 267)"//cmessage
            goto 9999
          End if

        End if

!   Aerosol optical depth diagnostics
!   (loop on wavelength)

!
! AOD Sulphate. Stash: sf_calc(284,2)
!
      If (LW_diag(j_lw)%L_aod_sulphate) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(284+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_sulphate(1,1,k),                         & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,284+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 284)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Dust. Stash: sf_calc(285,2)
!
      If (LW_diag(j_lw)%L_aod_dust) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(285+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_dust(1,1,k),                             & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,285+i_off,                                      & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 285)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Seasalt. Stash: sf_calc(286,2)
!
      If (LW_diag(j_lw)%L_aod_seasalt) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(286+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(j_lw)%aod_seasalt(1,1,k),                          & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,286+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 286)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Soot. Stash: sf_calc(287,2)
!
      If (LW_diag(j_lw)%L_aod_soot) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(287+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_soot(1,1,k),                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,287+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 287)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Biomass. Stash: sf_calc(288,2)
!
      If (LW_diag(j_lw)%L_aod_biomass) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(288+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_biomass(1,1,k),                          & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,288+i_off,                                      & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 288)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
      
!
! AOD Biogenic. Stash: sf_calc(289,2)
!
      If (LW_diag(j_lw)%L_aod_biogenic) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(289+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(j_lw)%aod_biogenic(1,1,k),                         &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,289+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 289)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Fossil-fuel organic carbon. Stash: sf_calc(295,2)
!
      If (LW_diag(j_lw)%L_aod_ocff) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(295+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(j_lw)%aod_ocff(1,1,k),                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,295+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 295)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! AOD Delta aerosol. Stash: sf_calc(296,2)
!
      If (LW_diag(j_lw)%L_aod_delta) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(296+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(j_lw)%aod_delta(1,1,k),                            &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,296+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 296)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
! 
! AOD Nitrate. Stash: sf_calc(297,2) 
! 
      If (LW_diag(j_lw)%L_aod_nitrate) Then 
        Do k = 1, n_aod_wavel 
         
! DEPENDS ON: copydiag   
          Call copydiag (STASHwork(si(297+i_off,2,im_index)             & 
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_nitrate(1,1,k),                          & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &  
     &       atmos_im,2,297+i_off,                                      & 
     &       icode,cmessage) 
 
          If (icode  >   0) then 
            cmessage=": error in copydiag( item 297)"//cmessage 
            goto 9999 
          End if 
 
        Enddo ! k 
      End if 
! 
! Total AOD in Radiation. Stash: sf_calc(298,2) 
! 
      If (LW_diag(j_lw)%L_aod_total_radn) Then 
        Do k = 1, n_aod_wavel 
         
! DEPENDS ON: copydiag   
          Call copydiag (STASHwork(si(298+i_off,2,im_index)             & 
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_total_radn(1,1,k),                       & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &  
     &       atmos_im,2,298+i_off,                                      & 
     &       icode,cmessage) 
 
          If (icode  >   0) then 
            cmessage=": error in copydiag( item 298)"//cmessage 
            goto 9999 
          End if 
 
        Enddo ! k 
      End if
! 
! Angstrom Exponent of the Total AOD in Radiation: Stash: sf_calc(299,2) 
! 
      If (LW_diag(j_lw)%L_angst_total_radn) Then 
        Do k = 1, n_aod_wavel 
         
! DEPENDS ON: copydiag   
          Call copydiag (STASHwork(si(299+i_off,2,im_index)             & 
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%angst_total_radn(1,1,k),                     & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &  
     &       atmos_im,2,299+i_off,                                      & 
     &       icode,cmessage) 
 
          If (icode  >   0) then 
            cmessage=": error in copydiag( item 299)"//cmessage 
            goto 9999 
          End if 
 
        Enddo ! k 
      End if
! 
! UKCA AOD Aitken soluble. Stash: sf_calc(300,2) 
! 
      item = 300 + i_off
      IF (LW_diag(j_lw)%L_aod_ukca_ait_sol) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(j_lw)%aod_ukca_ait_sol(1,1,k),                     & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 300)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF  
! 
! UKCA AOD accum. soluble. Stash: sf_calc(301,2) 
! 
      item = 301 + i_off
      IF (LW_diag(j_lw)%L_aod_ukca_acc_sol) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(j_lw)%aod_ukca_acc_sol(1,1,k),                     & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 301)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF 
  
! 
! UKCA AOD coarse soluble. Stash: sf_calc(302,2) 
! 
      item = 302 + i_off
      IF (LW_diag(j_lw)%L_aod_ukca_cor_sol) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(j_lw)%aod_ukca_cor_sol(1,1,k),                     & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 302)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF 
  
! 
! UKCA AOD Aitken insoluble. Stash: sf_calc(303,2) 
! 
      item = 303 + i_off
      IF (LW_diag(j_lw)%L_aod_ukca_ait_ins) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(j_lw)%aod_ukca_ait_ins(1,1,k),                     & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 303)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF 
  
! 
! UKCA AOD accum. insoluble. Stash: sf_calc(304,2) 
! 
      item = 304 + i_off
      IF (LW_diag(j_lw)%L_aod_ukca_acc_ins) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(j_lw)%aod_ukca_acc_ins(1,1,k),                     & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 304)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF 
  
! 
! UKCA AOD coarse insoluble. Stash: sf_calc(305,2) 
! 
      item = 305 + i_off
      IF (LW_diag(j_lw)%L_aod_ukca_cor_ins) THEN 
        DO k = 1, n_aod_wavel 
  
! DEPENDS ON: copydiag  
          CALL copydiag (STASHwork(si(item,2,im_index)                  & 
             +(row_length*rows*(k-1))),                                 & 
             LW_diag(j_lw)%aod_ukca_cor_ins(1,1,k),                     & 
             row_length,rows,0,0,0,0, at_extremity,                     &  
             atmos_im,2,item,                                           & 
             icode,cmessage) 
 
          IF (icode  >   0) THEN 
            cmessage=": error in copydiag( item 305)"//cmessage 
            GOTO 9999 
          END IF 
 
        END DO ! k 
 
      END IF 
 
!
! Prognositc Sulphate AOD. Stash: sf_calc(421,2)
!
      If (LW_diag(j_lw)%L_aod_prog_sulphate) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(421+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_prog_sulphate(1,1,k),                    & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,421+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 421)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! Prognostic Dust AOD. Stash: sf_calc(422,2)
!
      If (LW_diag(j_lw)%L_aod_prog_dust) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(422+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_prog_dust(1,1,k),                        & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,422+i_off,                                      & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 422)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! Diagnosed Seasalt AOD. Stash: sf_calc(423,2)
!
      If (LW_diag(j_lw)%L_aod_prog_seasalt) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(423+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(j_lw)%aod_prog_seasalt(1,1,k),                     & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,423+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 423)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! Prognostic Soot AOD. Stash: sf_calc(424,2)
!
      If (LW_diag(j_lw)%L_aod_prog_soot) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(424+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_prog_soot(1,1,k),                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,424+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 424)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

!
! Prognostic Biomass AOD. Stash: sf_calc(425,2)
!
      If (LW_diag(j_lw)%L_aod_prog_biomass) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(425+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_prog_biomass(1,1,k),                     & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,2,425+i_off,                                      & 
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 425)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if
      

!
! AOD Fossil-fuel organic carbon. Stash: sf_calc(426,2)
!
      If (LW_diag(j_lw)%L_aod_prog_ocff) Then
        Do k = 1, n_aod_wavel
        
! DEPENDS ON: copydiag  
          Call copydiag (STASHwork(si(426+i_off,2,im_index)             &
     &       +(row_length*rows*(k-1))),                                 &
     &       LW_diag(j_lw)%aod_prog_ocff(1,1,k),                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,2,426+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage=": error in copydiag( item 426)"//cmessage
            goto 9999
          End if

        Enddo ! k
      End if

! 
! Prognostic Nitrate AOD. Stash: sf_calc(427,2)
! 
      If (LW_diag(j_lw)%L_aod_prog_nitrate) Then 
        Do k = 1, n_aod_wavel 
         
! DEPENDS ON: copydiag   
          Call copydiag (STASHwork(si(427+i_off,2,im_index)             & 
     &       +(row_length*rows*(k-1))),                                 & 
     &       LW_diag(j_lw)%aod_prog_nitrate(1,1,k),                     & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &  
     &       atmos_im,2,427+i_off,                                      & 
     &       icode,cmessage) 
 
          If (icode  >   0) then 
            cmessage=": error in copydiag( item 427)"//cmessage 
            goto 9999 
          End if 
 
        Enddo ! k 
      End if 

  End if ! .not.l_rad_perturb

 9999 continue  ! exit point on error
      If(icode /= 0) Then

        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      IF (lhook) CALL dr_hook('DIAGNOSTICS_LW',zhook_out,zhook_handle)
      RETURN
    End Subroutine diagnostics_lw
