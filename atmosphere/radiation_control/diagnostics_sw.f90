! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************
!
!+ Subroutine diagnostics_sw
!
! Purpose:
!   Calculates diagnostics and outputs them.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------

      Subroutine diagnostics_sw(                                        &
     &      row_length, rows, model_levels                              &
      ,     wet_model_levels, cloud_levels, ntiles                      &
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
     &,     surfsw, itoasw, surfsw_cor, toasw_cor                       &
     &,     surfdir_cor, surfdif_cor                                    &
     &,     SWsea, flux_below_690nm_surf                                &
     &,     photosynth_act_rad, flxdirparsurf                           &
     &,     f_orog                                                      &
     &,     slope_aspect, slope_angle                                   &
     &,     sw_net_land,sw_net_sice                                     &
     &,     T_incr_diagnostic                                           &
     &,     n_channel                                                   &
     &,     sea_salt_film, sea_salt_jet                                 &
     &,     salt_dim1, salt_dim2, salt_dim3                             &
     &,     j_sw                                                        &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &      STASHwork )

      Use sw_diag_mod, Only: SW_diag
      USE spec_sw_lw, Only: sw_spectrum
      Use rad_input_mod, Only: l_rad_perturb
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE solinc_data, ONLY: horiz_ang, n_horiz_ang, l_skyview

      USE ereport_mod, ONLY : ereport
      USE Submodel_Mod
      IMPLICIT NONE
      
!
! Arguments with Intent IN. ie: Input variables.
!

      Logical :: at_extremity(4)  
!       Indicates if this processor is at north,
!       south, east or west of the processor grid

      Integer :: row_length       ! number of points on a row
      Integer :: rows             ! number of rows in a theta field
      Integer :: n_rows           ! number of rows in a v field
      Integer :: model_levels     ! number of model levels
      Integer :: cloud_levels     ! number of cloudy levels
      INTEGER :: ntiles           ! number of land surface tiles 
      Integer :: wet_model_levels ! number of model levels where moisture
                                  ! variables are held
      Integer :: number_format    ! switch controlling number format diagnostics
                                  ! are written out in. See PP_WRITE for details.
      Integer :: model_domain     ! indicator as to model type, ie global, lam
      Integer :: n_channel        ! Number of satellite channels used
      Integer :: salt_dim1        !
      Integer :: salt_dim2        ! Dimensions for sea-salt aerosol diagnostics.
      Integer :: salt_dim3        !

      
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

      Real :: timestep

      Integer, Intent(IN) :: i_off
!       Offset to diagnostics in multiple calls to radiation
      Integer, Intent(IN) :: j_sw
!       Call to SW radiation

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

      Real :: itoasw (row_length, rows)
      Real :: surfsw (row_length, rows)
      Real :: surfsw_cor(row_length,rows)
      Real :: toasw_cor(row_length,rows)
      Real :: surfdir_cor(row_length,rows)
      Real :: surfdif_cor(row_length,rows)
      Real :: SWsea(row_length, rows)   ! Net short-wave absorbed by planet
      Real :: flux_below_690nm_surf(row_length, rows)
      Real :: sw_net_land(row_length, rows)
      Real :: sw_net_sice(row_length, rows)
      Real :: photosynth_act_rad(row_length, rows)
      Real :: flxdirparsurf(row_length, rows)
      Real :: T_incr_diagnostic(row_length,rows,model_levels)
      Real :: sea_salt_film(salt_dim1, salt_dim2, salt_dim3)
      Real :: sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)
!
! Orography variables
!
      Real :: f_orog(row_length,rows)      ! Extra SW surf flux
      Real :: slope_aspect(row_length,rows)! Gridbox mean slope aspect
      Real :: slope_angle(row_length,rows) ! Gridbox mean slope angle

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
! Diagnostic variables
!
      Real ::  STASHwork(*)    ! STASH workspace
!
! Local array & variables
!
      
      Real :: T_plus_T_inc(row_length, rows, model_levels)
      Real :: heating_rate(row_length,rows,model_levels)
      Real :: work_3d(row_length,rows,model_levels)

      Integer :: i, j, k
      Integer :: icode           ! Return code  =0 Normal exit  >1 Error
      Integer :: ptr_stash       ! Pointer to position in STASH

      CHARACTER(LEN=80)  cmessage
      CHARACTER(LEN=*) RoutineName
      Parameter ( RoutineName='diagnostics_sw')

      Integer :: im_index                ! internal model index
      Integer :: item                    ! STASH item
      Integer, Parameter :: sect=1       ! STASH SW radiation section

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('DIAGNOSTICS_SW',zhook_in,zhook_handle)
      icode    = 0                       ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1
! First treat diagnostics that are not contained in the SW_Diag 
! structure. These can only be output once and are valid after the 
! last call to radiation where j_sw==1
! ----------------------------------------------------------------------

  If (j_sw == 1) Then
!
! Temperature
!
     If (sf_calc(004,1)) Then

         Do k = 1, model_levels
            Do j = 1, rows
               Do i = 1, row_length
                  T_plus_T_inc(i,j,k) = T_n(i,j,k) + T_inc(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(004,1,im_index)),                &
     &       T_plus_T_inc,                                              &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,004,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,004,                                            &
     &       icode,cmessage)

         If (icode >  0) Then
            cmessage="Error in copydiag_3d( item 004)"
            goto 9999
         End if

      End if

! Horizon angles
    IF (l_skyview) THEN
      DO i = 1, n_horiz_ang
        item = 100+i
        IF (sf_calc(item,sect)) THEN
! DEPENDS ON: copydiag
          CALL copydiag (STASHwork(si(item,sect,im_index)),             &
            horiz_ang(:,:,i),                                           &
            row_length,rows,0,0,0,0, at_extremity,                      &
            atmos_im,sect,item,icode,cmessage)
          IF (icode > 0) THEN
            WRITE(cmessage,'(A,I3,A)') 'Error in copydiag(item ',item,')'
            GOTO 9999
          END IF
        END IF
      END DO
    END IF

!
! Temperature Increment: increment diagnostics= modified - previous
!
      item = 161
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
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            & 
     &       work_3d,                                                   &
     &       row_length,rows,wet_model_levels,0,0,0,0, at_extremity,    &
     &       stlist(1,stindex(1,item,sect,im_index)),len_stlist,        & 
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,sect,item,                                        &
     &       icode,cmessage)

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
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
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
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
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
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
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
          cmessage=": error in copydiag_3d(item 198)"//cmessage
        END IF

      END IF !  sf_calc(item,sect)

!
! Surfsw
!
      If (sf_calc(201,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(201,1,im_index)),surfsw,           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,201,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if

!
! Incoming SW at TOA
!
      If (sf_calc(207,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(207,1,im_index)),itoasw,           &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,207,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 207)"
            goto 9999
         End if

      End if

!
! surfsw_cor: surface SW corrected for solar zenith angle
!
      If (sf_calc(202,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(202,1,im_index)),surfsw_cor,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,202,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 202)"
            goto 9999
         End if

      End if

!
! toasw_cor: outgoing SW corrected for solar zenith angle
!
      If (sf_calc(205,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(205,1,im_index)),toasw_cor,        &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,205,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 205)"
            goto 9999
         End if

      End if

!
! surfdir_cor: direct surface SW corrected for solar zenith angle
!
      If (sf_calc(215,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(215,1,im_index)),surfdir_cor,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,215,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 215)"
            goto 9999
         End if

      End if

!
! surfdif_cor: diffuse surface SW corrected for solar zenith angle
!
      If (sf_calc(216,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(216,1,im_index)),surfdif_cor,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,216,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 216)"
            goto 9999
         End if

      End if

!
! SWSea
!
      If (sf_calc(203,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(203,1,im_index)),SWsea,            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,203,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 203)"
            goto 9999
         End if

      End if

!
! Flux Below 690nm at Surface
!
      If (sf_calc(204,1)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(204,1,im_index)),                   &
     &        flux_below_690nm_surf,                                    & 
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,204,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 204)"
            goto 9999
         End if

      End if

!
! Surface Net Land
!
      If (sf_calc(257,1)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(257,1,im_index)),                   &
     &        sw_net_land,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,257,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 257)"
            goto 9999
         End if

      End if

!
! SW Net Sice
!
      If (sf_calc(258,1)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(258,1,im_index)),                   &
     &        sw_net_sice,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,1,258,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 258)"
            goto 9999
         End if

      End if

!
! SW Heating:  SW heating =
! sw radiation temperature increment per timestep / timestep
!
      If (icode <= 0 .and. sf_calc(232,1)) Then

        Do k=1,model_levels
         Do j=1,rows
          Do i=1,row_length
            heating_rate(i,j,k) =  T_incr_diagnostic(i,j,k) /           &
     &                                       timestep
          Enddo ! i
         Enddo ! j
        Enddo ! k

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(232,1,im_index)),                &
     &       heating_rate,                                              &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,232,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,232,                                            &
     &       icode,cmessage)

        If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 232)"//cmessage
        End if

      Endif

!
! Photosynth_act_rad (total PAR flux at surface)
! 
      If (sf_calc(290,1)) Then

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(290,1,im_index)),                  &
     &       photosynth_act_rad,                                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,290,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
             cmessage="Error in copydiag( item 290)"
             goto 9999
          End if

      End if

!
! Flux_direct_par (direct component of PAR flux at surface)
!
      if (sf_calc(291,1)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(291,1,im_index)),                   &
     &       flxdirparsurf,                                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,291,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 291)"
            goto 9999
         End if

      End if

!
! Slope Aspect
!      
      If (sf_calc(293,1)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(293,1,im_index)),                   &
     &       slope_aspect,                                              &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,293,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 293)"
            goto 9999
         End if

      End if

!
! Slope Angle
!
      If (sf_calc(294,1)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(294,1,im_index)),                   &
     &       slope_angle,                                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,294,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 294)"
            goto 9999
         End if
      End if

!      
! F_orog
!
      If (sf_calc(296,1)) Then
      
! DEPENDS ON: copydiag      
         Call copydiag(STASHwork(si(296,1,im_index)),                   &
     &       f_orog,                                                    &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,296,                                            & 
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 296)"
            goto 9999
         End if

      End if

!
! Sea Salt film
!      
      If (sf_calc(247,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(247,1,im_index)),                &
     &       sea_salt_film,                                             &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,247,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,247,                                            & 
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 247)"
            goto 9999
         End if

      End if

!
! Sea Salt Jet
!
      If (sf_calc(248,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(248,1,im_index)),                &
     &       sea_salt_jet,                                              &
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,248,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,248,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 248)"
            goto 9999
         End if

      End if      

  End if ! j_sw == 1


! ----------------------------------------------------------------------
! Section 2
! Now treat diagnostics contained in structure SW_DIAG. These may be 
! available for each of the radiation calls dependent on the scheme.
! ----------------------------------------------------------------------

  If (l_rad_perturb.and.(j_sw == 1)) Then

! ----------------------------------------------------------------------
! Section 2.1
! For the incremental time-stepping scheme many of the diagnostics are
! not calculated on the "cloud only" radiation calls. For these the 
! diagnostics from the last full radiation call are used.
! ----------------------------------------------------------------------

!   Clear-sky upward flux on levels. Stash: sf_calc(219,1)
    IF (sf_calc(219,1)) THEN
      CALL copydiag_3d(STASHwork(si(219,1,im_index)),                   &
           SW_diag(1)%flux_up_clear,                                    &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,219,1,im_index)),len_stlist,              &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,219,                                              &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 219)"
        GOTO 9999
      END IF
    END IF

!   Clear-sky downward flux on levels. Stash: sf_calc(220,1)
    IF (sf_calc(220,1)) THEN
      CALL copydiag_3d(STASHwork(si(220,1,im_index)),                   &
           SW_diag(1)%flux_down_clear,                                  &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,220,1,im_index)),len_stlist,              &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,220,                                              &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 220)"
        GOTO 9999
      END IF
    END IF

!
! Outwards Solar Clear Flux at TOA. Stash: Sf_calc(209,1)
! 
        If (sf_calc(209,1)) Then

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(209,1,im_index)),                  &
     &       SW_diag(1)%solar_out_clear,                                &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,209,                                            &
     &       icode,cmessage)                     

           If (icode  >   0) then
             cmessage="Error in copydiag( item 209)"
             goto 9999
           End if

        End if

!
! Surface Down Clear: Stash: sf_calc(210,1)
!
        If (sf_calc(210,1)) Then

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(210,1,im_index)),                  &
     &       SW_diag(1)%surf_down_clr,                                  &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,210,                                            &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 210)"
             goto 9999
           End if

        End if

!
! Surface up Clear. Stash: sf_calc(211,1)
!
        If (sf_calc(211,1)) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(211,1,im_index)),                 &
     &       SW_diag(1)%surf_up_clr,                                    &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,211,                                            &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 211)"
             goto 9999
           End if

        End if

! Direct UV-Flux. Stash: sf_calc(212,1)
    IF (sf_calc(212,1)) THEN
      CALL copydiag_3d(STASHwork(si(212,1,im_index)),                   &
           SW_diag(1)%uvflux_direct,                                    &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        &
           stlist(1,stindex(1,212,1,im_index)),len_stlist,              &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,212,                                              &
           icode,cmessage)
      IF (icode > 0) THEN
        cmessage="Error in copydiag( item 212)"
        GOTO 9999
      END IF
    END IF

! UV Flux Up. Stash: sf_calc(213,1)
    IF (sf_calc(213,1)) THEN
      CALL copydiag_3d(STASHwork(si(213,1,im_index)),                   & 
           SW_diag(1)%uvflux_up,                                        &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        &
           stlist(1,stindex(1,213,1,im_index)),len_stlist,              & 
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,213,                                              &
           icode,cmessage)
      IF (icode > 0) THEN
        cmessage="Error in copydiag( item 213)"
        GOTO 9999
      END IF
    END IF

! Net UV Flux. Stash: sf_calc(214,1)
    IF (sf_calc(214,1)) THEN
      CALL copydiag_3d(STASHwork(si(214,1,im_index)),                   &
           SW_diag(1)%uvflux_net,                                       &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,214,1,im_index)),len_stlist,              &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,214,                                              &
           icode,cmessage)
      IF (icode > 0) THEN
        cmessage="Error in copydiag( item 214)"
        GOTO 9999
      END IF
    END IF

! Surface Down UV Flux. Stash: sf_calc(288,1)
    IF (sf_calc(288,1)) THEN
      CALL copydiag (STASHwork(si(288,1,im_index)),                     &
           SW_diag(1)%surf_uv,                                          &
           row_length,rows,0,0,0,0, at_extremity,                       &
           atmos_im,1,288,                                              &
           icode,cmessage)
      IF (icode > 0) THEN
        cmessage="Error in copydiag( item 288)"
        GOTO 9999
      END IF
    END IF

! Surface Down clear-sky UV Flux. Stash: sf_calc(289,1)
    IF (sf_calc(289,1)) THEN
      CALL copydiag (STASHwork(si(289,1,im_index)),                     &
           SW_diag(1)%surf_uv_clr,                                      &
           row_length,rows,0,0,0,0, at_extremity,                       &
           atmos_im,1,289,                                              &
           icode,cmessage)
      IF (icode > 0) THEN
        cmessage="Error in copydiag( item 289)"
        GOTO 9999
      END IF
    END IF

!
! Clear Sky Heating Rates. Stash: sf_calc(233,1)
!      
        If (sf_calc(233,1)) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(233,1,im_index)),               &
     &       SW_diag(1)%clear_hr,                                       & 
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,233,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,233,                                            &
             icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 233)"
             goto 9999
           End if

        End if

!
! Cloud Extinction. Stash: sf_calc(262,1)
!
        If (sf_calc(262,1)) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(262,1,im_index)),               &
     &       SW_diag(1)%cloud_extinction,                               &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,262,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,262,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 262)"
            goto 9999
          End if

        End if

!
! Cloud Weight Extinction. Stash: sf_calc(263,1)
!
        If (sf_calc(263,1)) Then
 
 ! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(263,1,im_index)),               & 
     &       SW_diag(1)%cloud_weight_extinction,                        &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,263,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,263,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
             cmessage="Error in copydiag( item 263)"
             goto 9999
          End if

        End if

!
! Large-Scale Cloud Extinction. Stash: sf_calc(264,1)
!
        If (sf_calc(264,1)) Then
 
! DEPENDS ON: copydiag_3d
           Call copydiag_3d(STASHwork(si(264,1,im_index)),              &
     &       SW_diag(1)%ls_cloud_extinction,                            &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,264,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,264,                                            &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 264)"
             goto 9999
           End if
      
        End if

!
! Large-Scale Cloud Weight Extinction. Stash: sf_calc(265,1)
!
        If (sf_calc(265,1)) Then

! DEPENDS ON: copydiag_3d
           Call copydiag_3d(STASHwork(si(265,1,im_index)),              &
     &       SW_diag(1)%ls_cloud_weight_extinction,                     &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,265,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,265,                                            &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 265)"
             goto 9999
           End if

        End if

!
! Convective Cloud Extinction. Stash: sf_calc(266,1)
!
        If (sf_calc(266,1)) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(266,1,im_index)),               &
     &       SW_diag(1)%cnv_cloud_extinction,                           &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,266,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,266,                                            &
     &       icode,cmessage)

          If (icode  >   0) then
             cmessage="Error in copydiag( item 266)"
             goto 9999
          End if

        End if

!
! Convective Cloud weight Extinction. Stash: sf_calc(267,1)
!
        If (sf_calc(267,1)) Then

! DEPENDS ON: copydiag_3d
           Call copydiag_3d(STASHwork(si(267,1,im_index)),              &
     &       SW_diag(1)%cnv_cloud_weight_extinction,                    &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,267,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,267,                                            &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 267)"
             goto 9999
           End if

        End if

!
! Re. Strat. Stash: sf_calc(221,1)
!
      If (sf_calc(221,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(221,1,im_index)),                &
     &       SW_diag(1)%re_strat,                                       &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,221,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,221,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 221)"
            goto 9999
         End if

      End if

!
! Wgt. Strat. Stash: sf_calc(223,1)
!
      If (sf_calc(223,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223,1,im_index)),                &
     &       SW_diag(1)%wgt_strat,                                      &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,223,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,223,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 223)"
            goto 9999
         End if

      End if

!
! LWP. Strat. Stash: sf_calc(224,i)
!
      If (sf_calc(224,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(224,1,im_index)),                &
     &       SW_diag(1)%lwp_strat,                                      &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,224,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,224,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 224)"
            goto 9999
         End if

      End if

!
! Re. Conv. Stash: sf_calc(225,i)
!
      If (sf_calc(225,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(225,1,im_index)),                &
     &       SW_diag(1)%re_conv,                                        &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,225,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           & 
     &       atmos_im,1,225,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 225)"
            goto 9999
         End if

      End if

!
! Wgt. Conv. Stash: sf_calc(226,1)
!
      If (sf_calc(226,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(226,1,im_index)),                & 
     &       SW_diag(1)%wgt_conv,                                       &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,226,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,226,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 226)"
            goto 9999
         End if

      End if

!
! Ntot. Diag. Stash: sf_calc(241,1)
!
      If (sf_calc(241,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(241,1,im_index)),                &
     &       SW_diag(1)%ntot_diag,                                      &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,241,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,241,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 241)"
            goto 9999
         End if

      End if

!
! Strat. LWC Diag. Stash: sf_calc(242,1)
!
      If (sf_calc(242,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(242,1,im_index)),                &
     &       SW_diag(1)%strat_lwc_diag,                                 &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,242,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,242,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 242)"
            goto 9999
         End if

      End if

!
! SO4 Cloud Condensation Nuclei. Stash: sf_calc(243,1)
!
      If (sf_calc(243,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(243,1,im_index)),                &
     &       SW_diag(1)%so4_ccn_diag,                                   & 
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,243,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,243,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 243)"
            goto 9999
         End if

      End if

!
! Cond. Samp. Wgt. Stash(244,1)
!
      If (sf_calc(244,1)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(244,1,im_index)),                &
     &       SW_diag(1)%cond_samp_wgt,                                  &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,244,1,im_index)),len_stlist,            &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,244,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 244)"
            goto 9999
         End if

      End if

!
!Weighted Re. Stash: sf_calc(245,1)
! 
      If (sf_calc(245,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(245,1,im_index)),                  & 
     &       SW_diag(1)%weighted_re,                                    &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,245,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 245)"
            goto 9999
         End if

      End if

!
! Sum Weighted Re. Stash: sf_calc(246,1)
!
      If (sf_calc(246,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(246,1,im_index)),                  &
     &       SW_diag(1)%sum_weight_re,                                  &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,246,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 246)"
            goto 9999
         End if

      End if

!
! Weighted Warm Re. Stash: sf_calc(254,1)
!
      If (sf_calc(254,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(254,1,im_index)),                  &
     &       SW_diag(1)%weighted_warm_re,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,1,254,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 254)"
            goto 9999
         End if

      End if

!
! Sum Weighted Warm Re. Stash: sf_calc(255,1)
!
      If (sf_calc(255,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(255,1,im_index)),                  &
     &       SW_diag(1)%sum_weight_warm_re,                             & 
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,1,255,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 255)"
            goto 9999
         End if

      End if

!
! Nc. Diag. Stash: sf_calc(280,1)
!
      If (sf_calc(280,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(280,1,im_index)),                  &
     &       SW_diag(1)%Nc_diag,                                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,280,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 280)"
            goto 9999
         End if

      End if

!
! Nc. Weight. Stash: sf_calc(281,1)
!
      If (sf_calc(281,1)) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(281,1,im_index)),                  &
     &       SW_diag(1)%Nc_weight,                                      &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,281,                                            &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 281)"
            goto 9999
         End if

      End if

!
! Solar Bearing. Stash: sf_calc(292,1)
!
      If (sf_calc(292,1)) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(292,1,im_index)),                   &
     &       SW_diag(1)%sol_bearing,                                    &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,292,                                            &
     &       icode,cmessage)              

         If (icode  >   0) then
            cmessage="Error in copydiag( item 292)"
            goto 9999
         End if

      End if

!
! Direct Surface Albedo on SW bands. Stash 268,1
!
     IF (sf_calc(268,1)) THEN
       DO k = 1, sw_spectrum(1)%n_band
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si(268,1,im_index)                    &
             +(row_length*rows*(k-1))),                                 & 
             SW_diag(1)%direct_albedo(1,1,k),                           &
             row_length,rows,0,0,0,0, at_extremity,                     &
             atmos_im,1,268,                                            &
             icode,cmessage)

          IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 268)"
            goto 9999
          END IF
        END DO  
      END IF

!
! Diffuse Surface Albedo on SW bands. Stash 269,1
!
     IF (sf_calc(269,1)) THEN
       DO k = 1, sw_spectrum(1)%n_band
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si(269,1,im_index)                    &
             +(row_length*rows*(k-1))),                                 & 
             SW_diag(1)%diffuse_albedo(1,1,k),                          &
             row_length,rows,0,0,0,0, at_extremity,                     &
             atmos_im,1,269,                                            &
             icode,cmessage)

          IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 269)"
            goto 9999
          END IF
        END DO  
      END IF

  Else ! l_rad_perturb.and.(j_sw==1)

! ----------------------------------------------------------------------
! Section 2.2
! The following diagnostics are available on "cloud only" radiation
! calls for the incremental time-stepping scheme, and all calls for
! other schemes.
! ----------------------------------------------------------------------


!   Upward flux on levels. Stash: sf_calc(217,1)
    IF (SW_diag(j_sw)%l_flux_up) THEN
      CALL copydiag_3d(STASHwork(si(217+i_off,1,im_index)),             &
           SW_diag(j_sw)%flux_up,                                       &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,217+i_off,1,im_index)),len_stlist,        &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,217+i_off,                                        &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 217)"
        GOTO 9999
      END IF
    END IF

!   Downward flux on levels. Stash: sf_calc(218,1)
    IF (SW_diag(j_sw)%l_flux_down) THEN
      CALL copydiag_3d(STASHwork(si(218+i_off,1,im_index)),             &
           SW_diag(j_sw)%flux_down,                                     &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,218+i_off,1,im_index)),len_stlist,        &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,218+i_off,                                        &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 218)"
        GOTO 9999
      END IF
    END IF

!
! Outwards Solar Flux at TOA. Stash: Sf_calc(208,1)
!           
        If (SW_diag(j_sw)%L_solar_out_toa) Then

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(208+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%solar_out_toa,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,208+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
          End if

        End if

!
! Surface Down Flux. Stash: sf_calc(235,1)
!
        If (SW_diag(j_sw)%l_surface_down_flux) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(235+i_off,1,im_index)),           &
     &       SW_diag(j_sw)%surface_down_flux,                           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,235+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 235)"
             goto 9999
           End if

        End if

!
! Direct Flux. Stash: sf_calc(230,1)
!
      If (SW_diag(j_sw)%L_flux_direct) Then
      
         Call copydiag_3d(STASHwork(si(230+i_off,1,im_index)),          &
     &        SW_diag(j_sw)%flux_direct,                                &
     &        row_length,rows,model_levels+1,0,0,0,0, at_extremity,     & 
     &        stlist(1,stindex(1,230+i_off,1,im_index)),len_stlist,     &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,230+i_off,                                     &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 230)"
            goto 9999
         End if

      End if
!
! Diffuse Flux. Stash: sf_calc(231,1)
!
      If (SW_diag(j_sw)%L_flux_diffuse) Then

         Call copydiag_3d(STASHwork(si(231+i_off,1,im_index)),          &
     &        SW_diag(j_sw)%flux_diffuse,                               &
     &        row_length,rows,model_levels+1,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,231+i_off,1,im_index)),len_stlist,     &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,1,231+i_off,                                     &
     &        icode,cmessage)

         If (icode  >  0) then
            cmessage="Error in copydiag( item 231)"
            goto 9999
         End if

      End if    
        
!
! Net SW Flux at Tropopause: Stash: sf_calc(237,1)
!
        If (SW_diag(j_sw)%L_net_flux_trop) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(237+i_off,1,im_index)),           &
     &       SW_diag(j_sw)%net_flux_trop,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,237+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
              cmessage="Error in copydiag( item 237)"
              goto 9999
           End if

        End if

!
! SW Up Flux at Tropopause: Stash: sf_calc(238,i)
!
        If (SW_diag(j_sw)%l_up_flux_trop) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(238+i_off,1,im_index)),           &
     &       SW_diag(j_sw)%up_flux_trop,                                &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,238+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 238)"
            goto 9999
          End if

        End if

!
! Radiances at TOA: Stash: sf_calc(297,1)
!
        If (SW_diag(j_sw)%L_toa_radiance) Then

           Do i=1, n_channel

             ptr_stash = si(297+i_off,1,im_index)                       &
     &                 + (i-1) * row_length * rows
     
! DEPENDS ON: copydiag     
             Call copydiag (STASHwork(ptr_stash),                       &
     &           SW_diag(j_sw)%toa_radiance(1, 1, i),                   &
     &           row_length,rows,0,0,0,0, at_extremity,                 &
     &           atmos_im,1,297+i_off,                                  &
     &           icode,cmessage)

             If (icode  >   0) then
                cmessage="Error in copydiag( item 297)"
                goto 9999
             End if
           End do

        End if

!
! Fl_solid_below_690nm_surf. Stash: sf_calc(259,1)
!
      If (SW_diag(j_sw)%L_FlxSolBelow690nmSurf) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(259+i_off,1,im_index)),             &
     &       SW_diag(j_sw)%FlxSolBelow690nmSurf,                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,259+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 259)"
            goto 9999
         End if

      End if

!
! Fl_sea_below_690nm_surf. Stash: sf_calc(260,i)
! 
      If (SW_diag(j_sw)%L_FlxSeaBelow690nmSurf) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(260+i_off,1,im_index)),             &
     &       SW_diag(j_sw)%FlxSeaBelow690nmSurf,                        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,260+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 260)"
            goto 9999
         End if

      End if

!
! Orographic Correction. Stash: sf_calc(295,1)
!
      If (SW_diag(j_sw)%L_orog_corr) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(295+i_off,1,im_index)),             &
     &       SW_diag(j_sw)%orog_corr,                                   & 
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,295+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 295)"
            goto 9999
         End if

      End if

!
! VIS Albedo scaling to obs, on land surface tiles. Stash 270,1
!
     IF (SW_diag(j_sw)%l_vis_albedo_sc) THEN
       DO k = 1, ntiles
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si(270+i_off,1,im_index)              &
             +(row_length*rows*(k-1))),                                 & 
             SW_diag(j_sw)%vis_albedo_sc(1,1,k),                        &
             row_length,rows,0,0,0,0, at_extremity,                     &
             atmos_im,1,270+i_off,                                      &
             icode,cmessage)

          IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 270)"
            goto 9999
          END IF
        END DO  
      END IF

!
! NIR Albedo scaling to obs, on land surface tiles. Stash 271,1
!
     IF (SW_diag(j_sw)%l_nir_albedo_sc) THEN
       DO k = 1, ntiles
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si(271+i_off,1,im_index)              &
             +(row_length*rows*(k-1))),                                 & 
             SW_diag(j_sw)%nir_albedo_sc(1,1,k),                        &
             row_length,rows,0,0,0,0, at_extremity,                     &
             atmos_im,1,271+i_off,                                      &
             icode,cmessage)

          IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 271)"
            goto 9999
          END IF
        END DO  
      END IF

  End if  ! l_rad_perturb.and.(j_sw==1) - end of Else block



  If (.not.l_rad_perturb) Then

! ----------------------------------------------------------------------
! Section 2.3
! The following diagnostics have been output in section 2.1 for the
! incremental time-stepping scheme. They are now output for other
! schemes on any radiation call.
! ----------------------------------------------------------------------


!   Clear-sky upward flux on levels. Stash: sf_calc(219,1)
    IF (SW_diag(j_sw)%l_flux_up_clear) THEN
      CALL copydiag_3d(STASHwork(si(219+i_off,1,im_index)),             &
           SW_diag(j_sw)%flux_up_clear,                                 &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,219+i_off,1,im_index)),len_stlist,        &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,219+i_off,                                        &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 219)"
        GOTO 9999
      END IF
    END IF

!   Clear-sky downward flux on levels. Stash: sf_calc(220,1)
    IF (SW_diag(j_sw)%l_flux_down_clear) THEN
      CALL copydiag_3d(STASHwork(si(220+i_off,1,im_index)),             &
           SW_diag(j_sw)%flux_down_clear,                               &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,220+i_off,1,im_index)),len_stlist,        &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,220+i_off,                                        &
           icode,cmessage)
      IF (icode  >  0) THEN
        cmessage="Error in copydiag( item 220)"
        GOTO 9999
      END IF
    END IF

!
! Outwards Solar Clear Flux at TOA. Stash: Sf_calc(209,1)
! 
        If (SW_diag(j_sw)%L_solar_out_clear) Then

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(209+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%solar_out_clear,                             &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,209+i_off,                                      &
     &       icode,cmessage)                     

           If (icode  >   0) then
             cmessage="Error in copydiag( item 209)"
             goto 9999
           End if

        End if

!
! Surface Down Clear: Stash: sf_calc(210,1)
!
        If (SW_diag(j_sw)%L_surf_down_clr) Then

! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(210+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%surf_down_clr,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,210+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 210)"
             goto 9999
           End if

        End if

!
! Surface up Clear. Stash: sf_calc(211,1)
!
        If (SW_diag(j_sw)%L_surf_up_clr) Then

! DEPENDS ON: copydiag
          Call copydiag (STASHwork(si(211+i_off,1,im_index)),           &
     &       SW_diag(j_sw)%surf_up_clr,                                 &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,211+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 211)"
             goto 9999
           End if

        End if

! Direct UV-Flux. Stash: sf_calc(212,1)
    IF (SW_diag(j_sw)%l_uvflux_direct) THEN
      CALL copydiag_3d(STASHwork(si(212+i_off,1,im_index)),             &
           SW_diag(j_sw)%uvflux_direct,                                 &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        &
           stlist(1,stindex(1,212+i_off,1,im_index)),len_stlist,        &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,212+i_off,                                        &
           icode,cmessage)
      IF (icode > 0) THEN
        cmessage="Error in copydiag( item 212)"
        GOTO 9999
      END IF
    END IF

! UV Flux Up. Stash: sf_calc(213,1)
    IF (SW_diag(j_sw)%l_uvflux_up) THEN
      CALL copydiag_3d(STASHwork(si(213+i_off,1,im_index)),             & 
           SW_diag(j_sw)%uvflux_up,                                     &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        &
           stlist(1,stindex(1,213+i_off,1,im_index)),len_stlist,        & 
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,213+i_off,                                        &
           icode,cmessage)
      IF (icode > 0) THEN
        cmessage="Error in copydiag( item 213)"
        GOTO 9999
      END IF
    END IF

! Net UV Flux. Stash: sf_calc(214,1)
    IF (SW_diag(j_sw)%l_uvflux_net) THEN
      CALL copydiag_3d(STASHwork(si(214+i_off,1,im_index)),             &
           SW_diag(j_sw)%uvflux_net,                                    &
           row_length,rows,model_levels+1,0,0,0,0, at_extremity,        & 
           stlist(1,stindex(1,214+i_off,1,im_index)),len_stlist,        &
           stash_levels,num_stash_levels+1,                             &
           atmos_im,1,214+i_off,                                        &
           icode,cmessage)
      IF (icode > 0) THEN
        cmessage="Error in copydiag( item 214)"
        GOTO 9999
      END IF
    END IF

! Surface Down UV Flux. Stash: sf_calc(288,1)
    IF (SW_diag(j_sw)%l_surf_uv) THEN
      CALL copydiag (STASHwork(si(288+i_off,1,im_index)),               &
           SW_diag(j_sw)%surf_uv,                                       &
           row_length,rows,0,0,0,0, at_extremity,                       &
           atmos_im,1,288+i_off,                                        &
           icode,cmessage)
      IF (icode > 0) THEN
        cmessage="Error in copydiag( item 288)"
        GOTO 9999
      END IF
    END IF

! Surface Down clear-sky UV Flux. Stash: sf_calc(289,1)
    IF (SW_diag(j_sw)%l_surf_uv_clr) THEN
      CALL copydiag (STASHwork(si(289+i_off,1,im_index)),               &
           SW_diag(j_sw)%surf_uv_clr,                                   &
           row_length,rows,0,0,0,0, at_extremity,                       &
           atmos_im,1,289+i_off,                                        &
           icode,cmessage)
      IF (icode > 0) THEN
        cmessage="Error in copydiag( item 289)"
        GOTO 9999
      END IF
    END IF

!
! Clear Sky Heating Rates. Stash: sf_calc(233,1)
!      
        If (SW_diag(j_sw)%L_clear_hr) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(233+i_off,1,im_index)),         &
     &       SW_diag(j_sw)%clear_hr,                                    & 
     &       row_length,rows,model_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,233+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,233+i_off,                                      &
             icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 233)"
             goto 9999
           End if

        End if

!
! Cloud Extinction. Stash: sf_calc(262,1)
!
        If (SW_diag(j_sw)%L_cloud_extinction) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(262+i_off,1,im_index)),         &
     &       SW_diag(j_sw)%cloud_extinction,                            &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,262+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,262+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
            cmessage="Error in copydiag( item 262)"
            goto 9999
          End if

        End if

!
! Cloud Weight Extinction. Stash: sf_calc(263,1)
!
        If (SW_diag(j_sw)%L_cloud_weight_extinction) Then
 
 ! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(263+i_off,1,im_index)),         & 
     &       SW_diag(j_sw)%cloud_weight_extinction,                     &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,263+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,263+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
             cmessage="Error in copydiag( item 263)"
             goto 9999
          End if

        End if

!
! Large-Scale Cloud Extinction. Stash: sf_calc(264,1)
!
        If (SW_diag(j_sw)%L_ls_cloud_extinction) Then
 
! DEPENDS ON: copydiag_3d
           Call copydiag_3d(STASHwork(si(264+i_off,1,im_index)),        &
     &       SW_diag(j_sw)%ls_cloud_extinction,                         &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,264+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,264+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 264)"
             goto 9999
           End if
      
        End if

!
! Large-Scale Cloud Weight Extinction. Stash: sf_calc(265,1)
!
        If (SW_diag(j_sw)%L_ls_cloud_weight_extinction) Then

! DEPENDS ON: copydiag_3d
           Call copydiag_3d(STASHwork(si(265+i_off,1,im_index)),        &
     &       SW_diag(j_sw)%ls_cloud_weight_extinction,                  &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,265+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,265+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 265)"
             goto 9999
           End if

        End if

!
! Convective Cloud Extinction. Stash: sf_calc(266,1)
!
        If (SW_diag(j_sw)%L_cnv_cloud_extinction) Then

! DEPENDS ON: copydiag_3d
          Call copydiag_3d(STASHwork(si(266+i_off,1,im_index)),         &
     &       SW_diag(j_sw)%cnv_cloud_extinction,                        &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,266+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,266+i_off,                                      &
     &       icode,cmessage)

          If (icode  >   0) then
             cmessage="Error in copydiag( item 266)"
             goto 9999
          End if

        End if

!
! Convective Cloud weight Extinction. Stash: sf_calc(267,1)
!
        If (SW_diag(j_sw)%L_cnv_cloud_weight_extinction) Then

! DEPENDS ON: copydiag_3d
           Call copydiag_3d(STASHwork(si(267+i_off,1,im_index)),        &
     &       SW_diag(j_sw)%cnv_cloud_weight_extinction,                 &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,267+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,267+i_off,                                      &
     &       icode,cmessage)

           If (icode  >   0) then
             cmessage="Error in copydiag( item 267)"
             goto 9999
           End if

        End if

!
! Re. Strat. Stash: sf_calc(221,1)
!
      If (SW_diag(j_sw)%re_strat_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(221+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%re_strat,                                    &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,221+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,221+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 221)"
            goto 9999
         End if

      End if

!
! Wgt. Strat. Stash: sf_calc(223,1)
!
      If (SW_diag(j_sw)%wgt_strat_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%wgt_strat,                                   &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,223+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,223+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 223)"
            goto 9999
         End if

      End if

!
! LWP. Strat. Stash: sf_calc(224,i)
!
      If (SW_diag(j_sw)%lwp_strat_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(224+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%lwp_strat,                                   &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,224+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,224+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 224)"
            goto 9999
         End if

      End if

!
! Re. Conv. Stash: sf_calc(225,i)
!
      If (SW_diag(j_sw)%re_conv_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(225+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%re_conv,                                     &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,225+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           & 
     &       atmos_im,1,225+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 225)"
            goto 9999
         End if

      End if

!
! Wgt. Conv. Stash: sf_calc(226,1)
!
      If (SW_diag(j_sw)%wgt_conv_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(226+i_off,1,im_index)),          & 
     &       SW_diag(j_sw)%wgt_conv,                                    &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,226+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,226+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 226)"
            goto 9999
         End if

      End if

!
! Ntot. Diag. Stash: sf_calc(241,1)
!
      If (SW_diag(j_sw)%ntot_diag_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(241+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%ntot_diag,                                   &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,241+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,241+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 241)"
            goto 9999
         End if

      End if

!
! Strat. LWC Diag. Stash: sf_calc(242,1)
!
      If (SW_diag(j_sw)%strat_lwc_diag_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(242+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%strat_lwc_diag,                              &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,242+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,242+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 242)"
            goto 9999
         End if

      End if

!
! SO4 Cloud Condensation Nuclei. Stash: sf_calc(243,1)
!
      If (SW_diag(j_sw)%so4_ccn_diag_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(243+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%so4_ccn_diag,                                & 
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,243+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,243+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 243)"
            goto 9999
         End if

      End if

!
! Cond. Samp. Wgt. Stash(244,1)
!
      If (SW_diag(j_sw)%cond_samp_wgt_flag) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(244+i_off,1,im_index)),          &
     &       SW_diag(j_sw)%cond_samp_wgt,                               &
     &       row_length,rows,cloud_levels,0,0,0,0, at_extremity,        &
     &       stlist(1,stindex(1,244+i_off,1,im_index)),len_stlist,      &
     &       stash_levels,num_stash_levels+1,                           &
     &       atmos_im,1,244+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 244)"
            goto 9999
         End if

      End if

!
!Weighted Re. Stash: sf_calc(245,1)
! 
      If (SW_diag(j_sw)%weighted_re_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(245+i_off,1,im_index)),            & 
     &       SW_diag(j_sw)%weighted_re,                                 &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,245+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 245)"
            goto 9999
         End if

      End if

!
! Sum Weighted Re. Stash: sf_calc(246,1)
!
      If (SW_diag(j_sw)%sum_weight_re_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(246+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%sum_weight_re,                               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,246+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 246)"
            goto 9999
         End if

      End if

!
! Weighted Warm Re. Stash: sf_calc(254,1)
!
      If (SW_diag(j_sw)%wgtd_warm_re_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(254+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%weighted_warm_re,                            &
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,1,254+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 254)"
            goto 9999
         End if

      End if

!
! Sum Weighted Warm Re. Stash: sf_calc(255,1)
!
      If (SW_diag(j_sw)%sum_wgt_warm_re_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(255+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%sum_weight_warm_re,                          & 
     &       row_length,rows,0,0,0,0, at_extremity,                     & 
     &       atmos_im,1,255+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 255)"
            goto 9999
         End if

      End if

!
! Nc. Diag. Stash: sf_calc(280,1)
!
      If (SW_diag(j_sw)%Nc_diag_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(280+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%Nc_diag,                                     &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,280+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 280)"
            goto 9999
         End if

      End if

!
! Nc. Weight. Stash: sf_calc(281,1)
!
      If (SW_diag(j_sw)%Nc_weight_flag) Then

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(281+i_off,1,im_index)),            &
     &       SW_diag(j_sw)%Nc_weight,                                   &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,281+i_off,                                      &
     &       icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 281)"
            goto 9999
         End if

      End if

!
! Solar Bearing. Stash: sf_calc(292,1)
!
      If (SW_diag(j_sw)%L_sol_bearing) Then

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(292+i_off,1,im_index)),             &
     &       SW_diag(j_sw)%sol_bearing,                                 &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,1,292+i_off,                                      &
     &       icode,cmessage)              

         If (icode  >   0) then
            cmessage="Error in copydiag( item 292)"
            goto 9999
         End if

      End if

!
! Direct Surface Albedo on SW bands. Stash 268,1
!
     IF (SW_diag(j_sw)%l_direct_albedo) THEN
       DO k = 1, sw_spectrum(j_sw)%n_band
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si(268+i_off,1,im_index)              &
             +(row_length*rows*(k-1))),                                 & 
             SW_diag(j_sw)%direct_albedo(1,1,k),                        &
             row_length,rows,0,0,0,0, at_extremity,                     &
             atmos_im,1,268+i_off,                                      &
             icode,cmessage)

          IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 268)"
            goto 9999
          END IF
        END DO  
      END IF

!
! Diffuse Surface Albedo on SW bands. Stash 269,1
!
     IF (SW_diag(j_sw)%l_diffuse_albedo) THEN
       DO k = 1, sw_spectrum(j_sw)%n_band
! DEPENDS ON: copydiag
          CALL copydiag(STASHwork(si(269+i_off,1,im_index)              &
             +(row_length*rows*(k-1))),                                 & 
             SW_diag(j_sw)%diffuse_albedo(1,1,k),                       &
             row_length,rows,0,0,0,0, at_extremity,                     &
             atmos_im,1,269+i_off,                                      &
             icode,cmessage)

          IF (icode  >   0) THEN
            cmessage="Error in copydiag( item 269)"
            goto 9999
          END IF
        END DO  
      END IF

  End if ! .not.l_rad_perturb

 9999 continue  ! exit point on error
      If(icode /= 0) Then

        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      IF (lhook) CALL dr_hook('DIAGNOSTICS_SW',zhook_out,zhook_handle)
      RETURN
    End Subroutine diagnostics_sw
