! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_gwd
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: gravity wave drag

      Subroutine diagnostics_gwd(                                       &
     &                       row_length, rows, model_levels             &
     &,                      n_rows, u_rows, v_rows                     &
     &,                      off_x, off_y                               &
     &,                      at_extremity                               &
     &,                      u, v, R_u, R_v, T_inc                      &
     &,                      u_incr_diagnostic,v_incr_diagnostic        &
     &,                      t_incr_diagnostic                          &
     &,                      exner_theta_levels                         &
     &,                      stress_ud     ,  stress_vd                 &
     &,                      stress_ud_satn,  stress_vd_satn            &
     &,                      stress_ud_wake,  stress_vd_wake            &
     &,                      du_dt_satn    ,  dv_dt_satn                &
     &,                      du_dt_wake   ,  dv_dt_wake                 &
     &,                      u_s_d, v_s_d, nsq_s_d                      &
     &,                      num_lim_d, num_fac_d                       &
     &,                      fr_d, bld_d ,bldt_d                        &
     &, tausx_d, tausy_d, taus_scale_d                                  &
     &, sd_orog_land, land_sea_mask, land_points                        &
     &, orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land         &
     &, orog_slope_d, orog_anis_d, orog_dir_d                           &
     &, GWSPEC_EFLUX, GWSPEC_SFLUX, GWSPEC_WFLUX                        &
     &, GWSPEC_NFLUX, GWSPEC_EWACC, GWSPEC_NSACC,                       &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &     STASHwork                                                    &
     &     )


! Purpose:
!          Calculates diagnostics and outputs them.
!
! Method:
!   Item by item copying of diagnostic arrays made available from gwd
! scheme into the corresponding STASHwork (master array for output)
! diagnostic space. All sizes and address pointers have been
! calculated by the STASH initialisation routines previously. Each
! diagnostic is protected by a STASHflag logical switch which is
! re-set every timestep.
!   Note that diagnostics for du_dt,dv_dt are available on all
! model rho levels. Stress diagnostics are now available on all model
! theta_levels - running from level 0 at the ground to
! level=model_levels at the model lid.
! New single level diagnostics (6214-6222) are output at theta points
! in the horizontal, i.e. they are output at the points on which the
! GWD scheme operates. 6233 and 6234 are the same.
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
! Definitions of prognostic variable array sizes
      USE atm_fields_bounds_mod, ONLY:                                   &
          udims, vdims, wdims, tdims, pdims,                             &
          udims_s, vdims_s, wdims_s, tdims_s, pdims_s                     
      USE mask_compression, ONLY : expand_from_mask
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE Submodel_Mod

      IMPLICIT NONE
!
!     Fixed starting theta-level (e.g. for STASH diagnostic arrays)
      Integer, Parameter :: tkfix1start    = 1

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, number_format                                                   &
                         ! switch controlling number format diagnostics
                         ! are written out in. See PP_WRITE for details.
     &, model_domain                                                    &
                         ! indicator as to model type, ie global, lam
     &, land_points               ! Number of land points

      Integer                                                           &
     &  off_x                                                           &
                            !IN. size of small halo in x direction
     &, off_y               !IN. size of small halo in y direction

      INTEGER ::             &
     &   u_rows              & ! rows on u grid
     &,  v_rows                ! rows on v grid

      Integer                                                           &
        index

      Real                                                              &
        exner  ! Exner pressure surface value


! Primary Arrays used in all models
      REAL       ::                           &!intent(inout):
     & r_u(udims_s%i_start:udims_s%i_end,    & !u wind increment
     &     udims_s%j_start:udims_s%j_end,    &
     &     udims_s%k_start:udims_s%k_end)    &
     &,r_v(vdims_s%i_start:vdims_s%i_end,    & !v wind increment
     &     vdims_s%j_start:vdims_s%j_end,    &
     &     vdims_s%k_start:vdims_s%k_end)    &
     &,T_inc(tdims%i_start:tdims%i_end,      & !temperature increment
     &       tdims%j_start:tdims%j_end,      &
     &       1:tdims%k_end) 
   
       REAL ::                              &!block for intent(inout)   
     & u(udims_s%i_start:udims_s%i_end,    & ! primary u field (ms**-1)
     &   udims_s%j_start:udims_s%j_end,    &
     &   udims_s%k_start:udims_s%k_end)    &
     &,v(vdims_s%i_start:vdims_s%i_end,    & ! primary v field (ms**-1)
     &   vdims_s%j_start:vdims_s%j_end,    &
     &   vdims_s%k_start:vdims_s%k_end)    

   REAL ::                                            &
     & exner_theta_levels                            &
     &             (tdims_s%i_start:tdims_s%i_end,   & ! Exner on theta level
     &               tdims_s%j_start:tdims_s%j_end,  & !
     &               tdims_s%k_start:tdims_s%k_end)   

      Real                                                              &
     &  stress_ud      (row_length,  rows,0:model_levels)               &
     &, stress_vd      (row_length,n_rows,0:model_levels)               &
     &, stress_ud_satn (row_length,  rows,0:model_levels)               &
     &, stress_vd_satn (row_length,n_rows,0:model_levels)               &
     &, stress_ud_wake (row_length,  rows,0:model_levels)               &
     &, stress_vd_wake (row_length,n_rows,0:model_levels)               &
     &, du_dt_satn     (row_length,  rows,  model_levels)               &
     &, dv_dt_satn     (row_length,n_rows,  model_levels)               &
     &, du_dt_wake     (row_length, rows,   model_levels)               &
     &, dv_dt_wake     (row_length,n_rows,  model_levels)               &
     &, GWSPEC_EFLUX(udims%i_start:udims%i_end,                         &
     &               udims%j_start:udims%j_end,                         &
     &                 tkfix1start:tdims%k_end)  & ! IN Fp in E azimuth
     &, GWSPEC_SFLUX(vdims%i_start:vdims%i_end,                         &
     &               vdims%j_start:vdims%j_end,                         &
     &                 tkfix1start:tdims%k_end)  & ! IN Fp in S azimuth
     &, GWSPEC_WFLUX(udims%i_start:udims%i_end,                         &
     &               udims%j_start:udims%j_end,                         &
     &                 tkfix1start:tdims%k_end)  & ! IN Fp in W azimuth
     &, GWSPEC_NFLUX(vdims%i_start:vdims%i_end,                         &
     &               vdims%j_start:vdims%j_end,                         &
     &                 tkfix1start:tdims%k_end)  & ! IN Fp in N azimuth
     &, GWSPEC_EWACC(udims%i_start:udims%i_end,                         &
     &               udims%j_start:udims%j_end,                         &
     &               udims%k_start:udims%k_end)  & ! IN Accel of U wind
     &, GWSPEC_NSACC(vdims%i_start:vdims%i_end,                         &
     &               vdims%j_start:vdims%j_end,                         &
     &               vdims%k_start:vdims%k_end)  & ! IN Accel of V wind
! Note: space is only allocated for _incr_diagnostic arrays in calling
!       routine if coresponding STASHflags are set .true.
     &, u_incr_diagnostic(udims%i_start:udims%i_end,                    & 
     &                    udims%j_start:udims%j_end,                    &
     &                    udims%k_start:udims%k_end)                    &
     &, v_incr_diagnostic(vdims%i_start:vdims%i_end,                    & 
     &                    vdims%j_start:vdims%j_end,                    &
     &                    vdims%k_start:vdims%k_end)                    &
     &, t_incr_diagnostic(tdims%i_start:tdims%i_end,                    & 
     &                    tdims%j_start:tdims%j_end,                    &
     &                    tkfix1start:tdims%k_end)                      &
     &, num_lim_d(row_length, rows)                                     &
     &, num_fac_d(row_length, rows)                                     &
     &, u_s_d  (row_length, rows)                                       &
     &, v_s_d  (row_length, rows)                                       &
     &, nsq_s_d  (row_length, rows)                                     &
     &, fr_d   (row_length, rows)                                       &
     &, bld_d  (row_length, rows)                                       &
     &, bldt_d (row_length, rows)                                       &
     &, tausx_d (row_length, u_rows)                                    &
     &, tausy_d (row_length, v_rows)                                    &
     &, taus_scale_d (row_length, rows)                                 &
     &, orog_slope_d  (row_length,rows)                                 &
     &, orog_anis_d (row_length,rows)                                   &
     &, orog_dir_d  (row_length,rows)                                   &
     &, work_array(row_length, rows)                                    &
     &, sd_orog_land(land_points)                                       &
                                  ! orographic std dev on land points
     &, orog_grad_xx_land(land_points)                                  &
                                       ! sigma_xx on land points
     &, orog_grad_xy_land(land_points)                                  &
                                       ! sigma_xy on land points
     &, orog_grad_yy_land(land_points) ! sigma_yy on land points

      LOGICAL, INTENT(IN) ::                                            &
     &  land_sea_mask (land_points)  ! land_sea_mask

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
      Real                                                              &
     & STASHwork(*)     ! STASH workspace

! Local variables
      Integer                                                           &
     &  i, j, k, l                                                      &
                      ! loop indices
     &,    icode                                                        &
                                ! Return code  =0 Normal exit  >1 Error
     & ,item                                                            &
                        ! STASH item
     & ,sect            ! STASH section
      Parameter( sect = 6 ) ! for gravity wave drag

      CHARACTER(LEN=80)  cmessage
      CHARACTER(LEN=*) RoutineName
      Parameter ( RoutineName='diagnostics_gwd')

      Integer                                                           &
     &  im_index        ! internal model index

      Integer                                                           &
        gwspec_eflux_p_levels,                                          &
        stress_ud_p_levels,                                             &
        du_dt_satn_p_levels,                                            &
        gwspec_wflux_p_levels,                                          &
        gwspec_ewacc_p_levels

      Real                                                              &
     &  interp_data_3(row_length*rows*model_levels)                     &
     &, interp_data_3_u(row_length*u_rows*model_levels)                 &
     &, interp_data_3_v(row_length*v_rows*model_levels)                 &
     &, sd_orog     (row_length,rows)                                   &
                                       ! orog std dev full field
     &, orog_grad_xx(row_length,rows)                                   &
                                       !     sigma_xx full field
     &, orog_grad_xy(row_length,rows)                                   &
                                       !     sigma_xy full field
     &, orog_grad_yy(row_length,rows)  !     sigma_yy full field


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ND pp codes temporarily kept around :
!
!      PP_code_r                =   1
!      PP_code_p                =   8
!      PP_code_eta              =  10   !Note: currently no PP code
!!                             for eta, so set to 10 which is sigma.
!      PP_code_T                =  16
!      PP_code_msl              = 128
!      PP_code_surf             = 129
!      PP_code_du_dt_satn       =  68
!      PP_code_dv_dt_satn       =  69
!      PP_code_gwd_stress_u     = 161
!      PP_code_gwd_stress_v     = 162
!      PP_code_u_tend           = 975
!      PP_code_v_tend           = 976
!
      IF (lhook) CALL dr_hook('DIAGNOSTICS_GWD',zhook_in,zhook_handle)
      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 1.1 Diagnostics for GWD
! ----------------------------------------------------------------------

!-----------------------------------------------------------------------
! u,v increment diagnostics= modified - previous
!-----------------------------------------------------------------------

      item = 185  ! u increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k=1,model_levels
          Do j=udims%j_start,udims%j_end
           Do i=udims%i_start,udims%i_end
              u_incr_diagnostic(i,j,k) = R_u(i,j,k) -                   &
     &                                      u_incr_diagnostic(i,j,k)
           End Do !i
          End Do  !j
         End Do   !k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        u_incr_diagnostic,                                        &
     &        row_length,u_rows,model_levels,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 185)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 186  ! v increment
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k=1,model_levels
          Do j=vdims%j_start,vdims%j_end
           Do i=vdims%i_start,vdims%i_end
            v_incr_diagnostic(i,j,k) = R_v(i,j,k) -                     &
     &                                      v_incr_diagnostic(i,j,k)
           End Do !i
          End Do  !j
         End Do   !k

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        v_incr_diagnostic,                                        &
     &        row_length,v_rows,model_levels,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 186)"//cmessage
         End if

      Endif  !  sf(item,sect)

!-----------------------------------------------------------------------
! temperature increment diagnostic= modified - previous
!-----------------------------------------------------------------------

      item = 181  ! T increment
      IF (icode <= 0 .AND. sf(item,sect)) THEN

         DO k=tkfix1start,model_levels
          DO j=tdims%j_start,tdims%j_end
           DO i=tdims%i_start,tdims%i_end
              t_incr_diagnostic(i,j,k) = T_inc(i,j,k) -                 &
                                            t_incr_diagnostic(i,j,k)
           END DO !i
          END DO  !j
         END DO   !k

! DEPENDS ON: copydiag_3d
         CALL copydiag_3d (stashwork(si(item,sect,im_index)),           &
              t_incr_diagnostic,                                        &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         IF (icode >  0) THEN
            cmessage=": error in copydiag_3d(item 181)"//cmessage
         END IF

      END IF !  sf(item,sect)

! ----------------------------------------------------------------------
!  U component of stress
! ----------------------------------------------------------------------
! Item 201 stress_ud

      If (sf(201,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(201,6,im_index)),stress_ud,      &
     &        row_length,u_rows,model_levels+1,0,0,0,0,                 &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,201,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,201,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_ud)"
            goto 9999
         Endif

      End if

      IF (sf(241,6)) THEN
         stress_ud_p_levels = stash_levels(1,                           &
            -stlist(10,stindex(1,241,6,im_index)))

! Interpolate stress_ud on model levels/u points on the c grid onto
! pressure levels specified in the UMUI domain profile for stash 6115
! and uv positions (B grid). Copy interpolated diagnostic into STASHwork  
! for stash processing.
! DEPENDS ON: interp_2_press_at_uv_pos_b_grid
         CALL Interp_2_Press_At_UV_Pos_B_Grid(                          &
              im_index, 241,6, stress_ud,                               &
              row_length, rows, n_rows, stress_ud_p_levels,             &
              model_levels, off_x,off_y,                                &
              exner_theta_levels(:,:,1:tdims%k_end),                    &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                                          STASHwork(si(241,6,im_index)))
      ENDIF



! ----------------------------------------------------------------------
!  V component of stress
! ----------------------------------------------------------------------
! Item 202 stress_vd

      If (sf(202,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(202,6,im_index)),stress_vd,      &
     &        row_length,v_rows,model_levels+1,0,0,0,0,                 &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,202,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,202,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_vd)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  Orographic standard deviation
! ----------------------------------------------------------------------
! Item 203 Orographic standard deviation

      IF (sf(203,6)) THEN

         CALL expand_from_mask(sd_orog,sd_orog_land, land_sea_mask,        &
             row_length*rows,land_points)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(203,6,im_index)),                   &
     &        sd_orog,                                                  &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,6,203,                                           &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag(sd_orog)"
            goto 9999
         END IF

      END IF

! ----------------------------------------------------------------------
!  Orographic sigma_xx field
! ----------------------------------------------------------------------
! Item 204 Orographic sigma_xx field

      IF (sf(204,6)) THEN


         CALL expand_from_mask(orog_grad_xx,orog_grad_xx_land,          &
     &                land_sea_mask,row_length*rows,land_points)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(204,6,im_index)),                   &
     &        orog_grad_xx,                                             &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,6,204,                                           &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag(orog_grad_xx)"
            goto 9999
         END IF

      END IF

! ----------------------------------------------------------------------
!  Orographic sigma_xy field
! ----------------------------------------------------------------------
! Item 205 Orographic sigma_xy field

      IF (sf(205,6)) THEN

         CALL expand_from_mask(orog_grad_xy,orog_grad_xy_land,             &
             land_sea_mask,row_length*rows,land_points)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(205,6,im_index)),                   &
     &        orog_grad_xy,                                             &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,6,205,                                           &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag(orog_grad_xy)"
            goto 9999
         END IF

      END IF

! ----------------------------------------------------------------------
!  Orographic sigma_yy field
! ----------------------------------------------------------------------
! Item 206 Orographic sigma_yy field

      IF (sf(206,6)) THEN

         CALL expand_from_mask(orog_grad_yy,orog_grad_yy_land,             &
             land_sea_mask,row_length*rows,land_points)

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(206,6,im_index)),                   &
     &        orog_grad_yy,                                             &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,6,206,                                           &
     &        icode,cmessage)

         IF (icode  >   0) THEN
            cmessage="Error in copydiag(orog_grad_yy)"
            goto 9999
         END IF

      END IF

! ----------------------------------------------------------------------
!  U component of satn stress
! ----------------------------------------------------------------------
! Item 223 stress_ud_satn

      If (sf(223,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223,6,im_index)),stress_ud_satn, &
     &        row_length,u_rows,model_levels+1,0,0,0,0,                 &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,223,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,223,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_ud_satn)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  V component of satn stress
! ----------------------------------------------------------------------
! Item 224 stress_vd_satn

      If (sf(224,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(224,6,im_index)),stress_vd_satn, &
     &        row_length,v_rows,model_levels+1,0,0,0,0,                 &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,224,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,224,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_vd_satn)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  U component  of wake stress
! ----------------------------------------------------------------------
! Item 227 stress_ud_wake

      If (sf(227,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(227,6,im_index)),stress_ud_wake, &
     &        row_length,u_rows,model_levels+1,0,0,0,0,                 &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,227,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,227,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_ud_wake)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  V component of wake stress
! ----------------------------------------------------------------------
! Item 228 stress_vd_wake

      If (sf(228,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(228,6,im_index)),stress_vd_wake, &
     &        row_length,v_rows,model_levels+1,0,0,0,0,                 &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,228,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,228,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(stress_vd_wake)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
! du/dt saturation component
! ----------------------------------------------------------------------
! Item 207 du_dt_satn

      If (sf(207,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(207,6,im_index)),du_dt_satn,     &
     &        row_length,u_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,207,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,207,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(du_dt_satn)"
            goto 9999
         Endif

      End if


      IF (sf(247,6)) THEN
         du_dt_satn_p_levels = stash_levels(1,                          &
            -stlist(10,stindex(1,247,6,im_index)))

! Interpolate du_dt_satn on model levels/u points on the c grid onto
! pressure levels specified in the UMUI domain profile for stash 6115
! and uv positions (B grid). Copy interpolated diagnostic into STASHwork  
! for stash processing.
! DEPENDS ON: interp_2_press_at_uv_pos_b_grid
         CALL Interp_2_Press_At_UV_Pos_B_Grid(                          &
              im_index, 247,6, du_dt_satn,                              &
              row_length, rows, n_rows, du_dt_satn_p_levels,            &
              model_levels, off_x,off_y,                                &
              exner_theta_levels(:,:,1:tdims%k_end),                    &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                                          STASHwork(si(247,6,im_index)))
      ENDIF



! ----------------------------------------------------------------------
!  dv/dt saturation component
! ----------------------------------------------------------------------
! Item 208 dv_dt_satn

      If (sf(208,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(208,6,im_index)),dv_dt_satn,     &
     &        row_length,v_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,208,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,208,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(dv_dt_satn)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  du/dt wake component
! ----------------------------------------------------------------------
! Item 231 du_dt_wake

      If (sf(231,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(231,6,im_index)),du_dt_wake,     &
     &        row_length,u_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,231,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,231,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(du_dt_wake)"
            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
!  dv/dt wake component
! ----------------------------------------------------------------------
! Item 232 dv_dt_wake

      If (sf(232,6)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(232,6,im_index)),dv_dt_wake,     &
     &        row_length,v_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,232,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,232,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(dv_dt_wake)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  u latest.
! ----------------------------------------------------------------------
! Item 002 u

      If (sf(002,6)) Then

         Do k = 1, model_levels
            Do j = 1, u_rows
               Do i = 1, row_length
                 l = i + (j-1)*row_length + (k-1)*u_rows*row_length
                 interp_data_3_u(l) = u(i,j,k) + R_u(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(002,6,im_index)),interp_data_3_u, &
     &        row_length,u_rows,model_levels,0,0,0,0,                    &
     &        at_extremity,                                              &
     &        stlist(1,stindex(1,002,6,im_index)),len_stlist,            &
     &        stash_levels,num_stash_levels+1,                           &
     &        atmos_im,6,002,                                            &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(u+R_u)"
            goto 9999
         Endif

      End if



! ----------------------------------------------------------------------
!  v latest
! ----------------------------------------------------------------------
! Item 003 v

      If (sf(003,6)) Then

         Do k = 1, model_levels
            Do j = 1, v_rows
               Do i = 1, row_length
                 l = i + (j-1)*row_length + (k-1)*v_rows*row_length
                 interp_data_3_v(l) = v(i,j,k) + R_v(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(003,6,im_index)),interp_data_3_v, &
     &        row_length,v_rows,model_levels,0,0,0,0,                    &
     &        at_extremity,                                              &
     &        stlist(1,stindex(1,003,6,im_index)),len_stlist,            &
     &        stash_levels,num_stash_levels+1,                           &
     &        atmos_im,6,003,                                            &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(v+R_v)"
            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
!  U_S diag
! ----------------------------------------------------------------------
! Item 214 U_S_D

      If (sf(214,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(214,6,im_index)),                   &
     &        U_S_D,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,214,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(u_s)"
            goto 9999
         Endif

      End if
! ----------------------------------------------------------------------
!  V_S diag
! ----------------------------------------------------------------------
! Item 215 V_S_d

      If (sf(215,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(215,6,im_index)),                   &
     &        V_S_D,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,215,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(v_s)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  NSQ_S diag
! ----------------------------------------------------------------------
! Item 216 nsq_s_d

      If (sf(216,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(216,6,im_index)),                   &
     &        NSQ_S_D,                                                  &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,216,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(nsq_s)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  Froude number diag
! ----------------------------------------------------------------------
! Item 217 FR_d

      If (sf(217,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(217,6,im_index)),                   &
     &        FR_D,                                                     &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,217,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(Fr)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  Blocked Layer Depth diag
! ----------------------------------------------------------------------
! Item 218 BLD_d

      If (sf(218,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(218,6,im_index)),                   &
     &        BLD_D,                                                    &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,218,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(BLD)"
            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  % of time of blocked layer
! ----------------------------------------------------------------------
! Item 222 BLDT_D

      If (sf(222,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(222,6,im_index)),                   &
     &        BLDT_D,                                                   &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,222,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(BLDT)"

            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  % of time numerical limiter invoked (4A only)
! ----------------------------------------------------------------------
! Item 233 NUM_LIM_D

      If (sf(233,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(233,6,im_index)),                   &
     &        num_lim_d,                                                &
     &    row_length,rows,0,0,0,0, at_extremity,                        &
     &        atmos_im,6,233,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(num_lim)"

            goto 9999
         Endif

      End if



! ----------------------------------------------------------------------
!  % redn. of flow-blocking stress after numerical limiter invoked (4A only)
! ----------------------------------------------------------------------
! Item 234 NUM_FAC_D

      If (sf(234,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(234,6,im_index)),                   &
     &        num_fac_d,                                                &
     &    row_length,rows,0,0,0,0, at_extremity,                        &
     &        atmos_im,6,234,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(num_fac)"

            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
!  Orographic X-component of the total surface stress
!  Same as 6201, level 1 output, but here the pp headers are correctly
!  labelled
! ----------------------------------------------------------------------
! Item 235 TAUSX_D

      If (sf(235,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(235,6,im_index)),                   &
     &        tausx_d,                                                  &
     &        row_length,u_rows,0,0,0,0, at_extremity,                  &
     &        atmos_im,6,235,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(tausx)"

            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
!  Orographic Y-component of the total surface stress
!  Same as 6202, level 1 output, but here the pp headers are correctly
!  labelled
! ----------------------------------------------------------------------
! Item 236 TAUSY_D

      If (sf(236,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(236,6,im_index)),                   &
     &        tausy_d,                                                  &
     &        row_length,v_rows,0,0,0,0, at_extremity,                  &
     &        atmos_im,6,236,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(tausy)"

            goto 9999
         Endif

      End if


! ----------------------------------------------------------------------
!  Orographic surface stress Froude number dependent scaling factor (4A only)
! ----------------------------------------------------------------------
! Item 237 TAUS_SCALE_D

      If (sf(237,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(237,6,im_index)),                   &
     &        taus_scale_d,                                             &
     &    row_length,rows,0,0,0,0, at_extremity,                        &
     &        atmos_im,6,237,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(taus_scale)"

            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  Slope of sub-grid scale orography (5A scheme only)
! ----------------------------------------------------------------------
! Item 248 orog_slope_d

      If (sf(248,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(248,6,im_index)),                   &
     &        orog_slope_d,                                             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,248,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(OROG_SLOPE_D)"

            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  Anisotropy of sub-grid scale orography (5A scheme only)
! ----------------------------------------------------------------------
! Item 249 orog_anis_d

      If (sf(249,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(249,6,im_index)),                   &
     &        orog_anis_d,                                              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,249,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(OROG_ANIS_D)"

            goto 9999
         Endif

      End if

! ----------------------------------------------------------------------
!  Orientation of sub-grid scale orography (5A scheme only)
! ----------------------------------------------------------------------
! Item 250 orog_dir_d

      If (sf(250,6)) Then
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(250,6,im_index)),                   &
     &        orog_dir_d,                                               &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,6,250,                                           &
     &        icode,cmessage)


         If (icode  >   0) Then
            cmessage="Error in copydiag(OROG_DIR_D)"

            goto 9999
         Endif

      End if

!---------------------------------------------------------------------
! Section 2.
! SPECTRAL (NON-OROGRAPHIC) GRAVITY WAVE SCHEME DIAGNOSTICS
!
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!  USSP Eastward flux component
! --------------------------------------------------------------------
! Item 101 GWSPEC_EFLUX(i,j,k)
      If (sf(101,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(101,6,im_index)),GWSPEC_EFLUX,   &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,101,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,101,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_EFLUX)"
            goto 9999
         Endif
      Endif

      IF (sf(111,6)) THEN
         gwspec_eflux_p_levels = stash_levels(1,                        &
            -stlist(10,stindex(1,111,6,im_index)))

! Interpolate GWSPEC_EFLUX on model levels/u points on the c grid onto
! pressure levels specified in the UMUI domain profile for stash 6115
! and uv positions (B grid). Copy interpolated diagnostic into STASHwork  
! for stash processing.
! DEPENDS ON: interp_2_press_at_uv_pos_b_grid
         CALL Interp_2_Press_At_UV_Pos_B_Grid(                          &
              im_index, 111,6, GWSPEC_EFLUX,                            &
              row_length, rows, n_rows, gwspec_eflux_p_levels,          &
              model_levels, off_x,off_y,                                &
              exner_theta_levels(:,:,1:tdims%k_end),                    &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                                          STASHwork(si(111,6,im_index)))
      ENDIF
! --------------------------------------------------------------------
!  USSP Southward flux component
! --------------------------------------------------------------------
! Item 102 GWSPEC_SFLUX(i,j,k)
      If (sf(102,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(102,6,im_index)),GWSPEC_SFLUX,   &
     &        row_length,n_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,102,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,102,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_SFLUX)"
            goto 9999
         Endif
      Endif
! --------------------------------------------------------------------
!  USSP Westward flux component
! --------------------------------------------------------------------
! Item 103 GWSPEC_WFLUX(i,j,k)
      If (sf(103,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(103,6,im_index)),GWSPEC_WFLUX,   &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,103,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,103,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_WFLUX)"
            goto 9999
         Endif
      Endif

      IF (sf(113,6)) THEN
         gwspec_wflux_p_levels = stash_levels(1,                        &
            -stlist(10,stindex(1,113,6,im_index)))

! Interpolate GWSPEC_WFLUX on model levels/u points on the c grid onto
! pressure levels specified in the UMUI domain profile for stash 6115
! and uv positions (B grid). Copy interpolated diagnostic into STASHwork  
! for stash processing.
! DEPENDS ON: interp_2_press_at_uv_pos_b_grid
         CALL Interp_2_Press_At_UV_Pos_B_Grid(                          &
              im_index, 113,6, GWSPEC_WFLUX,                            &
              row_length, rows, n_rows, gwspec_wflux_p_levels,          &
              model_levels, off_x,off_y,                                &
              exner_theta_levels(:,:,1:tdims%k_end),                    &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                                          STASHwork(si(113,6,im_index)))
      ENDIF



! --------------------------------------------------------------------
!  USSP Northward flux component
! --------------------------------------------------------------------
! Item 104 GWSPEC_NFLUX(i,j,k)
      If (sf(104,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(104,6,im_index)),GWSPEC_NFLUX,   &
     &        row_length,n_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,104,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,104,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_NFLUX)"
            goto 9999
         Endif
      Endif
! --------------------------------------------------------------------
!  USSP EW acceleration component
! --------------------------------------------------------------------
! Item 105 GWSPEC_EWACC(i,j,k)
      If (sf(105,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(105,6,im_index)),GWSPEC_EWACC,   &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,105,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,105,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_EWACC)"
            goto 9999
         Endif
      Endif

      IF (sf(115,6)) THEN
         gwspec_ewacc_p_levels = stash_levels(1,                        &
            -stlist(10,stindex(1,115,6,im_index)))

! Interpolate GWSPEC_EWACC on model levels/u points on the c grid onto
! pressure levels specified in the UMUI domain profile for stash 6115
! and uv positions (B grid). Copy interpolated diagnostic into STASHwork  
! for stash processing.
! DEPENDS ON: interp_2_press_at_uv_pos_b_grid
         CALL Interp_2_Press_At_UV_Pos_B_Grid(                          &
              im_index, 115,6, GWSPEC_EWACC,                            &
              row_length, rows, n_rows, gwspec_ewacc_p_levels,          &
              model_levels, off_x,off_y,                                &
              exner_theta_levels(:,:,1:tdims%k_end),                    &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                                          STASHwork(si(115,6,im_index)))
      ENDIF



! --------------------------------------------------------------------
!  USSP NS acceleration component
! --------------------------------------------------------------------
! Item 106 GWSPEC_NSACC(i,j,k)
      If (sf(106,6)) Then
! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(106,6,im_index)),GWSPEC_NSACC,   &
     &        row_length,n_rows,model_levels,0,0,0,0,                   &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,106,6,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,6,106,                                           &
     &        icode,cmessage)
         If (icode  >   0) Then
            cmessage="Error in copydiag_3d(GWSPEC_NSACC)"
            goto 9999
         Endif
      Endif

 9999 continue                  !  single point error handling.
      If(icode /= 0) Then

        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      IF (lhook) CALL dr_hook('DIAGNOSTICS_GWD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE diagnostics_gwd
