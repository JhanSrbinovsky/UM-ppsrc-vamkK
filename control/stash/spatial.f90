! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: SPATIAL --------------------------------------------------
!LL
!LL  Purpose: Performs general spatial processing on an input field to
!LL           produce an output field or scalar within STASH.  Lower-
!LL           level routines are called to perform the various spatial
!LL           processing options.
!LL
!LL
!LL  Programming standard: UM Doc Paper 3
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!LL-----------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

      SUBROUTINE SPATIAL(fieldin,vx,vy,vz,GR,st_grid,                   &
     &                   fld_type,halo_type,                            &
     &                   halo_x,halo_y,                                 &
     &                   lcyclic,lmasswt,                               &
     &      n_cols_out,n_rows_out,                                      &
     &      this_index_lev,level_list,index_lev,no_of_levels,           &
     &      no_of_levels_masswt,                                        &
     &      p,pstar,                                                    &
     &      cos_v_latitude,cos_theta_latitude,land,sea,                 &
     &      row_length,rows,n_rows,no_rows,model_levels,                &
     &      fieldout,lenout,                                            &
     &      control,control_size,rmdi,                                  &
     &      icode,cmessage)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE sterr_mod, ONLY: st_upper_less_lower, st_not_supported,       &
                           st_no_data,st_nd, st_bad_array_param,        &
                           st_bad_address, st_unknown,                  &
                           st_bad_wraparound, st_illegal_weight,        &
                           unknown_weight, unknown_mask,                &
                           unknown_processing, nonsense
      IMPLICIT NONE
!

      INTEGER                                                           &
     &    vx,vy,vz,                                                     &
                                            ! IN size of fieldin
     &    lenout,                                                       &
                                            ! IN size of fieldout
     &    GR,                                                           &
                                            ! IN ppxref gridtype code
     &    st_grid,                                                      &
                                            ! IN STASH gridtype code
     &    fld_type,                                                     &
                                            ! IN field type (u/v/p)
     &    halo_type,                                                    &
                                            ! IN halo type
     &    halo_x,                                                       &
                                            ! IN EW halo of input
     &    halo_y,                                                       &
                                            ! IN NS halo of input
     &    n_rows_out,                                                   &
                                            ! OUT no. of output rows
     &    n_cols_out,                                                   &
                                            ! OUT no. of output cols
     &    this_index_lev,                                               &
                                            ! IN level index, this field
     &    row_length,rows,n_rows,                                       &
                                            ! IN horiz. sizes (C grid)
     &    no_rows,                                                      &
                                            ! IN number of rows used
     &    model_levels,                                                 &
                                            ! IN vertical size
     &    control_size,                                                 &
                                            ! IN size of control record
     &    control(control_size),                                        &
                                            ! IN control record
     &    icode,                                                        &
                                            ! OUT error code 0 if ok
     &    no_of_levels,                                                 &
                                            ! IN no of levels
     &    no_of_levels_masswt,                                          &
                                ! IN levels for mass weighting array
                                ! lmasswt F: =1; T: =no_of_levels
     &    index_lev(no_of_levels),                                      &
                                            ! IN index to levels
     &    level_list(no_of_levels)          ! IN model level list
      REAL                                                              &
     &    fieldin(vx,vy,vz),                                            &
                             ! IN fieldin which is to be acted on
     &    p(1-offx:row_length+offx,1-offy:rows+offy,model_levels),      &
                                            ! IN pressure (rho levels)
     &    pstar(row_length,rows),                                       &
                                            ! IN surface pressure
     &    cos_v_latitude(row_length,n_rows),                            &
                                            ! IN v-grid area fn
     &    cos_theta_latitude(row_length,rows),                          &
                                                ! IN T-grid area fn
     &    fieldout(lenout),                                             &
                                                ! OUT output field
     &    rmdi                                  ! IN  missing data indic
      LOGICAL                                                           &
     &    lcyclic,                                                      &
                                                ! IN .true. if cyclic EW
     &    lmasswt,                                                      &
                                                ! IN  TRUE if masswts OK
     &    land(row_length,rows),                                        &
                                                ! IN land mask
     &    sea(row_length,rows)                  ! IN sea mask
      CHARACTER(LEN=*) cmessage                    ! OUT error message

!*----------------------------------------------------------------------
!
! STPARAM
!
!  Purpose: Meaningful PARAMETER names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that first the address of the item is declare in
!           both long and short form. example is;
!             integer st_item_code,s_item  !Item number (declaration)
!             parameter(st_item_code=3,s_item=3)
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D70
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!--------------------------------------------------------------

      ! Internal model number address
      INTEGER,PARAMETER:: st_model_code = 28
      INTEGER,PARAMETER:: s_modl        = 28

      ! Section Number address
      INTEGER,PARAMETER:: st_sect_no_code = 2
      INTEGER,PARAMETER:: s_sect          = 2
      INTEGER,PARAMETER:: st_sect_code    = 2

      INTEGER,PARAMETER:: st_item_code=1,s_item=1 ! Item number address

      ! Processing Code address
      INTEGER,PARAMETER:: st_proc_no_code=3,s_proc=3

      ! subsidiary codes for st_proc_no_code now
      INTEGER,PARAMETER:: st_replace_code=1
      INTEGER,PARAMETER:: st_accum_code=2
      INTEGER,PARAMETER:: st_time_mean_code=3
      INTEGER,PARAMETER:: st_time_series_code=4
      INTEGER,PARAMETER:: st_max_code=5
      INTEGER,PARAMETER:: st_min_code=6
      INTEGER,PARAMETER:: st_append_traj_code=7
      INTEGER,PARAMETER:: st_time_series_mean=8
      INTEGER,PARAMETER:: st_variance_code=9

      ! Frequency (Input & output) addres
      INTEGER,PARAMETER:: st_freq_code=4,s_freq=4

      ! Offset for sampling
      INTEGER,PARAMETER:: st_offset_code=30,s_offs=30

      ! start timestep address
      INTEGER,PARAMETER:: st_start_time_code=5,s_times=5

      ! end timestep address
      INTEGER,PARAMETER:: st_end_time_code=6,s_timee=6

      ! period in timesteps address
      INTEGER,PARAMETER:: st_period_code=7,s_period=7

      ! infinite end/period value
      INTEGER,PARAMETER:: st_infinite_time=-1

      INTEGER,PARAMETER:: st_end_of_list=-1 !end-of-list marker in times

      ! grid point stuff
      ! gridpoint info address
      INTEGER,PARAMETER:: st_gridpoint_code=8,s_grid=8

      ! now subsid grid point stuff
      ! no masking done
      INTEGER,PARAMETER:: stash_null_mask_code=1,s_nomask=1

      ! land mask conds
      INTEGER,PARAMETER:: stash_land_mask_code=2,s_lndms=2

      ! sea mask code
      INTEGER,PARAMETER:: stash_sea_mask_code=3,s_seams =3

      ! processing options

      ! size of block for gridpoint code
      INTEGER,PARAMETER:: block_size=10

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: extract_base=block_size*0

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: extract_top=block_size*1

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_base=block_size*1

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_top=block_size*2

      ! max code for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_base=block_size*2

      ! base codes for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_top=block_size*3

      ! max code for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_base=block_size*3

      ! base codes for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_top=block_size*4

      ! max code for field mean subroutine
      INTEGER,PARAMETER:: field_mean_base=block_size*4

      ! base codes for field mean subroutine
      INTEGER,PARAMETER:: field_mean_top=block_size*5

      ! max code for global mean subroutine
      INTEGER,PARAMETER:: global_mean_base=block_size*5

      ! base codes for global mean subroutine
      INTEGER,PARAMETER:: global_mean_top=block_size*6

      ! Weighting

      ! weighting info address
      INTEGER,PARAMETER:: st_weight_code=9,s_weight=9

      INTEGER,PARAMETER:: stash_weight_null_code  =0,s_noweight  =0
      INTEGER,PARAMETER:: stash_weight_area_code  =1,s_areaweight=1
      INTEGER,PARAMETER:: stash_weight_volume_code=2,s_volweight =2
      INTEGER,PARAMETER:: stash_weight_mass_code  =3,s_massweight=3

      ! Domain definition

      ! row addresses
      INTEGER,PARAMETER:: st_north_code=12,s_north=12
      INTEGER,PARAMETER:: st_south_code=13,s_south=13
      INTEGER,PARAMETER:: st_west_code =14,s_west =14
      INTEGER,PARAMETER:: st_east_code =15,s_east =15

      ! Levels

      ! input bottom level address
      INTEGER,PARAMETER:: st_input_bottom=10,s_bottom =10

      ! special code
      INTEGER,PARAMETER:: st_special_code=100,s_special=100

      ! input top level address
      INTEGER,PARAMETER:: st_input_top=11,s_top=11

      ! output bottom level address
      INTEGER,PARAMETER:: st_output_bottom=21,s_outbot=21

      ! output top level address
      INTEGER,PARAMETER:: st_output_top=22,s_outtop=22

      INTEGER,PARAMETER:: st_model_level_code=1,s_model=1

      ! code for pressure leve
      INTEGER,PARAMETER:: st_pressure_level_code=2,s_press=2

      ! code for height levels
      INTEGER,PARAMETER:: st_height_level_code=3,s_height=3

      ! input code addres
      INTEGER,PARAMETER:: st_input_code=16,s_input=16

      ! input length of diagnostic address
      INTEGER,PARAMETER:: st_input_length=17,s_length=17

      ! output code address
      INTEGER,PARAMETER:: st_output_code=18,s_output=18

      ! Pointer to D1 addressing information
      ! Pos of item in D1 for relevant submodel
      INTEGER,PARAMETER:: st_position_in_d1=29,st_d1pos=29

      ! Output destination options

      INTEGER,PARAMETER:: st_dump=1
      INTEGER,PARAMETER:: st_secondary=2

      ! output length of diagnostic address
      INTEGER,PARAMETER:: st_output_length=19,s_outlen=19
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

      ! ptr to dump lookup header address
      INTEGER,PARAMETER:: st_lookup_ptr=23

      ! ptr into stash_series where control data address
      INTEGER,PARAMETER:: st_series_ptr=24

      ! subsid stuff for time series
      INTEGER,PARAMETER:: series_grid_type=1
      INTEGER,PARAMETER:: series_grid_code=0
      INTEGER,PARAMETER:: series_long_code=1
      INTEGER,PARAMETER:: series_size=2
      INTEGER,PARAMETER:: series_proc_code=3
      INTEGER,PARAMETER:: series_north=4
      INTEGER,PARAMETER:: series_south=5
      INTEGER,PARAMETER:: series_west=6
      INTEGER,PARAMETER:: series_east=7
      INTEGER,PARAMETER:: series_list_start=8
      INTEGER,PARAMETER:: series_list_end=9
      INTEGER,PARAMETER:: record_size=9

      ! Miscellaneous parameters

      ! system/user tag field in stlist address
      INTEGER,PARAMETER:: st_macrotag=25

      ! Pseudo-level list pointers

      ! pseudo-levels input list address
      INTEGER,PARAMETER:: st_pseudo_in=26

      ! pseudo-levels output list address
      INTEGER,PARAMETER:: st_pseudo_out=27

      ! Internal horizontal gridtype codes common to all diagnostics

      INTEGER,PARAMETER:: st_tp_grid =1 ! T-p grid
      INTEGER,PARAMETER:: st_uv_grid =2 ! u-v grid
      INTEGER,PARAMETER:: st_cu_grid =3 ! C-grid u point
      INTEGER,PARAMETER:: st_cv_grid =4 ! C-grid v point
      INTEGER,PARAMETER:: st_zt_grid =5 ! Zonal T-grid
      INTEGER,PARAMETER:: st_zu_grid =6 ! Zonal u-grid
      INTEGER,PARAMETER:: st_mt_grid =7 ! Meridional T-grid
      INTEGER,PARAMETER:: st_mu_grid =8 ! Meridional u-grid
      INTEGER,PARAMETER:: st_riv_grid= 23    ! river_routing grid
      INTEGER,PARAMETER:: st_scalar  =9 ! Scalar (ie. single value)
      INTEGER,PARAMETER:: st_wam_all= 60    ! Wam Field on Full Grid
      INTEGER,PARAMETER:: st_wam_sea= 62    ! Wam Field on Sea Points

! STPARAM end
!L
!L external routines
!L
      EXTERNAL stextc ! extracts the field
      EXTERNAL stcolm ! computes the column mean
      EXTERNAL stzonm ! computes the zonal mean
      EXTERNAL stmerm ! computes the meridional mean
      EXTERNAL stglom ! computes the global mean
      EXTERNAL stfieldm ! computes the field mean
!L
!L local variables
!L
      LOGICAL lwrap                ! TRUE if output field wraparound EW
      LOGICAL lmasswt_strict       ! copy of lmasswt - but set to false
!                                  ! if mass weighting is not requested
      INTEGER xstart,ystart        ! lower left hand corner coords
      INTEGER xend,yend            ! upper right hand corner coords
      INTEGER processing_code      ! what kind of mean  will be done
      INTEGER what_level           ! what type of input level
      INTEGER what_mask            ! what mask is used
      INTEGER what_weight          ! what weighting is used

      INTEGER i,j,k                                                     &
                                   ! loop counters
     &,model_level                                                      &
                                   ! model level
     &,this_index_lev_masswt       ! level index for mass weighting
                                   ! (=1 if no mass weighting or no
                                   !  model level weights available)

      INTEGER                                                           &
! global versions of the extracted area domain limits
     &  global_xstart,global_xend,global_ystart,global_yend

! workspace arrays containining weighting factors and masks.
      REAL                                                              &
     &  area_weight(1-halo_x:row_length+halo_x+1,                       &
     &              1-halo_y:no_rows+halo_y)                            &
     &, pstar_weight(1-halo_x:row_length+halo_x+1,                      &
     &               1-halo_y:no_rows+halo_y,                           &
     &               no_of_levels_masswt)                               &
     &, pstar_interp(row_length,no_rows)
      LOGICAL                                                           &
     &  mask(1-halo_x:row_length+halo_x+1,                              &
     &       1-halo_y:no_rows+halo_y)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!L----------------------------------------------------------------------
!L 1. Set up local variables
!L
      IF (lhook) CALL dr_hook('SPATIAL',zhook_in,zhook_handle)
      xstart=control(st_west_code)
      xend=control(st_east_code)
      ystart=control(st_south_code)  ! NOTE: Grid is assumed to be
      yend=control(st_north_code)    !       oriented south-to-north

      global_xstart=xstart
      global_ystart=ystart
      global_xend=xend
      global_yend=yend

! and calculate what the local subdomain limits are:
! DEPENDS ON: global_to_local_subdomain
        CALL GLOBAL_TO_LOCAL_SUBDOMAIN( .TRUE.,.TRUE.,                  &
     &                                  GR,halo_type,mype,              &
     &                                  global_ystart,global_xend,      &
     &                                  global_yend,global_xstart,      &
     &                                  ystart,xend,yend,xstart)

! Check if wraparound field
      IF (xstart >  xend) THEN

        IF (lcyclic) THEN
          xend=xend+blsize(1,fld_type)
! subtract two halos as we don't wish to include halos at the end
! and start of the row within the wrap around domain
          lwrap=.TRUE.
        ELSE
          icode=st_bad_wraparound    ! wraparound illegal unless cyclic
          GOTO 999
        ENDIF

      ELSE
        lwrap=.FALSE.
      ENDIF

      IF (global_xstart  >   global_xend) THEN
        IF (lcyclic) THEN
          global_xend=global_xend+glsize(1,fld_type)
        ELSE
          icode=st_bad_wraparound  ! wraparound illegal unless cyclic
          GOTO 999
        ENDIF
      ENDIF

      processing_code=control(st_gridpoint_code)
      what_level=control(st_input_bottom)
      what_mask=mod(processing_code,block_size)
      what_weight=control(st_weight_code)
!L
!L 1.1 Prevent masking or weighting if input field is not 2D in extent
!L     - weighting and masking is assumed to have been done outside.
!L
      IF ( (.NOT.(st_grid == st_tp_grid.OR.st_grid == st_uv_grid.OR.    &
     &            st_grid == st_cu_grid.OR.st_grid == st_cv_grid.OR.    &
     &             st_grid == st_riv_grid))                             &
     &      .AND.(what_mask   /= stash_null_mask_code   .OR.            &
     &            what_weight /= stash_weight_null_code) ) THEN
       icode=st_not_supported
       cmessage='SPATIAL : Masking/weighting unsupported - non 2D field'
       GOTO 999
      ENDIF

! Check for supported weighting and masking options

      IF (.NOT. ((what_weight  ==  stash_weight_null_code) .OR.         &
     &           (what_weight  ==  stash_weight_area_code) .OR.         &
     &           (what_weight  ==  stash_weight_volume_code) .OR.       &
     &           (what_weight  ==  stash_weight_mass_code) ) ) THEN
        cmessage='SPATIAL : Unrecognized weighting option'
        icode=unknown_weight
        GOTO 999
      ENDIF

      IF (.NOT. ((what_mask  ==  stash_null_mask_code) .OR.             &
     &           (what_mask  ==  stash_land_mask_code) .OR.             &
     &           (what_mask  ==  stash_sea_mask_code ) ) ) THEN
        cmessage='SPATIAL : Unrecognized masking option'
        icode=unknown_mask
        GOTO 999
      ENDIF

      IF (what_weight  ==  stash_weight_volume_code) THEN
        cmessage='SPATIAL : Volume-weighting not supported'
        icode=st_illegal_weight
        GOTO 999
      ENDIF

! Set lmasswt_strict - copy of lmasswt, but set to false is mass
! weighting not requested

      lmasswt_strict=                                                   &
     &  (lmasswt .AND. (what_weight  ==  stash_weight_mass_code))

! Precalculate weighting and masking arrays
! I've used IF tests inside the loops, but since the logical
! expressions are invariant wrt i and j, the compiler will
! move them outside the DO loops. It makes the code a lot shorter!

! NOTE that neither area-weights or mass-weights are normalised so
! that the interpretation of weighted diagnostics is non-obvious. Also
! the PP lookup header has no switch for indicating whether or not such
! processing has taken place. More careful treatment of horizontal
! interpolation is required for stricter accuracy.


! area weighting
      IF (what_weight  ==  stash_weight_null_code) THEN
          ! Ensure initialisation of weight array including halos
          area_weight (:,:)   = 1.0
      ELSE ! some form of area weighting will be required
          DO j=1,no_rows
          DO i=1,row_length
             IF (st_grid  ==  st_cv_grid) THEN
                 area_weight(i,j)=cos_v_latitude(i,j)
             ELSE
! NOTE that this is will not be accurate for C-u grid variables for
! LAMs, since cos_theta_latitude will differ between theta,u positions
                 area_weight(i,j)=cos_theta_latitude(i,j)
             ENDIF
         ENDDO ! i
         ENDDO ! j

      ENDIF    ! what_weight


! mass weighting
      IF ((what_weight  ==  stash_weight_null_code) .OR.                &
     &    (what_weight  ==  stash_weight_area_code)) THEN
! No mass weighting is required
          this_index_lev_masswt = 1

          ! Ensure initialisation of weight array including halos
          pstar_weight(:,:,this_index_lev_masswt) = 1.0
      ELSE

! Mass weighting requested

! Ensure that halos are initialised
        DO j=no_rows-1,no_rows
        DO i=1,row_length
           pstar_interp(i,j) =1.0
        ENDDO !i
        ENDDO !j

! Interpolate pstar to correct horizontal staggering
        IF     (st_grid  ==  st_cu_grid) THEN
! NOT YET CORRECT: pstar has no halos! So set to nearby value
!          CALL P_TO_U(pstar,row_length,rows,1,0,0,pstar_interp)
           DO j=1,no_rows
           DO i=1,row_length
              pstar_interp(i,j) =pstar(i,j)
           ENDDO !i
           ENDDO !j

        ELSEIF (st_grid  ==  st_cv_grid) THEN
! NOT YET CORRECT: pstar has no halos! So set to nearby value
!          CALL P_TO_V(pstar,row_length,rows,1,0,0,pstar_interp)
           DO j=1,no_rows
           DO i=1,row_length
              pstar_interp(i,j) =pstar(i,j)
           ENDDO !i
           ENDDO !j

        ELSE
          DO j=1,no_rows
            DO i=1,row_length
              pstar_interp(i,j)=pstar(i,j)
            ENDDO
          ENDDO
        ENDIF

        IF(lmasswt) THEN  ! model level mass weights available

          this_index_lev_masswt = this_index_lev

          DO k=1,no_of_levels_masswt
            model_level = level_list(k)
            IF(model_level == model_levels) THEN  ! top level

              DO j=1,no_rows
              DO i=1,row_length
                pstar_weight(i,j,k) = p(i,j,model_level)
              ENDDO !i
              ENDDO !j

            ELSE      ! not top level

              DO j=1,no_rows
              DO i=1,row_length
! Only accurate for variables on theta levels
                pstar_weight(i,j,k) = p(i,j,model_level) -              &
     &                                  p(i,j,model_level+1)
              ENDDO !i
              ENDDO !j
            ENDIF

          ENDDO ! k

        ELSE              ! no model level mass weights available:
                          ! weight by pstar

            this_index_lev_masswt = 1

            DO j=1,no_rows
            DO i=1,row_length
               pstar_weight(i,j,this_index_lev_masswt) =                &
     &                                               pstar_interp(i,j)
            ENDDO !i
            ENDDO !j
! Horizontal interpolation may be required - this should be done here:

        ENDIF   ! lmasswt
      ENDIF

! masking

      DO j=1,no_rows
        DO i=1,row_length
          IF (what_mask  ==  stash_land_mask_code) THEN
            mask(i,j)=land(i,j)
          ELSEIF (what_mask  ==  stash_sea_mask_code) THEN
            mask(i,j)=sea(i,j)
          ELSE
            mask(i,j)=.TRUE.
          ENDIF
        ENDDO
      ENDDO

! Update the halo areas of the weighting/masking arrays

!  [Note that for lams at UM5.2 and before, valid values for rim
!   boundaries would be needed using FILL_EXTERNAL_HALOS calls for
!   area_weight and pstar_weight arrays, but this should now be
!   superseded by initialising full arrays]

      ! Update halos only if halos present (standard diagnostics have no
      ! halos, and if some weighting required (probably defunct
      ! functionality)
      IF( (halo_x  /=  0 .OR. halo_y  /=  0) .AND.                      &
     &    (what_weight  /=  stash_weight_null_code .OR.                 &
     &     what_mask    /=  stash_null_mask_code )     ) THEN

! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS( AREA_WEIGHT(1-halo_x:row_length+halo_x,:),      &
     &     ROW_LENGTH,NO_ROWS,1,halo_x,halo_y,                          &
     &                fld_type_p,.FALSE.)
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS( PSTAR_WEIGHT(1-halo_x:row_length+halo_x,:,:),   &
     &     ROW_LENGTH,NO_ROWS,                                          &
     &                no_of_levels_masswt,halo_x,halo_y,                &
     &                fld_type_p,.FALSE.)
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS( MASK(1-halo_x:row_length+halo_x,:),             &
     &     ROW_LENGTH,NO_ROWS,1,halo_x,halo_y,                          &
     &                fld_type_p,.FALSE.)
      ENDIF
!L----------------------------------------------------------------------
!L 2. Call service routine to perform required processing
!L
!L 2.1 Extract sub field (single level at a time)
!L
      IF (processing_code <  extract_top.and.                           &
     &    processing_code >= extract_base) THEN
        n_rows_out=(yend+1)-ystart
        n_cols_out=(xend+1)-xstart

        IF (                                                            &
     &   (xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.    &
     &   (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN

! DEPENDS ON: stextc
        CALL STEXTC(fieldin,vx,vy,fld_type,halo_type,                   &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)

        ELSE  ! just set values to non NaN
          DO i=1,lenout
            fieldout(i)=0.0
          ENDDO
        ENDIF

!L
!L 2.2 Calculate column mean (over multiple levels indexed by index_lev)
!L
      ELSEIF (processing_code <  vert_mean_top.and.                     &
     &        processing_code >  vert_mean_base) THEN
        n_rows_out=yend+1-ystart
        n_cols_out=xend+1-xstart

        IF (                                                            &
     &   (xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.    &
     &   (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN

! DEPENDS ON: stcolm
        CALL STCOLM(fieldin,vx,vy,vz,fld_type,halo_type,                &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              fieldout,index_lev,no_of_levels,                    &
     &              pstar_weight,                                       &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)

        ELSE  ! just set values to non NaN
          DO i=1,lenout
            fieldout(i)=0.0
          ENDDO
        ENDIF

!L
!L 2.3 Calculate zonal mean (single level at a time)
!L
      ELSEIF (processing_code <  zonal_mean_top.and.                    &
     &        processing_code >  zonal_mean_base) THEN
        n_rows_out=yend+1-ystart
        n_cols_out=1
! DEPENDS ON: stzonm
        CALL STZONM(fieldin,vx,vy,fld_type,gr,halo_type,                &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.4 Calculate meridional mean (single level at a time)
!L
      ELSEIF (processing_code <  merid_mean_top.and.                    &
     &        processing_code >  merid_mean_base) THEN
        n_rows_out=1
        n_cols_out=xend+1-xstart
! DEPENDS ON: stmerm
        CALL STMERM(fieldin,vx,vy,fld_type,gr,halo_type,                &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.5 Calculate field mean (single level at a time)
!L
      ELSEIF (processing_code <  field_mean_top.and.                    &
     &        processing_code >  field_mean_base) THEN
        n_rows_out=1
        n_cols_out=1
! DEPENDS ON: stfieldm
        CALL STFIELDM(fieldin,vx,vy,fld_type,gr,halo_type,              &
     &                lwrap,lmasswt_strict,                             &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.6 Calculate global mean (over multiple levels)
!L
      ELSEIF (processing_code <  global_mean_top.and.                   &
     &        processing_code >  global_mean_base) THEN
        n_rows_out=1
        n_cols_out=1
! DEPENDS ON: stglom
        CALL STGLOM(fieldin,vx,vy,vz,fld_type,gr,halo_type,             &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,index_lev,no_of_levels,                    &
     &              pstar_weight,                                       &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.7 Invalid processing option
!L
      ELSE
        icode=unknown_processing
        write(cmessage,111)'unknown processing option',                 &
     &    processing_code
      ENDIF
!L
  999 CONTINUE
111   format('SPATIAL : >>FATAL ERROR <<',a40,i5)
!
      IF (lhook) CALL dr_hook('SPATIAL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SPATIAL
