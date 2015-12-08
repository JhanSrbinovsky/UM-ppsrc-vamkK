! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: MULTI_SPATIAL --------------------------------------------
!
! Purpose: Control routine for spatial processing when extracting a
!          timeseries within STASH.  Calls SPATIAL to extract global
!          mean samples from subdomains pointed to by the mother
!          STASHlist record using weighting and masking codes from
!          the STASHlist record within each subdomain.  The list of
!          subdomains is held as part of the STASH control file.
!          They may be in terms of gridpoints, or latitude/longitude
!          relative to the grid coordinates.  All timeseries samples
!          at a given step are appended to the output field.
!          On the first timestep it fills the entire output
!           vector to missing data (apart from values for the
!           first timestep and computes extra data for the time-
!           series -- tis prevents time meaning routines
!           failing due to uninitialised data.
!           however as a result of this the output vector length
!           will change from timestep to timestep
!
! External documentation:
!   Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                system (STASH)
!
! Interface and arguments: ------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: STASH
      SUBROUTINE MULTI_SPATIAL(fieldin,vx,vy,vz,gr,st_grid,             &
     &                         fld_type,halo_type,                      &
     &                         halo_x,halo_y,                           &
     &                         lcyclic,                                 &
     &     lmasswt,horiz_size,num_vert_levels,                          &
     &     no_of_levels_masswt,                                         &
     &     p,pstar,                                                     &
     &     cos_v_latitude,cos_theta_latitude,land,sea,                  &
     &     row_length,rows,n_rows,model_levels,                         &
     &     fieldout,lenout,amdi,                                        &
     &     control,control_size,                                        &
     &     stash_series,series_entry_size,no_records,                   &
     &     num_stash_levels,index_lev,level_list,start_ts,extraw,       &
     &     n_rows_out,n_cols_out,                                       &
     &     real_hd,len_realhd,int_hd,len_inthd,                         &
     &     icode,errmssg)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE gcom_mod, ONLY: gc_none
      USE UM_ParVars
      USE sterr_mod, ONLY: st_upper_less_lower, st_not_supported,       &
                           st_no_data,st_nd, st_bad_array_param,        &
                           st_bad_address, st_unknown,                  &
                           st_bad_wraparound, st_illegal_weight,        &
                           unknown_weight, unknown_mask,                &
                           unknown_processing, nonsense
      IMPLICIT NONE


      INTEGER vx,vy,vz          ! IN size of fieldin
      INTEGER gr                ! IN ppxref gridtype code
      INTEGER st_grid           ! IN STASH internal gridtype code
      INTEGER                                                           &
     & fld_type                                                         &
                               ! IN field type (u/v/p) of input field
     &,halo_type                                                        &
                               ! IN halo type of input field
     &,halo_x                                                           &
                               ! IN halo size East-West
     &,halo_y                  ! IN halo size North-South
      LOGICAL lcyclic           ! IN TRUE if cyclic EW BCs
      LOGICAL lmasswt           ! IN TRUE if level-dep mass-wts OK
      INTEGER row_length,rows,n_rows,model_levels ! IN size parameters
      REAL fieldin(vx*vy*vz)    ! IN data field
      INTEGER horiz_size        ! OUT no. of points in horizontal slice
      INTEGER num_vert_levels   ! OUT no of horizontal slices.
      INTEGER  no_of_levels_masswt ! IN: levels for mass weighting array
      INTEGER num_stash_levels  ! IN size of vertical levels list.
      INTEGER index_lev(num_stash_levels) ! IN offsets for each horiz fi
      INTEGER level_list(num_stash_levels) ! IN level for each horiz. fi

      REAL                                                              &
     & p(1-offx:row_length+offx,1-offy:rows+offy,model_levels)          &
                                              ! IN pressure (rho levels)
     &,pstar(row_length,rows)                                           &
                                              ! IN surface pressure
     &,cos_v_latitude(row_length,n_rows)                                &
                                              ! IN v-grid area fn
     &,cos_theta_latitude(row_length,rows)    ! IN T-grid area fn

      LOGICAL land(row_length,rows)           ! IN land mask
      LOGICAL  sea(row_length,rows)           ! IN sea mask
      INTEGER lenout               ! IN max size of output field
      REAL fieldout(lenout)        ! OUT output field
      REAL amdi                    ! IN missing data indicator
      INTEGER len_realhd           ! IN size of real header
      INTEGER len_inthd            ! IN size of integer header
      REAL real_hd(len_realhd)     ! IN real header
      INTEGER int_hd(len_inthd)    ! IN integer header
      INTEGER control_size         ! IN size of control array
      INTEGER control(control_size)! IN control array (mostly not used)
      INTEGER series_entry_size    ! IN no of entries in each record.
      INTEGER no_records           ! IN no of records to process.
      INTEGER extraw               ! OUT no of words required by extra d
      INTEGER n_rows_out,n_cols_out! OUT data-set size and extent
      LOGICAL start_ts             ! IN true if first time-series timest
      INTEGER stash_series(series_entry_size,no_records) ! IN
! IN control data for calls to spatial
      INTEGER icode                       ! OUT error code
      CHARACTER(LEN=80) errmssg              ! OUT error message
!*----------------------------------------------------------------------
! Parameters
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

!*L
! Subroutines called
!
      EXTERNAL SPATIAL
      EXTERNAL STASH_COMP_GRID,EXTRA_MAKE_VECTOR,EXTRA_TS_INFO
!*
! Local variables
!
      INTEGER fake_record(control_size) ! fake_record for SPATIAL call
      INTEGER pp_ptr ! ptr to pp_field for where output from spatial goe
      INTEGER i,j                       ! loop variable
      INTEGER data_size                 ! size of data.
      INTEGER index_size
      INTEGER this_index_lev   ! index of input level
      INTEGER kl ! loop count for levels
      INTEGER stash_list_start ! the start address in index_levs for lev
      INTEGER stash_list_end ! the end address in index_levs for levels
      INTEGER what_process ! what kind of processing is requested
      INTEGER what_mask ! what mask is wanted.
      INTEGER extra_start ! start address in fieldout for extra data
      REAL bdx,bzx,bdy,bzy ! grid descriptors

      INTEGER                                                           &
     &  proc_start_x                                                    &
                        ! processor co-ords of processor at origin
     &, proc_start_y                                                    &
                        ! (SW corner) of subarea
     &, start_pe                                                        &
                        ! processor id of this processor
     &, dummy1,dummy2                                                   &
                        ! unused return variables from
!                         GLOBAL_TO_LOCAL_RC
     &, info            ! GCOM return variable

      REAL                                                              &
     &  global_mean ! global mean returned by call to SPATIAL

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      COMMON /MPP_STATIC_VAR/ global_mean  !must be memory aligned
!                                           to allow fast shmem


!L----------------------------------------------------------------------
!L 1. Error checking
!L
!  Check we are in fact doing a time series. Error if not.
      IF (lhook) CALL dr_hook('MULTI_SPATIAL',zhook_in,zhook_handle)
      IF ((control(st_proc_no_code) /= st_time_series_code).and.        &
     &   (control(st_proc_no_code) /= st_time_series_mean)) THEN
        icode=st_unknown
        write(errmssg,99)control(st_proc_no_code)
99      format(3x,'MULTI_SP : unexpected call to extract timeseries',i5)
        goto 999            ! jump to return
      ENDIF
!L----------------------------------------------------------------------
!L 2. Workout what kind of processing we are doing and what mask is used
!L
       what_mask=mod(control(st_gridpoint_code),block_size)
       what_process=control(st_gridpoint_code)-what_mask
       extraw=0
!L note that the first word in fieldout is assumed to be
!L the word where the first spatial domain mean for this timeseries
!L  will be stored
!L
!L Next we compute the grid discriptors -- as used by extra data
       IF (start_ts) THEN ! compute grid descrip
         extraw=6*(no_records+1)
         extra_start=control(st_output_length)-extraw+1
! DEPENDS ON: stash_comp_grid
         CALL STASH_COMP_GRID(bzx,bzy,bdx,bdy,0,st_grid,                &
     &     1,1,real_hd,len_realhd,int_hd,len_inthd,extract_base+1,      &
     &     icode,errmssg)
       ENDIF
!L 3. Set up pseudo STASH record to be passed to SPATIAL on each call
!L    to extract a sample from the input field.
!L
      fake_record(st_input_bottom)=control(st_input_bottom)
      fake_record(st_input_top)=control(st_input_top)
      fake_record(st_weight_code)=control(st_weight_code)
! doing a field mean on this sub-domain with mask specified by control.
      pp_ptr=1
!L----------------------------------------------------------------------
!L 4. Loop over samples and extract global mean within subdomain for
!L    each, appending to output field
!L
      DO i=1,no_records ! loop over sub domains
        data_size=stash_series(series_size,i)
!L 4.1 Do preliminary verifications on stash_series
!L 4.1.1 Gridtype code
        IF (stash_series(series_grid_type,i) /= series_grid_code) THEN
!
! NB: Latitude/longitude range conversion to gridpoint range needs to
!     be added
!
          icode=st_not_supported
          errmssg='MULTI_SP : only support grid point processing'
          goto 999
        ENDIF
!L ---------------------------------------------------------------------
!L 5. Set up the fake record domain info depending on what kind of
!L    "primary" processing is requested.
!L    As far as stash is concerned everything looks like a global mean
!L    here and it is just a question of setting up the fake record
!L    correctly.
!L
        fake_record(st_gridpoint_code)=                                 &
     &    (stash_series(series_proc_code,i)/block_size)*block_size      &
     &      +what_mask
        IF (what_process == extract_base) THEN ! an extract
          fake_record(st_north_code)=stash_series(series_north,i)
          fake_record(st_south_code)=stash_series(series_south,i)
          fake_record(st_west_code)= stash_series(series_west,i)
          fake_record(st_east_code)= stash_series(series_east,i)
          stash_list_start=stash_series(series_list_start,i)
          stash_list_end=stash_series(series_list_end,i)
          fake_record(st_input_bottom)=stash_list_start
          fake_record(st_input_top)=stash_list_end
        ELSEIF (what_process == zonal_mean_base) THEN ! a zonal_mean
          fake_record(st_north_code)=stash_series(series_north,i)
          fake_record(st_south_code)=stash_series(series_south,i)
          fake_record(st_west_code)= control(st_west_code)
          fake_record(st_east_code)= control(st_east_code)
          stash_list_start=stash_series(series_list_start,i)
          stash_list_end=stash_series(series_list_end,i)
          fake_record(st_input_bottom)=stash_list_start
          fake_record(st_input_top)=stash_list_end
        ELSEIF (what_process == merid_mean_base) THEN ! a merid_mean
          fake_record(st_north_code)= control(st_north_code)
          fake_record(st_south_code)= control(st_south_code)
          fake_record(st_east_code)=stash_series(series_east,i)
          fake_record(st_west_code)=stash_series(series_west,i)
          stash_list_start=stash_series(series_list_start,i)
          stash_list_end=stash_series(series_list_end,i)
          fake_record(st_input_bottom)=stash_list_start
          fake_record(st_input_top)=stash_list_end
        ELSEIF (what_process == vert_mean_base) THEN ! a vert_mean
          fake_record(st_north_code)=stash_series(series_north,i)
          fake_record(st_south_code)=stash_series(series_south,i)
          fake_record(st_east_code)=stash_series(series_east,i)
          fake_record(st_west_code)=stash_series(series_west,i)
          stash_list_start=1
          stash_list_end=num_stash_levels
          fake_record(st_input_bottom)=stash_list_start
          fake_record(st_input_top)=stash_list_end
        ELSEIF (what_process == field_mean_base) THEN ! a field_mean
          fake_record(st_north_code)=control(st_north_code)
          fake_record(st_south_code)=control(st_south_code)
          fake_record(st_east_code)=control(st_east_code)
          fake_record(st_west_code)=control(st_west_code)
          stash_list_start=1
          stash_list_end=num_stash_levels
          fake_record(st_input_bottom)=stash_list_start
          fake_record(st_input_top)=stash_list_end
        ELSE ! error code...
          icode=unknown_processing
          write(errmssg,111) 'unknown processing option',what_process
          goto 999 ! jump to error return
        ENDIF

!L Check record (south > north and west < east)
        IF (fake_record(st_north_code) <                                &
     &    fake_record(st_south_code))then
          write(errmssg,101)fake_record(st_north_code),                 &
     &       fake_record(st_south_code),i
          icode=st_bad_array_param
          goto 999 ! error exit
        ENDIF
        IF (fake_record(st_west_code) >                                 &
     &    fake_record(st_east_code))then
          write(errmssg,102)fake_record(st_west_code),                  &
     &       fake_record(st_east_code),i
          icode=st_bad_array_param
          goto 999 ! error exit
        ENDIF

! Determine which pe holds the first point of the subdomain, ie the
! bottom left hand corner [=SW], that SPATIAL will process.
! This is the pe used to gather all points of the subdomain and
! calculate their global mean in STGLOM before sending these to pe0
! for storage and output.

! DEPENDS ON: global_to_local_rc
        CALL GLOBAL_TO_LOCAL_RC(gr,halo_type,                           &
     &    fake_record(st_west_code),fake_record(st_south_code),         &
     &    proc_start_x, proc_start_y,                                   &
     &    dummy1,dummy2)

        start_pe=proc_start_x+proc_start_y*nproc_x

!
! NB: At present timeseries samples are global (ie. 3D) means, so
!     there is no levels loop outside the call to SPATIAL here -
!     this may be extended at some point to allow multi-level
!     timeseries sampling inside a levels loop
!
!     n_cols_out and n_rows_out are recalculated within SPATIAL but are
!     now appropriate for an individual timeseries sample, not the whole
!     field.  They are reset for the whole field after subdomain loop.
!
        lcyclic=.false.
        this_index_lev = 1     ! no (none) levels loop
! DEPENDS ON: spatial
       CALL SPATIAL(fieldin,vx,vy,vz,gr,st_grid,                        &
     &       fld_type,halo_type,                                        &
     &       halo_x,halo_y,                                             &
     &       lcyclic,lmasswt,                                           &
     &       n_cols_out,n_rows_out,this_index_lev,                      &
     &       level_list(stash_list_start),                              &
     &       index_lev(stash_list_start),                               &
     &       (stash_list_end+1-stash_list_start),                       &
     &       no_of_levels_masswt,                                       &
     &       p,pstar,                                                   &
     &       cos_v_latitude,cos_theta_latitude,land,sea,                &
     &       row_length,rows,n_rows,blsize(2,fld_type),model_levels,    &
     &       global_mean,1,                                             &
     &       fake_record,control_size,amdi,                             &
     &       icode,errmssg)
        IF (icode /= 0) goto 999 ! got some error so jump to return

! Must move the global_mean data to PE 0 which stores all timeseries
! data
! (NB. This assumes that the output from SPATIAL is just a
!      single number)

        IF (mype  ==  start_pe) THEN
          info=GC_NONE
          CALL GC_RSEND(100,1,0,info,fieldout(pp_ptr),global_mean)
        ENDIF

        IF (mype  ==  0) THEN
          info=GC_NONE
          CALL GC_RRECV(100,1,start_pe,info,                            &
     &                  fieldout(pp_ptr),global_mean)
        ELSE
          fieldout(pp_ptr)=0.0
        ENDIF


        pp_ptr=pp_ptr+(n_cols_out*n_rows_out)  ! increment the pp_ptr
!
! NB: n_cols_out and n_rows_out should both be 1 as timeseries samples
!       are currently designed to be scalar quantities only.
!
!L check on n_cols_out and n_rows_out
        IF (n_cols_out /= 1) THEN
          errmssg='MULTI_SP : n_cols_out <> 1'
          icode=st_not_supported
          goto 999
        ENDIF
        IF (n_rows_out /= 1) THEN
          errmssg='MULTI_SP : n_rows_out <> 1'
          icode=st_not_supported
          goto 999
        ENDIF

        IF (mype  ==  0) THEN

        IF (start_ts) THEN ! put the descriptive info for this record
! DEPENDS ON: extra_make_vector
          CALL EXTRA_MAKE_VECTOR(fake_record,control_size,              &
     &      i,no_records,fieldout(extra_start),extraw,bzx,bzy,bdx,bdy)
        ENDIF

        ENDIF

      ENDDO   ! end the loop over sub-domains
!
      horiz_size=pp_ptr-1
      num_vert_levels=1
!L --------------------------------------------------------------------
!L 7. If this is the first time in a time-series then
!L     put the codes describing the extra data into the extra data fld
!L     In addition set pphoriz out to the total length
!L       as well as setting the input vetor to missing
!L       where no values are set
!L----------------------------------------------------------------------

       n_cols_out=no_records
       n_rows_out=control(st_period_code)/control(st_freq_code)
       horiz_size=n_cols_out
       IF (start_ts) THEN  ! on start timestep we have entire vector
         horiz_size=n_cols_out*n_rows_out+extraw

        IF (mype  ==  0) THEN

! DEPENDS ON: extra_ts_info
         CALL EXTRA_TS_INFO(fieldout(extra_start),extraw,no_records)
         do i=no_records+1,extra_start-1
           fieldout(i)=amdi
         enddo

        ELSE
          DO i=no_records+1,lenout
            fieldout(i)=0.0
          ENDDO
        ENDIF

       ENDIF
!
999   CONTINUE ! jump here for error exit
!
111   FORMAT('MULTI_SP :  >>>FATAL ERROR <<',a40,i5,i5)
101   FORMAT('MULTI_SP : north < south',2i5,' in record ',i5)
102   FORMAT('MULTI_SP : west > east',2i5,'in record ',i5 )
      IF (lhook) CALL dr_hook('MULTI_SPATIAL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE MULTI_SPATIAL
