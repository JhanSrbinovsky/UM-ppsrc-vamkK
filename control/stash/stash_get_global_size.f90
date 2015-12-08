! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates the "global" size of STASHed data.
!
! Subroutine Interface:
      SUBROUTINE STASH_GET_GLOBAL_SIZE(                                 &
     &  GLOBAL_NORTH_IN , GLOBAL_EAST_IN ,                              &
     &  GLOBAL_SOUTH_IN , GLOBAL_WEST_IN ,                              &
     &  LEVELS_IN ,                                                     &
     &  GRIDPOINT_CODE , PROCESSING_CODE ,                              &
     &  GLOBAL_FIELD_SIZE ,                                             &
     &  ICODE , CMESSAGE)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      IMPLICIT NONE

! Description:
! Calculates the global (ie. total size on disk) size of a
! STASH request.
!
! Method:
! Using the PROCESSING_CODE to indicate the type of STASH request,
! the GRIDPOINT_CODE to indicate the grid type of the data,
! and the subdomain limits, the total size of the data is calculated.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: stash
!
! Subroutine arguments:

      INTEGER                                                           &
     &  GLOBAL_NORTH_IN                                                 &
                          ! IN: specification of subdomain boundaries
     &, GLOBAL_EAST_IN                                                  &
                          ! IN: ""
     &, GLOBAL_SOUTH_IN                                                 &
                          ! IN: ""
     &, GLOBAL_WEST_IN                                                  &
                          ! IN: ""
     &, LEVELS_IN                                                       &
                          ! IN: number of levels
     &, GRIDPOINT_CODE                                                  &
                          ! IN: indicates the output grid type
     &, PROCESSING_CODE   ! IN: indicates the type of STASH processing

      INTEGER                                                           &
     &  GLOBAL_FIELD_SIZE ! OUT: size of STASH data on disk

      INTEGER                                                           &
     &  ICODE             ! OUT: Return code (0=OK)

      CHARACTER(LEN=80)                                                      &
     &  CMESSAGE          ! OUT: Error message

! Parameters and common blocks
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


! Local variables

      INTEGER                                                           &
! copies of input arguments, which get modified according the
! type of output grid
     &  global_north,global_east,global_south,global_west               &
     &, levels

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ------------------------------------------------------------------

      IF (lhook) CALL dr_hook('STASH_GET_GLOBAL_SIZE',zhook_in,zhook_handle)
      global_north = GLOBAL_NORTH_IN
      global_east  = GLOBAL_EAST_IN
      global_south = GLOBAL_SOUTH_IN
      global_west  = GLOBAL_WEST_IN
      levels       = LEVELS_IN

! Fix wrap-arounds s.t. east > west

      IF (global_west  >   global_east)                                 &
     &  global_east=global_east+glsize(1,fld_type_p)

! Full field or subdomain output:

      IF ((PROCESSING_CODE  ==  st_replace_code) .OR.                   &
     &    (PROCESSING_CODE  ==  st_accum_code) .OR.                     &
     &    (PROCESSING_CODE  ==  st_time_mean_code) .OR.                 &
     &    (PROCESSING_CODE  ==  st_max_code) .OR.                       &
     &    (PROCESSING_CODE  ==  st_min_code)) THEN

        IF ((GRIDPOINT_CODE  >=  vert_mean_base) .AND.                  &
                                                         ! vertical
     &      (GRIDPOINT_CODE  <   zonal_mean_base)) THEN ! mean
          levels=1

        ELSEIF ((GRIDPOINT_CODE  >=  zonal_mean_base) .AND.             &
                                                             ! zonal
     &          (GRIDPOINT_CODE  <   merid_mean_base)) THEN  ! mean
          global_east=global_west

        ELSEIF ((GRIDPOINT_CODE  >=  merid_mean_base) .AND.             &
                                                            ! merid.
     &          (GRIDPOINT_CODE  <   field_mean_base)) THEN ! mean
          global_south=global_north

        ELSEIF ((GRIDPOINT_CODE  >=  field_mean_base)  .AND.            &
                                                             ! field
     &          (GRIDPOINT_CODE  <   global_mean_base)) THEN ! fmean
          global_east=global_west
          global_south=global_north

        ELSEIF (GRIDPOINT_CODE  >=  global_mean_base) THEN
          global_east=global_west
          global_south=global_north
          levels=1

        ELSEIF (GRIDPOINT_CODE  >   global_mean_top) THEN
          ICODE=1
          WRITE(6,*) 'Grid type ',GRIDPOINT_CODE,                       &
     &               ' not yet supported by MPP STASH'
          CMESSAGE='Unsupported grid type'
          GOTO 9999
        ENDIF

        GLOBAL_FIELD_SIZE=(global_east-global_west+1)*                  &
     &                    (global_north-global_south+1)*                &
     &                    levels

      ELSE
        ICODE=2
        WRITE(6,*) 'Processing code ',PROCESSING_CODE,                  &
     &             ' not yet supported by MPP STASH'
        CMESSAGE='Unsupported processing code'
        GOTO 9999
      ENDIF

 9999 CONTINUE

      IF (lhook) CALL dr_hook('STASH_GET_GLOBAL_SIZE',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE STASH_GET_GLOBAL_SIZE

