! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SET_PSEUDO_LIST -----------------------------------------
! 
! Purpose : To set up a list of pseudo levels at which a diagnostic
!           is required, using information in the STASH list.
! 
! Service routine
! 
! Copy of Subroutine SET_LEVELS_LIST (Deck SETLST1) taken and
! adapted for pseudo levels.
! 
! Programming Standard : Unified Model Documentation paper number 3
!                      : Version no 4, dated 5/2/92
! 
! Logical components covered : D3
! 
! System task : P0
! 
! Documentation: U.M. Documentation paper number P0,C4
!
! -----------------------------------------------------------------
!
!*L Arguments
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

      SUBROUTINE SET_PSEUDO_LIST                                        &
     &      (N_LEVELS,LEN_STLIST,STLIST,PSEUDO_LIST,                    &
     &      STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,ICODE,CMESSAGE)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER                                                           &
     &   N_LEVELS                                                       &
                             ! IN Number of possible pseudo levels
     &,  LEN_STLIST                                                     &
                             ! IN Dimension of STLIST
     &,  STLIST(LEN_STLIST)                                             &
                             ! IN STASH list
     &,  NUM_STASH_PSEUDO                                               &
                             ! IN Dimension for STASH_PSEUDO_LEVELS
     &,  STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,*)                      &
                                                   ! IN Pseudo levels
     &,  ICODE               ! OUT Return code

      LOGICAL                                                           &
     &   PSEUDO_LIST(N_LEVELS) ! OUT List of pseudo levels required.

      CHARACTER(LEN=80) CMESSAGE ! Error message

!* ---------------------------------------------------------------------

!*L   Workspace Usage :- None

!*L   External Subroutines called :- None

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
!* ---------------------------------------------------------------------
!  Local variables

      INTEGER                                                           &
     &      JLEV                                                        &
                       ! Loop counter over levels
     &,     LEVEL_NO                                                    &
                       ! Level no in pseudo list
     &,     LIST_NO    ! Pseudo level list number

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!* ---------------------------------------------------------------------

!L Initialise pseudo levels list to false

      IF (lhook) CALL dr_hook('SET_PSEUDO_LIST',zhook_in,zhook_handle)
      DO JLEV=1,N_LEVELS
        PSEUDO_LIST(JLEV)= .FALSE.
      END DO

!L Get pseudo list number

      LIST_NO = STLIST(ST_PSEUDO_IN)

!L Check that Pseudo list number is valid (should be GE 0)

      IF (LIST_NO <  0) THEN

!       Illegal control data

        ICODE=1
        CMESSAGE = 'SET_PSEUDO_LIST: Illegal control data'
        WRITE(6,*) 'SET_PSEUDO_LIST: Illegal control data'
        WRITE(6,*) 'ST_PSEUDO_IN         = ',ST_PSEUDO_IN
        WRITE(6,*) 'STLIST(ST_PSEUDO_IN) = ',STLIST(ST_PSEUDO_IN)
        WRITE(6,*) 'Section and item numbers ',STLIST(2),STLIST(1)
        GO TO 999  !  Return

      ENDIF

!  Set logical array list to identify pseudo levels required.

      IF (LIST_NO >  0) THEN

        DO JLEV=2,STASH_PSEUDO_LEVELS(1,LIST_NO)+1
          LEVEL_NO = STASH_PSEUDO_LEVELS(JLEV,LIST_NO)
          IF (LEVEL_NO >= 1 .AND. LEVEL_NO <= N_LEVELS) THEN

!           Level is within range
            PSEUDO_LIST(LEVEL_NO) =.TRUE.

          ELSE

!           Level is out of range
            ICODE=2
            CMESSAGE=  ' SET_PSEUDO_LIST : level out of range'
            WRITE(6,*) ' SET_PSEUDO_LIST : level out of range'
            WRITE(6,*) ' pseudo list no = ',LIST_NO
            WRITE(6,*) ' level = ',LEVEL_NO
            WRITE(6,*) ' Section, Item = ',STLIST(2),STLIST(1)
            GO TO 999   !  Return

          END IF
        END DO

      END IF

999 CONTINUE

      IF (lhook) CALL dr_hook('SET_PSEUDO_LIST',zhook_out,zhook_handle)
      RETURN

      END SUBROUTINE SET_PSEUDO_LIST
