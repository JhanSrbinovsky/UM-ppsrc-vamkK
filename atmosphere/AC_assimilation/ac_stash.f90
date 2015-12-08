! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE AC_STASH -----------------------------------------------
!LL
!LL  Purpose : Stash processing for AC Scheme
!LL
!LL  Global and Limited area
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: AC Assimilation
MODULE ac_stash_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE AC_STASH (STASH_ITEM_NO,LEVEL_NO,N_LEVELS,             &
     &                     GROUP_NO,N_GROUPS,TIMESTEP_NO,               &
     &                     STINDEX,STLIST,LEN_STLIST,SI,SF,             &
     &                     STASHWORK,STASH_LEVELS,NUM_STASH_LEVELS,     &
     &                     STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,        &
     &                     FIELD,LEN_FLD,LABEL,                         &
     &                     ICODE,CMESSAGE)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER                                                           &
     &   STASH_ITEM_NO                                                  &
                              ! Stash Item Number
     &,  LEVEL_NO                                                       &
                              ! Level number
     &,  N_LEVELS                                                       &
                              ! No of levels
     &,  GROUP_NO                                                       &
                              ! Group number
     &,  N_GROUPS                                                       &
                              ! No of groups
     &,  TIMESTEP_NO                                                    &
                              ! Timestep number
     &,  LEN_FLD                                                        &
                              ! Length of field to be stashed
     &,  LEN_STLIST                                                     &
                              ! Dimension of STLIST
     &,  STINDEX(2,*)                                                   &
                              ! Start and no of items in STLIST
     &,  STLIST(LEN_STLIST,*)                                           &
                              ! Stash List of items to be output
     &,  SI(*)                                                          &
                              ! Address of Item in STASHWORK
     &,  NUM_STASH_LEVELS                                               &
                              ! Number of levels lists
     &,  NUM_STASH_PSEUDO                                               &
                              ! Number of pseudo lists
     &,  STASH_LEVELS(NUM_STASH_LEVELS+1,*)                             &
                                                   ! Levels lists
     &,  STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,*) ! Pseudo lists

      REAL                                                              &
     &   STASHWORK(*)                                                   &
                              ! Work array for stashed data
     &,  FIELD(LEN_FLD)       ! Field to be stashed

      LOGICAL SF(*)           ! Stash Flags

      CHARACTER(LEN=*) LABEL     !  Label to indicate field being stashed

      INTEGER ICODE           !  Return code
      CHARACTER(LEN=256) CMESSAGE  !  Error message

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

!*L   Dynamic allocated arrays

      LOGICAL LEVELS_LIST(N_LEVELS)  ! Expanded levels list
      LOGICAL PSEUDO_LIST(N_GROUPS)  ! Expanded pseudo list

!*L   External Subroutines called

      EXTERNAL SET_LEVELS_LIST, SET_PSEUDO_LIST

!     Local variables

      LOGICAL                                                           &
     &   L_SINGLE_LEV                                                   &
     &,  L_LEVELS_LIST                                                  &
     &,  L_PSEUDO_LIST                                                  &
     &,  LSTASH

      INTEGER                                                           &
     &   J,JLEV,JGRP                                                    &
                          !  Loop counters over levels/groups
     &,  J0                                                             &
                          !  Pointer in STASHWORK
     &,  IPOS                                                           &
                          !  Position in STLIST for this STASH_ITEM_NO
     &,  N_LEVELS_LIST                                                  &
                          !  No of levels in levels list
     &,  FLD_NO                                                         &
                          !  Field number in STASHWORK
     &,  GRP_NO           !  Group number in STASHWORK

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     Check that LEVEL_NO le N_LEVELS

      IF (lhook) CALL dr_hook('AC_STASH',zhook_in,zhook_handle)
      IF (LEVEL_NO >  N_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = ' AC_STASH : LEVEL_NO gt N_LEVELS ?'
        WRITE (6,*) 'AC_STASH : LEVEL_NO must be LE to N_LEVELS'
        WRITE (6,*) 'Stash Item No = ',STASH_ITEM_NO
        WRITE (6,*) 'LEVEL_NO = ',LEVEL_NO,' N_LEVELS = ',N_LEVELS
        GO TO 999   !  Return
      ENDIF

!     Check that GROUP_NO le N_GROUPS
      IF (GROUP_NO >  N_GROUPS) THEN
        ICODE    = 1
        CMESSAGE = ' AC_STASH : GROUP_NO gt N_GROUPS ?'
        WRITE (6,*) 'AC_STASH : GROUP_NO must be LE to N_GROUPS'
        WRITE (6,*) 'Stash Item No = ',STASH_ITEM_NO
        WRITE (6,*) 'GROUP_NO = ',GROUP_NO,' N_GROUPS = ',N_GROUPS
        GO TO 999   !  Return
      ENDIF

      DO JLEV = 1,N_LEVELS
        LEVELS_LIST(JLEV) = .FALSE.
      ENDDO
      DO JGRP = 1,N_GROUPS
        PSEUDO_LIST(JGRP) = .FALSE.
      ENDDO

!     Get position in STLIST
      IPOS = STINDEX(1,STASH_ITEM_NO)

!     Determine if levels list used (Entry 10 in STLIST = negative)
      L_LEVELS_LIST = STLIST(ST_INPUT_BOTTOM,IPOS) <  0

!     Determine if single level (Entry 10 in STLIST = 100)
      L_SINGLE_LEV  = STLIST(ST_INPUT_BOTTOM,IPOS) == 100

!     Determine if pseudo list used (Entry 26 in STLIST = positive)
      L_PSEUDO_LIST = STLIST(ST_PSEUDO_IN,IPOS) >  0

      N_LEVELS_LIST = 1

      IF (L_LEVELS_LIST) THEN

        N_LEVELS_LIST = STASH_LEVELS(1,-STLIST(ST_INPUT_BOTTOM,IPOS))

!----   Get levels required for this field
!----   (Sets up LEVELS_LIST from STASH_LEVELS)
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST (N_LEVELS,LEN_STLIST,STLIST(1,IPOS),       &
     &       LEVELS_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,               &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) GO TO 999   !  Return

      ELSEIF (.NOT. L_SINGLE_LEV) THEN

        ICODE = 1
        CMESSAGE = 'AC_STASH ; No levels list ?'
        GO TO 999   !  Return

      ELSE
      ENDIF

      IF (L_PSEUDO_LIST) THEN

!----   Get pseudo levels list required for this field
!----   (Sets up PSEUDO_LIST from STASH_PSEUDO_LEVELS)
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST (N_GROUPS,LEN_STLIST,STLIST(1,IPOS),       &
     &       PSEUDO_LIST,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,          &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) GO TO 999   !  Return

      ENDIF

!     Determine if this field is to stashed
      LSTASH = .FALSE.
      IF ( L_SINGLE_LEV .OR.                                            &
     &    (L_LEVELS_LIST .AND. LEVELS_LIST(LEVEL_NO)) ) THEN
        IF (L_PSEUDO_LIST) THEN
          IF (PSEUDO_LIST(GROUP_NO)) LSTASH = .TRUE.
        ELSE
          LSTASH = .TRUE.
        ENDIF
      ENDIF

      IF (LSTASH) THEN

!----   Determine position in STASHWORK for this field
        FLD_NO = 0
        IF (L_SINGLE_LEV) THEN
          FLD_NO = 1
!         Could have pseudo levels - look at later
        ELSE
          DO JLEV=1,LEVEL_NO
            IF (LEVELS_LIST(JLEV)) THEN
              FLD_NO = FLD_NO+1
            ENDIF
          ENDDO
        ENDIF

        GRP_NO = 0
        IF (L_PSEUDO_LIST) THEN
          DO JLEV=1,GROUP_NO
            IF (PSEUDO_LIST(JLEV)) THEN
              GRP_NO = GRP_NO+1
            ENDIF
          ENDDO
        ELSE
          GRP_NO = 1
        ENDIF

!----   Check FLD_NO
        IF (FLD_NO == 0) THEN
          ICODE    = 1
          CMESSAGE = ' AC_STASH : FLD_NO = 0 ?'
          WRITE (6,*) 'AC_STASH : FLD_NO must be GT than 0'
          WRITE (6,*) 'Stash Item No = ',STASH_ITEM_NO
          WRITE (6,*) 'Level,Group,FLD_NO = ',LEVEL_NO,GROUP_NO,FLD_NO
          GO TO 999   !  Return
        ENDIF

!----   Check GRP_NO
        IF (GRP_NO == 0) THEN
          ICODE    = 1
          CMESSAGE = ' AC_STASH : GRP_NO = 0 ?'
          WRITE (6,*) 'AC_STASH : GRP_NO must be GT than 0'
          WRITE (6,*) 'Stash Item No = ',STASH_ITEM_NO
          WRITE (6,*) 'Level,Group,FLD_NO = ',LEVEL_NO,GROUP_NO,FLD_NO
          GO TO 999   !  Return
        ENDIF

!----   Set up pointer for this field in STASHWORK
        J0 = SI(STASH_ITEM_NO) - 1                                      &
     &     + (GRP_NO-1)*(N_LEVELS_LIST*LEN_FLD)                         &
     &     + (FLD_NO-1)*LEN_FLD

!----   Copy field into work space
        DO J=1,LEN_FLD
          STASHWORK(J0 + J) = FIELD (J)
        ENDDO


      ENDIF

 999  CONTINUE

      IF (lhook) CALL dr_hook('AC_STASH',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE AC_STASH
END MODULE ac_stash_mod
