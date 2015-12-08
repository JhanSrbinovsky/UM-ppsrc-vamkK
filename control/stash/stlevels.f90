! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine STLEVELS -----------------------------------------------
!
! Purpose: Generate a level index from STASHrecord and level_lists
!          and number of levels tailored to a particular diagnostic.
!          Also set levels and pseudo-levels information for encoding
!          PPheader details.  (This subroutine based on a merger
!          between GEN_INDEX and PP_COMPUTE_LEVEL).
!                 New subroutine STLEVELS is based on GEN_INDEX and
!                 PP_COMPUTE_LEVEL with merged functionality.
!          A general note as levels list is an integer
!          real values are multiplied by a 1000.0.
!          When computing the real value of the level for the
!          pp header it is necessary to divide by a 1000.0.
!          Levels that are affected by this are theta, pressure and
!          height.
!
! Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
! Logical components covered : C4?
!
! Project task: C4
!
! External documentation : UMDP no C4
!
! Interface and arguments: ------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: STASH
      SUBROUTINE STLEVELS(stash_control,stash_control_size,             &
     &     stash_levels,num_stash_levels,num_level_lists,               &
     &     stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,       &
     &     max_stash_levs,num_levs_in,num_levs_out,index_size,          &
     &     index_lev,level_list,                                        &
     &     lbvcl,level,pseudo_level,                                    &
     &     icode,cmessage)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE sterr_mod, ONLY: st_upper_less_lower, st_not_supported,       &
                           st_no_data,st_nd, st_bad_array_param,        &
                           st_bad_address, st_unknown,                  &
                           st_bad_wraparound, st_illegal_weight,        &
                           unknown_weight, unknown_mask,                &
                           unknown_processing, nonsense
      USE cppxref_mod, ONLY: ppx_lbvc_height, ppx_lbvc_pressure,        &
                             ppx_lbvc_theta, ppx_lbvc_PV
      IMPLICIT NONE
!
      INTEGER                                                           &
     &       stash_control_size,                                        &
                                 ! IN size of stash control record
     &       stash_control(stash_control_size),                         &
                                               ! IN  stash control
     &       num_stash_levels,                                          &
                                 ! IN max. no of hts for a levels list
     &       num_level_lists,                                           &
                                 ! IN max. no of level lists
     &       stash_levels(num_stash_levels+1,num_level_lists),          &
                                                               ! IN
!                                !    lookup table for level lists
     &       num_stash_pseudo,num_pseudo_lists,                         &
                                               ! IN dims of pseudo_levs
     &       stash_pseudo_levels(num_stash_pseudo+1,num_pseudo_lists),  &
!                                ! IN lookup table for pseudo-lev lists
     &       max_stash_levs,                                            &
                                 ! IN max. no of output levels
     &       num_levs_in,                                               &
                                 ! OUT no of levels in input data
     &       num_levs_out,                                              &
                                 ! OUT no of levels in output data
     &       index_size,                                                &
                                 ! OUT no of levels in levels index
     &       index_lev(max_stash_levs),                                 &
                                        ! OUT index of output level
!                                               relative to input level
     &       level_list(max_stash_levs),                                &
                                         ! OUT value of model level
     &       pseudo_level(max_stash_levs),                              &
                                           ! OUT Value of pseudo levels
     &       lbvcl,                                                     &
                                 ! IN  vertical coordinate PP code
     &       icode               ! OUT error code
      REAL                                                              &
     &       level(max_stash_levs)  ! OUT Value of output levels (real)
      CHARACTER(LEN=*)                                                     &
     &       cmessage            ! OUT error message
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
!
! Local variables
!
      INTEGER                                                           &
     &       index_pseudo_lev(max_stash_levs),                          &
                                               ! Pseudo-level 1D index
     &       num_pseudo_in,num_pseudo_out,                              &
                                               ! Number of pseudo levs
     &       k2,ml,kl,                                                  &
                                       ! loop counts
     &       NI,NO,                                                     &
                                       ! Number In/Out
     &       indx1,                                                     &
                                       ! index count
     &       ilev,                                                      &
                                       ! Integer level/pseudo-level
     &       what_mean,what_proc       ! Meaning and processing code

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
! First compute the index for physical levels
!
      IF (lhook) CALL dr_hook('STLEVELS',zhook_in,zhook_handle)
      IF(STASH_CONTROL(st_input_bottom) <  0) THEN ! Input LEVELS list
        NI=-STASH_CONTROL(st_input_bottom)
        NUM_LEVS_IN=STASH_LEVELS(1,NI)
        IF(STASH_CONTROL(st_output_bottom) <  0) THEN ! LEVELS LIST out
          NO=-STASH_CONTROL(st_output_bottom)
          NUM_LEVS_OUT=STASH_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_LEVS_OUT
            ilev=STASH_LEVELS(ML+1,NO)    !  Level required
            DO KL=1,NUM_LEVS_IN
              IF(STASH_LEVELS(KL+1,NI) == ilev) THEN
                INDX1=INDX1+1
                INDEX_LEV(INDX1)=KL   ! Relative position of Input to Ou
                level_list(indx1)=ilev
                GOTO 400
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output level ',ilev,                   &
     &                          ' not found in input levels list'
            GOTO 999
 400        CONTINUE
            ENDDO
        ELSE           !  Output as a Level range
          NUM_LEVS_OUT=STASH_CONTROL(st_output_top)-                    &
     &                 STASH_CONTROL(st_output_bottom)+1
          ilev=STASH_CONTROL(st_output_bottom) !1st output model level
          DO KL=1,NUM_LEVS_IN
            IF(STASH_LEVELS(KL+1,NI) == ilev) THEN
              INDEX_LEV(1)=KL ! Relative posn of Input to the 1st level
              level_list(1)=ilev
              GOTO 401
            ENDIF
          ENDDO
          ICODE=nonsense
          WRITE(CMESSAGE,101) 'Output bottom model level ',ilev,        &
     &                        ' not found in input levels list'
          GOTO 999
 401      CONTINUE
          DO KL=2,NUM_LEVS_OUT
            INDEX_LEV(KL)=INDEX_LEV(KL-1)+1
            level_list(kl)=level_list(kl-1)+1
          ENDDO
        ENDIF
      ELSEIF(STASH_CONTROL(st_input_bottom) == 100) THEN !Special level
          NUM_LEVS_IN=1
          NUM_LEVS_OUT=1
          INDEX_LEV(1)=1
          level_list(1)=1 ! could be worth setting to some nonsense no.
      ELSE     !  Input is Model level range
        NUM_LEVS_IN=STASH_CONTROL(st_input_top)-                        &
     &              STASH_CONTROL(st_input_bottom)+1
        IF(STASH_CONTROL(st_output_bottom) <  0) THEN ! LEVELS LIST out
          NO=-STASH_CONTROL(st_output_bottom)
          NUM_LEVS_OUT=STASH_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_LEVS_OUT
            ilev=STASH_LEVELS(ML+1,NO)    ! Output level reqd
            DO KL=1,NUM_LEVS_IN
              IF((STASH_CONTROL(st_input_bottom)+KL-1) == ilev) THEN
                INDX1=INDX1+1
                INDEX_LEV(INDX1)=KL   ! Relative posn of output to inpt
                level_list(INDX1)=ilev
                GOTO 402
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output model level ',ilev,             &
     &                          ' not in input model level range'
            GOTO 999
 402        CONTINUE
          ENDDO
        ELSE     !   Output as model level range
! Do some consistency checks here to ensure valid processing request
! output bottom should be greater or equal to input bottom
          IF (stash_control(st_output_bottom) <                         &
     &       stash_control(st_input_bottom)) THEN
            icode=nonsense
            write(cmessage,103)'bad level spec, bot input>output',      &
     &       stash_control(st_input_bottom),                            &
     &       stash_control(st_output_bottom)
            goto 999 ! jump to error
          ELSEIF (stash_control(st_output_top) >                        &
     &         stash_control(st_input_top)) THEN
            icode=nonsense
            write(cmessage,103)'bad level spec, top input<output',      &
     &        stash_control(st_input_top),                              &
     &        stash_control(st_output_top)
              goto 999 ! jump to error
          ENDIF
          NUM_LEVS_OUT=STASH_CONTROL(st_output_top)-                    &
     &                 STASH_CONTROL(st_output_bottom)+1
          INDEX_LEV(1)=STASH_CONTROL(st_output_bottom)-                 &
     &                 STASH_CONTROL(st_input_bottom)+1
          level_list(1)=stash_control(st_output_bottom)
          DO kl=2,NUM_LEVS_OUT
            INDEX_LEV(kl)=INDEX_LEV(kl-1)+1
            level_list(kl)=level_list(kl-1)+1
          ENDDO
        ENDIF
      ENDIF
      index_size=num_levs_out
      IF (num_levs_out >  num_levs_in) THEN   ! things very badly wrong
        icode=nonsense
        write(cmessage,103)'asking for num_levs_out>num_levs_in',       &
     &   num_levs_out,num_levs_in
        goto 999 ! jump to return
      ENDIF
!
! Next, compute actual (physical) levels for encoding PPheaders
!
      IF (STASH_CONTROL(st_output_bottom) <  0) THEN ! Levels List ?
        NO=-STASH_CONTROL(st_output_bottom)     ! Index of Levels list

          ! Remove scaling (by factor 1000) of vertical level coord
          ! for certain types of STASH output [originally needed to
          ! store in an intermediary integer array]
          IF( LBVCL  ==  ppx_lbvc_height   .OR.                         &
                                                !  height levels
     &        LBVCL  ==  ppx_lbvc_pressure .OR.                         &
                                                ! pressure levels
     &        LBVCL  ==  ppx_lbvc_theta    .OR.                         &
                                                ! theta levels
     &        LBVCL  ==  ppx_lbvc_PV ) THEN     ! potential vorticity


          DO ML=1,NUM_LEVS_OUT
            LEVEL(ML)=REAL(STASH_LEVELS(ML+1,NO))*0.001+1.0E-10
          ENDDO
        ELSE
          DO ML=1,NUM_LEVS_OUT
            LEVEL(ML)=REAL(STASH_LEVELS(ML+1,NO))
          ENDDO
        ENDIF
      ELSEIF (STASH_CONTROL(st_output_bottom) == st_special_code) THEN
       ! Special level.
       ! The LEVEL array is not used by the model except to construct pp
       ! header items at output. The value of -1.0 is set as a flag for
       ! special levels so that routine PP_HEAD will insert the lbvc
       ! item in STASHmaster record.
        DO ML=1,NUM_LEVS_OUT
          LEVEL(ML)=-1.0
        ENDDO
      ELSE
        DO ML=1,NUM_LEVS_OUT
          LEVEL(ML)=REAL(STASH_CONTROL(st_output_bottom)+ML-1)
        ENDDO
      ENDIF
!
!
! Now reset the number of output levels to 1 if vertical compression is
! to be done in SPATIAL.  NB: index_lev and level_list need to be filled
! with values corresponding to the full range of levels processed.
!
      what_proc=STASH_CONTROL(st_proc_no_code)
      what_mean=(STASH_CONTROL(st_gridpoint_code)/block_size)*block_size
      IF(what_mean == vert_mean_base .OR. what_mean == global_mean_base &
     &   .OR. what_proc == st_time_series_code                          &
     &   .OR. what_proc == st_time_series_mean                          &
     &   .OR. what_proc == st_append_traj_code) num_levs_out=1
!
! Next compute the index for pseudo levels, if there are any
!
      IF(STASH_CONTROL(st_pseudo_in) >  0) THEN ! Input PSEUDO_LEVELS
        NI=STASH_CONTROL(st_pseudo_in)
        num_pseudo_in=STASH_PSEUDO_LEVELS(1,NI)
        IF(STASH_CONTROL(st_pseudo_out) >  0) THEN ! Output PSEUDO_LEVS
          NO=STASH_CONTROL(st_pseudo_out)
          num_pseudo_out=STASH_PSEUDO_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_PSEUDO_OUT
            ilev=STASH_PSEUDO_LEVELS(ML+1,NO)   !  Level required
            DO KL=1,NUM_PSEUDO_IN
              IF(STASH_PSEUDO_LEVELS(KL+1,NI) == ilev) THEN
                INDX1=INDX1+1
                INDEX_PSEUDO_LEV(INDX1)=KL
                pseudo_level(indx1)=ilev
                GOTO 500
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output pseudo level ',ilev,            &
     &                          ' not found in input levels list'
            GOTO 999
 500        CONTINUE
          ENDDO
        ELSE  ! Illegal combination
          ICODE=nonsense
          WRITE(CMESSAGE,101) 'Input pseudo level list ',NI,            &
     &         ' has illegal output pseudo levels list'
          GOTO 999
        ENDIF
      ELSE  ! Only levels lists are supported for pseudo levels
        num_pseudo_out=0
      ENDIF
!
! Next expand the separate indexes and physical levels arrays into
! combined arrays if necessary, taking care not to overwrite earlier
! parts of the arrays.  If no pseudo-levels, set pseudo-level to 0.
!
      IF (num_pseudo_out >  0) THEN
        DO K2=num_pseudo_out,1,-1
          DO ML=1,num_levs_out
            INDEX_LEV(ML+(K2-1)*num_levs_out)=                          &
     &        (INDEX_PSEUDO_LEV(K2)-1)*num_levs_in+INDEX_LEV(ML)
            level(ML+(K2-1)*num_levs_out)=level(ML)
          ENDDO
          DO ML=num_levs_out,1,-1
            pseudo_level(ML+(K2-1)*num_levs_out)=pseudo_level(K2)
          ENDDO
        ENDDO
        num_levs_out=num_levs_out*num_pseudo_out
      ELSE
        DO ML=1,num_levs_out
          pseudo_level(ML)=0
        ENDDO
      ENDIF
!
999   CONTINUE ! jump here for error return
 101  FORMAT('STLEVELS : ',a,i6,a)
 103  FORMAT('STLEVELS : >> FATAL ERROR <<',a,2i5)
      IF (lhook) CALL dr_hook('STLEVELS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE STLEVELS

