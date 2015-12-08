! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine DIAGDESC -----------------------------------------------
!
!  Purpose: Prints a formatted diagnostic description using the name
!           of a diagnostic plus it's PPXREF and STASH record.  Gives
!           a hardcopy record of the diagnostics included in a run.
!
!  Programming standard: UM Doc Paper 3
!
!  External documentation:
!  Unified Model Doc Paper C4 - Storage handling and diagnostic
!                               system (STASH)
! --------------------------------------------------------------------
!
! Interface and arguments: ------------------------------------------
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Misc

SUBROUTINE DIAGDESC(seqno,name,stlist,ppxref,                           &
              stash_levels,num_stash_levels,num_level_lists,            &
              stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,    &
              sttabl,nsttims,nsttabl,                                   &
              stash_series,stash_series_rec_len,stash_series_len,       &
              stash_series_index,stash_ser_index_size)
!
USE UM_ParVars, ONLY: mype
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE cppxref_mod
USE Submodel_Mod
USE nlstcall_mod, ONLY : pp_pack_code

USE chsunits_mod, ONLY : nunits

IMPLICIT NONE
!

CHARACTER(LEN=36), INTENT(IN)  ::  name        ! IN  diagnostic name
INTEGER, INTENT(IN)         ::                                          &
                                seqno,                                  &
                                            ! IN  sequence number
                                stlist(*),                              &
                                            ! IN  STASHlist record
                                ppxref(*)   ! IN  PPXREF record

! STASH levels list information
INTEGER, INTENT(IN)         ::  num_stash_levels                        &
                                            ! IN Max levels in a list
                               ,num_level_lists                         &
                                            ! IN Number of lists
                    ,stash_levels(num_stash_levels+1,num_level_lists)
                    
! STASH pseudo-levels list information      
INTEGER, INTENT(IN)         ::  num_stash_pseudo                        &
                                            ! IN Max ps-levs in a list
                               ,num_pseudo_lists                        &
                                            ! IN No of ps-lev lists
            ,stash_pseudo_levels(num_stash_pseudo+1,num_pseudo_lists)
            
! STASH time list information
INTEGER, INTENT(IN)         ::  nsttims                                 &
                                            ! IN Max times in a list
                               ,nsttabl                                 &
                                            ! IN Number of lists
                               ,sttabl(nsttims,nsttabl)
                               
! STASH timeseries information
INTEGER, INTENT(IN)         ::  stash_series_len                        &
                                            ! IN Total no of records
                               ,stash_series_rec_len                    &
                                            ! IN Length of each record
                 ,stash_series(stash_series_rec_len,stash_series_len)   &
                                            ! IN array of records
                               ,stash_ser_index_size                    &
                                            ! IN No of index records
                          ,stash_series_index(2,stash_ser_index_size)
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
!LL  Comdeck: CCONTROL -------------------------------------------------
!LL
!LL  Purpose: COMMON block for top level switches and 2nd level switches
!LL           needed by the top level (C0) and 2nd level routines, but
!LL           not held in the history COMMON block.
!LL
!LLEND ---------------------------------------------------------------

!#include "cntlall.h"
! cntlgen.h was replaced by control/top_level/nlstgen_mod.F90
! #include "cntlgen.h"


!
! Local variables
!
CHARACTER(LEN=80)  :: ch                ! Working character string variable

CHARACTER(LEN=8)   :: submodel,datatype,grid_type
CHARACTER(LEN=9)   :: level_type,freq,weight
CHARACTER(LEN=7)   :: packing,mask
CHARACTER(LEN=15)  :: time_proc,pseudo,levels
CHARACTER(LEN=6)   :: from, to,period
CHARACTER(LEN=10)  :: source
CHARACTER(LEN=17)  :: dest
CHARACTER(LEN=18)  :: spatial
CHARACTER(LEN=27)  :: horiz

INTEGER         :: i1,i2,k           ! Array indices
INTEGER         :: j                 ! Code value
INTEGER         :: time_list,lev_list                                   &
                                     ! pointers to time and levels lists
                  ,plev_list                                            &
                                     ! pointer  to pseudo-level list
                  ,tser_list         ! pointer  to time series record list
INTEGER         :: ntimes            ! no of times in a time list
INTEGER         :: packing_profile   ! packing profile for output PPfield




! Header
CHARACTER(LEN=*), PARAMETER ::                                             &
 stars='   ********************************************************'    &
,starp='   **                                                    **'    &
,list ='   **    LIST OF USER-DEFINED DIAGNOSTICS IN THIS RUN    **'    &
,notes='   ** NOTES:                                             **'    &  
,note1='   **   Time processing details are in timesteps, where  **'    &      
,note2='   **     ... represents "for ever".                     **'    &
,note3='   **   Spatial processing domain is in gridpoints.      **'    &  
,headend=                                                               &
      '=================================================================&
     &=========================================================='     

CHARACTER(LEN=*), PARAMETER ::                                             &
    diag1=' #No Diagnostic Description-------------- Submodel Item Section '//&
    'PPfcode Datatype Gridtype Leveltype MetO8lv MetO8fc Packacc',     &
    diag2=' Time-processing -From- --To-- Frequency Period --Source-- ---De'//&
    'stination---                                               ',     &
    diag3=' Spatial-processing -Levels-domain- -Pseudo-levels- -----Horizon'//&
    'tal-domain----- Weighting Masking                          '



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!
!----------------------------------------------------------------------
! 0. Write header if sequence no indicates first item
!
      
IF (lhook) CALL dr_hook('DIAGDESC',zhook_in,zhook_handle)

IF (seqno == 1 .AND. mype == 0) THEN
  WRITE(200,'(15(a,/),/,a)')                                            &
       stars,stars,starp,list,starp,stars,stars,starp,                  &
       notes,note1,note2,note3,starp,stars,stars,headend
END IF

!----------------------------------------------------------------------
! 1. For each diagnostic processing request in the STASHlist,
!    print the diagnostic name followed by a summary of the processing
!    information on 3 lines.
!
! 1.0 If diagnostic is not required for output, exit routine
!
IF (stlist(st_proc_no_code) /= 0) THEN
!
! 1.1 Line 1.
!
! #No               -- supports up to 999 diagnostic fields, i3.

! Name              

! Submodel
  j=stlist(st_sect_no_code)
  IF      (ppxref(ppx_model_number) == atmos_im) THEN
    submodel=' ATMOS  '
  ELSE
    WRITE(6,'(A)')' Error in DIAGDESC: Unknown model'
    submodel=' UNKNOWN'
  ENDIF

! Item

! Section

! PPfcode

! Datatype
  j=ppxref(ppx_data_type)
  IF (j == 1.OR.j == 4) THEN
    datatype='  REAL  '
  ELSE IF (j == 2.OR.j == 5) THEN
    datatype='INTEGER '
  ELSE IF (j == 3) THEN
    datatype='LOGICAL '
  ELSE
    datatype='UNKNOWN '
  END IF
    
! Gridtype
  j=ppxref(ppx_grid_type)
  IF (j == ppx_atm_nonstd.OR.j == ppx_ocn_nonstd) THEN
    grid_type=' NONSTD '
  ELSE IF ((j >  ppx_atm_nonstd.AND.j <= ppx_atm_tsea) .OR.         &
            j == ppx_atm_compressed.OR.j == ppx_atm_ozone) THEN
    grid_type =' P-GRID '
  ELSE IF (j >= ppx_atm_uall.AND.j <= ppx_atm_usea) THEN
    grid_type=' UV-GRID'
  ELSE IF (j == ppx_atm_cuall.OR.j == ppx_ocn_cuall) THEN
    grid_type=' CU-GRID'
  ELSE IF (j == ppx_atm_cvall.OR.j == ppx_ocn_cvall) THEN
    grid_type=' CV-GRID'
  ELSE IF (j == ppx_atm_tzonal) THEN
    grid_type=' PZ-GRID'
  ELSE IF (j == ppx_atm_uzonal) THEN
    grid_type=' UZ-GRID'
  ELSE IF (j == ppx_atm_tmerid) THEN
    grid_type=' PM-GRID'
  ELSE IF (j == ppx_atm_umerid) THEN
    grid_type=' UM-GRID'
  ELSE IF (j == ppx_atm_rim.OR.j == ppx_ocn_rim) THEN
    grid_type='   RIM  '
  ELSE IF (j == ppx_ocn_tcomp.OR.j == ppx_ocn_tall.OR.              &
           j == ppx_ocn_tfield) THEN
    grid_type=' T-GRID '
  ELSE IF (j == ppx_ocn_tzonal) THEN
    grid_type=' TZ-GRID'
  ELSE IF (j == ppx_ocn_uzonal) THEN
    grid_type=' UZ-GRID'
  ELSE IF (j == ppx_ocn_tmerid) THEN
    grid_type=' TM-GRID'
  ELSE IF (j == ppx_ocn_umerid) THEN
    grid_type=' UM-GRID'
  ELSE IF (j == ppx_ocn_ucomp.OR.j == ppx_ocn_uall.OR.              &
           j == ppx_ocn_ufield) THEN
    grid_type=' UV-GRID'
  ELSE IF (j == ppx_atm_scalar.OR.j == ppx_ocn_scalar) THEN
    grid_type=' SCALAR '
  ELSE IF (j == ppx_wam_all.OR.j == ppx_wam_sea) THEN
    grid_type=' WAVE   '
  ELSE IF (j == ppx_wam_rim) THEN
    grid_type=' RIM    '
  ELSE
    grid_type=' UNKNOWN'
  END IF

! Leveltype
  j=ppxref(ppx_lv_code)
  IF (j == ppx_full_level) THEN
    level_type='FULLLEVEL'
  ELSE IF (j == ppx_half_level) THEN
    level_type='HALFLEVEL'
  ELSE
    level_type='STD-LEVEL'
  END IF

! Meto8LV

! Meto8FC

! PackAcc
  j=stlist(st_output_code)
  IF (j == 1) THEN
    IF (stlist(st_macrotag) >= 1000) THEN
      packing_profile=pp_pack_code(27)
    ELSE
      packing_profile=0
    ENDIF
  ELSE IF(j == 2) THEN
    packing_profile=0
  ELSE IF(j <  0) THEN
    packing_profile=pp_pack_code(-j)
  ELSE
    packing_profile=0
  END IF
  IF (packing_profile == 0) THEN
    packing='       '
  ELSE
    j=ppxref(ppx_packing_acc+packing_profile-1)
    WRITE(packing,'(i7)') j
  END IF

!
! 1.2 Line 2.
!
! Time-processing
  j=stlist(st_proc_no_code)
  tser_list=0
  IF (j == st_replace_code) THEN
    time_proc='   EXTRACT     '
  ELSE IF (j == st_accum_code) THEN
    time_proc=' ACCUMULATION  '
  ELSE IF (j == st_time_mean_code) THEN
    time_proc='  TIME MEAN    '
  ELSE IF (j == st_time_series_code) THEN
    WRITE(time_proc,'(''  TIME SERIES  '')')
    tser_list=stlist(st_series_ptr)
  ELSE IF (j == st_max_code) THEN
    time_proc='MAX OVER PERIOD'
  ELSE IF (j == st_min_code) THEN
    time_proc='MIN OVER PERIOD'
  ELSE IF (j == st_append_traj_code) THEN
    time_proc='  TRAJECTORY   '
  ELSE IF (j == st_variance_code) THEN
    time_proc=' TIME VARIANCE '
  ELSE IF (j == st_time_series_mean) THEN
    time_proc='MEAN TIMESERIES'
  ELSE
    time_proc='  UNKNOWN      '
  END IF

! -From-
  IF (stlist(st_freq_code) <  0) THEN
    from='      '
  ELSE
    j=stlist(st_start_time_code)
    WRITE(from,'(i6)') j
  END IF

! --To--
  IF (stlist(st_freq_code) <  0) THEN
    to='      '
  ELSE
    j=stlist(st_end_time_code)
    IF (j == st_infinite_time) THEN
      to='  ... '
    ELSE
      WRITE(to,'(i6)') j
    END IF
  END IF

! Frequency
  j=stlist(st_freq_code)
  IF (j <  0) THEN
    j=-j
    WRITE(freq,'(''TIME LIST'')')
    time_list=j
  ELSE
    WRITE(freq,'(i9)') j
    time_list=0
  END IF
      
! Period
  IF (stlist(st_freq_code) <  0) THEN
    period='      '
  ELSE
    j=stlist(st_period_code)
    IF (stlist(st_proc_no_code) == st_replace_code) THEN
      period='      '
    ELSE IF (j == st_infinite_time) THEN
      period='  ... '
    ELSE
      WRITE(period,'(i6)') j
    END IF
  END IF

! __Source__
  j=stlist(st_input_code)
  IF (j == 0) THEN
    source='PROGNOSTIC'
  ELSE IF(j == 1) THEN
    source='  STWORK  '
  ELSE IF(j <  0) THEN
    j=-j
    WRITE(source,'(''DUMP #'',i4)') j
  ELSE
    source=' UNKNOWN  '
  END IF
      
! ___Destination___
  j=stlist(st_output_code)
  IF (j == 1) THEN
    IF (stlist(st_macrotag) >= 1000) THEN
      dest='MEAN PP VIA DUMP'
    ELSE IF (stlist(st_macrotag) >  0) THEN
      WRITE(dest,'(''DUMP WITH TAG '',i3)') stlist(st_macrotag)
    ELSE
      dest='      DUMP       '
    END IF
  ELSE IF(j == 2) THEN
    dest='   SECONDARY     '
  ELSE IF(j <  0) THEN
    j=-j
    IF (j == 27) THEN
      dest='MEAN PP (DIRECT) '
    ELSE
      WRITE(dest,'(''   PP UNIT'',i3)') j
    END IF
  ELSE
    dest='  UNKNOWN  '
  END IF


!
! 1.3 Line 3.
!
! Spatial-Processing
  j=stlist(st_gridpoint_code)
  IF (j >= extract_base.AND.j <  extract_top) THEN
    spatial='    FULL FIELD    '
  ELSE IF (j >= vert_mean_base.AND.j <  vert_mean_top) THEN
    spatial='  VERTICAL MEAN   '
  ELSE IF (j >= zonal_mean_base.AND.j <  zonal_mean_top) THEN
    spatial='   ZONAL MEAN     '
  ELSE IF (j >= merid_mean_base.AND.j <  merid_mean_top) THEN
    spatial=' MERIDIONAL MEAN  '
  ELSE IF (j >= field_mean_base.AND.j <  field_mean_top) THEN
    spatial=' FIELD MEAN - 2D  '
  ELSE IF (j >= global_mean_base.AND.j <  global_mean_top) THEN
    spatial=' GLOBAL MEAN - 3D '
  ELSE
    spatial='  ** UNKNOWN **   '
  END IF

! Levels-domain
  j=stlist(st_output_bottom)
  lev_list=0
  IF (j == st_special_code) THEN
    levels ='STANDARD LEV '
  ELSE IF (j >  0) THEN
    WRITE(levels,'(''LEVELS '',i3,''-'',i3)') j,stlist(st_output_top)
  ELSE IF (j <  0) THEN
    j=-j
    WRITE(levels,'('' LEVELS LIST '')')
    lev_list=j
  END IF

! Pseudo-levels
  j=stlist(st_pseudo_out)
  plev_list=0
  IF (j >  0) THEN
    pseudo='PSEUDO-LEV LIST'
    plev_list=j
  ELSE
    pseudo='     NONE      '
  END IF

! Horizontal-domain..... suports up 9999 in each direction.

  WRITE(horiz,'(''ROW:'',i4,''-'',i4,'' COL:'',i4,''-'',i4)')           &
       stlist(st_south_code),stlist(st_north_code),                     &
       stlist(st_west_code),stlist(st_east_code)

! Weighting
  j=stlist(st_weight_code)
  IF (j == stash_weight_null_code) THEN
    weight='  NONE   '
  ELSE IF (j == stash_weight_area_code) THEN
    weight='  AREA   '
  ELSE IF (j == stash_weight_volume_code) THEN
    weight=' VOLUME  '
  ELSE IF (j == stash_weight_mass_code) THEN
    weight='  MASS   '
  END IF

! Masking
  j=mod(stlist(st_gridpoint_code),block_size)
  IF (j == stash_null_mask_code) THEN
    mask=' NONE  '
  ELSE IF (j == stash_land_mask_code) THEN
    mask=' LAND  '
  ELSE IF (j == stash_sea_mask_code) THEN
    mask='  SEA  '
  ELSE
    mask='UNKNOWN'
  END IF

!
! 1.4 Print the main part of the summary
!

IF(mype == 0) THEN
  WRITE(200,'(a)') diag1
  WRITE(200,'(1x,i3,1x,a,1x,a,1x,i4,1x,i7,1x,i7,'//                     &
              '1x,a,1x,a,1x,a,1x,i7,1x,i7,1x,a)')                       &
       seqno,name,submodel,stlist(st_item_code),stlist(st_sect_no_code),&
       ppxref(ppx_field_code),datatype,grid_type,level_type,            &
       ppxref(ppx_meto8_levelcode),ppxref(ppx_meto8_fieldcode),packing
  WRITE(200,'(a)') diag2       
  WRITE(200,'(7(1x,a))') &
         time_proc,from,to,freq,period,source,dest  
  WRITE(200,'(a)') diag3   
  WRITE(200,'(6(1x,a))') spatial,levels,pseudo,horiz,weight,mask     
END IF ! my = 0

!
! 1.5 Print associated time and levels lists if appropriate
!
! 1.5.1 Time list
!
  IF (time_list /= 0 .AND. mype == 0) THEN
    DO j=1,nsttims
      IF (sttabl(j,time_list) == st_end_of_list) THEN
        ntimes=j-1
        EXIT
      END IF
    END DO
        
    WRITE(200,'(a,i3,a)')                                               &
        ' ***** TIME LIST ***** ',ntimes,' times are as follows:-'

    WRITE(200,'(10(1x,i7))') sttabl(1:ntimes,time_list)       
        
  END IF
!
! 1.5.2 Levels list
!
  IF (lev_list /= 0 .AND. mype == 0) THEN
  
    WRITE(200,'(a,i3,a)')                                               &
        ' ***** LEVELS LIST ***** ',stash_levels(1,lev_list),           &
        ' levels are as follows:-'
    
    k=1+stash_levels(1,lev_list)
        
    WRITE(200,'(10(1x,i7))') stash_levels(2:k,lev_list)          
        
  END IF
  
!
! 1.5.3 Pseudo-levels list
!
  IF (plev_list /= 0 .AND. mype == 0) THEN
    WRITE(200,'(a,i3,a)')                                               &
        ' ***** PSEUDO-LEVELS LIST ***** ',                             &
         stash_pseudo_levels(1,plev_list),                              &
        ' pseudo-levels are as follows:-'

    k=1+stash_pseudo_levels(1,plev_list)
    WRITE(200,'(10(1x,i7))')  stash_pseudo_levels(2:k,plev_list)     
        
  END IF

!
! 1.5.4 Time series subdomain record list
!
  IF (tser_list /= 0 .AND. mype == 0) THEN
    i1=stash_series_index(1,tser_list)
    i2=stash_series_index(2,tser_list)
    WRITE(200,'('' ***** TIME SERIES ***** '',i3,  '//&
        ' '' subdomain records are as follows:-''/ '//&
        ' '' Record      North/South       West/ East     Bottom/  Top'')')&
        i2
    DO j=1,i2
      WRITE(200,'(3x,i4,1x,3(5x,i5,1x,i5,1x))')                         &
 &        j,(stash_series(3+k,i1+j-1),k=1,6)
    END DO
  END IF
!
! 1.5.5 Print final ruler line
!
IF (mype == 0) THEN
  WRITE(200,'(a)') headend
END IF


END IF
IF (lhook) CALL dr_hook('DIAGDESC',zhook_out,zhook_handle)
RETURN
END SUBROUTINE DIAGDESC
