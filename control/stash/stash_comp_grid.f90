! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE STASH_COMP_GRID----------------------
!LL
!LL  Compute grid descriptors
!LL  out_bdx,out_bdy,out_bzx and out_bzy from inputs
!LL  1)  samples
!LL  2)  grid type code
!LL  3)  ocean (REMOVED)
!LL  4)  real and integer headers
!LL  5)  w_mostcol and n_mostrow
!LL  6) processing code
!LL
!LL  Programming standard: U M DOC  Paper NO. 3
!LL
!LL  External documentation  C4
!LL
!LLEND-------------------------------------------------------------

!
!*L  INTERFACE and ARGUMENTS:------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH
      SUBROUTINE STASH_COMP_GRID(                                       &
     &  out_bzx,out_bzy,out_bdx,out_bdy,                                &
     &  samples,st_grid,                                                &
     &  w_mostcol,n_mostrow,                                            &
     &  realhd,len_realhd,inthd,len_inthd,gr,                           &
     &  ICODE,CMESSAGE)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE lookup_addresses

      IMPLICIT NONE


!
      INTEGER samples                ! IN no of samples in period (times
!
      INTEGER                                                           &
     &  ICODE                                                           &
                          !OUT   Return code from the routine
     &, LEN_REALHD                                                      &
                          !IN    Length of the Real Constants
     &, LEN_INTHD                                                       &
                          !IN    Length of integer constants
     &, INTHD(LEN_INTHD)  !IN    Integer constants
!
      INTEGER                                                           &
     &  st_grid                                                         &
                          !IN    STASH horizontal grid type
     &, N_MOSTROW                                                       &
                          !IN    The most nrthrly row.
     &, W_MOSTCOL                                                       &
                          !IN    The most westerly column
     &, gr                !IN    The type of processing done
!
      REAL                                                              &
     &  REALHD(LEN_REALHD)  !IN  Real header
      REAL                                                              &
     &  out_bzx,out_bdx,out_bzy,out_bdy  ! OUT grid descriptors
      CHARACTER(LEN=*)                                                     &
     &  cmessage           ! OUT error messages
!*
!
!L Local Variables
!
      INTEGER mean_code

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!*---------------------------------------------------------------------
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
! Description: This include file contains information needed when
!              generating variable horizontal grid data in the
!              STASH extra data vector. Introduced UM 5.4 - R. Hill
!===================================================================
      LOGICAL :: X_VAR_GRID ! Whether variable grid in E-W direction
      LOGICAL :: Y_VAR_GRID ! and/or in S-N direction

      INTEGER :: VAR_GRID_TYPE ! 0 = none
                               ! 1 = T grid
                               ! 2 = U/V grid

      ! Grid boundaries for T and U,V
      REAL :: X_BOUNDARY(ROW_LENGTH_MAX+1,2)
      REAL :: Y_BOUNDARY(ROWS_MAX+1,2)

      ! Grid Points for T and U,V
      REAL :: X_GRID(ROW_LENGTH_MAX,2)
      REAL :: Y_GRID(ROWS_MAX,2)

      COMMON /OVARGRID/ X_VAR_GRID,Y_VAR_GRID                           &
     & ,X_BOUNDARY,Y_BOUNDARY,X_GRID,Y_GRID,VAR_GRID_TYPE

      ! The following parameters correspond to the extra data
      ! vector descriptors expected, for e.g., in PV-WAVE
      ! plotting routines (e.g. decode_extra.pro). There are
      ! numerous other areas of code where these integer
      ! descriptors must be handled (e.g. FIELDCOS, PPI2H, FTT)
      ! So it is not a trivial matter to introduce new code
      ! descriptors. Furthermore, ieee -32 will destroy these
      ! integers so PP data must always be processed via the
      ! long winded route: QXFIELDCOS -> PUTIBM -> FTT/PPI2H.
      ! (This thoroughly unsatisfactory state of affairs may
      ! be correctable with developments to ieee and convpp).
      INTEGER,PARAMETER :: x_coord_vector=1
                                     ! Indicates that an extra
                                     ! data vector gives LBNPT
                                     ! x-coordinate values
      INTEGER,PARAMETER :: y_coord_vector=2
                                     ! Indicates that an extra
                                     ! data vector gives LBROW
                                     ! Y-coordinate values
      INTEGER,PARAMETER :: x_lbnd_vector=12
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! x-boundary values
      INTEGER,PARAMETER :: x_ubnd_vector=13
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! x-boundary values
      INTEGER,PARAMETER :: y_lbnd_vector=14
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! y-boundary values
      INTEGER,PARAMETER :: y_ubnd_vector=15
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! y-boundary values

!LL   Construct PP header
      IF (lhook) CALL dr_hook('STASH_COMP_GRID',zhook_in,zhook_handle)
      IF (samples >  0) THEN   ! Indicates a timeseries/trajectory
        OUT_BZX=0.0
        OUT_BDX=0.0
        OUT_BZY=0.0
        OUT_BDY=0.0
      ELSE
          !   set OUT_BZY,OUT_BZX,OUT_BDY,OUT_BDX for
          IF(st_grid == st_riv_grid)THEN
            OUT_BDY = 180.0/180
            OUT_BZY = REALHD(3) - OUT_BDY*0.5
            OUT_BDX = 360.0/360
            OUT_BZX = REALHD(4) - OUT_BDX*0.5
          ELSE
            IF(st_grid == st_uv_grid.OR.st_grid == st_cv_grid) THEN
              OUT_BZY=REALHD(3)-REALHD(2)/2.0 ! UV pts
            ELSE
              OUT_BZY=REALHD(3)-REALHD(2) ! Zeroth Lat OUT_BZY
            ENDIF
!
            IF(st_grid == st_uv_grid.OR.st_grid == st_cu_grid) THEN
              OUT_BZX=REALHD(4)-REALHD(1)/2.0 !UV points
            ELSE
              OUT_BZX=REALHD(4)-REALHD(1) ! Zeroth Long OUT_BZX
            ENDIF
            OUT_BDX=REALHD(1) ! Long intvl OUT_BDX
            OUT_BDY=REALHD(2) ! Lat intvl OUT_BDY

          ENDIF
!
! Add on offset for fields not starting from the origin
!
      ! If this is variable horizontal grid information
      ! then the horizontal start points are NOT equally spaced
      ! we must get our actual start point some other way!
      ! This is all theoretical since this routine currently
      ! gets called with hard coded values of N_MOSTROW
      ! and W_MOSTCOL. Furthermore, the use of N_MOSTROW
      ! is entirely wrong in this routine!
      ! It should be the southern-most row for this to do
      ! anything useful!
      IF (VAR_GRID_TYPE >  0) THEN
         IF (X_VAR_GRID) OUT_BZX = X_BOUNDARY(W_MOSTCOL,VAR_GRID_TYPE)
         IF (Y_VAR_GRID) OUT_BZY = Y_BOUNDARY(N_MOSTROW,VAR_GRID_TYPE)
      ELSE
        OUT_BZY=OUT_BZY                                                 &
     &       +(N_MOSTROW-1)*OUT_BDY
        OUT_BZX=OUT_BZX                                                 &
     &     +(W_MOSTCOL-1)*OUT_BDX

      ENDIF
        IF(OUT_BZX >= 360.0)                                            &
     &     OUT_BZX=OUT_BZX-360.0
!
! If horizontal averaging has been applied to the output field,
! set OUT_BDX and/or OUT_BDY to the full domain extent
!
        mean_code=(GR/block_size)*block_size
        IF (mean_code == zonal_mean_base .OR.                           &
     &      mean_code == field_mean_base .OR.                           &
     &      mean_code == global_mean_base) THEN
          OUT_BDX=REAL(INTHD(6))*REALHD(1)
        ENDIF
        IF (mean_code == merid_mean_base .OR.                           &
     &      mean_code == field_mean_base .OR.                           &
     &      mean_code == global_mean_base) THEN
          OUT_BDY=REAL(INTHD(7))*REALHD(2)
        ENDIF
      ENDIF
!
      IF (lhook) CALL dr_hook('STASH_COMP_GRID',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE STASH_COMP_GRID
