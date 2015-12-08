! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Routine: STMERM ---------------------------------------------------
!
!    Purpose: Calculate weighted meridional mean within a region
!             specified by lower left hand and upper right hand corner.
!             (STASH service routine).
!
!
!
!    Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!
!    Project task: D7
!
!    External documentation:
!      Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                   system (STASH)
!
!    Interface and arguments: ------------------------------------------


!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: STASH
SUBROUTINE stmerm(fieldin,vx,vy,fld_type,gr,halo_type,            &
                  lwrap,lmasswt,                                  &
                  xstart,ystart,xend,yend,                        &
                  global_xstart,global_ystart,                    &
                  global_xend,global_yend,                        &
                  fieldout,                                       &
                  pstar_weight,                                   &
                  area_weight,mask,                               &
                  level_code,mask_code,weight_code,rmdi,          &
                  icode,cmessage)


USE um_parvars
USE global_2d_sums_mod, ONLY: global_2d_sums
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE sterr_mod, ONLY: st_upper_less_lower, st_not_supported,       &
                     st_no_data,st_nd, st_bad_array_param,        &
                     st_bad_address, st_unknown,                  &
                     st_bad_wraparound, st_illegal_weight,        &
                     unknown_weight, unknown_mask,                &
                     unknown_processing, nonsense
IMPLICIT NONE

INTEGER, INTENT(IN) ::  vx                ! input field x size
INTEGER, INTENT(IN) ::  vy                ! input field y size
INTEGER, INTENT(IN) ::  fld_type          ! field type (u/v/p)
INTEGER, INTENT(IN) ::  gr                ! input fld grid
INTEGER, INTENT(IN) ::  halo_type         ! halo type
INTEGER, INTENT(IN) ::  xstart            ! lower LH corner
INTEGER, INTENT(IN) ::  ystart            ! lower LH corner
INTEGER, INTENT(IN) ::  xend              ! upper RH corener
INTEGER, INTENT(IN) ::  yend              ! upper RH corner
INTEGER, INTENT(IN) ::  global_xstart     ! global version
INTEGER, INTENT(IN) ::  global_ystart     ! global version
INTEGER, INTENT(IN) ::  global_xend       ! global version
INTEGER, INTENT(IN) ::  global_yend       ! global version
INTEGER, INTENT(IN) ::  level_code        ! input level code
INTEGER, INTENT(IN) ::  mask_code         ! masking code
INTEGER, INTENT(IN) ::  weight_code       ! weighting code
INTEGER, INTENT(INOUT) :: icode           ! error return code

CHARACTER(LEN=*),INTENT(INOUT) :: cmessage  ! error return msg

LOGICAL, INTENT(IN) :: lwrap              ! TRUE if wraparound
LOGICAL, INTENT(IN) :: lmasswt            ! TRUE if masswts OK
LOGICAL, INTENT(IN) :: mask(vx+1,vy)      ! mask array

REAL, INTENT(IN) ::  rmdi                   ! missing data indicator
REAL, INTENT(IN) ::  fieldin(vx,vy)         ! input field
REAL, INTENT(IN) ::  pstar_weight(vx+1,vy)  ! mass weight factor
REAL, INTENT(IN) ::  area_weight(vx+1,vy)   ! area weight factor

REAL, INTENT(OUT) :: fieldout(xstart:xend)  ! output field

!*----------------------------------------------------------------------

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

INTEGER :: i              ! Looper
INTEGER :: ii             ! Looper
INTEGER :: j              ! Looper



! size of subarea on this processor, not including halo areas
INTEGER :: local_sum_xstart
INTEGER :: local_sum_xend
INTEGER :: local_sum_ystart
INTEGER :: local_sum_yend

! number of columns of meridional data to sum on this processor
! (xend-xstart+1)
INTEGER :: partial_sum_data_sizex
INTEGER :: partial_sum_data_sizey

! Sums and partial sums
REAL, ALLOCATABLE :: partial_sumtop_2d(:,:)
REAL, ALLOCATABLE :: partial_sumbot_2d(:,:)
REAL, ALLOCATABLE :: sumtop(:)
REAL, ALLOCATABLE :: sumbot(:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('STMERM',zhook_in,zhook_handle)

! pstar_weight and area_weight arrays contain appropriate
! weighting factors, interpolated to the correct grid, for
! mass weighting and area weighting respectively. If either type
! of weighting is not required, the relevant array is set to 1.0
! The mask array contains appropriate masking


! 1.0 Find the bounds of the actual data required in the summation
!    (ie. excluding the halos, contained within
!    xstart,xend,ystart,yend.

! DEPENDS ON: global_to_local_subdomain
CALL global_to_local_subdomain(.FALSE.,.FALSE.,                   &
  gr,halo_type,mype,                                              &
  global_ystart,global_xend,                                      &
  global_yend,global_xstart,                                      &
  local_sum_ystart,local_sum_xend,                                &
  local_sum_yend,local_sum_xstart)

IF (local_sum_xstart  >   local_sum_xend)                         &
  local_sum_xend=local_sum_xend+vx-                               &
                 2*halosize(1,halo_type)

! 1.1 And the number of partial sums on this processor
!     (If there is a wrap around, with the subarea starting
!      and ending on this processor, the halo cols are included
!      in this number - although no sums will actually be
!      carried out there)

partial_sum_data_sizex=                                           &
  local_sum_xend-local_sum_xstart+1
partial_sum_data_sizey=                                           &
  local_sum_yend-local_sum_ystart+1


! 1.2 Initialise the sum arrays
!     Note that partial arrays are transposed to make global sums
!     possible and straightforward later...

ALLOCATE (partial_sumtop_2d(local_sum_ystart:local_sum_yend, &
                      local_sum_xstart:local_sum_xend))
ALLOCATE (partial_sumbot_2d(local_sum_ystart:local_sum_yend, &
                      local_sum_xstart:local_sum_xend))
ALLOCATE (sumtop(local_sum_xstart:local_sum_xend))
ALLOCATE (sumbot(local_sum_xstart:local_sum_xend))
partial_sumtop_2d(:,:) = 0.0
partial_sumbot_2d(:,:) = 0.0
sumtop(:) = 0.0
sumbot(:) = 0.0

! 2.0 Calculate the partial sums

! Only do the calculations if some of the subdomain exists on this
! processor

IF ( (local_sum_xstart  /=  st_no_data) .AND.                     &
     (local_sum_xend    /=  st_no_data) .AND.                     &
     (local_sum_ystart  /=  st_no_data) .AND.                     &
     (local_sum_yend    /=  st_no_data)) THEN

! 2.2 Do the actual sum

  DO i=local_sum_xstart,local_sum_xend

      IF ( lwrap .AND.                                            &
          (i  >   (lasize(1,fld_type,halo_type)-                  &
                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
        ii=i-blsize(1,fld_type)
    ELSE
      ii=i
    END IF

    DO j=local_sum_ystart,local_sum_yend
      IF (mask(ii,j)) THEN

          partial_sumbot_2d(j,i)=                  &
            pstar_weight(ii,j)*area_weight(ii,j)
          partial_sumtop_2d(j,i)=                  &
            fieldin(ii,j)*pstar_weight(ii,j)*area_weight(ii,j)
      END IF ! if this point is to be processed
    END DO ! j : loop over rows
  END DO ! i : loop over columns
END IF ! if subdomain covers this processor

! 3.0 Sums up the partial sums down each column to make a full sum

! So a sum down the processor column if the subdomain covers any
! processor(s) along the column

IF ((local_sum_xstart  /=  st_no_data) .AND.                      &
    (local_sum_xend  /=  st_no_data)) THEN

  CALL global_2d_sums(partial_sumbot_2d, partial_sum_data_sizey,  &
                      1, 0, 0, partial_sum_data_sizex,            &
                      sumbot, gc_proc_col_group)
  CALL global_2d_sums(partial_sumtop_2d, partial_sum_data_sizey,  &
                      1, 0, 0, partial_sum_data_sizex,            &
                      sumtop, gc_proc_col_group)

! So now the partial_* arrays actually contain the full sums
! along the column

! 3.1 And put the mean meridional values into the fieldout array

! Only processors in the subdomain area need to record the
! results
  IF ((local_sum_ystart  /=  st_no_data) .AND.                    &
      (local_sum_yend  /=  st_no_data)) THEN

    DO i=local_sum_xstart,local_sum_xend
      ii = (local_sum_xstart - xstart) + i
      IF (sumbot(i)  ==  0.0) THEN
        fieldout(ii)=rmdi
      ELSE
        fieldout(ii)=sumtop(i)/sumbot(i)
      END IF
    END DO

  END IF ! is this processor in the subdomain

END IF ! does the subdomain intersect with this processor column

IF (lhook) CALL dr_hook('STMERM',zhook_out,zhook_handle)
RETURN
END SUBROUTINE stmerm
