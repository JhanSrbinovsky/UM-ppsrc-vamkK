! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Gathers STASHed data from many processors to one processor

! Subroutine interface:
SUBROUTINE stash_gather_field (                                   &
   local_field , global_field ,                                    &
   local_size, global_size, levels,                                &
   global_north , global_east_in , global_south , global_west,     &
   gridtype_code ,halo_type,                                       &
   gather_pe,                                                      &
   data_extracted,                                                 &
   packing, im_ident, packing_type,                          &
   num_out,                                                        &
   comp_accrcy, loc_rmdi,                                          &
   icode, cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE gcom_mod
USE UM_ParVars
USE cppxref_mod, ONLY:                                            &
    ppx_atm_tzonal, ppx_atm_uzonal,                               &
    ppx_ocn_tzonal, ppx_ocn_uzonal
IMPLICIT NONE

! Description:
! Takes a decomposed, STASH processed field and gathers
! it to a single processor, ready for I/O,

! Method:
! See in-line documentation

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine arguments:


INTEGER, INTENT(IN) ::                                            &
  local_size                                                      &
                  ! IN: size of level of LOCAL_FIELD
, global_size                                                     &
                  ! IN: size of level of GLOBAL_FIELD
, levels                                                          &
                  ! IN: number of levels
, global_north                                                    &
                  ! IN: specification of subdomain boundaries
, global_east_in                                                  &
                  ! IN: ""
, global_south                                                    &
                  ! IN: ""
, global_west                                                     &
                  ! IN: ""
, gridtype_code                                                   &
                  ! IN: indicates the type of grid output
, halo_type                                                       &
                  ! IN: type of halo on this field
, gather_pe       ! IN: the PE to gather the global field to

INTEGER, INTENT(OUT) ::                                           &
  icode           ! OUT: return code, 0=OK

! Optional Arguments to handle the COEX packing if necessary

LOGICAL, INTENT(IN), OPTIONAL ::                                  &
  packing
                  ! IN: Set .true. if packing of the input
                  !     field is to be packed!

INTEGER, INTENT(IN), OPTIONAL ::                                  &
  im_ident        ! IN: Internal model identifier

INTEGER, INTENT(INOUT), OPTIONAL ::                               &
  packing_type    ! IN/OUT: This flag is zero on input,
                  !         then stash packing is selected,
                  !         and the routine returns the
                  !         packing flag.
                  !
                  !         If the variable is set to 1 on input
                  !         then 32-bit packing for dumpfiles
                  !         is selected

INTEGER, INTENT(OUT), OPTIONAL ::                                 &
  num_out         ! OUT: Number of 32-bit IBM words in the Packed
                  !      field for WDGOS packing

INTEGER, INTENT(IN), OPTIONAL ::                                  &
  comp_accrcy     ! IN: Packing Accuracy in Power of 2

REAL, INTENT(IN), OPTIONAL ::                                     &
  loc_rmdi        ! IN: Missing data indicator

! Remaining Non-Optional Arguments

LOGICAL, INTENT(IN) ::                                            &
  data_extracted  ! IN: TRUE if the data in LOCAL_FIELD has
                  !     already been extracted, or FALSE if
                  !     the extraction must be done here.

REAL, INTENT(IN) ::                                               &
  local_field(local_size,levels)
                  ! IN : local data

REAL, INTENT(OUT) ::                                              &
  global_field(global_size,levels)
                  ! OUT : (PE GATHER_PE only) - full gathered
                  !       field

CHARACTER(LEN=*), INTENT(OUT) ::                                     &
  cmessage        ! OUT: Error message if ICODE  /=  0

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
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Local variables

INTEGER    ::                                                     &
  global_east                                                     &
                ! copy of GLOBAL_EAST_IN with wrap around s.t.
!                     ! GLOBAL_EAST > GLOBAL_ROW_LEN
, global_x                                                        &
                ! size of global data EW
, global_y                                                        &
                ! size of global data NS
, fld_type                                                        &
                ! indicates if field is on P or U grid
, level                                                           &
                ! loop index for loop over levels
, i                                                               &
                ! loop counter
, proc_topleft_x,proc_topleft_y                                   &
                                  ! processors at corners of
, proc_botright_x,proc_botright_y                                 &
                                  ! the subarea
, dummy1,dummy2                                                   &
                ! ignored return arguments
, procx,procy                                                     &
                ! loop indexes for loops over processors
, eff_procx                                                       &
                ! real x co-ord of processor column procx
, procid                                                          &
                ! processor id of (procx,procy)
, local_xstart,local_xend                                         &
                           ! boundaries of subdomain for
, local_ystart,local_yend                                         &
                           ! processor procid
, local_start_row                                                 &
                    ! first row to send from procid
, local_start_col                                                 &
                    ! first column to send from procid
, sendsize_x                                                      &
                    ! number of points on each row to send
, nrows_to_send                                                   &
                    ! number of rows to send from procid
, local_row_length                                                &
                    ! size of sending array EW
, global_start_row                                                &
                    ! first row to receive at on GATHER_PE
, global_start_col                                                &
                    ! first col. to recieve on GATHER_PE
, global_row_length                                               &
                    ! size of receiving array EW
, flag,info         ! GCOM arguments

! Copies of arguments / variables used to decide if we can use the
! send/receive maps used in the last call

! Save all the variables that may be used in the next call

INTEGER, SAVE ::                                                  &
  old_local_size , old_global_size                                &
, old_global_north , old_global_east_in                           &
, old_global_south , old_global_west                              &
, old_gridtype_code , old_gather_pe                               &
, old_current_decomp_type, old_halo_type

! variables defining send and receive maps to be passed to
! GC_RALL_TO_ALL, defining the data transposition
INTEGER, SAVE :: send_map(7,2), n_sends, n_recvs 
INTEGER, SAVE, ALLOCATABLE :: receive_map(:,:)


LOGICAL  ::                                                       &
  wrap                                                            &
               ! if the subdomain wraps over 0 degree meridion
, wrap_special                                                    &
               ! if there is a wrap around, which starts and
!                      ends on the same processor
, zonal_data                                                      &
               ! if this is a zonal data grid
, fullfield                                                       &
               ! if this is a full field - NOT a subarea
, l_packing    ! if packing of data is required



! Set all the old_* variables to a number indicating they've
! not been used yet

DATA                                                              &
  old_local_size , old_global_size                                &
, old_global_north , old_global_east_in                           &
, old_global_south , old_global_west                              &
, old_gridtype_code , old_gather_pe                               &
, old_current_decomp_type, old_halo_type                          &
  / -1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /

! Functions

INTEGER :: get_fld_type

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! ------------------------------------------------------------------

IF (lhook) CALL dr_hook('STASH_GATHER_FIELD',zhook_in,zhook_handle)

IF (.NOT.ALLOCATED(receive_map)) &
  ALLOCATE(receive_map(7,2*nproc_max))

icode=0
IF (PRESENT(packing)) THEN
  l_packing = packing
ELSE
  l_packing = .FALSE.
END IF

! DEPENDS ON: get_fld_type
 fld_type=get_fld_type(gridtype_code)
! See if there is wrap around over meridion, and if so make
! sure that GLOBAL_EAST is > glsize(1)

global_east=global_east_in
IF (global_east  >   glsize(1,fld_type)) THEN
  wrap=.TRUE.
ELSE IF (global_east  <   global_west) THEN
  wrap=.TRUE.
  global_east=global_east_in+glsize(1,fld_type)
ELSE
  wrap=.FALSE.
END IF

IF ((gridtype_code  ==  ppx_atm_tzonal) .OR.                      &
                                             ! Atmos T zonal
   ( gridtype_code  ==  ppx_atm_uzonal) .OR.                      &
                                             ! Atmos U zonal
   ( gridtype_code  ==  ppx_ocn_tzonal) .OR.                      &
                                             ! Ocean T zonal
   ( gridtype_code  ==  ppx_ocn_uzonal))                          &
                                             ! Atmos U zonal
  THEN

! This is a zonal field

  zonal_data=.TRUE.
  global_x=1

  IF ((gridtype_code  ==  ppx_atm_tzonal) .OR.                    &
                                               ! Atmos T zonal
      ( gridtype_code  ==  ppx_ocn_tzonal))                       &
                                               ! Ocean T zonal
  THEN
    fld_type=fld_type_p
  ELSE
    fld_type=fld_type_u
  END IF
ELSE

! This is a normal field

  zonal_data=.FALSE.
  global_x=glsize(1,fld_type)
END IF

global_y=glsize(2,fld_type)

! Set up logical indicating if this is a full field, or just
! a subdomain

IF (zonal_data) THEN

  fullfield= ( ( global_south  ==  1) .AND.                       &
             ( global_north  ==  global_y))

ELSE

  fullfield = (( global_west  ==  1) .AND.                        &
               ( global_east  ==  global_x) .AND.                 &
               ( global_south  ==  1) .AND.                       &
               ( global_north  ==  global_y))

END IF

! Dealing with fields not in model grid

IF((global_x == 0).OR.(global_x == imdi)) THEN
  WRITE(6,*)'local_size=',local_size
  WRITE(6,*)'global_size=',global_size
  DO level=1,levels
    DO i=1,global_size
      global_field(i,level)=local_field(i,level)
    END DO
  END DO
ELSE

! If this a fullfield, we can simply use the standard
! GATHER_FIELD routine

IF (fullfield) THEN

  IF (zonal_data) THEN

! DEPENDS ON: gather_zonal_field
    CALL gather_zonal_field( local_field,global_field,            &
                          lasize(2,fld_type,halo_type),global_y,  &
                            levels,gridtype_code,fld_type,        &
                            halo_type,gather_pe)

  ELSE

    DO level=1,levels

      IF (l_packing) THEN
! DEPENDS ON: gather_pack_field
        CALL gather_pack_field(                                   &
                        local_field(1,level),                     &
                        global_field(1,level),                    &
                        lasize(1,fld_type,halo_type),             &
                        lasize(2,fld_type,halo_type),             &
                        global_x,global_y,                        &
                        fld_type,halo_type,                       &
                        gather_pe,gc_all_proc_group,              &
                        packing, im_ident, packing_type,    &
                        num_out,                                  &
                        comp_accrcy, loc_rmdi)
      ELSE
! DEPENDS ON: gather_field
        CALL gather_field( local_field(1,level),                  &
                           global_field(1,level),                 &
                           lasize(1,fld_type,halo_type),          &
                           lasize(2,fld_type,halo_type),          &
                           global_x, global_y,                    &
                           fld_type, halo_type,                   &
                           gather_pe, gc_all_proc_group,          &
                           icode, cmessage)
      END IF

      IF (icode  /=  0) THEN
        WRITE(6,*)                                                &
         'STASH_GATHER_FIELD : Failed during GATHER_FIELD'
        WRITE(6,*) 'Error code : ',icode
        WRITE(6,*) 'Message : ',cmessage

        icode=2
        cmessage='Failed to gather field'
        GO TO 9999
      END IF

     END DO

   END IF
 ELSE
! for subdomains, life is not so easy - we must explicitly
! calculate our own send and receive maps, and use GCG_RALLTOALLE
! to shift the data around.

! If the same arguments are used as were used in the last call
! to this routine, we can just use the previously calculated
! send and receive maps, otherwise we need to calculate new maps

  IF (.NOT. (                                                     &
    (local_size  ==  old_local_size) .AND.                        &
    (global_size  ==  old_global_size) .AND.                      &
    (global_north  ==  old_global_north) .AND.                    &
    (global_east_in  ==  old_global_east_in) .AND.                &
    (global_south  ==  old_global_south) .AND.                    &
    (global_west  ==  old_global_west) .AND.                      &
    (gridtype_code  ==  old_gridtype_code) .AND.                  &
    (gather_pe  ==  old_gather_pe) .AND.                          &
    (halo_type  ==  old_halo_type) .AND.                          &
    (current_decomp_type  ==  old_current_decomp_type ))) THEN

    old_local_size=local_size
    old_global_size=global_size
    old_global_north=global_north
    old_global_east_in=global_east_in
    old_global_south=global_south
    old_global_west=global_west
    old_gridtype_code=gridtype_code
    old_gather_pe=gather_pe
    old_current_decomp_type=current_decomp_type
    old_halo_type=halo_type

! Find out what the boundaries of the subdomain area

! DEPENDS ON: global_to_local_rc
    CALL global_to_local_rc(gridtype_code,halo_type,              &
                            global_west,global_north,             &
                            proc_topleft_x,proc_topleft_y,        &
                            dummy1,dummy2)
! DEPENDS ON: global_to_local_rc
    CALL global_to_local_rc(gridtype_code,halo_type,              &
                            global_east,global_south,             &
                            proc_botright_x,proc_botright_y,      &
                            dummy1,dummy2)

! Ensure that the processor x co-ords are such that the botright_x is
! always greater than (or equal to) top_left_x.
    IF (wrap) proc_botright_x=gridsize(1)+proc_botright_x

! wrap_special is set to true if there is a wrap around which starts
! and ends on the same processor. This case requires extra work as
! the processor in question
    IF (wrap .AND. (proc_topleft_x+gridsize(1)  ==                &
                    proc_botright_x)) THEN
      wrap_special=.TRUE.
    ELSE
      wrap_special=.FALSE.
    END IF

    n_sends=0
    n_recvs=0

    DO procy=proc_botright_y,proc_topleft_y
      DO procx=proc_topleft_x,proc_botright_x

        eff_procx=MOD(procx,gridsize(1))
        procid=eff_procx+procy*gridsize(1)

! DEPENDS ON: global_to_local_subdomain
        CALL global_to_local_subdomain(                           &
          .TRUE.,.TRUE.,                                          &
          gridtype_code,halo_type,procid,                         &
          global_south,global_east,                               &
          global_north,global_west,                               &
          local_ystart,local_xend,                                &
          local_yend  ,local_xstart)

! Calculate the shape of the arrays, and where to start sending/
! receiving data, and how many rows to send

        IF (data_extracted) THEN
          local_start_row=1
        ELSE
          local_start_row=local_ystart
        END IF
        nrows_to_send=local_yend-local_ystart+1

       global_start_row=g_datastart(2,procid)+local_ystart-       &
                         halosize(2,halo_type) - global_south
        global_row_length=global_east-global_west+1

! Calculate the following variables:
! local_row_length : the X dimension size of the local array
! local_send_offx  : the offset into each row to start sending from
! sendsize_x       : the number of points on each row to send
! The calculation of these numbers is different for processors
! at the start and end of a wrap_special case
! Note that when DATA_EXTRACTED is true, then local_field has no
! halos.

        IF (wrap_special .AND. procx  ==  proc_topleft_x) THEN
          IF (data_extracted) THEN
local_row_length = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
                   + local_xend - local_xstart + 1
sendsize_x       = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
                   - local_xstart + 1
local_start_col  = 1
          ELSE
local_row_length = g_lasize(1,fld_type,halo_type,procid)
sendsize_x       = g_lasize(1,fld_type,halo_type,procid)          &
                   - local_xstart
local_start_col  = local_xstart

          END IF
          global_start_col=1

        ELSE IF (wrap_special .AND. procx  ==  proc_botright_x)    &
        THEN
          IF (data_extracted) THEN
local_row_length = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
                   + local_xend - local_xstart + 1
local_start_col  = local_row_length - local_xend + 1
sendsize_x       = local_xend
          ELSE
local_row_length = g_lasize(1,fld_type,halo_type,procid)
local_start_col  = offx + 1
sendsize_x       = local_xend - offx
          END IF
          global_start_col=global_row_length-sendsize_x+1

        ELSE
          IF (data_extracted) THEN
            local_row_length=local_xend-local_xstart+1
            local_start_col=1
          ELSE
    local_row_length=g_lasize(1,fld_type,halo_type,procid)
            local_start_col=local_xstart
          END IF
          sendsize_x=local_xend-local_xstart+1
          global_start_col=local_xstart-halosize(1,halo_type)+    &
                           g_datastart(1,procid)-global_west
        END IF

        IF (global_start_col  <   0) THEN
! Wrapped around field, but this processor is not start or end
! processor
  global_start_col=global_start_col+glsize(1,fld_type)
        END IF

! Now we can set up the send and receive map entries for the data on
! this processor

        IF (mype  ==  procid) THEN  ! I need to send some data


          n_sends=n_sends+1

          send_map(s_destination_pe,n_sends) = gather_pe
          send_map(s_base_address_in_send_array,n_sends) =        &
            (local_start_row-1)*local_row_length +                &
            local_start_col
          send_map(s_number_of_elements_in_item,n_sends) =        &
            nrows_to_send
          send_map(s_stride_in_send_array,n_sends) =              &
            local_row_length
          send_map(s_element_length,n_sends) = sendsize_x
          send_map(s_base_address_in_recv_array,n_sends) =        &
            (global_start_row-1)*global_row_length +              &
            global_start_col
          send_map(s_stride_in_recv_array,n_sends) =              &
            global_row_length

        END IF ! if I'm sending data

        IF (mype  ==  gather_pe) THEN  ! I need to receive data

          n_recvs=n_recvs+1

          receive_map(r_source_pe,n_recvs) = procid
          receive_map(r_base_address_in_recv_array,n_recvs) =     &
            (global_start_row-1)*global_row_length +              &
            global_start_col
          receive_map(r_number_of_elements_in_item,n_recvs) =     &
            nrows_to_send
          receive_map(r_stride_in_recv_array,n_recvs) =           &
            global_row_length
          receive_map(r_element_length,n_recvs) = sendsize_x
          receive_map(r_base_address_in_send_array,n_recvs) =     &
            (local_start_row-1)*local_row_length +                &
            local_start_col
          receive_map(r_stride_in_send_array,n_recvs) =           &
            local_row_length

        END IF ! if I'm receiving data

      END DO ! procx : loop along processor row

    END DO ! procy : loop down processor column

  END IF ! if I need to recalculate my send/receive maps

! Send / receive the data using GCG_RALLTOALLE


  DO level=1,levels

    flag=0  ! This is currently ignored at GCG v1.1
    info=gc_none

    CALL gcg_ralltoalle(                                          &
      local_field(1,level)  ,                                     &
      send_map    , n_sends  ,local_size  ,                       &
      global_field(1,level) ,                                     &
      receive_map , n_recvs , global_size ,                       &
      gc_all_proc_group , flag, info)

  END DO

  END IF ! if this is a full or extracted field

END IF


 9999 CONTINUE

IF (lhook) CALL dr_hook('STASH_GATHER_FIELD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE stash_gather_field

