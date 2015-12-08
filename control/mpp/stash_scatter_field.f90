! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Scatters STASHed data from one processor to many processors

! Subroutine interface:
SUBROUTINE stash_scatter_field (                                  &
  local_field , global_field ,                                    &
  local_size, global_size, levels,                                &
  global_north , global_east_in , global_south , global_west,     &
  gridtype_code, halo_type,                                       &
  scatter_pe,                                                     &
  icode, cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE gcom_mod
USE UM_ParVars
USE cppxref_mod, ONLY:                                            &
    ppx_atm_tzonal, ppx_ocn_uzonal, ppx_ocn_tzonal
IMPLICIT NONE


! Description:
! Takes a decomposed, STASH processed field and gathers
! it to a single processor, ready for I/O,

! Method:
! See in-line documentation

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine arguments:

INTEGER, INTENT(IN) :: local_size      ! IN: size of level of LOCAL_FIELD
INTEGER, INTENT(IN) :: global_size     ! IN: size of level of GLOBAL_FIELD
INTEGER, INTENT(IN) :: levels          ! IN: number of levels
INTEGER, INTENT(IN) :: global_north    ! IN: specification of subdomain boundaries
INTEGER, INTENT(IN) :: global_east_in  ! IN: ""
INTEGER, INTENT(IN) :: global_south    ! IN: ""
INTEGER, INTENT(IN) :: global_west     ! IN: ""
INTEGER, INTENT(IN) :: gridtype_code   ! IN: indicates the type of grid output
INTEGER, INTENT(IN) :: halo_type       ! IN: type of halo on this field
INTEGER, INTENT(IN) :: scatter_pe      ! IN: the PE to scatter global field from

INTEGER, INTENT(OUT) :: icode          ! OUT: return code, 0=OK

REAL, INTENT(OUT) :: local_field(local_size,levels) 
                                       ! OUT : local scattered data
REAL, INTENT(OUT) :: global_field(global_size,levels)
                                       ! IN : (PE SCATTER_PE only) - full field

CHARACTER(LEN=80) :: cmessage          ! OUT: Error message if ICODE  /=  0

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

INTEGER  ::  global_east        ! copy of GLOBAL_EAST_IN with wrap around s.t.
                                !         GLOBAL_EAST > GLOBAL_ROW_LEN
INTEGER  ::  global_x           ! size of global data EW
INTEGER  ::  global_y           ! size of global data NS
INTEGER  ::  fld_type           ! indicates if field is on P or U grid
INTEGER  ::  level              ! loop index for loop over levels
INTEGER  ::  i                  ! loop counter
INTEGER  ::  proc_topleft_x,proc_topleft_y    ! processors at corners of
INTEGER  ::  proc_botright_x,proc_botright_y  ! the subarea
INTEGER  ::  dummy1,dummy2      ! ignored return arguments
INTEGER  ::  procx,procy        ! loop indexes for loops over processors
INTEGER  ::  eff_procx          ! real x co-ord of processor column procx
INTEGER  ::  procid             ! processor id of (procx,procy)
INTEGER  ::  local_xstart,local_xend          ! boundaries of subdomain for
INTEGER  ::  local_ystart,local_yend          ! processor procid
INTEGER  ::  local_start_row    ! first row to receive on procid
INTEGER  ::  local_start_col    ! first column to receive on procid
INTEGER  ::  sendsize_x         ! number of points on each row to send to procid
INTEGER  ::  nrows_to_send      ! number of rows to send to procid
INTEGER  ::  local_row_length   ! size of receiving array EW
INTEGER  ::  global_start_row   ! first row to send on SCATTER_PE
INTEGER  ::  global_start_col   ! first col. to send on SCATTER_PE
INTEGER  ::  global_row_length  ! size of sending array EW
INTEGER  ::  flag,info          ! GCOM arguments

LOGICAL  ::  l_vec              ! Indicates if a field is a vector quantity
                                ! (SWAPBOUNDS argument)
! Copies of arguments / variables used to decide if we can use the
! send/receive maps used in the last call

! Save all the variables that may be used in the next call

INTEGER, SAVE ::                                                  &
  old_local_size , old_global_size                                &
, old_global_north , old_global_east_in                           &
, old_global_south , old_global_west                              &
, old_gridtype_code , old_scatter_pe                              &
, old_halo_type,old_current_decomp_type

INTEGER, SAVE ::                                                  &
! variables defining send and receive maps to be passed to
! GC_RALL_TO_ALL, defining the data transposition
  receive_map(7,2)                                                &
, n_sends,n_recvs  ! number of sends and receives
INTEGER, SAVE, ALLOCATABLE ::   send_map(:,:)


LOGICAL :: wrap          ! if the subdomain wraps over 0 degree meridion
LOGICAL :: wrap_special  ! if there is a wrap around, which starts and
                         ! ends on the same processor
LOGICAL :: zonal_data    ! if this is a zonal data grid
LOGICAL :: fullfield     ! if this is a full field - NOT a subarea


! Set all the old_* variables to a number indicating they've
! not been used yet

DATA                                                              &
  old_local_size , old_global_size                                &
, old_global_north , old_global_east_in                           &
, old_global_south , old_global_west                              &
, old_gridtype_code , old_scatter_pe                              &
, old_halo_type,old_current_decomp_type                           &
  / -1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /

! Functions

INTEGER :: get_fld_type

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! ------------------------------------------------------------------

IF (lhook) CALL dr_hook('STASH_SCATTER_FIELD',zhook_in,zhook_handle)

IF (.NOT.ALLOCATED(send_map)) &
  ALLOCATE(send_map(7,2*nproc_max))

! DEPENDS ON: get_fld_type
fld_type=get_fld_type(gridtype_code)

l_vec=((fld_type  ==  fld_type_u) .OR.                            &
       (fld_type  ==  fld_type_v))

! See if there is wrap around over meridion, and if so make
! sure that GLOBAL_EAST is > glsize(1)

global_east=global_east_in
IF (fld_type  ==  fld_type_unknown) THEN
  WRITE(6,*)                                                      &
  'STASH_SCATTER_FIELD cannot process field with ppx gridtype ',  &
  gridtype_code
  icode=1
  cmessage='STASH_SCATTER_FIELD : Incompatible GRIDTYPE code'
  GO TO 9999
END IF

IF (global_east  >   glsize(1,fld_type)) THEN
  wrap=.TRUE.
ELSE IF (global_east  <   global_west) THEN
  wrap=.TRUE.
  global_east=global_east_in+glsize(1,fld_type)
ELSE
  wrap=.FALSE.
END IF


IF ((gridtype_code  ==  ppx_atm_tzonal) .OR.                      &
    (gridtype_code  ==  ppx_ocn_uzonal) .OR.                      &
    (gridtype_code  ==  ppx_ocn_tzonal)) THEN

  zonal_data=.TRUE.
  global_x=1
ELSE

! This is a normal field

  zonal_data=.FALSE.
  global_x=glsize(1,fld_type)

END IF

global_y=glsize(2,fld_type)

! Set up logical indicating if this is a full field, or just
! a subdomain

IF (zonal_data) THEN

  fullfield= ( ( global_north  ==  global_y) .AND.                &
             ( global_south  ==  1))

ELSE

  fullfield = (( global_west  ==  1) .AND.                        &
               ( global_east  ==  global_x) .AND.                 &
               ( global_north  ==  global_y) .AND.                &
               ( global_south  ==  1))

END IF

! Dealing with fields not in model grid

IF((global_x == 0).OR.(global_x == imdi)) THEN
  DO level=1,levels
    DO i=1,global_size
      local_field(i,level)=global_field(i,level)
    END DO
  END DO
ELSE

! If this is a fullfield, we can simply use the standard
! SCATTER_FIELD routine

IF (fullfield) THEN

  IF (zonal_data) THEN

! DEPENDS ON: scatter_zonal_field
    CALL scatter_zonal_field( local_field,global_field,           &
                              lasize(2,fld_type,halo_type),       &
                              global_y,                           &
                              levels,gridtype_code,fld_type,      &
                              halo_type,                          &
                              scatter_pe)

  ELSE

    DO level=1,levels

! DEPENDS ON: scatter_field
      CALL scatter_field( local_field(1,level) ,                  &
                          global_field(1,level),                  &
                          lasize(1,fld_type,halo_type),           &
                          lasize(2,fld_type,halo_type),           &
                          global_x,global_y,                      &
                          fld_type,halo_type,                     &
                          scatter_pe,gc_all_proc_group,           &
                          icode,cmessage)

      IF (icode  /=  0) THEN
        WRITE(6,*) 'STASH_SCATTER_FIELD : SCATTER_FIELD failed'
        WRITE(6,*) 'Return code was : ',icode
        WRITE(6,*) 'Error message was : ',cmessage
        cmessage='SCATTER_FIELD failed'
        GO TO 9999
      END IF

    END DO
   END IF
 ELSE

! Check for WAM Wave fields - no code exists now

  IF(fld_type  ==  fld_type_comp_wave .OR.                        &
     fld_type  ==  fld_type_full_wave) THEN
    WRITE(6,*)                                                    &
    'STASH_SCATTER_FIELD: no code exists to generate sub-domains',&
    ' for the wam wave model'
    WRITE(6,*) 'Unable to process this field.'
    cmessage='STASH_SCATTER_FIELD:  Unable to process this field.'
    icode=100
    GO TO 9999
  END IF

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
    (halo_type  ==  old_halo_type) .AND.                          &
    (scatter_pe  ==  old_scatter_pe) .AND.                        &
    (current_decomp_type  ==  old_current_decomp_type ))) THEN

    old_local_size=local_size
    old_global_size=global_size
    old_global_north=global_north
    old_global_east_in=global_east_in
    old_global_south=global_south
    old_global_west=global_west
    old_gridtype_code=gridtype_code
    old_halo_type=halo_type
    old_scatter_pe=scatter_pe
    old_current_decomp_type=current_decomp_type

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

        local_start_row=1
        nrows_to_send=local_yend-local_ystart+1

        global_start_row=g_datastart(2,procid)+local_ystart -     &
                         halosize(2,halo_type) -                  &
                         global_south
        global_row_length=global_east-global_west+1

! Calculate the following variables:
! local_row_length : the X dimension size of the local array
! local_send_offx  : the offset into each row to start sending from
! sendsize_x       : the number of points on each row to send
! The calculation of these numbers is different for processors
! at the start and end of a wrap_special case

        IF (wrap_special .AND. procx  ==  proc_topleft_x) THEN
          local_row_length=g_blsize(1,fld_type,procid) +          &
                           local_xend - local_xstart + 1
          local_start_col=1
          sendsize_x=g_lasize(1,fld_type,halo_type,procid) -      &
                     local_xstart
          global_start_col=1

        ELSE IF (wrap_special .AND. procx  ==  proc_botright_x)    &
        THEN
          local_row_length=g_blsize(1,fld_type,procid) +          &
                           local_xend - local_xstart + 1
          local_start_col=local_row_length - local_xend +         &
                          halosize(1,halo_type) + 1
          sendsize_x=local_xend - halosize(1,halo_type)
          global_start_col=global_row_length-sendsize_x+1

        ELSE
          local_row_length=local_xend-local_xstart+1
          local_start_col=1
          sendsize_x=local_xend-local_xstart+1
          global_start_col=local_xstart -                         &
                           (halosize(1,halo_type) + 1 ) +         &
                           g_datastart(1,procid)-global_west+1
        END IF

        IF (global_start_col  <   0) THEN
! Wrapped around field, but this processor is not start or end
! processor
          global_start_col=global_start_col+glsize(1,fld_type)
        END IF

! Now we can set up the send and receive map entries for the data on
! this processor

        IF (mype  ==  procid) THEN  ! I need to receive some data

            n_recvs=n_recvs+1

          receive_map(r_source_pe,n_recvs) = scatter_pe
          receive_map(r_base_address_in_recv_array,n_recvs) =     &
            (local_start_row-1)*local_row_length +                &
            local_start_col
          receive_map(r_number_of_elements_in_item,n_recvs) =     &
            nrows_to_send
          receive_map(r_stride_in_recv_array,n_recvs) =           &
            local_row_length
          receive_map(r_element_length,n_recvs) = sendsize_x
          receive_map(r_base_address_in_send_array,n_recvs) =     &
            (global_start_row-1)*global_row_length +              &
            global_start_col
          receive_map(r_stride_in_send_array,n_recvs) =           &
            global_row_length

        END IF ! if I'm receiving data

        IF (mype  ==  scatter_pe) THEN ! I need to send data

          n_sends=n_sends+1

          send_map(s_destination_pe,n_sends) = procid
          send_map(s_base_address_in_send_array,n_sends) =        &
            (global_start_row-1)*global_row_length +              &
            global_start_col
          send_map(s_number_of_elements_in_item,n_sends) =        &
            nrows_to_send
          send_map(s_stride_in_send_array,n_sends) =              &
            global_row_length
          send_map(s_element_length,n_sends) = sendsize_x
          send_map(s_base_address_in_recv_array,n_sends) =        &
            (local_start_row-1)*local_row_length +                &
            local_start_col
          send_map(s_stride_in_recv_array,n_sends) =              &
            local_row_length

        END IF ! if I'm sending data

      END DO ! procx : loop along processor row

    END DO ! procy : loop down processor column

  END IF ! if I need to recalculate my send/receive maps

! Send / receive the data using GCG_RALLTOALLE

  DO level=1,levels

    flag=0  ! This is currently ignored at GCG v1.1
    info=gc_none

    CALL gcg_ralltoalle(                                          &
      global_field(1,level)  ,                                    &
      send_map    , n_sends  ,global_size  ,                      &
      local_field(1,level) ,                                      &
      receive_map , n_recvs , local_size ,                        &
      gc_all_proc_group , flag, info)

  END DO

END IF ! if this is a full or extracted field

END IF

 9999 CONTINUE

IF (lhook) CALL dr_hook('STASH_SCATTER_FIELD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE stash_scatter_field

