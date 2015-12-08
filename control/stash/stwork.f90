! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine STWORK -------------------------------------------------
!
!    Purpose: Processes all the STASHlist entries for a particular
!   item and section number after a timestep.  Raw input diagnostics to
!   stwork come either from D1, the main data array, or stash_work, a
!   work array dimensioned by the control level and passed in.  Each
!   field is spatially processed, temporally processed and then output.
!   The output destination can either be to an address in D1 or to a PP
!   fieldsfile on a given unit.  In either case a PP-type LOOKUP header
!   is generated to describe the contents of the output field.  Now
!   only handles atmosphere diagnostics.

!   Method:
!   Input information describing all STASH processing requests for a
!   single model variable (STASH item,STASH section,internal model).
!   Data input is either D1 array (primary variables) or stwork array
!   (diagnostic variables).

! 0. Initialisation. Generate flags describing the grid location
!    for this STASH variable.
! Start loop over STASH requests for this STASH variable......
! 1.   Determine input and output lengths and levels, and relative
!      positions, dependent on type of processing:
!      simple extraction, spatial, time-series.
!      Determine which type of processing required.
! 2.   If spatial processing:
! 2.1  Timeseries. Extracts data and processes via multi_spatial
!      routine.
! 2.2  Standard spatial processing: loop over output levels (unless
!      multi-level processing required) and process via spatial.
! 3.   Not spatial processing.
! 4.   Output.
! 4.1  Determine output unit, PP lookup location and output lengths.
! 4.2  Pack output data if requested, via:
!      pp_file - coex (WGDOS packing), or
!      grib_file - grib packing.
!      Set up PP lookup values in pp_head
!      Write out data to PP file.
! 4.3  Write PP lookup field to PP file.
! 4.4  Output to D1 array (either to model dump or secondary space).
!      Time-processing (eg time means, accumulations) through temporal
!      routine if requested.
! End of loop over STASH requests..................
! 9.   Error exits.

!
!    Programming standard: UM Doc Paper 3
!
!    External documentation : UMDP no C4
!

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

SUBROUTINE stwork (                                                  &
     d1,len_tot,stash_work,stash_work_len,lenout,                    &
     global_lenout,                                                  &
     is,im,ilstart,ilend,step,steps_per_period,secs_per_period,      &
     previous_time,                                                  &
     stlist,len_stlist,totitems,si,nsects,nitems,                    &
     stash_levels,num_stash_levels,num_level_lists,                  &
     stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,          &
     max_stash_levs,sttabl,nsttims,                                  &
     nsttabl,stash_series,stash_series_len,                          &
     stash_series_rec_len,stash_series_index,stash_ser_index_size,   &
     pp_pack_code,model_ft_unit,ft_steps,ft_firststep,               &
     fixhd,inthd,realhd,len_fixhd,len_inthd,len_realhd,              &
     levdepc,len1_levdepc,len2_levdepc,                              &
     lookup,rlookup,len1_lookup,len2_lookup,                         &
     pp_len2_lookup,nunits,pp_len2_look,                             &
     lcyclic,                                                   &
     t_rows,u_rows,row_length,t_field,u_field,t_levels,              &
     river_rows, river_row_length,                                   &
     fcst_prd,run_indic_op,elf,ft_lastfield,                         &
     sm_ident,im_ident,dump_pack,                                    &
     stsuparrlen, stsuparr, istsuparr, sa_idx, sa_idxlen,            &
     icode,cmessage)

  USE yomhook, ONLY: lhook, dr_hook

  USE parkind1, ONLY: jprb, jpim

  USE filenamelength_mod, ONLY: filenamelength

  USE io, ONLY : is_unit_open, ioFileTypeUM

  USE model_file, ONLY :          &
     get_file_address,            &
     mf_data_address,             &
     get_mf_information,          &
     mf_num_preallocated_headers, &
     attachlookups,               &
     mf_data_missing,             &
     model_file_open,             &
     model_file_managed

  USE ios_stash

  USE stwork_aux

  USE ios_model_geometry, ONLY :   &
     maxfielddomain

  USE io_configuration_mod, ONLY : &
     io_field_padding

  USE missing_data_mod, ONLY: rmdi, imdi

  USE ereport_mod, ONLY : ereport
  USE um_parvars
  USE decomp_db
  USE sterr_mod, ONLY: st_upper_less_lower, st_not_supported,       &
                       st_no_data,st_nd, st_bad_array_param,        &
                       st_bad_address, st_unknown,                  &
                       st_bad_wraparound, st_illegal_weight,        &
                       unknown_weight, unknown_mask,                &
                       unknown_processing, nonsense
  USE lookup_addresses
  USE Submodel_Mod

  USE cppxref_mod

  IMPLICIT NONE

  INTEGER                            &
     totitems,                       &!IN   Max no of items in STASHlist
     nsects,                         &!IN   Max no of sections
     nitems,                         &!IN   Max no of items in a section
     len_tot,                        &!IN   Length of real data array D1
     len_fixhd,                      &!IN   Length of fixed constants
     len_inthd,                      &!IN   Length of integer constants
     len_realhd,                     &!IN   Length of real constants
     len1_levdepc,                   &!IN   First dimension of levdepc
     len2_levdepc,                   &!IN   Second dimension of levdepc
     num_stash_levels,               &!IN   Dimension of stash levels
     num_level_lists,                &!IN   Dimension of stash levels
     num_stash_pseudo,               &!IN   Maximum num pseudo-levels in a list
     num_pseudo_lists,               &!IN   Number of pseudo-level lists
     max_stash_levs,                 &!IN   Max no of output levels for any diag
     len1_lookup,                    &!IN   First dimension of lookup/ipplook
     len2_lookup,                    &!IN   Second dimension of lookup
     pp_len2_lookup,                 &!IN   Largest poss. value in pp_len2_look
     nunits,                         &!IN   Max i/o ft unit no
     pp_len2_look(20:nunits),        &!IN   Individual pp_len2_looks per unit
     pp_pack_code(20:nunits),        &!IN   Packing code per unit
     ft_lastfield(20:nunits),        &!IN   Current write posn in each pp file
     ft_steps(20:nunits),            &!IN   File reinitialisation freq per unit
     ft_firststep(20:nunits),        &!IN   First step file initialised
     nsttims,                        &!IN   Number of times against to test
     nsttabl,                        &!IN   Number of STASH timetables
     num_words,                      &!IN   Number of 64 Bit words to hold DATA
     sm_ident,                       &!IN   Submodel identifier
     im_ident,                       &!IN   Internal model identifier
     dump_pack,                      &!IN   Packing Indicator for Dump
     sa_idxlen,                      &!IN   Superarray index length
     sa_idx(sa_idxlen),              &!IN   Superarray index
     stsuparrlen,                    &!IN   Superarray index length
     istsuparr(stsuparrlen)           !IN   Integer superarray

  CHARACTER(LEN=filenamelength) ::   &
     model_ft_unit(nunits)            !IN   Current table of file associations

  INTEGER, INTENT(IN) ::             &
     fixhd(len_fixhd),               &!IN   Array of fixed constants
     inthd(len_inthd),               &!IN   Array of integer constants
     ilstart,                        &!IN   Start of loop over entries
     ilend,                          &!IN   End of loop over entries
     is,                             &!IN   Section numbers
     im,                             &!IN   Item number
     step,                           &!IN   Model step number
     steps_per_period,               &!IN   No of steps in defining period
     secs_per_period,                &!IN   No of secs in period
                                      !     (define timestep)
     previous_time(7)                 !IN   Time at start of current step

  INTEGER                               &
     lookup(len1_lookup,len2_lookup),   &! Integer lookup headers
     rlookup(len1_lookup,len2_lookup),  &! Real version of lookup
     icode                               !OUT   Return code from routine

  INTEGER, INTENT(IN) :: &
     lenout,            &!IN     Length of largest workfield needed
     global_lenout,     &!IN     Output length of largest field
     t_field,           &!IN     No of temp/press points
     u_field,           &!IN     No of u,v points
     row_length,        &!IN     No of points per row
     u_rows,            &!IN     No of u,v rows
     t_rows,            &!IN     No of press/temp rows
     t_levels,          &!IN     No of model press/temp levels
     stash_work_len,    &!IN     Length of stash_work
     len_stlist,        &!IN     No of entries in stashlist
     stlist(len_stlist,totitems),                        &!IN stashlist
     si(nitems,0:nsects,n_internal_model),               &!IN STASH in address
     sttabl(nsttims,nsttabl),                            &!IN  STASH time tables
     stash_levels(num_stash_levels+1,num_level_lists),   &
     stash_pseudo_levels(num_stash_pseudo+1,num_pseudo_lists), &
     run_indic_op,           &! IN     Operational Run indicator (ITAB)
! STASH timeseries information
     stash_series_len,       &! IN no of STASH series records
     stash_series_rec_len,   &! IN length of each record
     stash_series(stash_series_rec_len,stash_series_len), &
! IN individual sample records
     stash_ser_index_size,   &! IN no. of index blocks
     stash_series_index(2,stash_ser_index_size)
! IN index block (1=start, 2=no of records)

  INTEGER                    &
     exppxi,                 &! Function to extract ppxref info
     im_index,               &! Internal model index number
     n1                       ! Packing Indicator for Lookup(21)

  REAL fcst_prd  !IN     Forecast period

  INTEGER, INTENT(IN)  :: river_rows         ! River routeing
  INTEGER, INTENT(IN)  :: river_row_length   ! dimensions

  CHARACTER(LEN=80) cmessage ! OUT MESSAGE FROM ROUTINE


  LOGICAL              &
     lcyclic,         &!IN TRUE if cyclic EW BCs
     elf,             &!IN TRUE if the input grid is rotated Equatorial
     packing_hold

  REAL                                     &
     d1(len_tot),                          &!IN  Real data array
     realhd(len_realhd),                   &!IN  Array of real constants
     levdepc(len1_levdepc*len2_levdepc+1), &!IN  Level dep constants
     stash_work(stash_work_len),           &!IN  Input work array to STASH
     stsuparr(stsuparrlen)                  !IN  Real superarray

! Common blocks and PARAMETERs

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
! Start i_stgfld

! Description:
!   This file contains an interface to STASH_GATHER_FIELD and
!   must be included whenever this routines is used so as to
!   get declarations of optional arguments correct.
!
      INTERFACE
        SUBROUTINE STASH_GATHER_FIELD (                                 &
     &    LOCAL_FIELD , GLOBAL_FIELD ,                                  &
     &    LOCAL_SIZE, GLOBAL_SIZE, LEVELS,                              &
     &    GLOBAL_NORTH , GLOBAL_EAST_IN , GLOBAL_SOUTH , GLOBAL_WEST,   &
     &    GRIDTYPE_CODE ,HALO_TYPE,                                     &
     &    GATHER_PE,                                                    &
     &    DATA_EXTRACTED,                                               &
     &    PACKING, IM_IDENT, PACKING_TYPE,                        &
     &    NUM_OUT,                                                      &
     &    COMP_ACCRCY, loc_RMDI,                                        &
     &    ICODE, CMESSAGE)

        INTEGER, INTENT(IN) ::                                          &
     &    LOCAL_SIZE                                                    &
                          ! IN: size of level of LOCAL_FIELD
     &  , GLOBAL_SIZE                                                   &
                          ! IN: size of level of GLOBAL_FIELD
     &  , LEVELS                                                        &
                          ! IN: number of levels
     &  , GLOBAL_NORTH                                                  &
                          ! IN: specification of subdomain boundaries
     &  , GLOBAL_EAST_IN                                                &
                          ! IN: ""
     &  , GLOBAL_SOUTH                                                  &
                          ! IN: ""
     &  , GLOBAL_WEST                                                   &
                          ! IN: ""
     &  , GRIDTYPE_CODE                                                 &
                          ! IN: indicates the type of grid output
     &  , HALO_TYPE                                                     &
                          ! IN: type of halo on this field
     &  , GATHER_PE       ! IN: the PE to gather the global field to

        INTEGER, INTENT(OUT) ::                                         &
     &    ICODE           ! OUT: return code, 0=OK
!
! Optional Arguments to handle the COEX packing if necessary
!
        LOGICAL, INTENT(IN), OPTIONAL ::   &
     &    PACKING
                          ! IN: Set .true. if packing of the input
                          !     field is to be packed

        INTEGER, INTENT(IN), OPTIONAL ::                                &
     &    IM_IDENT        ! IN: Internal model identifier

        INTEGER, INTENT(INOUT), OPTIONAL ::                             &
     &    PACKING_TYPE    ! IN/OUT: This flag is zero on input,
                          !         then stash packing is selected,
                          !         and the routine returns the
                          !         packing flag.
                          !
                          !         If the variable is set to 1 on input
                          !         then 32-bit packing for dumpfiles
                          !         is selected

        INTEGER, INTENT(OUT), OPTIONAL ::                               &
     &    NUM_OUT         ! OUT: Number of 32-bit IBM words in the
                          !      Packed field for WDGOS packing

        INTEGER, INTENT(IN), OPTIONAL ::                                &
     &    COMP_ACCRCY     ! IN: Packing Accuracy in Power of 2

        REAL, INTENT(IN), OPTIONAL ::                                   &
     &    loc_RMDI        ! IN: Missing data indicator
!
! Remaining Non-Optional Arguments
!
        LOGICAL, INTENT(IN) ::                                          &
     &    DATA_EXTRACTED  ! IN: TRUE if the data in LOCAL_FIELD has
                          !     already been extracted, or FALSE if
                          !     the extraction must be done here.

        REAL, INTENT(IN) ::                                             &
     &    LOCAL_FIELD(LOCAL_SIZE,LEVELS)
                          ! IN : local data

        REAL, INTENT(OUT) ::                                            &
     &    GLOBAL_FIELD(GLOBAL_SIZE,LEVELS)
                          ! OUT : (PE GATHER_PE only) - full gathered
                          !       field

        CHARACTER(LEN=*), INTENT(OUT) ::                                   &
     &    CMESSAGE        ! OUT: Error message if ICODE .NE. 0

        END SUBROUTINE STASH_GATHER_FIELD
      END INTERFACE
! End i_stgfld
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

  INTEGER extra_var_data ! Size of extra grid data
 
! External function:
  INTEGER get_fld_type
  
! Local variables

  REAL                      &
     ppfield(lenout),       &! Main internal work array
     level(max_stash_levs), &! The levels of the data as real nos
     sample_prd,            &! Sampling period in hours for means, etc
     a_io                    ! The output code from the unit command.

! I/O buffer workspace arrays - used for holding the data to be
! written to disk.

  REAL :: buf(global_lenout)

  INTEGER, POINTER :: ipplook(:,:)=>NULL()
  INTEGER                        &
     n_rows_out,                 &! No of rows used for a diagnostic
     n_cols_out,                 &! No of cols used pphoriz=n_rows*n_cols
     srow_in,srow_out,           &! North, South, East and West
     nrow_in,nrow_out,           &! subdomain limits in the horizontal sense
     wcol_in,wcol_out,           &! corresponding to the subarea before being
     ecol_in,ecol_out,           &! processed (IN) and after processing (OUT)
     lev_in_addr,                &! The num of pts skipped in the Input
     gr,                         &! Grid point code
     lbproc_comp(14),            &! Array of 1/0 denoting lbproc components
     ppunit,                     &! PPinit number (also used in pp_file)
     lenbuf,                     &! PPhoriz_out rnd to 512 words (used PP)
     comp_accrcy,                &! Packing accuracy in power of 2
     pphoriz_out,                &! No of points in the output field
     pphoriz_in,                 &! No of points in the input field
     iwa,                        &! Record location
     len_buf_words,              &! Number of 64 Bit words (rounded to 512)
     num_levs_out,               &! Number of output levels
     num_levs_in,                &! Number of input levels
     this_index_lev,             &! index of level in output array for multi-
                                  ! level processing in SPATIAL. Note:
                                  ! loop over output levels:-
                                  ! this_index_level=1,num_levs_out
                                  ! level_list(this_ ) is model level or
                                  ! pressure level
                                  ! index_lev(this_ ) is index to level in
                                  ! input array
     index_lev(max_stash_levs),  &! Index used to relate input and
                                  ! output levels
     level_list(max_stash_levs), &! model level for each output level
     pseudo_level(max_stash_levs),&! pseudo-level at each output level
     lv,                         &! LV code for field from STASHmaster
     samples,                    &! no of samples (timeseries/trajec)
     icurrll_dump_ptr,           &! pointer to mother record lookup
     start_time(7),              &! start time for PPheader in
                                  ! timeseries/time mean/accumulation
     no_records,                 &! no of records processed by multi_spatial
     record_start,               &! the start record for multi_spatial
     jp_1,                       &! Jpointer to first level of p in D1
     jpstar,                     &! Jpointer to pstar in D1
     len_field,                  &! Holds the amount of data available
                                  ! to pp_file
     packing_type_hold,          &! Holds the packing type aftr packing by
                                  ! (STASH_)gather_field
     num_out                      ! Number of 32-bit words from coex

  REAL                    &
     dummy3d(1,1,1), dummy2d(1,1)

  INTEGER :: end_time(7)
  INTEGER             &
     jl,ii,il,jj,it,  &
     ntab,            &! Number of the STASH times table
     ikount,          &! Local  counter
     points,          &! No of points in a field
     kl,              &! Local  counter
     ml,              &! Local  counter
     len_io,          &! The length of data transferred.
     ilprev,          &! The counter of the first of a pair of stlist
     ilcurr,          &! The current value of il
     lbvcl,           &! Vertical coordinate code in lookup table
     icurrll,         &! Current position in the PP Lookup table
     i,               &! loop variable
     j,               &! Level indicator used in call to grib_file
     packing_type,    &! 0 No packing, 1 for wgdos, 3 for grib
     grib_packing      ! ppxref profile number used to determine
                       ! grib packing method

  INTEGER vx,vy,vz,   &! Sizes of arrays.
     st_grid,         &! Horizontal grid type, (eg. T-p or u-v)
     len_ppname
  INTEGER input_code   ! Value of input code in stashlist
  INTEGER addr         ! Address of STASH variable in either dump or
  INTEGER addr_out     ! Address of spatially processed field
  INTEGER elap_time    ! No of timesteps elapsed in period.
  INTEGER series_ptr   ! The address in stash_series where domin inf
  INTEGER index_size   ! The number of levels in the index
  INTEGER base_level   ! Base model level needed for mass weighting
  INTEGER base_level0  ! Ref base levels in levs loop
  INTEGER what_proc    ! What kind of processing will be done
  INTEGER what_mean    ! What kind of meaning will be done
  INTEGER output_code  ! Output destination code from stlist
  INTEGER expected_len ! expected length of output field
  LOGICAL             &
     s_f,             &! TRUE for items requiring processing
     llproc,          &! TRUE if spatial processing required
     lnullproc,       &! TRUE if null processing indicated
     lfullfield,      &! TRUE if output field on full horiz domain
     lmasswt,         &! TRUE if level-by-level mass-weights exist
     start_step,      &! TRUE at start of time period (temporal)
     packing,         &! TRUE if packing required
     grib_out,        &! TRUE if output to be in grib code.
     rotate            ! TRUE if input data to be rotated
  LOGICAL :: found_end

  INTEGER expected_extra ! expected length of extra data
  INTEGER extraw         ! number of extra words this timestep
  INTEGER extraw_hdr     ! number of extra words for the header
  INTEGER data_type_code ! ppx_data_type code copied from PPX file
  INTEGER rotatecode     ! code for rotated grid
  INTEGER nt_dim         ! Number of tracers

  INTEGER              & ! local versions of the subdomain limits
     local_nrow_out,local_srow_out,local_ecol_out,local_wcol_out,&
     local_nrow_in,local_srow_in,local_ecol_in,local_wcol_in,    &
! global versions of the total size of output
     global_n_rows_out,global_n_cols_out, global_pphoriz_out,    &
     info                ! return variable from gcom
  INTEGER, PARAMETER  :: pe_zero = 0

  INTEGER                  &
     grid_type_code,       &! grid type of field being processed
     fld_type,             &! field type: u-,v- or p- location on C grid
     no_of_levels_masswt,  &! no. levels for mass weights array
     halo_type              ! halo type of the field being processed

  INTEGER :: ios_q_slot ! Opaque handle for referencing async operations
  INTEGER :: ios_packing_flag
  INTEGER :: ios_fullfield_flag
  INTEGER             ::  istat  ! gcom return flag

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('STWORK',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!  0. Initialise: set constants relating to input grid type and size
!
!  0.1  Set up internal grid type st_grid and input field size
!       according to the master GR code for the diagnostic
!
  icurrll_dump_ptr=-1

! Get internal model index
  im_index = internal_model_index(im_ident)
  nt_dim = (sa_idx(2) - sa_idx(1))/2

! [care needed! stsuparr passed as real but integer value here.]
  jp_1   = stsuparr(sa_idx(7))
  jpstar = stsuparr(sa_idx(8))


! Get STASHmaster gridtype code
! DEPENDS ON: exppxi
  gr = exppxi( im_ident, is, im, ppx_grid_type,                     &
  icode, cmessage)
  grid_type_code=gr

! Get field type, ie u or v or p location in Arakawa C grid staggering
! DEPENDS ON: get_fld_type
  fld_type=get_fld_type(grid_type_code) ! Function

! Find out what halo_type (and hence halo width) this field has

! DEPENDS ON: exppxi
  halo_type=exppxi( im_ident, is, im, ppx_halo_type,                &
  icode, cmessage)

! Determine the input length pphoriz_in and stash grid staggering for
! this field.

  IF (gr == ppx_atm_tall.OR.gr == ppx_atm_tland.OR.                 &
     gr == ppx_atm_tsea) THEN
! Atmosphere data on T-grid
    st_grid=st_tp_grid
    pphoriz_in=t_field
  ELSE IF (gr == ppx_atm_uall.OR.gr == ppx_atm_uland.OR.            &
     gr == ppx_atm_usea) THEN
! Atmosphere data on U-grid
    st_grid=st_uv_grid
    pphoriz_in=u_field
  ELSE IF (gr == ppx_atm_compressed) THEN
! Atmosphere data on T-grid (compressed)
    st_grid=st_tp_grid
    pphoriz_in=t_field
  ELSE IF (gr == ppx_atm_cuall) THEN
! Atmosphere data on C-grid (u-points)
    st_grid=st_cu_grid
    pphoriz_in=t_field
  ELSE IF (gr == ppx_atm_cvall) THEN
! Atmosphere data on C-grid (v-points)
    st_grid=st_cv_grid
    pphoriz_in=u_field
  ELSE IF (gr == ppx_atm_tzonal) THEN
! Atmosphere zonal data on T-grid
    st_grid=st_zt_grid
    pphoriz_in=t_rows
  ELSE IF (gr == ppx_atm_uzonal) THEN
! Atmosphere zonal data on u-grid
    st_grid=st_zu_grid
    pphoriz_in=u_rows
  ELSE IF (gr == ppx_atm_tmerid) THEN
! Atmosphere meridional data on T-grid
    st_grid=st_mt_grid
    pphoriz_in=row_length
  ELSE IF (gr == ppx_atm_umerid) THEN
! Atmosphere meridional data on u-grid
    st_grid=st_mu_grid
    pphoriz_in=row_length
  ELSE IF (gr == ppx_atm_scalar) THEN
! Atmosphere scalar
    st_grid=st_scalar
    pphoriz_in=1
  ELSE IF (gr == ppx_atm_river) THEN
! Atmosphere river routing
    st_grid=st_riv_grid
    pphoriz_in=river_rows * river_row_length
  ELSE
! Unknown grid type
    icode=1
    cmessage='STWORK   : Unknown grid type found in STASHmaster'
    GO TO 9999
  END IF
! The input length pphoriz_in has been calculated explicitly depending
! on grid_type where there is an implicit assumption that diagnostic
! fields have no halos.

! For input fields in section 0, ie primary fields held in D1, this is
! not necessarily the case and input lengths are extracted using an
! addressing service routine dependent on grid and halo codes:
  IF(is == 0 .OR. is == 33 .OR. is == 34) THEN ! section 0,33 or 34
! DEPENDS ON: addrln
    CALL addrln(grid_type_code,halo_type,pphoriz_in,local_data)
  END IF

!
!  0.2 Set up ROTATE to flag fields which are rotated (eg. ELF winds)
!      (this is used to set alternative fieldcodes in pphead)
!

! DEPENDS ON: exppxi
  rotatecode = exppxi( im_ident, is, im, ppx_rotate_code,           &
  icode, cmessage)

  IF (rotatecode  ==  ppx_elf_rotated .AND. elf)                    &
     THEN
    rotate=.TRUE.
  ELSE
    rotate=.FALSE.
  END IF


!  1. Loop over entries with this section/item number

  DO il=ilstart,ilend  !loop over num entries for each item/sec
    extraw=0 ! no extra data by default


!  1.1 Set up S_F which has to be set for each stashlist entry. The
!      stashflag is be set for a particular item/section the s_f for
!      the stashlist entry.

    s_f=.FALSE.
    IF (stlist(st_freq_code,il) == 1.AND.                            &
       step >= stlist(st_start_time_code,il).AND.                    &
       (step <= stlist(st_end_time_code,il).OR.                      &
       stlist(st_end_time_code,il) == st_infinite_time)) THEN
!     ... if required every step between start and end
      s_f=.TRUE.
    ELSE IF (stlist(st_freq_code,il) <  0) THEN
!     ... if required at specified times and this is one of them
      ntab=-stlist(st_freq_code,il)
      found_end=.FALSE.
      it=1
      DO WHILE (it <= nsttims .AND. .NOT. found_end)
        IF (sttabl(it,ntab) == st_end_of_list) THEN
          found_end=.TRUE.
        ELSE IF (step == sttabl(it,ntab)) THEN
          s_f=.TRUE.
        END IF
        it=it+1
      END DO
    ELSE IF (stlist(st_freq_code,il) >  0) THEN
      IF   (MOD((step-stlist(st_start_time_code,il)),                &
         stlist(st_freq_code,il)) == 0.AND.                          &
         step >= stlist(st_start_time_code,il).AND.                  &
         (step <= stlist(st_end_time_code,il).OR.                    &
         stlist(st_end_time_code,il) == st_infinite_time))           &
!        ... if required every N timesteps and this is one of them
         s_f=.TRUE.
    END IF

!  S_F now set - Start of IF (s_f) block .......
    IF (s_f) THEN

!  1.2 Find number of input and output levels and relative positions
!      and set up levels and pseudo-levels arrays for PPheaders.
!      Set indicator lmasswt if level-by-level mass weighting possible
!      - only currently available with atmosphere model full levels.

! Special case of mean timeseries leave ilcurr pointing to il

      ilcurr=il   ! The current stashlist entry il
      IF (stlist(st_input_code,il) <  0.AND.                          &
         stlist(st_proc_no_code,il) /= st_time_series_mean) THEN
        ilcurr=-stlist(st_input_code,il) ! points to prev entry
      END IF

! Get STASHmaster lbvc code
! DEPENDS ON: exppxi
      lbvcl = exppxi( im_ident, is, im, ppx_lbvc_code,                &
      icode, cmessage)

! Get STASHmaster lv code
! DEPENDS ON: exppxi
      lv = exppxi( im_ident, is, im, ppx_lv_code,                     &
      icode, cmessage)

! DEPENDS ON: stlevels
      CALL stlevels(stlist(1,ilcurr),len_stlist,                     &
         stash_levels,num_stash_levels,num_level_lists,              &
         stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,      &
         max_stash_levs,num_levs_in,num_levs_out,index_size,         &
         index_lev,level_list,lbvcl,level,pseudo_level,              &
         icode,cmessage)
      IF (icode >  0) GO TO 9999
      vz=num_levs_in

! Set switch to indicate whether level by level mass weights are
! available ( whether or not mass weighting requested)
      no_of_levels_masswt=1 ! default no. of levels for pre-
! calculating mass weights array
      IF (lv == ppx_full_level) THEN
        lmasswt=.TRUE.
! Level by level mass weights available and mass weighting required:
! Extra space only needed if mass weighting requested.
        IF(stlist(st_weight_code,il) == stash_weight_mass_code) THEN
          no_of_levels_masswt = index_size  ! req. no. model levels
        END IF
      ELSE
        lmasswt=.FALSE.
      END IF

!  1.3 Find the horizontal dimensions for the output grid
!      and the input field subdomain limits in which processing happens.

      wcol_in= stlist(st_west_code,ilcurr)    ! Input subdomain limits
      ecol_in= stlist(st_east_code,ilcurr)
      nrow_in= stlist(st_north_code,ilcurr)
      srow_in= stlist(st_south_code,ilcurr)

!  1.3.2 Other output types need to calculate lengths in detail
!        (ie. number of rows, columns, and horizontal field size)
!        according to processing options

! Calculate local versions of the subdomain limits and area

! DEPENDS ON: global_to_local_subdomain
      CALL global_to_local_subdomain(                                &
         .TRUE.,.TRUE.,                                              &
         grid_type_code,halo_type,mype,                              &
         srow_in,ecol_in,nrow_in,wcol_in,                            &
         local_srow_in,local_ecol_in,                                &
         local_nrow_in,local_wcol_in)

      what_proc=stlist(st_proc_no_code,ilcurr)
      what_mean=(stlist(st_gridpoint_code,ilcurr)/block_size)*       &
         block_size
      samples=0         ! Initialise value for non-timeseries

!  1.3.2.1 Time series or trajectory processing

      IF (what_proc == st_time_series_code.OR.                       &
         what_proc == st_append_traj_code.OR.                        &
         what_proc == st_time_series_mean) THEN

!  1.3.2.2 Compute number of samples in period for timeseries for
!          input to PP_HEAD.
!          No of output rows and cols are set to the no of points
!          in each time sample and number of time samples in the
!          period spanned by the output field, respectively.

        samples=stlist(st_period_code,ilcurr)/                       &
           stlist(st_freq_code,ilcurr)
        wcol_out=1
        ecol_out=samples
        nrow_out=stlist(st_output_length,ilcurr)/samples
        srow_out=1

        local_wcol_out=1
        local_ecol_out=samples
        local_nrow_out=stlist(st_output_length,ilcurr)/samples
        local_srow_out=1

!  1.3.2.3 Multi spatial processing of some other type (not supported)
      ELSE IF (stlist(st_series_ptr,ilcurr) /= 0) THEN
        icode=1323
        cmessage='STWORK   : Illegal timeseries processing selected'
        GO TO 9999
!  1.3.2.4 Primary record requesting an extract
      ELSE IF (what_mean == extract_base) THEN
        wcol_out= wcol_in
        ecol_out= ecol_in
        nrow_out= nrow_in
        srow_out= srow_in

        local_wcol_out= local_wcol_in
        local_ecol_out= local_ecol_in
        local_nrow_out= local_nrow_in
        local_srow_out= local_srow_in

!  1.3.2.5 Primary record requesting a vertical mean
      ELSE IF (what_mean == vert_mean_base) THEN
        wcol_out= wcol_in
        ecol_out= ecol_in
        nrow_out= nrow_in
        srow_out= srow_in

        local_wcol_out= local_wcol_in
        local_ecol_out= local_ecol_in
        local_nrow_out= local_nrow_in
        local_srow_out= local_srow_in

!  1.3.2.6 Primary record requesting a zonal mean
      ELSE IF (what_mean == zonal_mean_base) THEN
        wcol_out= 1
        ecol_out= 1
        nrow_out= nrow_in
        srow_out= srow_in

        local_wcol_out= 1
        local_ecol_out= 1
        local_nrow_out= local_nrow_in
        local_srow_out= local_srow_in

!  1.3.2.7 Primary record requesting a meridional mean
      ELSE IF (what_mean == merid_mean_base) THEN
        wcol_out= wcol_in
        ecol_out= ecol_in
        nrow_out= 1
        srow_out= 1

        local_wcol_out= local_wcol_in
        local_ecol_out= local_ecol_in
        local_nrow_out= 1
        local_srow_out= 1

!  1.3.2.8 Primary record requesting a global mean
      ELSE IF (what_mean == global_mean_base) THEN
        wcol_out= 1
        ecol_out= 1
        nrow_out= 1
        srow_out= 1

        local_wcol_out= 1
        local_ecol_out= 1
        local_nrow_out= 1
        local_srow_out= 1

!  1.3.2.9 Primary record requesting a field mean
      ELSE IF (what_mean == field_mean_base) THEN
        wcol_out= 1
        ecol_out= 1
        nrow_out= 1
        srow_out= 1

        local_wcol_out= 1
        local_ecol_out= 1
        local_nrow_out= 1
        local_srow_out= 1

!  1.3.2.10 Error trap for unknown request
      ELSE         ! Invalid option
        icode=st_unknown
        WRITE(cmessage,'(A,A,1x,i5)')            &
           'STWORK: FATAL ERROR: ',             &
           'unknown option in setup',what_mean
        GO TO 9999 ! jump to error return
      END IF

!  1.3.3 Compute expected length. This differs from total output length
!        when data is appended from multiple timesteps into the same
!        field, being output_length/number_of_appends in this case.

      IF (stlist(st_output_code,il) >= 0.AND.                         &
         (what_proc == st_time_series_code.OR.                        &
         what_proc == st_append_traj_code.OR.                         &
         what_proc == st_time_series_mean)) THEN
        series_ptr=stlist(st_series_ptr,il) !set up ptr to stashseries
        expected_extra=(stash_series_index(2,series_ptr)+1)*6
        extraw_hdr=expected_extra
        expected_len=((stlist(st_output_length,ilcurr)                &
           -expected_extra)*                                          &
           stlist(st_freq_code,ilcurr))/stlist(st_period_code,ilcurr)
      ELSE
        expected_len=stlist(st_output_length,ilcurr)
        expected_extra=0 ! no extra data for non timeseries stuff
        extraw_hdr=0
      END IF

!  1.3.6 Compute number of rows and columns and field size for output
!        - first adjust easternmost column if field wraps EW

      IF (wcol_in  >  ecol_in .AND.lcyclic)                           &
         ecol_in =ecol_in + glsize(1,fld_type)

      IF (wcol_out >  ecol_out.AND.lcyclic)                           &
         ecol_out=ecol_out + glsize(1,fld_type)

      IF (local_wcol_out  >   local_ecol_out)                         &
         local_ecol_out=local_ecol_out+blsize(1,fld_type)

      n_rows_out = local_nrow_out - local_srow_out + 1
      n_cols_out = local_ecol_out - local_wcol_out + 1
      global_n_rows_out =  nrow_out - srow_out + 1
      global_n_cols_out = ecol_out - wcol_out + 1

      pphoriz_out= n_rows_out*n_cols_out
      global_pphoriz_out=global_n_rows_out*global_n_cols_out

!  1.4 Check to see if any processing is required.
!      Set flag LLPROC if some spatial processing indicated.
!      NB: If input and output bottom levels differ (or the input and
!          output pseudo-levels lists differ), level-by-level
!          processing in the spatial loop IS required.
!          Multi-spatial processing is always required for timeseries.

      lfullfield=((st_grid == st_tp_grid .OR. st_grid == st_cu_grid  &
         .OR. st_grid == st_riv_grid)                                &
         .AND.  stlist(st_west_code,il) == 1.AND.                    &
         stlist(st_east_code,il) == glsize(1,fld_type).AND.          &
         stlist(st_south_code,il) == 1.AND.                          &
         stlist(st_north_code,il) == glsize(2,fld_type)).OR.         &
         ((st_grid == st_uv_grid .OR. st_grid == st_cv_grid)         &
         .AND.  stlist(st_west_code,il) == 1.AND.                    &
         stlist(st_east_code,il) == glsize(1,fld_type).AND.          &
         stlist(st_north_code,il) == glsize(2,fld_type).AND.         &
         stlist(st_south_code,il) == 1) .OR.                         &
         ((st_grid == st_zt_grid)                                    &
         .AND.  stlist(st_north_code,il) == glsize(2,fld_type).AND.  &
         stlist(st_south_code,il) == 1) .OR.                         &
         ((st_grid == st_zu_grid)                                    &
         .AND.  stlist(st_north_code,il) == glsize(2,fld_type).AND.  &
         stlist(st_south_code,il) == 1) .OR.                         &
         ((st_grid == st_mt_grid .OR. st_grid == st_mu_grid)         &
         .AND.  stlist(st_west_code,il) == 1.AND.                    &
         stlist(st_east_code,il) == glsize(1,fld_type)) .OR.         &
         (st_grid == st_scalar)

      lnullproc= lfullfield .AND.                                    &
         (stlist(st_input_bottom,il) ==                              &
         stlist(st_output_bottom,il)) .AND.                          &
         (stlist(st_pseudo_in,il) ==                                 &
         stlist(st_pseudo_out,il)) .AND.                             &
         (stlist(st_gridpoint_code,il) ==                            &
         (extract_base+stash_null_mask_code) .AND.                   &
         stlist(st_weight_code,il) ==                                &
         stash_weight_null_code )
      IF (stlist(st_series_ptr,il) >  0) THEN
        lnullproc=.FALSE.     ! Timeseries always requires processing
      END IF
!  LLPROC must be false for output from a prev stlist
!  or simple extraction of full field with no weighting
      IF ((stlist(st_input_code,il) < 0 .AND.                        &
         stlist(st_proc_no_code,il) /= st_time_series_mean) .OR.     &
         lnullproc) THEN
        llproc=.FALSE.
      ELSE
        llproc=.TRUE.
      END IF
     
!  1.5 Check that no spatial processing is requested if the input field
!      is of integer or logical type -- these types of fields can be
!      passed directly through STASH, for example for coupling purposes,
!      but no arithmetic is allowed at present.

! DEPENDS ON: exppxi
      data_type_code=exppxi(im_ident,is,im,ppx_data_type,            &
      icode,cmessage)
 
      IF ((       data_type_code  == ppx_type_int    .OR.            &
         data_type_code  == ppx_type_log)   .AND.                    &
         .NOT. lnullproc)  THEN
        icode=st_not_supported
        cmessage='STWORK  : Spatial processing of INT/LOGICAL illegal'
        GO TO 9999
      END IF
  
!  2. Perform spatial processing (loop over output levels)

! If this is a variable horizontal grid set flags
! to ensure that bdx and bdy are set properly
! in the pp header. Also work out the size of
! extra data which will be generated. This must be done
! here because the size is needed for the lookup headers
! which are created long before the actual ouput data
! is set up.
      var_grid_type = 0
      extra_var_data = 0

      IF (x_var_grid .OR. y_var_grid) THEN

! Each dimesion has 3 lots of extra grid data -
! grid coordinates, lower and upper box boundaries.
! May need to incorporate test to exclude certain
! grid types (eg LBCs) if they're not appropriate.
        IF (fld_type == fld_type_p) THEN
          var_grid_type = 1 ! It's a T field
          IF (x_var_grid) extra_var_data =                            &
             extra_var_data + (3*(global_n_cols_out+1))
          IF (y_var_grid) extra_var_data =                            &
             extra_var_data + (3*(global_n_rows_out+1))
        ELSE IF (fld_type == fld_type_u.OR.                           &
           fld_type == fld_type_v) THEN
          var_grid_type = 2 ! It's a U,V field
          IF (x_var_grid) extra_var_data =                            &
             extra_var_data + (3*(global_n_cols_out+1))
          IF (y_var_grid) extra_var_data =                            &
             extra_var_data + (3*(global_n_rows_out+1))
        END IF

        extraw_hdr = extraw_hdr + extra_var_data

      END IF

      IF (llproc) THEN   ! Processing is required
        input_code=stlist(st_input_code,il)
!  Make sure no volume processing asked for as not supported
!  this will need adding at some point
        IF (stlist(st_weight_code,il) == stash_weight_volume_code)    &
           THEN
          icode=st_not_supported
          cmessage='STWORK  : volume processing not supported'
          GO TO 9999
        END IF
!  Work out vx,vy (depends on kind of grid)
        vx=lasize(1,fld_type,halo_type)
        vy=lasize(2,fld_type,halo_type)

        IF ((st_grid == st_zt_grid) .OR.                              &
           (st_grid == st_zu_grid)) THEN
          vx=1
        ELSE IF ((st_grid == st_mt_grid) .OR.                         &
           (st_grid == st_mu_grid)) THEN
          vy=1
        ELSE IF (st_grid == st_scalar) THEN
          vx=1
          vy=1
        END IF

!  Work out if this is the first timestep in a timeseries.
!  This is required so that the extra data can be generated
 
        series_ptr=stlist(st_series_ptr,il)
        IF (series_ptr >  0) THEN ! multi spatial processing reqd.
!  Recompute expected sizes
          elap_time=step-stlist(st_start_time_code,il)
          elap_time=MOD(elap_time,stlist(st_period_code,il))
          start_step=(elap_time == 0)
          IF (start_step) THEN
            expected_len=stlist(st_output_length,ilcurr)
          ELSE
            expected_extra=0  ! reset to zero as no extra data
            expected_len=((stlist(st_output_length,ilcurr)           &
               -((stash_series_index(2,series_ptr)+1)*6))*           &
               stlist(st_freq_code,ilcurr))/stlist(st_period_code,   &
               ilcurr)
          END IF
          
!  2.1 Timeseries extraction section follows

          no_records=stash_series_index(2,series_ptr)
          record_start=stash_series_index(1,series_ptr)
          
!  2.1.1 Process a primary field from D1 (timeseries)

          IF (input_code == 0) THEN
            addr=si(im,is,im_index)
! DEPENDS ON: multi_spatial
            CALL multi_spatial(d1(addr),                             &
               vx,vy,vz,grid_type_code,st_grid,fld_type,halo_type,   &
               halosize(1,halo_type),halosize(2,halo_type),          &
               lcyclic,lmasswt,                                      &
               pphoriz_out,num_levs_out,                             &
               no_of_levels_masswt,                                  &
! Extra arguments for atmos sub-model
               d1(jp_1),d1(jpstar),                                  &
! pressure,pstar
               stsuparr(sa_idx(9)),                                  &
! cos_v_latitude
               stsuparr(sa_idx(10)),                                 &
! cos_theta_latitude
               stsuparr(sa_idx(11)),                                 &
! land mask
               stsuparr(sa_idx(12)),                                 &
! sea mask
               row_length,t_rows,u_rows,t_levels,                    &
               ppfield,lenout,                                       &
               rmdi,stlist(1,il),len_stlist,                         &
               stash_series(1,record_start),                         &
               stash_series_rec_len,no_records,                      &
               index_size,index_lev,level_list,                      &
               start_step,extraw,n_rows_out,n_cols_out,              &
               realhd,len_realhd,inthd,len_inthd,                    &
               icode,cmessage)
               
!  2.1.2 Process a field from STASHwork (timeseries)
           
          ELSE IF (input_code == 1) THEN
            addr=si(im,is,im_index)
! DEPENDS ON: multi_spatial
            CALL multi_spatial(stash_work(addr),                     &
               vx,vy,vz,grid_type_code,st_grid,fld_type,halo_type,   &
               halosize(1,halo_type),halosize(2,halo_type),          &
               lcyclic,lmasswt,                                      &
               pphoriz_out,num_levs_out,                             &
               no_of_levels_masswt,                                  &
!     Extra arguments for atmos sub-model
               d1(jp_1),d1(jpstar),                                  &
! pressure,pstar
               stsuparr(sa_idx(9)),                                  &
! cos_v_latitude
               stsuparr(sa_idx(10)),                                 &
! cos_theta_latitude
               stsuparr(sa_idx(11)),                                 &
! land mask
               stsuparr(sa_idx(12)),                                 &
! sea mask
               row_length,t_rows,u_rows,t_levels,                    &
               ppfield,lenout,                                       &
               rmdi,stlist(1,il),len_stlist,                         &
               stash_series(1,record_start),                         &
               stash_series_rec_len,no_records,                      &
               index_size,index_lev,level_list,                      &
               start_step,extraw,n_rows_out,n_cols_out,              &
               realhd,len_realhd,inthd,len_inthd,                    &
               icode,cmessage)
          ELSE IF (input_code <  0) THEN
            
!  2.1.3 Process a field from previously STASHed position in D1
!      (currently unsupported since diagnostic-of-diagnostic)
         
            IF (what_proc == st_time_series_mean) THEN
! special case of mean timeseries
!   Mother record
              ilprev=-stlist(st_input_code,il)
! address of mother record in D1
              addr=stlist(20,ilprev)
! DEPENDS ON: multi_spatial
              CALL multi_spatial(d1(addr),                             &
                 vx,vy,vz,grid_type_code,st_grid,fld_type,halo_type,   &
                 halosize(1,halo_type),halosize(2,halo_type),          &
                 lcyclic,lmasswt,                                      &
                 pphoriz_out,num_levs_out,                             &
                 no_of_levels_masswt,                                  &
!     Extra arguments for atmos sub-model
                 d1(jp_1),d1(jpstar),                                  &
! pressure,pstar
                 stsuparr(sa_idx(9)),                                  &
! cos_v_latitude
                 stsuparr(sa_idx(10)),                                 &
! cos_theta_latitude
                 stsuparr(sa_idx(11)),                                 &
! land mask
                 stsuparr(sa_idx(12)),                                 &
! sea mask                
                 row_length,t_rows,u_rows,t_levels,                    &
                 ppfield,lenout,                                       &
                 rmdi,stlist(1,il),len_stlist,                         &
                 stash_series(1,record_start),                         &
                 stash_series_rec_len,no_records,                      &
                 index_size,index_lev,level_list,                      &
                 start_step,extraw,n_rows_out,n_cols_out,              &
                 realhd,len_realhd,inthd,len_inthd,                    &
                 icode,cmessage)
            ELSE
              icode=st_not_supported
              cmessage='STWORK1  : diag-of-diagnostic unsupported'
              GO TO 9999 ! jump to error return
            END IF
          ELSE
            icode=st_unknown
            WRITE(cmessage,'(A,A,1x,i5)')                              &
               'STWORK: FATAL ERROR: ',                                &
               'unknown input option',input_code
            GO TO 9999
          END IF
          IF (icode /= 0) GO TO 9999 ! error exit
          IF (start_step.AND.(extraw /= expected_extra)) THEN
            icode=st_bad_array_param
            WRITE(cmessage,'(A,i8,1x,i8)')                             &
               'STWORK : Inconsistent length for extra data ',         &
               extraw,expected_extra
            GO TO 9999
          END IF

!          pphoriz_out has been computed by multi_spatial
!          it is the size of the output field
        ELSE  ! do "normal" spatial processing
!
!  2.2 Standard spatial processing section follows
!
!  If multi-level processing (ie. vertical, global mean) is performed
!  by spatial, the input field is passed in with the original start
!  address but if single-level processing is done by spatial, the
!  field is passed in with an address pointing to the single
!  level required.
          base_level0=stlist(st_input_bottom,il)
          what_proc=stlist(st_gridpoint_code,il)
          
          addr_out=1                   ! Initialise output address
          
          DO kl=1,num_levs_out         ! --- Start of levels loop ---
! Work out model level if model level range, otherwise set to 1
            IF (base_level0 <  0.OR.base_level0 == st_special_code)     &
               THEN
              base_level=1
            ELSE
              base_level=base_level0+index_lev(kl)-1
            END IF
! Pass level index into spatial (instead of base_level)
            this_index_lev = kl
            
!  2.2.1 Process a primary field from D1
          
            IF (input_code == 0) THEN
              IF ((what_proc <  vert_mean_top .AND.                     &
                 what_proc >  vert_mean_base) .OR.                      &
                 (what_proc <  global_mean_top .AND.                    &
                 what_proc >  global_mean_base)) THEN
                addr=si(im,is,im_index)
              ELSE
                addr=si(im,is,im_index)+                                &
                   (index_lev(kl)-1)*pphoriz_in
              END IF
              IF (addr <  1.OR.addr >  len_tot) THEN
                icode=st_bad_address
                cmessage='STWORK  : D1 address out of bounds'
                GO TO 9999
              END IF
! DEPENDS ON: spatial
              CALL spatial(d1(addr),vx,vy,vz,                           &
                 grid_type_code,st_grid,                                &
                 fld_type,halo_type,                                    &
                 halosize(1,halo_type),halosize(2,halo_type),           &
                 lcyclic,lmasswt,                                       &
                 n_cols_out,n_rows_out,this_index_lev,                  &
                 level_list,index_lev,index_size,                       &
                 no_of_levels_masswt,                                   &
!     Extra arguments for atmos sub-model
                 d1(jp_1),d1(jpstar),                                   &
! pressure,pstar
                 stsuparr(sa_idx(9)),                              &
! cos_v_latitude
                 stsuparr(sa_idx(10)),                             &
! cos_theta_latitude
                 stsuparr(sa_idx(11)),                             &
! land mask
                 stsuparr(sa_idx(12)),                             &
! sea mask
                 row_length,t_rows,u_rows,                         &
                 blsize(2,fld_type),t_levels,                      &
                 ppfield(addr_out),pphoriz_out,                    &
                 stlist(1,il),len_stlist,rmdi,                     &
                 icode,cmessage)
             
!  2.2.2 Process a field from STASHwork
       
            ELSE IF (input_code == 1) THEN
              IF ((what_proc <  vert_mean_top .AND.                &
                 what_proc >  vert_mean_base) .OR.                 &
                 (what_proc <  global_mean_top .AND.               &
                 what_proc >  global_mean_base)) THEN
                addr=si(im,is,im_index)
              ELSE
                addr=si(im,is,im_index)+                           &
                   (index_lev(kl)-1)*pphoriz_in
              END IF
              IF (addr <  1.OR.addr >  stash_work_len) THEN
                icode=st_bad_address
                cmessage='STWORK  : STASHWORK addr out of bounds'
                GO TO 9999
              END IF
! DEPENDS ON: spatial
              CALL spatial(stash_work(addr),vx,vy,vz,              &
                 grid_type_code,st_grid,                           &
                 fld_type,halo_type,                               &
                 halosize(1,halo_type),halosize(2,halo_type),      &
                 lcyclic,lmasswt,                                  &
                 n_cols_out,n_rows_out,this_index_lev,             &
                 level_list,index_lev,index_size,                  &
                 no_of_levels_masswt,                              &
!     Extra arguments for atmos sub-model
                 d1(jp_1),d1(jpstar),                              &
! pressure,pstar
                 stsuparr(sa_idx(9)),                              &
! cos_v_latitude
                 stsuparr(sa_idx(10)),                             &
! cos_theta_latitude
                 stsuparr(sa_idx(11)),                             &
! land mask
                 stsuparr(sa_idx(12)),                             &
! sea mask
                 row_length,t_rows,u_rows,                         &
                 blsize(2,fld_type),t_levels,                      &
                 ppfield(addr_out),pphoriz_out,                    &
                 stlist(1,il),len_stlist,rmdi,                     &
                 icode,cmessage)
            ELSE IF (input_code <  0) THEN
              
!  2.2.3 Process a field from previously STASHed position in D1
!        (currently unsupported since diagnostic-of-diagnostic)
            
              icode=st_not_supported
              cmessage='STWORK1  : diag-of-diagnostic unsupported'
              GO TO 9999
            ELSE
              icode=st_unknown
              WRITE(cmessage,'(A,A,1x,i5)')     &
                 'STWORK1 : >>FATAL ERROR <<', &
                 'unknown input option',input_code
              GO TO 9999
            END IF
            !
            IF (icode >  0) GO TO 9999 ! Trap error

!  Compute pphoriz_out
!  pphoriz_out is the size of the output vector
!  we should not be doing timeseries processing here.

!  NOTE: n_cols_out and n_rows_out should agree with values
!        calculated before, but are not checked for consistency.

            pphoriz_out=n_cols_out*n_rows_out
            addr_out=addr_out+pphoriz_out ! increment output address
          END DO                          ! --- End of levels loop ---
      
        END IF         ! End of multi-spatial/spatial IF block
       
        IF (icode >  0) GO TO 9999     ! Trap processing error

!  2.3 Set length of output field and check against expected length

!  Calculate size of global pphoriz_out - the size on disk

        IF (what_proc  ==  st_time_series_code .OR.                   &
           what_proc  ==  st_time_series_mean) THEN
          global_pphoriz_out=pphoriz_out
          global_n_rows_out=n_rows_out
          global_n_cols_out=n_cols_out
        ELSE
! DEPENDS ON: stash_get_global_size
          CALL stash_get_global_size(                                   &
             stlist(st_north_code,il) , stlist(st_east_code,il),        &
             stlist(st_south_code,il) , stlist(st_west_code,il),        &
             1,                                                         &
             stlist(st_gridpoint_code,il) , stlist(st_proc_no_code,il), &
             global_pphoriz_out,                                        &
             icode, cmessage)

          IF (icode  /=  0) GO TO 9999

        END IF


        IF (pphoriz_out*num_levs_out /= expected_len) THEN
          icode=st_bad_array_param
          WRITE(cmessage,'(A,i8,1x,i8)')                         &
             'STWORK   : Inconsistent length for output field ', &
             pphoriz_out*num_levs_out,expected_len
          GO TO 9999
        END IF
      ELSE

!  3. No SPATIAL processing - extract output field by direct copy

!  Determine input source

        input_code=stlist(st_input_code,il)
        
!  Other fields are simply copied
        
        IF (input_code == 0) THEN
! Simple extraction with no weighting from primary field in D1
! except for those needing special extraction on funny grids
          addr=si(im,is,im_index)
          DO jl=1,stlist(st_output_length,il)
            ppfield(jl)=d1(addr+jl-1)
          END DO
        ELSE IF (input_code == 1) THEN
! Simple extraction with no weighting from STASH_WORK
! except for those needing special extraction on funny grids
          addr=si(im,is,im_index)
          DO jl=1,stlist(st_output_length,il)
            ppfield(jl)=stash_work(addr+jl-1)
          END DO
        ELSE IF (input_code <  0) THEN
! Previously STASHed entry in D1
! for all sub-models as diagnostic D1 is always on a proper grid.

          addr=stlist(st_output_addr,-input_code)
          DO jl=1,stlist(st_output_length,il)
            ppfield(jl)=d1(addr+jl-1)
          END DO
        ELSE
! Illegal input code
          icode=st_unknown
          cmessage='STWORK   : Unknown input code encountered'
          GO TO 9999
        END IF
      END IF      ! End of LLPROC IF block    ************************

!  4. Output section.

!     The data is in ppfield with a length lenout.
!     The horizontal field size pphoriz_out and number of output levels
!     num_levs_out were calculated in section 1.
!     Output option depends on the stlist code.

!  4.0 Find mother STASHlist record if necessary.

!   Packing_type not set from PP_FILE in the ELSE IF part of block.
      packing_type = 0  ! Default is unpacked.
      IF(stlist(st_input_code,il) <  0) THEN ! Second of two stlist
        ilprev=-stlist(st_input_code,il)
      ELSE
        ilprev=il ! no daughter record
      END IF

!  4.0.1 Set up lbproc sub-components based on STASH processing info.

      DO jj=1,14
        lbproc_comp(jj)=0
      END DO
    
      IF(stlist(st_gridpoint_code,ilprev) >= zonal_mean_base)         &
         lbproc_comp(7)=1
   
      IF((stlist(st_gridpoint_code,ilprev) >= vert_mean_base .AND.    &
         stlist(st_gridpoint_code,ilprev) <  vert_mean_top) .OR.     &
         (stlist(st_gridpoint_code,ilprev) >= global_mean_base .AND. &
         stlist(st_gridpoint_code,ilprev) <  global_mean_top))       &
         lbproc_comp(12)=1
      
      IF((stlist(st_proc_no_code,ilprev) == st_accum_code) .OR.       &
         (stlist(st_proc_no_code,ilprev) == st_time_mean_code).OR.   &
         (stlist(st_proc_no_code,ilprev) == st_time_series_mean))    &
         lbproc_comp(8)=1
      
      IF(stlist(st_proc_no_code,ilprev) == st_min_code)               &
         lbproc_comp(13)=1
      
      IF(stlist(st_proc_no_code,ilprev) == st_max_code)               &
         lbproc_comp(14)=1
      
      output_code=stlist(st_output_code,il)
      
!  4.1 Output to ppfile
      
      IF (output_code <  0) THEN                  ! PP Output
!  Find appropriate dump header if a daughter record
        IF (il /= ilcurr) THEN
          icurrll_dump_ptr=stlist(st_lookup_ptr,ilcurr)
        END IF

!  4.1.0 Determine output PP unit and associated filename; OPEN file

!  If preattached files are used the file is left OPEN by ppctl following
!  the initial OPEN; if reinitialised files are used the unit must be
!  OPENed and CLOSEd explicitly every time it is used.

        ppunit=-output_code

        IF (ft_steps(ppunit) /= 0) THEN ! Filename generated by model
! Check if re-initialised file stream has been opened yet
          IF (step <  ft_firststep(ppunit)) THEN ! File stream opened?
            icode=1
            cmessage='STWORK  : Re-initialised file not yet created'
            WRITE(6,*)                                                  &
               'STWORK  : FATAL ERROR. Attempt to WRITE to ',           &
               're-initialised file stream before file first opened:'
            WRITE(6,*)                                                  &
               '        : check that output on UNIT ',ppunit,' is not', &
               ' requested before first initialisation of output file:'
            WRITE(6,*)                                                  &
               '        :  See UMUI window (Initialisation of PP file', &
               's) accessed from (Post Processing) from (Submodel ',    &
               'independent).'
            GO TO 9999
          END IF                                ! File stream opened ?

          len_ppname=LEN_TRIM(model_ft_unit(ppunit))
! To avoid unnecessary open/close/read/write's, only
! open this file once more, after the initialisation.
          IF(.NOT.model_file_managed(ppunit)) THEN
            CALL model_file_open &
                (ppunit,model_ft_unit(ppunit),len_ppname,1,1,icode,      &
                fileType=ioFileTypeUM)
          END IF
        END IF

!  4.1.2 Find the first available pp lookup record.

        icurrll=ft_lastfield(ppunit) ! Position of the last field
        icurrll=icurrll+1            ! Position of the next field

! Attach the lookup table

        iwa = imdi
        ipplook=>attachLookups(ppunit)

!  4.1.3 Find the first available position for the next data record(s)
!        by reading last pp lookup record.

        IF (icurrll == 1) THEN      ! First record
          iwa=get_file_address(ppunit,mf_data_address)!Loc of start of data
        ELSE
! Pointer to next available data location in output file
          IF (mype == pe_zero) THEN
            iwa = ipplook(lbegin, icurrll-1)+                         &
               ipplook(lbnrec, icurrll-1)
          END IF
        END IF                     ! Test on first record

!  4.1.4 If a daughter record is being processed then recover
!          size information from dump LOOKUP header referenced by
!          mother record.

        IF (il /= ilcurr) THEN
          extraw_hdr=lookup(lbext,icurrll_dump_ptr)

          global_pphoriz_out=lookup(lblrec,icurrll_dump_ptr)
          global_n_rows_out=lookup(lbrow,icurrll_dump_ptr)
          global_n_cols_out=lookup(lbnpt,icurrll_dump_ptr)
          IF (what_proc == st_time_series_mean) THEN 
            pphoriz_out=global_pphoriz_out
            n_rows_out=global_n_rows_out
            n_cols_out=global_n_cols_out
          END IF
          IF (what_proc  ==  st_time_series_code) THEN
            pphoriz_out=global_pphoriz_out
            n_rows_out=global_n_rows_out
            n_cols_out=global_n_cols_out
          END IF

        END IF

!  4.1.5 Check pp_pack_code for grib output. Set grib flag and reset
!        pp_pack_code to give packing profile.
        IF (pp_pack_code(ppunit) >= 100) THEN
          grib_out=.TRUE.
          pp_pack_code(ppunit)=pp_pack_code(ppunit)-100
          grib_packing=pp_pack_code(ppunit)
        ELSE
          grib_out=.FALSE.
        END IF

!  4.1.6 Set packing accuracy for output data field and buffer length
!        Multiple packing profiles are held in STASHmaster and chosen o
!        a per-unit basis through PP_PACK_CODE.  Profile 0 means
!        unpacked. If the field has any extra data switch off packing.

        IF (pp_pack_code(ppunit) == 0.OR.extraw_hdr /= 0) THEN
          packing=.FALSE.
          comp_accrcy=-99
        ELSE
          packing=.TRUE.
! DEPENDS ON: exppxi
          comp_accrcy= exppxi( im_ident, is, im,           &
             ppx_packing_acc+pp_pack_code(ppunit)-1,       &
          icode, cmessage)
        END IF
        IF (grib_out) THEN  ! reset packing code
          pp_pack_code(ppunit)=pp_pack_code(ppunit)+100
        END IF

        lenbuf=((global_pphoriz_out+io_field_padding-1)/io_field_padding)*&
           io_field_padding
! Output length before pack

!  4.2 Select routine to output data using logical grib.
!      If data to be output in grib code then call grib_file
!      If data to be output in pp   code then call pp_file

! If we are using the async stash approach then we prepare a new
! dispatch buffer. We will fill this with all output arising from
! the subsequent levels loop.

        IF (isusingasyncstash())THEN
          iwa=mf_data_missing
          num_words=mf_data_missing
        END IF

        DO  ii=1,num_levs_out           ! --- Start of levels loop ---
! Gather together distributed field to pe 0
! Distributed data is in pp_field, gathered data will be
! in the buf array

          IF ( (what_proc == st_time_series_code) .OR.                  &
             (what_proc == st_time_series_mean)  ) THEN
! If it's timeseries output - just copy on pe_zero
! (usually PE 0)

            IF (mype  ==  pe_zero) THEN
              DO i=1,pphoriz_out
                buf(i)=ppfield(i+(ii-1)*pphoriz_out)
              END DO
            END IF

          ELSE ! not timeseries output

            IF (isusingasyncstash())THEN
              ios_q_slot=getslotfornextlevel(ppunit)

! Prepare the arguments that are understood by the IO stash-server
              ios_packing_flag=ios_no_packing
              IF (packing)ios_packing_flag=ios_packing

! Do we have a full field?
              ios_fullfield_flag=ios_partial_field
              IF ( global_n_rows_out ==                  &
                 glsize(2,get_fld_type(grid_type_code)) &
                 .AND.                                  &
                 global_n_cols_out ==                   &
                 glsize(1,get_fld_type(grid_type_code)))&
                 ios_fullfield_flag=ios_full_field

! Add our local component of the output field to the ios_stash
! buffer we requested before the levels loop
! No communications happen yet we are simply acrueing things
! to send, and a metadata description of what to do with them

              CALL ios_stash_pack_pp_data ( &
                  IOS_Q_Slot,                   &! communications handle
                  ppfield(1+(ii-1)*pphoriz_out:1+ii*pphoriz_out-1) , &  
                                                 ! data in
                  pphoriz_out,                  &! len data in
                  get_fld_type(grid_type_code), &! field type
                  halo_type,                    &! halo type 
                  ios_stash_preprocess,         &! preprocess active
                  ios_packing_flag,             &! packing control
                  ios_fullfield_flag,           &! subdomain control
                  -1,                           &! pack type (not used)
                  comp_accrcy,                  &! coex accuracy
                  rmdi,                         &! 'data missing' value
                                                 ! for reals
                  io_field_padding,             &! field padding 
                  icurrll,                      &! the record number 
                  stlist(st_south_code,il),     &! subdomain boundaries
                  stlist(st_south_code,il)+global_n_rows_out-1, &
                  stlist(st_west_code,il),                      &
                  stlist(st_west_code,il)+global_n_cols_out-1,  &
                  is,                                           &
                  im,                                           &
                  ii                                            &
                  )
            ELSE ! We are not async stash so we do the conventional gather_field 
              packing_type=0
              num_words=0
              num_out=0
! DEPENDS ON: stash_gather_field
              CALL stash_gather_field (                                    &
                 ppfield(1+(ii-1)*pphoriz_out) , buf,                      &
                 pphoriz_out , global_pphoriz_out , 1,                     &
                 stlist(st_south_code,il)+global_n_rows_out-1 ,            &
                 stlist(st_west_code,il)+global_n_cols_out-1,              &
                 stlist(st_south_code,il),                                 &
                 stlist(st_west_code,il),                                  &
                 grid_type_code,halo_type,pe_zero,.TRUE.,            &
                 packing, im_ident, packing_type,                    &
                 num_out,                                                  &
                 comp_accrcy, rmdi,                                        &
                 icode=icode, cmessage=cmessage)

              IF (icode  >   0) THEN
                WRITE(6,*) 'Error occured in STASH while gathering ',       &
                   'data for output.'
                GO TO 9999
              END IF
            END IF ! isUsingAsyncStash()

          END IF ! is a time series output

!  Reset index for pp_head if there is a levels list of hybrid levels
          IF (stlist(st_output_bottom,il) <  0.AND.                         &
             lbvcl  ==  ppx_lbvc_hybrid) THEN
            jj = level_list(ii)
!  or a range of model levels, as they may not be consecutive,
          ELSE IF (stlist(st_output_bottom,il) >= 1.AND.                    &
             stlist(st_output_top,il) <= t_levels) THEN
            jj = level_list(ii)
!  otherwise use of level index is alright
          ELSE
            jj = ii
          END IF
      
!  Check that pp output file has sufficient headers pre-allocated
          IF (icurrll > &
             get_mf_information(ppunit,mf_num_preallocated_headers)) THEN
            icode=4
            WRITE(6,'(A)') 'ERROR detected in routine STWORK'
            WRITE(6,'(A,I5,A,I3)') ': no. of output fields (=',icurrll,')'//  &
               ' exceeds no. of reserved PP headers for unit ',ppunit
            WRITE(cmessage,'(A,I3)')                                         &
               'STWORK: Number of fields exceeds reserved headers for unit ' &
               ,ppunit
            GO TO 9999
          END IF     ! end  no. of pp fields check
         
          IF (grib_out) THEN
            IF (isusingasyncstash())                                        &
               CALL ereport('stwork',99,                                   &
               'grib data not supported for asyncStash')

! NOTE cannot pack data into grib before pphead correctly setup
   
            num_words = -99   ! ie unset before call to pp_head
            packing_type=3
          ELSE
!  Pack data into PP code.

! Set the normal default lengths for no packing
           
            packing_hold=packing
            packing_type_hold=packing_type
            len_buf_words=lenbuf
            len_field=global_pphoriz_out

! Check if we have already packed this data, and if so
! record the packing flag, and turn off packing

            IF (packing.AND.packing_type == 1) THEN

              packing=.FALSE.
              packing_type=0

              num_words=(num_out+1)/2 ! Round up to the nearest 64 Bit words.
              len_buf_words=((num_words+io_field_padding-1)/io_field_padding)*&
                 io_field_padding
              len_field=num_words

            END IF ! packing.AND.packing_type == 1

            IF ( isusingasyncstash()) THEN
              iwa=mf_data_missing
              num_words=mf_data_missing
            ELSE ! Output the gathered field via pp_file

! Restrict the call to output pe only
              IF (mype == pe_zero) THEN
! DEPENDS ON: pp_file
                CALL pp_file(buf,icurrll,                                    &
                   len_buf_words, num_words, rmdi, comp_accrcy,              &
                   len_field, ppunit, iwa,                                   &
                   global_n_cols_out,global_n_rows_out,                      &
                   packing,im_ident,                                         &
                   packing_type,pe_zero,extra_var_data,                &
                   srow_out,wcol_out,icode,cmessage)
              END IF ! mype == pe_zero
              
! Restore the packing flags
 
              IF (packing_type == 0) THEN
                packing=packing_hold
                packing_type=packing_type_hold
              END IF
! Make sure all processors get the return code

              CALL gc_ibcast(101,1,0,nproc,info,icode)

! Num_words is the no of 64 bit words required
              IF (icode >  0) THEN
                GO TO 9999
              END IF
            END IF ! isUsingAsyncStash()
          END IF ! grib out
! Add any extra data concerning variable grids to the extra data
! for lookup header

          IF (isusingasyncstash()) THEN
            len_buf_words=mf_data_missing
          ELSE
            len_buf_words=((num_words+io_field_padding-1)/io_field_padding)*  &
               io_field_padding ! No of words output
          END IF

!  4.2.1 Set STASH processing codes and sampling period for ppheader

          gr=stlist(st_gridpoint_code,ilprev)! Grid point code
! Any time-processed field has a (non-zero) sample_prd set -
! this will be translated by pp_head into an lbtim subcode ib of 2
          sample_prd=0.0
          IF (stlist(st_proc_no_code,ilprev) >  st_replace_code) THEN
            sample_prd=REAL(stlist(st_freq_code,ilprev)*secs_per_period) &
               /REAL(steps_per_period*3600)
          END IF

!  4.2.2 Verification time comes from fixhd(28), current time fixhd(21)
!        2 cases that require consideration here:

!       (1) this record is not a daughter record.
!           in which case, set start_step=.TRUE., verif time from fixhd
!           present time also from fixhd

!       (2) this record IS a daughter record.
!           in which case, will need to retreive info on start_time
!           from dump

          start_step=.TRUE.
! We need to protect the call
          IF (mype == pe_zero) THEN
            IF (il == ilcurr) THEN ! not daughter record
              start_time(1:7)=fixhd(28:34)
              end_time(1:7)=fixhd(21:27)
            ELSE
! Set up start_time from data in LOOKUP(lbyr,icurrll_dump_ptr)
              start_time(1)=lookup(lbyr,icurrll_dump_ptr)
              start_time(2)=lookup(lbmon,icurrll_dump_ptr)
              start_time(3)=lookup(lbdat,icurrll_dump_ptr)
              start_time(4)=lookup(lbhr,icurrll_dump_ptr)
              start_time(5)=lookup(lbmin,icurrll_dump_ptr)
              start_time(6)=lookup(lbsec,icurrll_dump_ptr)
              end_time(1:7)=fixhd(28:34)
            END IF  ! end if block over daughter/mother record
! DEPENDS ON: pp_head
            CALL pp_head(                                                &
               im_ident,fixhd,inthd,realhd,                              &
               len_inthd,len_realhd,                                     &
               im,is,gr,lfullfield,                                      &
               level(ii),pseudo_level(ii),                               &
               samples,start_step,start_time,end_time,len1_lookup,       &
               extraw_hdr,ipplook(1,icurrll),ipplook(1,icurrll),         &
               global_n_cols_out,num_words,len_buf_words,                &
               global_n_rows_out,nrow_in,srow_in,wcol_in,ecol_in,        &
               lbproc_comp,sample_prd,                                   &
               fcst_prd,comp_accrcy,packing_type,                        &
               st_grid,iwa,                                              &
! superarray 1-4: zseak_rho,Ck_rho,zseak_theta,Ck_theta
               stsuparr(sa_idx(1)),stsuparr(sa_idx(2)),                  &
               stsuparr(sa_idx(3)),stsuparr(sa_idx(4)),                  &
               t_levels,jj,                                              &
               rotate,elf,                                               &
               icode,cmessage)

            IF (icode >  0) GO TO 9999 ! An error has occured
          END IF ! mype == pe_zero

          IF (grib_out) THEN
!  Now safe to call grib coder as pphead correctly set apart from
!  length of data
!  Pack data into grib code

            IF (mype  ==  pe_zero) THEN
! DEPENDS ON: grib_file
              CALL grib_file(len1_lookup,pp_len2_lookup,         &
                 ipplook,ipplook,icurrll,                        &
                 buf,global_pphoriz_out,                         &
                 lenbuf,num_words,ppunit,iwa,grib_packing,       &
                 icode,cmessage)
            END IF  ! (IF mype == pe_zero)

! Make sure all processors get the return code

            CALL gc_ibcast(101,1,0,nproc,info,icode)

            IF (icode >  0)THEN
              cmessage='STWORK  : Error in GRIB_FILE'
              GO TO 9999
            END IF

          END IF       ! end of grib_out

          iwa = mf_data_missing
          IF (mype == pe_zero) THEN
            iwa = ipplook(lbegin,icurrll) + ipplook(lbnrec,icurrll)
          END IF

          icurrll=icurrll+1  ! Update the counter
          icurrll_dump_ptr=icurrll_dump_ptr+1
! strictly only needs doing if a daughter record
          
        END DO                           ! --- End of levels loop --

        ft_lastfield(ppunit)=icurrll-1  ! Position of the last field

      ELSE IF (output_code == st_dump.OR.output_code == st_secondary)  &
         THEN ! if(output_code...

!  4.4 Output to dump or secondary D1 space - this implies some
!      temporal processing possibly.  If destination is secondary D1
!      space, there will be no associated lookup header to update.

!      Length is calculated from STASHlist
!      NB: Full field length must be coded here, even for partial timeseries

        num_words=stlist(st_dump_output_length,il)/num_levs_out
        icurrll=stlist(st_lookup_ptr,il) ! Location of dump header

        DO  ii=1,num_levs_out          ! --- Start of levels loop ---
!  Reset index for PP_HEAD if there is a levels list of hybrid levels
          IF (stlist(st_output_bottom,il)  <   0 .AND.                &
             lbvcl == ppx_lbvc_hybrid) THEN
            jj = level_list(ii)
!  or a range of model levels, as they may not be consecutive,
          ELSE IF (stlist(st_output_bottom,il) >= 1.AND.              &
             stlist(st_output_top,il) <= t_levels) THEN
            jj = level_list(ii)
!  otherwise use of level index is alright
          ELSE
            jj = ii
          END IF
          addr=stlist(st_output_addr,il) ! start address
          IF (what_proc == st_time_series_code.OR.                    &
             what_proc == st_time_series_mean) THEN


!  4.4.1 Timeseries addresses are incremented according to timestep
 
            IF (stlist(st_freq_code,il) <  1) THEN
              icode=st_not_supported
              cmessage=                                               &
                 'STWORK  : STASHtime for timeseries not supported'
              GO TO 9999 ! got an error so jump to return
            END IF
            elap_time=step-stlist(st_start_time_code,il)
            elap_time=(MOD(elap_time,stlist(st_period_code,il)))/     &
               stlist(st_freq_code,il)
            addr=addr+(elap_time*pphoriz_out)
!  On the first time step of a timeseries processing
!  pphoriz_out is the length of the entire output vector
!  including extra data -- on other timesteps it is
!  the length of a single record (data for just one timestep)
          END IF

!  4.4.2 Temporal processing from ppfield array to D1

! DEPENDS ON: temporal
          CALL temporal(ppfield(1+(ii-1)*pphoriz_out),                  &
             d1(addr+(ii-1)*pphoriz_out),pphoriz_out,extraw,           &
             stlist(1,il),len_stlist,.FALSE.,step,                     &
             icode,cmessage,start_step,rmdi)

          IF (icode >  0) GO TO 9999
          
!  4.4.3 Set up lookup header if destination is main part of D1
      
          IF (output_code == st_dump) THEN

!  4.4.3 Set other information for input to pphead

            gr=stlist(st_gridpoint_code,ilprev)! Grid point code
! Any time-processed field has a (non-zero) sample_prd set -
! this will be translated by pp_head into an lbtim subcode ib of 2
            sample_prd=0.0
            IF (stlist(st_proc_no_code,ilprev) >  st_replace_code)    &
               THEN
              sample_prd=REAL(stlist(st_freq_code,ilprev)*            &
                 secs_per_period)/REAL(steps_per_period*3600)
            END IF

! Address of whole field is calculated from STASHlist
            iwa=stlist(st_dump_output_addr,il)+(ii-1)*num_words


!  4.4.4 Call pphead to set lookup header for field STASHed to D1.
!        Here pass previous_time as well as start_step from temporal
!        if start_step is TRUE then start time will be updated.
!        Value of end time is unimportant as that is handled properly
!        when data is written out to pp file.
!        Note that lbnrec is hardwired to 0 and so too is bacc.

! DEPENDS ON: pp_head
            CALL pp_head(                                            &
               im_ident,fixhd,inthd,realhd,                          &
               len_inthd,len_realhd,                                 &
               im,is,gr,lfullfield,                                  &
               level(ii),pseudo_level(ii),                           &
               samples,start_step,previous_time,fixhd(28),           &
               len1_lookup,                                          &
               extraw_hdr,lookup(1,icurrll),rlookup(1,icurrll),      &
               global_n_cols_out,num_words,0,                        &
               global_n_rows_out,nrow_in,srow_in,wcol_in,ecol_in,    &
               lbproc_comp,sample_prd,                               &
               fcst_prd,0,packing_type,                              &
               st_grid,iwa,                                          &
! superarray 1-4: zseak_rho,ck_rho,zseak_theta,ck_theta
               stsuparr(sa_idx(1)),stsuparr(sa_idx(2)),              &
               stsuparr(sa_idx(3)),stsuparr(sa_idx(4)),              &
               t_levels,jj,                                          &
               rotate,elf,                                           &
               icode,cmessage)

            IF (icode >  0) GO TO 9999 ! An error has occured

! Only (optionally) pack fields if no extra words of data
            IF (extraw_hdr  ==  0) THEN
! DEPENDS ON: exppxi
              lookup(lbpack,icurrll) =                               &
                 exppxi( im_ident, is, im, ppx_dump_packing,         &
              icode, cmessage)
              IF (dump_pack == 3 ) THEN
! Override packing indicator from PPXREF
                n1 = 0   !   No packing
                lookup(lbpack,icurrll) =                             &
                   (lookup(lbpack,icurrll)/10)*10 + n1
              END IF
            ELSE
              lookup(lbpack,icurrll)=0
            END IF
! Set data type (real/integer) from STASHmaster (-ve for timeseries)
            IF (stlist(st_series_ptr,ilprev) >  0 .OR.                &
               stlist(st_proc_no_code,ilprev) == st_time_series_mean) &
               THEN
! DEPENDS ON: exppxi
              lookup(data_type,icurrll) =                             &
                 -exppxi( im_ident, is, im, ppx_data_type,            &
              icode, cmessage)
            ELSE
! DEPENDS ON: exppxi
              lookup(data_type,icurrll) =                             &
                 exppxi( im_ident, is, im, ppx_data_type,             &
              icode, cmessage)
            END IF
            icurrll=icurrll+1 ! Update the counter for the next field

          END IF ! End of IF output_code=dump

        END DO ! --- End of levels loop ---

      ELSE
        icode=9
        cmessage='STWORK  : Illegal output destination in STLIST'
        GO TO 9999
      END IF      ! End of STLIST output destination IF block

    END IF      ! END OF S_F IF Block ---------------------------------
    
!  5. End of loop over STASHlist entries - Return to calling routine
    
  END DO

9999 CONTINUE
  IF (icode >  0) THEN
    WRITE(6,'(A,I5,A,I4,A,I3)')                               &
       'STWORK: Error when processing diagnostic section ',  &
       is,', item ',im,', code ',icode
    CALL ereport('STWORK',icode,cmessage)
  END IF

  IF (lhook) CALL dr_hook('STWORK',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE stwork
