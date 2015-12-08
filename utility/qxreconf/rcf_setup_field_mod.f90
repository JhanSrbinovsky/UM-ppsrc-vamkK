! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialises the field data type for a given dump and grid

MODULE Rcf_Setup_Field_mod

!  Subroutine Rcf_Setup_Field - sets up a field data type array
!
! Description:
!   Sets up the field array for a given grid and dump and "fills in"
!   all the relevant data parts.
!
! Method:
!   Sizes are calculated based on grid code.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

  SUBROUTINE Rcf_Setup_Field( field, hdr, grid, field_count, title,            &
      local_lsm_size )

    USE Rcf_Address_Length_Mod, ONLY :                                         &
        Rcf_Address_Length

    USE Rcf_set_interp_flags_mod, ONLY :                                       &
        interp_no_op

    USE Rcf_Field_Type_Mod, ONLY :                                             &
        field_type

    USE UM_ParVars, ONLY  :                                                    &
        mype,                                                                  &
        nproc

    USE Rcf_UMhead_Mod, ONLY :                                                 &
        UM_header_type

    USE Rcf_Grid_Type_Mod, ONLY :                                              &
        grid_type

    USE Ereport_mod, ONLY :                                                    &
        Ereport

    USE PrintStatus_mod, ONLY :                                                &
        PrintStatus,                                                           &
        PrStatus_Diag,                                                         &
        PrStatus_Oper

    USE Rcf_Exppx_Mod, ONLY :                                                  &
        Rcf_Exppx

    USE Rcf_Recon_Mod, ONLY :                                                  &
        Var_Recon,                                                             &
        l_interp_input_only

    USE Rcf_Level_Code_Mod, ONLY :                                             &
        Rcf_Level_Code

    USE Rcf_Headaddress_Mod, ONLY :                                            &
        FH_Dataset,                                                            &
        FH_Dataset_Ancil,                                                      &
        IC_PLevels

    USE lookup_addresses

    USE rcf_address_mod, ONLY :                                                &
        rcf_disct_lev

    USE cppxref_mod, ONLY :                                                    &
        ppx_lbvc_surface,     ppx_atm_lbc_theta,                               &
        ppx_atm_lbc_u,        ppx_atm_lbc_v,                                   &
        ppx_atm_tall,         ppx_atm_uall,                                    &
        ppx_ocn_tall,         ppx_ocn_uall,                                    &
        ppx_atm_usea,         ppx_atm_tsea,                                    &
        ppx_atm_cuall,        ppx_atm_cvall,                                   &
        ppx_atm_ozone,        ppx_atm_river,                                   &
        ppx_atm_compressed

    IMPLICIT NONE

! Arguments
    TYPE (um_header_type), INTENT (IN)  :: hdr      ! Dump header
    TYPE (grid_type), INTENT (IN)       :: grid     ! grid for fields
    TYPE (field_type), POINTER          :: field( : )
    INTEGER, INTENT (OUT)               :: field_count ! no. of fields
    CHARACTER (LEN=*), INTENT(IN)       :: title
    INTEGER, INTENT(IN), OPTIONAL       :: local_lsm_size
                                       ! size of local land point
                                       ! fields if known

! Comdecks
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

! Local data
    INTEGER            :: i, k        ! loopers
    INTEGER            :: ErrorStatus
    INTEGER            :: model
    INTEGER            :: sec_item
    INTEGER            :: item
    INTEGER            :: section
    INTEGER            :: start_level
    INTEGER            :: end_level
    INTEGER            :: levels
    INTEGER            :: land_sea
    INTEGER            :: mysize
    INTEGER            :: lsm_size
    INTEGER            :: discard_count
    LOGICAL            :: new_field
    CHARACTER (LEN=80) :: Cmessage
    CHARACTER (LEN=80) :: Phrase
    CHARACTER (LEN=*), PARAMETER :: RoutineName = 'setupfield'


! Note that the grid (will be either input or output) should
! correspond to the one referred to in the header lookups etc.

!-----------------------------------------------------------------
! Initialisation of lsm_size from optional value input
!-----------------------------------------------------------------
    IF (PRESENT( local_lsm_size ) ) THEN
      lsm_size = local_lsm_size
    ELSE
      lsm_size = imdi
    END IF

!--------------------------------------------------------------
! Tests for allowed behaviour/values
!-------------------------------------------------------------
    IF (ASSOCIATED(field)) THEN
      Cmessage = 'Field sizes already set'
      ErrorStatus = -20
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
      GO TO 9999
    END IF

! Since some users want to read a single level ancillary using ITEMS namelist
! only check for level consistency if not a ancillary or the lookup contains
! any non-single level data.
    IF ( hdr % FixHd(FH_Dataset)   /= FH_Dataset_Ancil .OR.                    &
        ANY(hdr % Lookup(lbvc, :) /= ppx_lbvc_surface) ) THEN
      IF ( grid % model_levels /= hdr % IntC(IC_PLevels) ) THEN
        Cmessage = 'Grid and headers disagree on number of levels'
        ErrorStatus = 10
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )
      END IF
    END IF

!--------------------------------------------------------------
! Need to count the total number of fields for allocation
! of memory.
!---------------------------------------------------------------
    field_count = 1
    DO i = 1, hdr % Len2Lookup
  ! Check field is valid for reconfiguring and skip if not valid.
      IF (.NOT. is_valid_field( title, hdr % lookup(:,i) ) ) THEN
        CYCLE
      END IF

      IF ( i /= 1 ) THEN
    ! If not the same field as before
        IF ( .NOT. (hdr % Lookup(item_code, i) ==                              &
            hdr % Lookup(item_code, i - 1) ) ) THEN
          field_count = field_count + 1
        END IF
      END IF
    END DO

! Allocate the space
    ALLOCATE( field( field_count ) )

!---------------------------------------------------------------
! Now initialise the fields
!---------------------------------------------------------------
    discard_count = 0
    field_count     = 1
    new_field = .TRUE.

    field( field_count ) % levels   = 1
    field( field_count ) % dump_pos = 1

    DO i = 1, hdr % Len2Lookup
  ! Check field is valid for reconfiguring and skip if not valid.
      IF (.NOT. is_valid_field( title, hdr % lookup(:,i) ) ) THEN
        discard_count=discard_count+1
        CYCLE
      END IF

      IF ( i /= 1 ) THEN
        IF (hdr % Lookup(item_code, i) ==                                      &
            hdr % Lookup(item_code, i - 1) )  THEN

          field( field_count ) % levels = field( field_count ) % levels + 1
          new_field = .FALSE.
        ELSE
          field_count = field_count + 1
          field( field_count ) % levels = 1
          field( field_count ) % dump_pos =                                    &
              field( field_count - 1) % dump_pos +                             &
              field(field_count - 1) % levels
          new_field = .TRUE.
        END IF
      END IF

  ! Default - nothing doing....
      field( field_count ) % interp =  interp_no_op

! Need STASHmaster record for some information.
      sec_item = hdr % Lookup(item_code, i)
      item     = MOD( hdr % Lookup(item_code, i), 1000)
      section  = (hdr % Lookup(item_code, i) - item) / 1000
      model    = hdr % Lookup(model_code , i)

      IF (title == 'Input grid') THEN
        field( field_count) % stashmaster => Rcf_Exppx( model, section, item,  &
            stash_in_arg = .TRUE. )
      ELSE IF (title == 'Auxillary File') THEN
        field( field_count) % stashmaster => Rcf_Exppx( model, section, item )
      ELSE IF (title == 'Output grid') THEN
        field( field_count) % stashmaster => Rcf_Exppx( model, section, item )
      ELSE
        cmessage = 'Unrecognised usage.  Assuming output STASHmaster for usage.'
        errorstatus = -10
        CALL ereport( routinename, errorstatus, cmessage )
        field( field_count) % stashmaster => Rcf_Exppx( model, section, item )
      END IF

! Need number of levels for lbc sizing
      IF ( field( field_count ) % stashmaster % grid_type ==                   &
          ppx_atm_lbc_theta .OR.                                               &
          field( field_count ) % stashmaster % grid_type ==                    &
          ppx_atm_lbc_u     .OR.                                               &
          field( field_count ) % stashmaster % grid_type ==                    &
          ppx_atm_lbc_v ) THEN

        CALL Rcf_Level_Code( field( field_count ) % stashmaster % lb_code,     &
            start_level, grid )
        CALL Rcf_Level_Code( field( field_count ) % stashmaster % lt_code,     &
            end_level,   grid )

        levels = (end_level - start_level) + 1
      END IF
      IF (rcf_disct_lev( field( field_count ) % stashmaster % lv_code )) THEN
        ! Set the bottom and top level information.
        CALL Rcf_Level_Code( field( field_count ) % stashmaster % lb_code,     &
          field( field_count ) % bottom_level, grid )
        CALL Rcf_Level_Code( field( field_count ) % stashmaster % lt_code,     &
          field( field_count ) % top_level,    grid )
      ELSE
        ! Otherwise label it as unset.
        field( field_count ) % bottom_level = -1
        field( field_count ) % top_level    = -1
      END IF

!------------------------------------------------------------------
! Sizes depend on the grid type - cases will need to be added as
! they are needed. Assume we don't mix our grids in a single
! call to this function (ie some fields B, some C grid etc)
!------------------------------------------------------------------
      SELECT CASE ( field( field_count ) % stashmaster % grid_type )
      CASE (ppx_atm_tall, ppx_atm_tsea, ppx_ocn_tall)      ! Theta points
        field( field_count ) % rows            = grid % loc_p_rows
        field( field_count ) % row_len         = grid % loc_p_row_length
        field( field_count ) % level_size      = grid % loc_p_rows *           &
            grid % loc_p_row_length

        field( field_count ) % glob_rows       = grid % glob_p_rows
        field( field_count ) % glob_row_len    = grid % glob_p_row_length
        field( field_count ) % glob_level_size = grid % glob_p_rows *          &
            grid % glob_p_row_length

      CASE (ppx_atm_cuall)      ! C grid u points
        field( field_count ) % rows            = grid % loc_u_rows
        field( field_count ) % row_len         = grid % loc_u_row_length
        field( field_count ) % level_size      = grid % loc_u_rows *           &
            grid % loc_u_row_length

        field( field_count ) % glob_rows       = grid % glob_u_rows
        field( field_count ) % glob_row_len    = grid % glob_u_row_length
        field( field_count ) % glob_level_size = grid % glob_u_rows *          &
            grid % glob_u_row_length

      CASE (ppx_atm_cvall)      ! C grid v points
        field( field_count ) % rows            = grid % loc_v_rows
        field( field_count ) % row_len         = grid % loc_v_row_length
        field( field_count ) % level_size      = grid % loc_v_rows *           &
            grid % loc_v_row_length

        field( field_count ) % glob_rows       = grid % glob_v_rows
        field( field_count ) % glob_row_len    = grid % glob_v_row_length
        field( field_count ) % glob_level_size = grid % glob_v_rows *          &
            grid % glob_v_row_length

      CASE (ppx_atm_ozone)      ! Ozone grid
        field( field_count ) % rows            = grid % loc_p_rows
        field( field_count ) % glob_rows       = grid % glob_p_rows

        IF ( hdr % Lookup(lbnpt,i) == 1) THEN
          field( field_count ) % row_len         = 1
          field( field_count ) % glob_row_len    = 1
        ELSE
          field( field_count ) % row_len         = grid % loc_p_row_length
          field( field_count ) % glob_row_len    =                             &
              grid % glob_p_row_length
        END IF

        field( field_count ) % level_size      =                               &
            field(field_count) % row_len * field(field_count) % rows
        field( field_count ) % glob_level_size =                               &
            field(field_count) % glob_row_len * field(field_count) % glob_rows

      CASE (ppx_atm_compressed)     ! Land compressed points
      ! Can only set global size here. Local level_size will need to
      ! be set when sizes are available.
        field( field_count ) % rows            = imdi
        field( field_count ) % row_len         = imdi
        field( field_count ) % level_size      = lsm_size


        field( field_count ) % glob_rows       = imdi
        field( field_count ) % glob_row_len    = imdi
        field( field_count ) % glob_level_size = hdr % Lookup(lblrec,i)

      CASE (ppx_atm_river)    ! River routing points
        field( field_count ) % rows            = grid % loc_r_rows
        field( field_count ) % row_len         = grid % loc_r_row_length
        field( field_count ) % level_size      = grid % loc_r_rows *           &
            grid % loc_r_row_length

        field( field_count ) % glob_rows       = grid % glob_r_rows
        field( field_count ) % glob_row_len    = grid % glob_r_row_length
        field( field_count ) % glob_level_size = grid % glob_r_rows *          &
            grid % glob_r_row_length

      CASE Default

        IF (section /= 0 .AND. section /= 33 .AND. section /= 34               &
            .AND. .NOT. var_recon .AND. .NOT. l_interp_input_only) THEN
          ! Diagnostic in input dump.
          ! Rcf does not process diagnostics in input dump, so set
          ! relevant field dimensions to imdi in field.  Only allow to
          ! reconfigure diagnostics for VAR or if we are only interpolating the
          ! input fields.
          field( field_count ) % glob_rows       = imdi
          field( field_count ) % glob_row_len    = imdi
          field( field_count ) % rows            = imdi
          field( field_count ) % row_len         = imdi
          field( field_count ) % glob_level_size = imdi
          field( field_count ) % level_size      = imdi

        ELSE
          ! If we are trying to interpolate a grid which is
          ! unknown in section 0/33/34, or especially for VAR and when
          ! interpolating input fields only, we need know about it.

          WRITE (6,*) 'Unsupported Grid-type ',                                &
              field( field_count ) % stashmaster % grid_type
          Cmessage = 'Grid type not yet catered for - size will be&
          &unset IN field data STRUCTURE'
          ErrorStatus = -10
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )

        END IF
      END SELECT
    END DO

!------------------------------------------------------------------
! This should have all values calculated correctly
! Only remains to print them out if required.
!------------------------------------------------------------------
    IF (discard_count == hdr % Len2Lookup) THEN
       WRITE(Cmessage,'(A,A15)') 'No instantaneous fields found in ',title
       ErrorStatus = 40
       CALL Ereport( RoutineName, ErrorStatus, Cmessage )  
    ELSE

      IF (mype == 0 .AND. PrintStatus >= PrStatus_Oper ) THEN
        WRITE(6,'(''  '',/,'' '',A/)')title

        i = 1
        DO k = 1, field_count

          Phrase = field(k) % stashmaster % name
          IF ( field(k) % stashmaster % grid_type == ppx_atm_compressed )      &
          THEN
            land_sea = 1
          ELSE
            land_sea = 0
          END IF

          i = i + field(k) % levels
          WRITE(6,'('' '',I4,I5,I8,I4,3I6,2x,A36)')                            &
            land_sea, field(k) % levels,                                       &
            field(k) % glob_level_size,                                        &
            field(k) % stashmaster % data_type,                                &
            field(k) % stashmaster % section,                                  &
            field(k) % stashmaster % item,                                     &
            field(k) % dump_pos, phrase
        END DO
      END IF
    END IF
    9999 CONTINUE
  CONTAINS
    LOGICAL FUNCTION is_valid_field(title, fld_lookup)

      USE rcf_ppx_info_mod, ONLY :                                             &
          stm_record_type

  ! Arguments:
  ! Identifier for type of field.
      CHARACTER (LEN=*), INTENT(IN)       :: title
  ! The 64 word header for the field.
      INTEGER, INTENT(IN)                 :: fld_lookup(:)

  ! Parameters:
      INTEGER, PARAMETER                  :: lbcode_reg_lat_lon     = 1
      INTEGER, PARAMETER                  :: lbcode_reg_lat_lon_rot = 101
      INTEGER, PARAMETER                  :: lbproc_no_processing   = 0

  ! Local variables:
      TYPE (stm_record_type), POINTER :: stm_entry

! Retrieve STASHmaster entry for field.
! If Nofindarg is false then rcf_exppx will abort if nothing is found.
! We use a test on the return pointer to decide whether field is valid.
      stm_entry => rcf_exppx(fld_lookup(model_code),                           &
          fld_lookup(item_code) / 1000,                                        &
          MOD( fld_lookup (item_code), 1000),                                  &
          nofindarg    = title == 'Input grid',                                &
          stash_in_arg = title == 'Input grid')

! Check grid code and lbproc since we should only be using lat/lon grid and
! processing instantaneous fields.  Reconfiguration for VAR might require
! interpolation of statistically processed fields.
      IF ( (fld_lookup(lbcode) /= lbcode_reg_lat_lon       .AND.               &
          fld_lookup(lbcode) /= lbcode_reg_lat_lon_rot)  .OR.                  &
          (fld_lookup(lbproc) /= lbproc_no_processing     .AND.                &
          .NOT. var_recon)) THEN
        IF (mype == 0 .AND. PrintStatus >= PrStatus_Diag ) THEN
          WRITE(6,'(A)')                                                       &
              "INFO: Detected incompatible field in is_valid_field."
          WRITE(6,'(A,3I7,A)') "INFO: Disabling (item_code,lbcode,lbproc): (", &
              fld_lookup(item_code), fld_lookup(lbcode), fld_lookup(lbproc), ")"
        END IF
        is_valid_field = .FALSE.
! If there is no stashmaster for input field then field is not valid.
! If not an input field then rcf_exppx would have aborted with no find.
      ELSE IF ( .NOT. ASSOCIATED(stm_entry) ) THEN
        IF (mype == 0 .AND. PrintStatus >= PrStatus_Diag ) THEN
          WRITE(6,'(A)')    "INFO: No STASHmaster available in is_valid_field."
          WRITE(6,'(A,I7)') "INFO: Disabling item_code ",fld_lookup(item_code)
        END IF

        is_valid_field = .FALSE.

      ELSE
    ! Default to true.
        is_valid_field = .TRUE.
      END IF

    END FUNCTION is_valid_field
  END SUBROUTINE Rcf_Setup_Field
END MODULE Rcf_Setup_Field_Mod
