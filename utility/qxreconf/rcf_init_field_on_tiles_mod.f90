! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialisation of a field on tiles from the same field, not tiled.

Module Rcf_Init_Field_On_Tiles_Mod

! Subroutine Rcf_Init_Field_On_Tiles
!
! Description:
!   Initialises a field on land points AND tiles
!   from the same field just on land points.
!
! Method:
!   Copies the input field into each of the ntiles pseudo-levels
!   of the output field.
!   (This routine was modelled on Rcf_Init_Canopy_Water.)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Reconfiguration

Contains

Subroutine Rcf_Init_Field_on_tiles( fields_in,  field_count_in,  hdr_in, &
                                  fld_tiles, field_stashcode )

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_prog_sec

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

USE UM_ParVars, Only : &
    mype

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

USE nlsizes_namelist_mod, ONLY : &
    ntiles

USE decomp_params, ONLY : &
    decomp_rcf_input

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,               &
    Output_Grid

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   &
    interp_copy

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Implicit None

! Arguments
Type( field_type ), Pointer          :: fields_in(:)
Type( field_type ), Intent( InOut )  :: fld_tiles
Type( um_header_type ), Intent( In ) :: hdr_in

Integer, Intent( In )                :: field_count_in
Integer, Intent( In )                :: field_stashcode

! Local variables
Type (field_type), Pointer           :: fld_tiles_in
Type (field_type)                    :: fld_tiles_out
Type (field_type)                    :: dummy  ! pretend orography
Integer                              :: pos    ! field position
Integer                              :: i      ! looper


!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,'(A,A,I6)') 'Initialising field on tiles from ', &
      'input single-layer field: stashcode = ',field_stashcode
End If

!----------------------------------------------------------------------
! Find single-layer field in input fields and read it in
!----------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, field_stashcode,         &
                 fields_in, field_count_in, pos)
fld_tiles_in => fields_in(pos)
Call Rcf_Alloc_Field( fld_tiles_in )
Call Rcf_Read_Field( fld_tiles_in, hdr_in, decomp_rcf_input )

!----------------------------------------------------------------------
! Initialise fld_tiles_out - and allocate space
!----------------------------------------------------------------------
Call Rcf_Field_Equals( fld_tiles_out, fld_tiles_in )
fld_tiles_out % rows            = fld_tiles % rows
fld_tiles_out % row_len         = fld_tiles % row_len
fld_tiles_out % level_size      = fld_tiles % level_size
fld_tiles_out % glob_rows       = fld_tiles % glob_rows
fld_tiles_out % glob_row_len    = fld_tiles % glob_row_len
fld_tiles_out % glob_level_size = fld_tiles % glob_level_size

If (h_int_active) Then
  fld_tiles_in % interp = interp_h_only
Else
  fld_tiles_in % interp = interp_copy
End If

Call Rcf_Alloc_Field( fld_tiles_out )

!----------------------------------------------------------------------
! Interpolate input field to dummy single level output field
! and then copy result to all ntiles pseudolevels.
!----------------------------------------------------------------------
Call Rcf_Interpolate( fld_tiles_in, fld_tiles_out, input_grid, &
                      output_grid, dummy, dummy )

Do i = 1, ntiles
  fld_tiles % Data(:,i) = fld_tiles_out % Data(:,1)
End Do

!----------------------------------------------------------------------
! Tidy up
!----------------------------------------------------------------------
Call Rcf_Dealloc_Field( fld_tiles_out )
Call Rcf_Dealloc_Field( fld_tiles_in )

End Subroutine Rcf_Init_Field_On_Tiles

End Module Rcf_Init_Field_On_Tiles_Mod
