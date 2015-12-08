! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Deriving 2d convective cloud amounts from 2d.

Module Rcf_Derv_2D_CCA_Mod

!  Subroutine Rcf_Derv_2D_CCA
!
! Description:
!    Top level routine obtaining variables (and some conversions)
!    for the Rcf_Calc_2D_CCA routine.
!    To avoid a 3D interpolation, calculation is done with
!    input grid variables and the result is interpolated.
!    Some correction of ccb and cct is then required.
!
! Method:
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Reconfiguration

Contains

Subroutine Rcf_Derv_2D_CCA( fields_in, fields_out, field_count_in, &
                            field_count_out, hdr_in, hdr_out, cca_2d)

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,           &
    Output_Grid

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE decomp_params, ONLY : &
    decomp_rcf_output,   &
    decomp_rcf_input

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   &
    interp_copy

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Calc_2d_CCA_Mod, Only : &
    Rcf_Calc_2d_CCA

USE UM_ParVars, Only : &
    mype

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_3d_cca,          &
    stashcode_ccb,             &
    stashcode_cct,             &
    stashcode_prog_sec

Implicit None

! Arguments
Type( field_type ), Pointer        :: fields_in(:)
Type( field_type ), Pointer        :: fields_out(:)
Type( field_type ), Intent(InOut)  :: cca_2d
Type( um_header_type ), Intent(In) :: hdr_in
Type( um_header_type ), Intent(In) :: hdr_out
Integer, Intent(In)                :: field_count_in
Integer, Intent(In)                :: field_count_out

! Local variables
Integer                       :: i           ! looper
Integer                       :: pos         ! field position

Type( field_type ), Pointer   :: cca_3d
Type( field_type ), Pointer   :: ccb
Type( field_type ), Pointer   :: cct
Type( field_type ), Pointer   :: ccb_in
Type( field_type ), Pointer   :: cct_in
Type( field_type )            :: dummy
Type( field_type )            :: cca_2d_in

!--------------------------------------------------------------
! Write out out action if appropriate
!--------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Initialising 2D CCA'
End If

!--------------------------------------------------------------
! Find and setup 3D CCA from input
!--------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_3d_cca,               &
                 fields_in, field_count_in, pos )
cca_3d => fields_in(pos)
Call Rcf_Alloc_Field( cca_3d )
Call Rcf_Read_Field( cca_3d, hdr_in, decomp_rcf_input )

!---------------------------------------------------------------
! Set up the cca_2d_in field
!---------------------------------------------------------------
Call rcf_field_equals( cca_2d_in, cca_2d )
cca_2d_in % rows            = cca_3d % rows
cca_2d_in % row_len         = cca_3d % row_len
cca_2d_in % level_size      = cca_3d % level_size
cca_2d_in % glob_rows       = cca_3d % glob_rows
cca_2d_in % glob_row_len    = cca_3d % glob_row_len
cca_2d_in % glob_level_size = cca_3d % glob_level_size

Call Rcf_Alloc_Field( cca_2d_in )
!---------------------------------------------------------------
! Read conv cloud base and top as these are required
!---------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_ccb,                  &
                 fields_in, field_count_in, pos )
ccb_in => fields_in(pos)
Call Rcf_Alloc_Field( ccb_in )
Call Rcf_Read_Field( ccb_in, hdr_in, decomp_rcf_input )

Call Rcf_Locate( stashcode_prog_sec, stashcode_cct,                  &
                 fields_in, field_count_in, pos )
cct_in => fields_in(pos)
Call Rcf_Alloc_Field( cct_in )
Call Rcf_Read_Field( cct_in, hdr_in, decomp_rcf_input )

!-----------------------------------------------------------------
! Have all fields we require - call the calculation routine
!-----------------------------------------------------------------
Call Rcf_Calc_2D_CCA( cca_3d, ccb_in, cct_in, cca_2d_in )

!--------------------------------------------------------------
! Get rid of the fields we no longer require
!--------------------------------------------------------------
Call Rcf_Dealloc_Field( cca_3d )
Call Rcf_Dealloc_Field( ccb_in )
Call Rcf_Dealloc_Field( cct_in )


!--------------------------------------------------------------
! We now need to interpolate cca_2d_in to the output grid
!--------------------------------------------------------------
If (h_int_active) Then
  cca_2d_in % interp = interp_h_only
Else
  cca_2d_in % interp = interp_copy
End If

Call Rcf_Interpolate( cca_2d_in, cca_2d, input_grid, output_grid, &
                      dummy, dummy)

!----------------------------------------------------------------
! Interpolation can mean that ccb, cct and cca_2d are not
! consistent. A quick `clean up' should prevent later crashes.
!----------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_ccb,                  &
                 fields_out, field_count_out, pos )
ccb => fields_out(pos)
Call Rcf_Alloc_Field( ccb )
Call Rcf_Read_Field( ccb, hdr_out, decomp_rcf_output )

Call Rcf_Locate( stashcode_prog_sec, stashcode_cct,                  &
                 fields_out, field_count_out, pos )
cct => fields_out(pos)
Call Rcf_Alloc_Field( cct )
Call Rcf_Read_Field( cct, hdr_out, decomp_rcf_output )

Do i = 1, cca_2d % level_size
  ! Make sure cca is +ive
  If (cca_2d % Data(i,1) < 0. ) Then
    cca_2d % Data(i,1) = 0.
  End If

  ! Make sure cca only exists where ccb and cct are both non-zero
  If ( ccb % Data_Int(i,1) == 0 .AND.     &
       cct % Data_Int(i,1) == 0 .AND.     &
       cca_2d % Data(i,1) > 0. ) Then
    cca_2d % Data(i,1) = 0.0
  End If
End Do

!-------------------------------------------------------------
! Write out the required fields and tidy up memory
!-------------------------------------------------------------
Call Rcf_Write_Field( ccb, hdr_out, decomp_rcf_output )
Call Rcf_Write_Field( cct, hdr_out, decomp_rcf_output )

Call Rcf_Dealloc_Field( cca_2d_in )
Call Rcf_Dealloc_Field( ccb )
Call Rcf_Dealloc_Field( cct )

Return
End Subroutine Rcf_Derv_2D_CCA
End Module Rcf_Derv_2D_CCA_Mod
