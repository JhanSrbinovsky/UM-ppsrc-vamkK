! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  magic numbers and code for setting up field interpolation flags

Module Rcf_Set_Interp_Flags_Mod

IMPLICIT NONE

!  Subroutine Rcf_Set_Interp_Flags - set field interpolation flags
!
! Description:
! This Module contains magic numbers and code for setting up
! field interpolation flags
!
! Method:
!  Based on the v_int_active and h_int_active logical switches.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


!------------------------------------------------------------------
! Magic Numbers
!------------------------------------------------------------------
Integer, Parameter     :: interp_all    = 1    ! h and v to be done
Integer, Parameter     :: interp_h_only = 2    ! Do h - copy for v
Integer, Parameter     :: interp_v_only = 3    ! Do v - copy for h
Integer, Parameter     :: interp_copy   = 4    ! copy for h and v
Integer, Parameter     :: interp_no_op  = 5    ! No interp. operations
Integer, Parameter     :: interp_done   = 6    ! interp. completed
Integer, Parameter     :: interp_copied = 7    ! copy completed

Contains

!--------------------------------------------------------------------
! Subroutine to set interp flags for intput fields based on
! h_int_active, v_int_active and internal rules
!--------------------------------------------------------------------
Subroutine Rcf_Set_Interp_Flags( fields_in, fields_out, field_count_in,&
                                 field_count_out, data_source )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Data_Source_Mod, Only : &
    Data_Source_Type,           &
    Input_Dump,                 &
    Other_Field

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active,         &
    v_int_active_soil

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active,                  &
    h_int_active_u,                &
    h_int_active_v,                &
    h_int_method,                  &
    nearest_neighbour

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_ccb,             &
    stashcode_cct,             &
    stashcode_rho,             &
    stashcode_exner,           &
    stashcode_land_frac,       &
    stashcode_p,               &
    stashcode_prog_sec

Use Rcf_locate_alt_field_mod, Only : &
    rcf_locate_alt_field

USE cppxref_mod, ONLY: &
    ppx_single_level,  &
    ppx_soil_level,    &
    ppx_atm_ozone,     &
    ppx_atm_cuall,     &
    ppx_atm_cvall

USE rcf_recon_mod, ONLY: &
    l_interp_input_only

Implicit None

! Arguments
Type( field_type ), Pointer       :: fields_in(:)
Type( field_type ), Pointer       :: fields_out(:)
Type( data_source_type ), Pointer :: data_source(:)
Integer, Intent(In)               :: field_count_in
Integer, Intent(In)               :: field_count_out

! Local variables
Integer                           :: pos
Integer                           :: i
Logical                           :: vertical
Logical                           :: horizontal

!------------------------------------------------------------------
! Find the output fields that should be sourced from the input dump
!------------------------------------------------------------------
Do i = 1, field_count_out
  If ( data_source(i) % Source == Input_Dump ) Then
    Call Rcf_Locate( fields_out(i) % stashmaster % section,          &
                     fields_out(i) % stashmaster % item,             &
                     fields_in, field_count_in, pos )
  Else If ( data_source(i) % Source == Other_Field ) Then
    Call Rcf_Locate_Alt_Field( fields_out(i),                        &
                               fields_in, field_count_in, pos)
  Else
    CYCLE
  End If



  ! Interpolation decisions based on standard interpolation flags
  vertical   = .FALSE.
  horizontal = .FALSE.
  If (v_int_active) Then
    vertical = .TRUE.
  End If

  If (h_int_active) Then
    horizontal = .TRUE.
  End If

!------------------------------------------------------------------
! Override standard decisions
!------------------------------------------------------------------
  ! Only do h interpolation if single level (not ccb or cct)
  If ( fields_in(pos) % stashmaster % lv_code ==              &
                                      ppx_single_level ) Then
    If ( .NOT.                                                        &
        ( fields_in(pos) % stashmaster % section ==                   &
          stashcode_prog_sec .AND.                                    &
          ( fields_in(pos) % stashmaster % item == stashcode_ccb .OR. &
            fields_in(pos) % stashmaster % item == stashcode_cct ) ) ) Then

      vertical = .FALSE.

    End If
  End If

  ! Soil levels need special treatment
  If ( fields_in(pos) % stashmaster % lv_code == ppx_soil_level) Then
    If ( v_int_active_soil ) Then

      vertical = .TRUE.

    Else

      vertical = .FALSE.

    End If
  End If

  ! We could have P grid the same but the U or V grid are different.  So lets
  ! allow this to occur to minimise interpolation (e.g. ND LAM -> EG LAM)
  IF ( fields_in(pos) % stashmaster % grid_type == ppx_atm_cuall) Then
    horizontal = h_int_active_u  
  END IF
  IF ( fields_in(pos) % stashmaster % grid_type == ppx_atm_cvall) Then
    horizontal = h_int_active_v
  END IF
  
  ! Treat ozone due to it using its own grid.
  If ( fields_in(pos) % stashmaster % grid_type == ppx_atm_ozone) Then

    If ( fields_in(pos) % levels /= fields_out(i) % levels )     Then

    ! If Ozone grid and levels differ then turn on vertical
    ! interpolation whatever the main interpolation flag
      vertical = .TRUE.
    End If

    If ( fields_in(pos) % glob_row_len /= fields_out(i) % glob_row_len ) Then
    ! We might be moving from full field to zonal (or vice versa).
      horizontal = .TRUE.
    End If
  End If


!------------------------------------------------------------------
! Set the field interp flag
!------------------------------------------------------------------
  If (vertical .AND. horizontal) Then
    fields_in(pos) % interp = interp_all

  Else If (vertical) Then
    fields_in(pos) % interp = interp_v_only

  Else If (horizontal) Then
    fields_in(pos) % interp = interp_h_only

  Else
    fields_in(pos) % interp = interp_copy

  End If

! For nearest neighbour we dont care about rebalancing atmosphere.
  IF (h_int_method /= nearest_neighbour .AND. .NOT. l_interp_input_only) THEN
        
  ! Special treatment for rho, exner and p as need calculation rather
  ! than interpolation if would have been interpolated
  If ( fields_in(pos) % stashmaster % section == stashcode_prog_sec .AND. & 
       (fields_in(pos) % stashmaster % item == stashcode_rho        .OR.  &
        fields_in(pos) % stashmaster % item == stashcode_exner      .OR.  &
        fields_in(pos) % stashmaster % item == stashcode_p    )           &
       .AND. ( horizontal .OR. vertical )             ) Then

    fields_in(pos) % interp = interp_no_op

  End If
  END IF

  If ( fields_in(pos) % stashmaster % section == stashcode_prog_sec .AND. &
       fields_in(pos) % stashmaster % item == stashcode_land_frac   .AND. &
       ( horizontal .OR. vertical )             ) Then
    If (vertical .AND. horizontal) Then
      fields_in(pos) % interp = interp_v_only
    Else If (vertical) Then
      fields_in(pos) % interp = interp_v_only
    Else
      fields_in(pos) % interp = interp_no_op
    End If
  End If

End Do

Return
End Subroutine Rcf_Set_Interp_Flags
End Module Rcf_Set_Interp_Flags_Mod
