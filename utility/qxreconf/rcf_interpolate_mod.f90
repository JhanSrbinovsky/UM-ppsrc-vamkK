! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ top level wrapper for interpolation

Module Rcf_interpolate_mod
IMPLICIT NONE

!  Subroutine Rcf_Interpolate - top level interpolation wrapper
!
! Description:
! This module contains a top-level wrapper subroutine for
! interpolation (both horizontal and vertical). It handles
! conversion between data-types (log,int and real) .
! It is worth noting that the orography fields that are fed in
! need not exist if only horizontal interpolation is done.
!
! Method:
!  Intermediate temporary fields are set up - heights are generated,
!  horizontal and then vertical interpolation is done - and then
!  data is reconverted as required.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_interpolate( field_in, field_out, grid_in, grid_out, &
                            interp_orography_in, orography_out )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Grid_Type_Mod, Only : &
    Grid_type

Use Rcf_horizontal_mod, Only : &
    Rcf_horizontal

Use Rcf_vertical_mod, Only : &
    Rcf_vertical

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,             &
    Rcf_Dealloc_Field

Use Rcf_generate_heights_mod, Only : &
    Rcf_generate_heights

USE decomp_params, ONLY : &
    decomp_rcf_output

Use Rcf_Vert_Cloud_Mod, Only : &
    Rcf_Vert_Cloud

Use Rcf_Set_Interp_Flags_Mod, Only :                         &
    interp_v_only,          interp_h_only,                   &
    interp_all,             interp_done,                     &
    interp_copied,          interp_copy

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_cct,             &
    stashcode_ccb,             &
    stashcode_lvl_bse_dp_sc,   &
    stashcode_lvl_top_dp_sc,   &
    stashcode_prog_sec

USE UM_Parvars, ONLY :         &
    change_decomposition

USE cppxref_mod, ONLY :             &
    ppx_type_int,                   &
    ppx_type_log,                   &
    ppx_rho_level,                  &
    ppx_single_level

IMPLICIT NONE

! Arguments
Type (field_type), Intent(InOut) :: field_in
Type (field_type), Intent(InOut) :: field_out
Type (grid_type), Intent(In)     :: grid_in
Type (grid_type), Intent(In)     :: grid_out

! This argument (pos dummy) is for the *FINAL OUTPUT* orography
Type (field_type), Intent(In)    :: orography_out

! This argument (pos dummy) is for the *INTERPOLATED INPUT* orography
Type (field_type), Intent(In)    :: interp_orography_in

! Local variables/parameters
Type (grid_type)                 :: grid_middle
Type (field_type)                :: field_middle
Character (Len=*), Parameter     :: RoutineName = 'Rcf_Interpolate'
Character (Len=80)               :: Cmessage
Integer                          :: ErrorStatus
Integer                          :: i
Integer                          :: j
Integer                          :: log_to_int
Real, Allocatable                :: heights_middle( :, : )
Real, Allocatable                :: heights_out( :, : )
Logical                          :: vertical_required
Logical                          :: horizontal_required
Logical                          :: interp_required


vertical_required   = (field_in % interp == interp_v_only .OR. &
                       field_in % interp == interp_all )
horizontal_required = (field_in % interp == interp_h_only .OR. &
                       field_in % interp == interp_all )
interp_required     = (vertical_required .OR. horizontal_required)

!------------------------------------------------------------
! Tests to see if orography is set should it be required
!------------------------------------------------------------
If ( vertical_required ) Then
  If ( .NOT. Associated( interp_orography_in % Data ) .OR. &
       .NOT. Associated( orography_out % Data ) ) Then
    Cmessage = 'Orography data not present when vertical '//&
               'interpolation required'
    ErrorStatus = 10
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If
End If

!------------------------------------------------------------
! Test to see if Data is associated where necessary
!------------------------------------------------------------
If ( (.NOT. Associated( field_in % Data  )       .AND. &
      .NOT. Associated( field_in % Data_Int  )   .AND. &
      .NOT. Associated( field_in % Data_Log  ) ) .OR.  &
     (.NOT. Associated( field_out % Data )       .AND. &
      .NOT. Associated( field_out % Data_Int )   .AND. &
      .NOT. Associated( field_out % Data_Log ) ) ) Then
  Cmessage = 'An interpolation field has not had Data space allocated!'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If


!------------------------------------------------------------------
! Do the interpolation in the relevant way.
!------------------------------------------------------------------
! Force convective cloud base and convective cloud top to have
! lv_code = ppx_rho_level for the duration of this routine to
! generate correct height field
! Note that have to check section as well as item number
! THIS IS DANGEROUS AS IT ALTERS THE INTERNAL STASHMASTER!!!!
If( (field_out % stashmaster % item    == stashcode_cct  .OR.  &
     field_out % stashmaster % item    == stashcode_ccb  .OR.  &
     field_out % stashmaster % item    == stashcode_lvl_bse_dp_sc  .OR.  &
     field_out % stashmaster % item    == stashcode_lvl_top_dp_sc) .AND. &
     field_out % stashmaster % section == stashcode_prog_sec ) Then
  field_in  % stashmaster % lv_code = ppx_rho_level
  field_out % stashmaster % lv_code = ppx_rho_level
End If

! middle grid/field have out horizont resolution and in vert.
grid_middle = grid_out     ! horizontal resolutions correct
grid_middle % model_levels = grid_in % model_levels
grid_middle % wet_levels   = grid_in % wet_levels
grid_middle % cloud_levels = grid_in % cloud_levels
grid_middle % st_levels    = grid_in % st_levels
grid_middle % sm_levels    = grid_in % sm_levels
grid_middle % bl_levels    = grid_in % bl_levels
grid_middle % ozone_levels = grid_in % ozone_levels
grid_middle % tr_levels    = grid_in % tr_levels
grid_middle % z_top_of_model = grid_in % z_top_of_model
grid_middle % first_constant_r_rho_level =                          &
                               grid_in % first_constant_r_rho_level
grid_middle % height_gen_method = grid_in % height_gen_method

grid_middle % eta_theta_levels => grid_in % eta_theta_levels
grid_middle % eta_rho_levels   => grid_in % eta_rho_levels
grid_middle % soil_depths      => grid_in % soil_depths

field_middle % levels          = field_in  % levels
field_middle % bottom_level    = field_in  % bottom_level
field_middle % top_level       = field_in  % top_level
field_middle % interp          = field_in  % interp
field_middle % rows            = field_out % rows
field_middle % row_len         = field_out % row_len
field_middle % level_size      = field_out % level_size
field_middle % glob_rows       = field_out % glob_rows
field_middle % glob_row_len    = field_out % glob_row_len
field_middle % glob_level_size = field_out % glob_level_size
field_middle % stashmaster     => field_in % stashmaster

! Allocate field_middle Data space
Call Rcf_Alloc_Field( field_middle )

! Allocate heights - need for vertical function call
Allocate( heights_middle( field_middle % level_size,            &
                      0 :  grid_middle % model_levels + 1) )


Allocate( heights_out( field_out % level_size,                  &
                        0 :  grid_out % model_levels + 1) )


! Only need to generate heights if vertical interpolation is actually
! done.  Check at start of routine is done so we can safely use orography data.
If (vertical_required) Then

  Call rcf_generate_heights( grid_middle, interp_orography_in,        &
                         field_middle % stashmaster % grid_type,      &
                         field_middle % stashmaster % lv_code,        &
                         heights_middle,                              &
                         field_middle % level_size )

  Call rcf_generate_heights( grid_out, orography_out,             &
                         field_out % stashmaster % grid_type,     &
                         field_out % stashmaster % lv_code,       &
                         heights_out,                             &
                         field_out % level_size )
End If

!------------------------------------------------------------------
! If data is integer or logical, we will need to convert it to
! Real to do any interpolation on it. Only need to do the conversion
! if some interpolation is actually required!
!------------------------------------------------------------------
If ( interp_required ) Then

  If (field_in % stashmaster % data_type == ppx_type_int ) Then
    Allocate( field_in % Data( field_in % level_size,       &
                               field_in % levels ) )
    Allocate( field_out % Data( field_out % level_size,     &
                                field_out % levels ) )
    Allocate( field_middle % Data( field_middle % level_size,  &
                                   field_middle % levels ) )

    field_in % Data(:,:) = Real( field_in % Data_Int(:,:) )

  Else If (field_in % stashmaster % data_type == ppx_type_log ) Then

    Allocate( field_in % Data( field_in % level_size,       &
                               field_in % levels ) )
    Allocate( field_out % Data( field_out % level_size,     &
                                field_out % levels ) )
    Allocate( field_middle % Data( field_middle % level_size, &
                                   field_middle % levels ) )
    Do j = 1, field_in % levels
      Do i = 1, field_in % level_size
        If ( field_in % Data_Log(i,j) ) Then
          field_in % Data(i,j) = 1.0
        Else
          field_in % Data(i,j) = 0.0
        End If
      End Do
    End Do
  End If
End If

!-------------------------------------------------------------------
! Now can interpolate
!-------------------------------------------------------------------
Call Rcf_horizontal( field_in, field_middle, grid_in, grid_middle )

! Ensure we are on the correct (ie output) decomposition before
! we do vertical interpolation
Call Change_Decomposition( decomp_rcf_output )

If( (field_out % stashmaster % item    == stashcode_cct  .OR.  &
     field_out % stashmaster % item    == stashcode_ccb  .OR.  &
     field_out % stashmaster % item    == stashcode_lvl_bse_dp_sc  .OR.  &
     field_out % stashmaster % item    == stashcode_lvl_top_dp_sc) .AND. &
     field_out % stashmaster % section == stashcode_prog_sec ) Then

  Call Rcf_Vert_Cloud( field_middle, field_out, grid_middle, grid_out, &
                       heights_middle, heights_out)

  ! convert lv_code in stashmaster back to the correct value
  field_in  % stashmaster % lv_code = ppx_single_level
  field_out % stashmaster % lv_code = ppx_single_level
Else
  Call Rcf_vertical( field_middle, field_out, grid_middle, grid_out,&
                     heights_middle, heights_out)
End If

! Deallocate heights
Deallocate( heights_middle )
Deallocate( heights_out )

!-------------------------------------------------------------------
! Deallocate field_middle Data
!-------------------------------------------------------------------
Call Rcf_Dealloc_Field( field_middle )

!-------------------------------------------------------------------
! Need to convert back to integer or logical if field is such. Only
! need to do this if interpolation active!
!-------------------------------------------------------------------
If ( interp_required ) Then
  If (field_out % stashmaster % data_type == ppx_type_int ) Then
    field_out % Data_Int(:,:) = Nint( field_out % Data(:,:) )

    Deallocate( field_out % Data )
    Deallocate( field_in % Data )
    ! middle data already deallocated by dealloc_field above

  Else If (field_out % stashmaster % data_type == ppx_type_log ) Then
    Do i = 1, field_out % level_size
      Do j = 1, field_out % levels
        log_to_int = Int( field_out % Data(i,j) + 0.5 )
        If ( log_to_int == 1 ) Then
          field_out % Data_Log(i,j) =  .TRUE.
        Else If (log_to_int == 0 ) Then
          field_out % Data_Log(i,j) =  .FALSE.
        Else
          Cmessage = 'Conversion from Logical to Integer: Illegal value'
          ErrorStatus = 30
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
      End Do
    End Do

    Deallocate( field_out % Data )
    Deallocate( field_in % Data )
  End If
End If

!--------------------------------------------------------------------
! Set the interp flag for the output field
!--------------------------------------------------------------------
If (interp_required) Then
  field_out % interp = interp_done

Else If (field_in % interp == interp_copy) Then
  field_out % interp = interp_copied

End If

End Subroutine Rcf_interpolate

End Module Rcf_interpolate_mod
