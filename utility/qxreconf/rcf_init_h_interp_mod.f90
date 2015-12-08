! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up horizontal interpolation weights

Module Rcf_init_h_interp_mod

IMPLICIT NONE

!  Subroutine Rcf_init_h_interp - initialises horizontal interp. weights
!
! Description:
!   Sets up interpolation weights based on choice of scheme.
!
! Method:
!   Allocates space for and sets up weights
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Init_H_Interp( grid_in, grid_out, hdr_in, hdr_out )

Use Rcf_Interp_Weights_Mod  ! Almost all of this used.

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_UMhead_Mod, Only :    &
    um_header_type

Use Ereport_mod, Only : &
    Ereport

Use Rcf_Lsm_Mod, Only : &
    cyclic

Use rcf_h_int_init_bl_mod, Only : &
    rcf_h_int_init_bl

USE Rcf_HeadAddress_Mod, ONLY:&
    FH_GridStagger_Endgame, FH_GridStagger

USE PrintStatus_mod, Only : &
    LTimer

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

! Arguments
Type (grid_type), Intent(In)      :: grid_in
Type (grid_type), Intent(In)      :: grid_out
Type (um_header_type), Intent(In) :: hdr_in
Type (um_header_type), Intent(In) :: hdr_out

! Local data
Integer                       :: interp_arr_sz
INTEGER                       :: interp_rows_sz
INTEGER                       :: max_rows_out, max_field_out
INTEGER                       :: max_rows_in
Integer                       :: ErrorStatus
Character (Len=*), Parameter  :: RoutineName = 'Init_H_Interp'
Character (Len=80)            :: Cmessage

External Timer

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

! Set the cyclic flag for coastal adjustment
If ( Hdr_Out % Fixhd(4) /= 3 .AND. Hdr_Out % FixHd(4) /= 103 ) Then
  CYCLIC=.TRUE.
Else
  CYCLIC=.FALSE.
End If

!----------------------------------------------------------------
! Test if we know which interpolation scheme we are doing!
!----------------------------------------------------------------
IF ( h_int_method == imdi ) THEN
  Cmessage = 'No interpolation method specified - cannot set weights'
  ErrorStatus = 10
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! To deal with Endgame need to find maximum grid size.
! IN EG v has more rows than p and thus glob_v_field can be larger.

interp_arr_sz = MAX(grid_out % glob_v_field,                      &
                    grid_out % glob_p_field)
interp_rows_sz = MAX(grid_out % glob_v_rows + 1,                  &
                     grid_out % glob_p_rows + 1)

Select Case ( h_int_method )
  CASE ( bilinear, nearest_neighbour )

    Allocate( bl_index_b_l( interp_arr_sz, idim ) )
    Allocate( bl_index_b_r( interp_arr_sz, idim ) )
    Allocate( coeff1( interp_arr_sz ) )
    Allocate( coeff2( interp_arr_sz ) )
    Allocate( coeff3( grid_in % glob_p_field ) )
    Allocate( coeff4( grid_in % glob_p_field ) )
    Allocate( weight_t_r( interp_arr_sz, idim ) )
    Allocate( weight_b_r( interp_arr_sz, idim ) )
    Allocate( weight_t_l( interp_arr_sz, idim ) )
    Allocate( weight_b_l( interp_arr_sz, idim ) )
   
    Call rcf_h_int_init_bl( grid_in, grid_out, hdr_in, hdr_out )

  Case ( area_weighted )
    
    Allocate( aw_index_targ_lhs( grid_out % glob_p_row_length + 1, idim ))
    Allocate( aw_index_targ_top( interp_rows_sz, idim ))
    Allocate( aw_colat_t( interp_rows_sz, idim ))
    Allocate( aw_long_l( grid_Out % glob_p_row_length + 1, idim ))
    Allocate( bl_index_b_l( interp_arr_sz, idim ) )
    Allocate( bl_index_b_r( interp_arr_sz, idim ) )
    Allocate( weight_t_r( interp_arr_sz, idim ) )
    Allocate( weight_b_r( interp_arr_sz, idim ) )
    Allocate( weight_t_l( interp_arr_sz, idim ) )
    Allocate( weight_b_l( interp_arr_sz, idim ) )

    IF (hdr_out % FixHd(FH_GridStagger) == FH_GridStagger_Endgame) THEN
      max_rows_out = grid_out % glob_v_rows
      max_field_out = grid_out % glob_v_field
    ELSE
      max_rows_out = grid_out % glob_p_rows
      max_field_out = grid_out % glob_p_field
    END IF

    IF (hdr_in % FixHd(FH_GridStagger) == FH_GridStagger_Endgame) THEN
      max_rows_in = grid_in % glob_v_rows
    ELSE
      max_rows_in = grid_in % glob_p_rows
    END IF

! DEPENDS ON: h_int_init_aw
    CALL h_int_init_aw( icof, idim, max_field_out,                          &
                 max_rows_in, max_rows_out,                                 &
                 grid_in % glob_p_row_length, grid_out % glob_p_row_length, &
                 grid_out % global, hdr_in % FixHd,                         &
                 hdr_out % FixHd, hdr_in % RealC, hdr_out % RealC,          &
                 aw_area_box, aw_index_targ_lhs, aw_index_targ_top,         &
                 bl_index_b_l, bl_index_b_r,                                &
                 aw_colat_t, aw_long_l, weight_t_r, weight_b_r,             &
                 weight_t_l, weight_b_l )

  Case Default
    Cmessage = 'Interpolation method not supported'
    ErrorStatus = 20
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

Return
End Subroutine Rcf_Init_H_Interp

End Module Rcf_init_h_interp_mod
