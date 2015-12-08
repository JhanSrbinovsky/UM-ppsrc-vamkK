! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Rotates winds for LAM dumps

Module Rcf_Rotate_Mod

!  Subroutine Rcf_Rotate - rotates winds
!
! Description:
! This routine reads u/v (and u_adv/v_adv if available) and
! rotates them to the standard orientation (Mode = ToStandard)
! or from the standard orientation (Mode = FromStandard).
!
! Method:
! The algorithm is a little different from UM <= 4.5,
! rather interpolating u/v onto theta points, then calculating
! the rotated winds before interpolating back.
! Uses same gather/scatter as horizontal interpolation code.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE w_eqtoll_mod, ONLY: w_eqtoll
USE w_lltoeq_mod, ONLY: w_lltoeq
IMPLICIT NONE

! magic numbers
Integer, Parameter        :: ToStandard   = 0
Integer, Parameter        :: FromStandard = 1

Contains

Subroutine Rcf_Rotate( fields, field_count, grid, hdr, decomp, mode )

Use Rcf_Gather_Field_Mod, Only : &
    Rcf_Gather_Field_Real

Use Rcf_Scatter_Field_Mod, Only : &
    Rcf_Scatter_Field_Real

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

USE UM_ParVars, Only :  &
    mype,               &
    nproc,              &
    gc_all_proc_group,  &
    current_decomp_type,&
    change_decomposition

Use Rcf_Interp_Weights_Mod, Only : &
    Coeff1,                    &
    Coeff2,                    &
    Coeff3,                    &
    Coeff4

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_read_field_mod, Only : &
    Rcf_read_field

Use Rcf_write_field_mod, Only : &
    Rcf_write_field

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use rcf_AC_interp_mod, Only : &
    rcf_p_to_v,               &
    rcf_p_to_u,               &
    rcf_u_to_p,               &
    rcf_v_to_p

USE PrintStatus_mod, Only : &
    LTimer

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_u,              stashcode_v,    &
    stashcode_u_adv,          stashcode_v_adv,&
    stashcode_prog_sec


Implicit None

! Arguments
Type( field_type ), Pointer        :: fields(:)
Type( grid_type ), Intent(In)      :: grid
Type( um_header_type ), Intent(In) :: hdr
Integer, Intent(In)                :: field_count
Integer, Intent(In)                :: decomp
Integer, Intent(In)                :: mode

! Local variables
Integer                            :: i
Integer                            :: j
Integer                            :: upos
Integer                            :: vpos
Integer                            :: pe
Integer                            :: rem
Integer                            :: div
Integer                            :: p_size
Integer                            :: outer
Integer                            :: decomp_original
Integer                            :: ErrorStatus

Real, Allocatable                  :: u_level(:,:)
Real, Allocatable                  :: v_level(:,:)
Real, Allocatable                  :: u_at_p (:,:)
Real, Allocatable                  :: v_at_p (:,:)
Real, Allocatable                  :: u_at_p_out (:,:)
Real, Allocatable                  :: v_at_p_out (:,:)

Type( field_type ), Pointer        :: u_field
Type( field_type ), Pointer        :: v_field

! DEPENDS ON: timer
If (LTimer) Call Timer( 'Rcf_Rotate', 3)

p_size = grid % glob_p_rows * grid % glob_p_row_length;

decomp_original = current_decomp_type
Call Change_Decomposition( decomp )

! We need to do this all twice - once for u/v and then for
! u_adv/u_adv
Do outer = 1, 2

!-----------------------------------------------------------------
! Locate the u and v fields, allocate space for the data and
! read in the field from file
!-----------------------------------------------------------------
  ! find the u-field (u/u_adv)
  If ( outer == 1) Then
    Call Rcf_Locate( stashcode_prog_sec, stashcode_u,                &
                     fields, field_count, upos, .TRUE. )
  Else     ! outer = 2
    Call Rcf_Locate( stashcode_prog_sec, stashcode_u_adv,            &
                     fields, field_count, upos, .TRUE.)
  End If

  If ( upos == 0 ) Then
    Cycle;
  End If

  ! find the v-field (v/v_adv)
  If ( outer == 1) Then
    Call Rcf_Locate( stashcode_prog_sec, stashcode_v,                &
                     fields, field_count, vpos, .TRUE. )
  Else     ! outer = 2
    Call Rcf_Locate( stashcode_prog_sec, stashcode_v_adv,            &
                     fields, field_count, vpos, .TRUE.)
  End If

  If ( vpos == 0 ) Then
    Cycle;
  End If

  ! Allocate space and read in fields
  u_field => fields( upos )
  Call Rcf_Alloc_Field( u_field )
  Call Rcf_Read_Field( u_field, hdr, decomp )

  v_field => fields( vpos )
  Call Rcf_Alloc_Field( v_field )
  Call Rcf_Read_Field( v_field, hdr, decomp )

!-----------------------------------------------------------------
! Gather the fields - one per pe - for the interpolation/rotation
!-----------------------------------------------------------------
  Allocate( u_level( grid % glob_u_row_length, grid % glob_u_rows) );
  Allocate( v_level( grid % glob_v_row_length, grid % glob_v_rows) );
  Allocate( u_at_p ( grid % glob_p_row_length, grid % glob_p_rows) );
  Allocate( v_at_p ( grid % glob_p_row_length, grid % glob_p_rows) );
  Allocate( u_at_p_out( grid % glob_p_row_length, grid % glob_p_rows) );
  Allocate( v_at_p_out( grid % glob_p_row_length, grid % glob_p_rows) );

  div = u_field % levels / nproc
  rem = Mod( u_field % levels, nproc )
  pe = 0

  Do i = 1, div
    Do j = ((i-1) * nproc) + 1, i * nproc
      ! Will gather level j on pe
      ! Cannot use generic subroutine as u_level is 2D

      Call Rcf_Gather_Field_Real( u_field % Data(:,j), u_level,    &
                                  u_field % row_len,               &
                                  u_field % rows,                  &
                                  u_field % glob_row_len,          &
                                  u_field % glob_rows, pe,         &
                                  gc_all_proc_group )

      Call Rcf_Gather_Field_Real( v_field % Data(:,j), v_level,    &
                                  v_field % row_len,               &
                                  v_field % rows,                  &
                                  v_field % glob_row_len,          &
                                  v_field % glob_rows, pe,         &
                                  gc_all_proc_group )

      pe = pe + 1
      If (pe == nproc) pe = 0
    End Do

    CALL Rcf_U_to_P( u_level, u_at_p, grid )
    CALL Rcf_V_to_P( v_level, v_at_p, grid )

    If ( Mode == ToStandard ) Then
      Call W_EQTOLL( Coeff3, Coeff4, u_at_p, v_at_p, u_at_p_out, &
                     v_at_p_out, p_size, .FALSE.)
    Else If ( Mode == FromStandard ) Then
      Call W_LLTOEQ( Coeff1, Coeff2, u_at_p, v_at_p, u_at_p_out, &
                     v_at_p_out, p_size, .FALSE.)
    End If

    CALL Rcf_P_to_U( u_at_p_out, u_level, grid )
    CALL Rcf_P_to_V( v_at_p_out, v_level, grid )

    Do j = ((i-1) * nproc) + 1, i * nproc

      Call Rcf_Scatter_Field_Real( u_field % Data(:,j), u_level,   &
                                   u_field % row_len,              &
                                   u_field % rows,                 &
                                   u_field % glob_row_len,         &
                                   u_field % glob_rows,pe,         &
                                   gc_all_proc_group )

      Call Rcf_Scatter_Field_Real( v_field % Data(:,j), v_level,   &
                                   v_field % row_len,              &
                                   v_field % rows,                 &
                                   v_field % glob_row_len,         &
                                   v_field % glob_rows,pe,         &
                                   gc_all_proc_group )

      pe = pe + 1

      If (pe == nproc) pe = 0
    End Do
  End Do

!-------------------------------------------------------------------
! There are rem levels left to process. Will do these now.
!-------------------------------------------------------------------
  pe = 0
  Do i = 1, rem
    j = nproc * div + i

    Call Rcf_Gather_Field_Real( u_field % Data(:,j), u_level,         &
                                u_field % row_len,                    &
                                u_field % rows,                       &
                                u_field % glob_row_len,               &
                                u_field % glob_rows, pe,              &
                                gc_all_proc_group )

    Call Rcf_Gather_Field_Real( v_field % Data(:,j), v_level,         &
                                v_field % row_len,                    &
                                v_field % rows,                       &
                                v_field % glob_row_len,               &
                                v_field % glob_rows, pe,              &
                                gc_all_proc_group )

    pe = pe + 1
  End Do

  If (mype < pe) Then
    CALL Rcf_U_to_P( u_level, u_at_p, grid )
    CALL Rcf_V_to_P( v_level, v_at_p, grid )

    If ( Mode == ToStandard ) Then
      Call W_EQTOLL( Coeff3, Coeff4, u_at_p, v_at_p, u_at_p_out, &
                     v_at_p_out, p_size, .FALSE. )
    Else If ( Mode == FromStandard ) Then
      Call W_LLTOEQ( Coeff1, Coeff2, u_at_p, v_at_p, u_at_p_out, &
                     v_at_p_out, p_size, .FALSE. )
    End If

    CALL Rcf_P_to_U( u_at_p_out, u_level, grid )
    CALL Rcf_P_to_V( v_at_p_out, v_level, grid )
  End If

  pe = 0
  Do i = 1, rem
    j = nproc * div + i

    Call Rcf_Scatter_Field_Real( u_field % Data(:,j), u_level,         &
                                 u_field % row_len,                    &
                                 u_field % rows,                       &
                                 u_field % glob_row_len,               &
                                 u_field % glob_rows,pe,               &
                                 gc_all_proc_group )

    Call Rcf_Scatter_Field_Real( v_field % Data(:,j), v_level,         &
                                 v_field % row_len,                    &
                                 v_field % rows,                       &
                                 v_field % glob_row_len,               &
                                 v_field % glob_rows,pe,               &
                                 gc_all_proc_group )

    pe = pe + 1

    If (pe == nproc) pe = 0
  End Do


  Deallocate( u_level );
  Deallocate( v_level );
  Deallocate( u_at_p  );
  Deallocate( v_at_p  );
  Deallocate( u_at_p_out );
  Deallocate( v_at_p_out );

!------------------------------------------------------------------
! Write the rotated data back to file and clear up
!------------------------------------------------------------------
  Call Rcf_Write_Field( u_field, hdr, decomp )
  Call Rcf_Write_Field( v_field, hdr, decomp )

  Call Rcf_DeAlloc_Field( u_field )
  Call Rcf_DeAlloc_Field( v_field )

End Do

!-----------------------------------------------------------------
! Back to the original decomposition
!-----------------------------------------------------------------
If ( decomp_original /= current_decomp_type ) Then
  Call Change_Decomposition( decomp_original )
End If

! DEPENDS ON: timer
If (LTimer) Call Timer( 'Rcf_Rotate', 4)

Return
End Subroutine Rcf_Rotate
End Module Rcf_Rotate_Mod
