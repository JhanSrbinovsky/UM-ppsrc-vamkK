! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Performs a sanity check on convective cloud variables

Module Rcf_Conv_Cld_Chk_Mod

!  Subroutine Rcf_Conv_Cld_Chk - sanity checks.
!
! Description:
!   Ensures that convective cloud top is at least as high as
!   convective cloud base!
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2   24/01/01   Remove 0 ccb when cct is +ive.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Conv_Cld_Chk( fields, field_count, grid, decomp, hdr )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_done

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_ccb,             &
    stashcode_cct,             &
    stashcode_prog_sec

USE UM_ParVars, Only : &
    nproc

Implicit None

! Arguments
Type( field_type ), Pointer        :: fields(:)
Type( grid_type ), Intent(In)      :: grid
Type( um_header_type ), Intent(In) :: hdr
Integer, Intent(In)                :: field_count
Integer, Intent(In)                :: decomp

! Local variables
Integer                            :: pos_cct
Integer                            :: pos_ccb
Integer                            :: i
Integer                            :: base_changed
Integer                            :: top_changed
Integer                            :: istat
Type( field_type ), Pointer        :: cct
Type( field_type ), Pointer        :: ccb

!-------------------------------------------------------------------
! Find and read conv_cloud_top and conv_cloud_base
!-------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_cct,                  &
                 fields, field_count, pos_cct, .TRUE. )
Call Rcf_Locate( stashcode_prog_sec, stashcode_ccb,                  &
                 fields, field_count, pos_ccb, .TRUE. )

If (pos_ccb /= 0 .AND. pos_cct /= 0 ) Then
  cct => fields( pos_cct )
  Call Rcf_Alloc_Field( cct )
  Call Rcf_Read_Field( cct, hdr, decomp )

  ccb => fields( pos_ccb )
  Call Rcf_Alloc_Field( ccb )
  Call Rcf_Read_Field( ccb, hdr, decomp )
!--------------------------------------------------------------------
! Compare all points
!--------------------------------------------------------------------
  If (ccb % interp == interp_done .OR. cct % interp == interp_done )Then
    base_changed = 0
    top_changed  = 0

    Do i = 1, cct % level_size
      If ( ccb % Data_Int(i,1) == 0 .AND. cct % Data_Int(i,1) > 0 ) Then
        ccb % Data_Int(i,1) = 1
        base_changed = 1
      End If

      If ( cct % Data_Int(i,1) == ccb % Data_Int(i,1) .AND.  &
           cct % Data_Int(i,1) /= 0) Then

        If ( cct % Data_Int(i,1) /= grid % wet_levels ) Then
          ! Convective cloud top and base the same
          ! increment top by one
          cct % Data_Int(i,1) = cct % Data_Int(i,1) + 1
          top_changed = 1

        Else
          ! Convective Cloud top and base the same
          ! decrement base by 1
          ccb % Data_Int(i,1) = ccb % Data_Int(i,1) - 1
          base_changed = 1
        End If
      End If
    End Do

!---------------------------------------------------------------------
! Synchronise `changed' flags
!---------------------------------------------------------------------
    Call GC_Imax( 1, nproc, istat, base_changed )
    Call GC_Imax( 1, nproc, istat, top_changed  )

!---------------------------------------------------------------------
! Write out changed fields
!---------------------------------------------------------------------
    If (top_changed == 1) Then
      Call Rcf_Write_Field( cct, hdr, decomp )
    End If

    If (base_changed == 1) Then
      Call Rcf_Write_Field( ccb, hdr, decomp )
    End If

  End If

  Call Rcf_Dealloc_Field( cct )
  Call Rcf_Dealloc_Field( ccb )

End If

Return
End Subroutine Rcf_Conv_Cld_Chk
End Module Rcf_Conv_Cld_Chk_Mod
