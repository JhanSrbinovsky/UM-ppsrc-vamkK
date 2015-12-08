! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Checks that snow amount fields contains no negative values.

Module Rcf_Snow_Amount_Chk_Mod

!  Subroutine Rcf_Snow_Amount_Chk
!
! Description:
!   Ensures that snow amount field contains no negative values.
!
! Method:
!   Reset any negative values to zero.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.4   30/08/02   Original code.  D. Robinson
!   6.4   15/12/06   Remove snow on sea points. P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Snow_Amount_Chk( fields, field_count, grid, decomp, hdr )

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

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_mean_snow,       &
    stashcode_prog_sec

USE UM_ParVars, Only : &
    nproc,                  &
    mype

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_LSM_Mod, Only : &
    local_lsm_out

Implicit None

! Arguments
Type( field_type ), Pointer        :: fields(:)
Type( grid_type ), Intent(In)      :: grid
Type( um_header_type ), Intent(In) :: hdr
Integer, Intent(In)                :: field_count
Integer, Intent(In)                :: decomp

! Local variables
Integer                            :: pos_snowa
Integer                            :: snowa_changed

Integer                            :: i
Integer                            :: istat
Integer                            :: count
Type( field_type ), Pointer        :: snow_amount

!-------------------------------------------------------------------
! Loacte snow amount in output dump
!-------------------------------------------------------------------

Call Rcf_Locate ( stashcode_prog_sec, stashcode_mean_snow, &
                  fields, field_count, pos_snowa, .TRUE. )

!--------------------------------------------------------------------
! Reset snow amount to zero where appropriate
!--------------------------------------------------------------------

snowa_changed = 0

If (pos_snowa /= 0 ) Then

  snow_amount => fields( pos_snowa )
  Call Rcf_Alloc_Field( snow_amount )
  Call Rcf_Read_Field( snow_amount, hdr, decomp )

  count = 0
  Do i = 1, snow_amount % level_size

! If snow amount is negative, reset it to 0
    If ( snow_amount % data (i,1) <  0.0 ) Then

      snow_amount % data (i,1) = 0.0
      snowa_changed = 1
      count = count + 1

    End If

! If snow amount is positive over sea, reset it to 0
    If ( snow_amount % data (i,1) > 0.0 .AND.  &
                    .NOT. local_lsm_out(i) ) Then

      snow_amount % data (i,1) = 0.0
      snowa_changed = 1
      count = count + 1

    End If

  End Do

  Call gc_isum (1, nproc, istat, count)
  If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
    Write (6,*) 'Snow Amount : No of inappropriate values reset to zero ', &
                 count
  End If

End If

!---------------------------------------------------------------------
! Synchronise `changed' flag
!---------------------------------------------------------------------
Call GC_Imax( 1, nproc, istat, snowa_changed )

!---------------------------------------------------------------------
! Write out changed field
!---------------------------------------------------------------------

If (snowa_changed == 1) Then
  Call Rcf_Write_Field( snow_amount, hdr, decomp )
End If

If (pos_snowa /= 0) Then
  Call Rcf_Dealloc_Field( snow_amount )
End If

Return
End Subroutine Rcf_Snow_Amount_Chk
End Module Rcf_Snow_Amount_Chk_Mod
