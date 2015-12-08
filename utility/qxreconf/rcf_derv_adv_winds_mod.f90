! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reinitialises advected winds as prognostic winds.

Module Rcf_Derv_Adv_Winds_Mod

!  Subroutine Rcf_Derv_Adv_Winds
!
! Description:
!   Derive advected wind components from actual wind components
!
! Method:
!   A basic copy with no calculation.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Derv_Adv_Winds( stash_item, fields_out,            &
                               field_count_out, hdr_out,          &
                               advect_wind )

  Use Rcf_Locate_Mod, Only : &
      Rcf_Locate

  Use Rcf_Alloc_Field_Mod, Only : &
      Rcf_Alloc_Field,            &
      Rcf_Dealloc_Field

  Use Rcf_Read_Field_Mod, Only : &
      Rcf_Read_Field

  Use Rcf_Stashcodes_Mod, Only : &
      stashcode_u_adv,     stashcode_v_adv,         &
      stashcode_u,         stashcode_v,             &
      stashcode_prog_sec

  Use Rcf_Field_Type_Mod, Only : &
      Field_Type

  Use Rcf_UMhead_Mod, Only : &
      um_header_type
  
  Use PrintStatus_Mod, Only : &
      PrintStatus,                &
      PrStatus_Normal

  Use UM_Parvars, Only : &
      mype

  USE decomp_params, ONLY : &
      decomp_rcf_output

  Implicit None

  ! Arguments
  Type( field_type ), Pointer       :: fields_out(:)
  Type( um_header_type), Intent(In) :: hdr_out
  Type( field_type ), Intent(InOut), Target :: advect_wind
  Integer, Intent(In)               :: STASH_Item
  Integer, Intent(In)               :: field_count_out

  ! Internal variables
  Type( field_type ), Pointer       ::  wind

  Integer                           ::  pos   ! position in array
  Integer                           ::  i,j,k ! loop index

  !----------------------------------------------------------------------
  ! Write out action if appropriate
  !----------------------------------------------------------------------
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
    If ( STASH_Item == stashcode_u_adv ) Then
      Write (6,*) 'Reinitialising Advected U Wind as U Wind'
    Else If ( STASH_Item == stashcode_v_adv ) Then
      Write (6,*) 'Reinitialising Advected V Wind as V Wind'
    End If
  End If

  !----------------------------------------------------------------------
  ! Find required fields in output dump and read them in where available
  !----------------------------------------------------------------------
  ! Prognostic U or V Wind; will abort if wind not found

  If ( STASH_Item == stashcode_u_adv ) Then
    Call Rcf_Locate(stashcode_prog_sec, stashcode_u, &
                    fields_out, field_count_out, pos)
  Else If ( STASH_Item == stashcode_v_adv ) Then
    Call Rcf_Locate(stashcode_prog_sec, stashcode_v, &
                    fields_out, field_count_out, pos)
  End If

  wind => fields_out(pos)
  Call Rcf_Alloc_Field( wind )
  Call Rcf_Read_Field( wind, hdr_out, decomp_rcf_output )

  !----------------------------------------------------------------------
  ! Loop through wind levels
  !----------------------------------------------------------------------

  Do i = 1,wind % levels
    advect_wind % data(:,i) = wind % data(:,i)
  End Do

  !----------------------------------------------------------------------
  ! Clear up dynamic memory used
  !----------------------------------------------------------------------

  Call Rcf_Dealloc_Field( wind )

  Return
End Subroutine Rcf_Derv_Adv_Winds
End Module Rcf_Derv_Adv_Winds_mod
