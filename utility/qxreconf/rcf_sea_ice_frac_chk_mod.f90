! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Checks that sea ice fraction has a minimum value.

Module Rcf_Sea_Ice_Frac_Chk_Mod

!  Subroutine Rcf_Sea_Ice_Frac_Chk
!
! Description:
!   Ensures that sea ice fraction field is reset to zero
!   when it goes below a threshold.
!
! Method:
!   Reset any values below 0.1 to 0.0
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.4   15/12/06   Original code.  P. Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Sea_Ice_Frac_Chk ( fields, field_count, decomp,        &
                                  hdr, data_source)

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

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
    stashcode_icefrac,         &
    stashcode_icethick,        &
    stashcode_prog_sec

Use Rcf_Data_Source_Mod, Only : &
    data_source_type,           &
    Field_Calcs

USE UM_ParVars, Only : &
    nproc,                  &
    mype

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Implicit None

! Arguments
Type( field_type ), Pointer        :: fields(:)
Type( um_header_type ), Intent(In) :: hdr
Integer, Intent(In)                :: field_count
Integer, Intent(In)                :: decomp
Type( data_source_type ), Pointer  :: data_source(:)

! Local variables
Integer                            :: pos_IceFrac
Integer                            :: pos_IceThick
Integer                            :: IceFrac_changed

Real, Parameter                    :: Threshold = 0.1

Integer                            :: i
Integer                            :: istat
Integer                            :: count
Type( field_type ), Pointer        :: Ice_Fraction

!-------------------------------------------------------------------
! Loacte Ice Fraction field in output dump
!-------------------------------------------------------------------

Call Rcf_Locate ( stashcode_prog_sec, stashcode_icefrac, &
                  fields, field_count, pos_IceFrac, .TRUE. )

!--------------------------------------------------------------------
! Reset sea ice fraction to 0 if it is below given threshold value.
!--------------------------------------------------------------------

IceFrac_changed = 0

If (pos_IceFrac /= 0 ) Then

  Ice_Fraction => fields( pos_IceFrac )
  Call Rcf_Alloc_Field( Ice_Fraction )
  Call Rcf_Read_Field( Ice_Fraction, hdr, decomp )

  count = 0
  Do i = 1, Ice_Fraction % level_size

    If ( Ice_Fraction % data (i,1) <  Threshold .AND.                 &
         Ice_Fraction % data (i,1) /= 0.0     ) Then

      Ice_Fraction % data (i,1) = 0.0
      IceFrac_changed = 1
      count = count + 1

    End If

  End Do

!---------------------------------------------------------------------
! Calculate sum of all changes on all PEs
!---------------------------------------------------------------------
  Call gc_isum (1, nproc, istat, count)
  If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
    write (6,*) 'Ice Fraction : No of small values reset to zero ', &
                 count
  End If

!---------------------------------------------------------------------
! Synchronise `changed' flag
!---------------------------------------------------------------------
  Call GC_Imax( 1, nproc, istat, IceFrac_changed )

!-------------------------------------------------------------------
! If there have been changes to IceFrac, find IceThick and reset
! datasource to recalculate IceThick from new IceFrac
!-------------------------------------------------------------------
  If (IceFrac_changed == 1) Then
    Call Rcf_Locate ( stashcode_prog_sec, stashcode_iceThick, &
                      fields, field_count, pos_IceThick, .TRUE.)

    If (pos_IceThick /= 0 ) Then
      If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
        write (6,*) 'Ice Fraction : Setting Ice Thickness to be',  &
                    ' recalculated'
      End If
      data_source( pos_IceThick ) % source      = Field_Calcs
    End If !(pos_IceThick /= 0 )

!---------------------------------------------------------------------
! Write out changed field
!---------------------------------------------------------------------

    Call Rcf_Write_Field( Ice_Fraction, hdr, decomp )
  End If


  Call Rcf_Dealloc_Field( Ice_Fraction )
End If !(pos_IceFrac /= 0 ) - found no icefrac, don't do anything.

Return
End Subroutine Rcf_Sea_Ice_Frac_Chk
End Module Rcf_Sea_Ice_Frac_Chk_Mod
