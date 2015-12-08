! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Derive sea-ice thickness from sea-ice fraction.

Module Rcf_Derv_Ice_Thick_Mod

!  Subroutine Rcf_Derv_Ice_Thick
!
! Description:
!   Derive Ice Thickness from the Ice Fraction field.
!
! Method:
!   Ice Thickness is derived from the Ice Fraction. The ice thickness
!   is set to 2 metres and 1 metre for Northern and Southern hemisphere
!   grid points respectively if the value of ice fraction is greater
!   than zero, ie. if sea-ice is present.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

Contains

Subroutine Rcf_Derv_Ice_Thick ( fields, field_count, hdr,    &
                                ice_thickness )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    grid_type,                &
    output_grid

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE decomp_params, ONLY : &
    decomp_rcf_output

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_icefrac,         &
    stashcode_icethick,        &
    stashcode_prog_sec

Use Rcf_HeadAddress_Mod, Only : &
    RC_LongSpacing,     RC_LatSpacing,     &
    RC_FirstLat,        RC_FirstLong,      &
    RC_PoleLat,         RC_PoleLong

USE UM_ParVars, Only : &
    g_datastart,        mype

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

USE eqtoll_mod, ONLY: eqtoll
Implicit None

! Arguments
Type( field_type ),     Pointer       :: fields(:)
Integer,                Intent(In)    :: field_count
Type( um_header_type ), Intent(In)    :: hdr
Type( field_type ),     Intent(InOut) :: ice_thickness

! Local variables
Real, Parameter                    :: NP_ice_thickness = 2.0
Real, Parameter                    :: SP_ice_thickness = 1.0
Integer                            :: pos_icef
Integer                            :: i,j,ij

Real                               :: latitude                         &
                                      (output_grid % loc_p_row_length, &
                                       output_grid % loc_p_rows)
Real                               :: longitude                        &
                                      (output_grid % loc_p_row_length, &
                                       output_grid % loc_p_rows)

Type( field_type ), Pointer        :: ice_fraction

!-------------------------------------------------------------------
! Locate Ice Fraction in output dump
!-------------------------------------------------------------------

Call Rcf_Locate ( stashcode_prog_sec, stashcode_icefrac,               &
                  fields, field_count, pos_icef )

!--------------------------------------------------------------------
! Allocate space for ice fraction and read in.
!--------------------------------------------------------------------

ice_fraction => fields( pos_icef )
Call Rcf_Alloc_Field( ice_fraction )
Call Rcf_Read_Field( ice_fraction, hdr, decomp_rcf_output )

!--------------------------------------------------------------------
! Set up a latitude/longitude field for the grid.
!--------------------------------------------------------------------

Do j = 1, output_grid % loc_p_rows
  Do i = 1, output_grid % loc_p_row_length

    latitude (i,j) = hdr % RealC(RC_FirstLat) +  &
                    (g_datastart(2,mype)+j-2) *  &
                     hdr % RealC(RC_LatSpacing)

    longitude(i,j) = hdr % RealC(RC_FirstLong) +  &
                    (g_datastart(1,mype)+i-2)  *  &
                     hdr % RealC(RC_LongSpacing)

  End Do
End Do

!--------------------------------------------------------------------
! For rotated grids, get true lats and longs.
!--------------------------------------------------------------------

If (output_grid % rotated) Then

! Use the same arrays to store true lats & longs

  Call EqToLL (latitude, longitude, latitude, longitude,           &
               hdr % RealC(RC_PoleLat), hdr % RealC(RC_Polelong),  &
               output_grid % loc_p_field)

End If

!--------------------------------------------------------------------
! Set up the ice thickness field.
!--------------------------------------------------------------------

ice_thickness % data (:,1) = 0.0

ij =0
Do j = 1, output_grid % loc_p_rows
  Do i = 1, output_grid % loc_p_row_length

    ij = ij + 1

    If ( ice_fraction % data (ij,1) >  0.0 ) Then

      If ( latitude (i,j) >  0.0 ) Then               ! In Northern Hem
        ice_thickness % data (ij,1) = NP_ice_thickness
      Else                                            ! In Southern Hem
        ice_thickness % data (ij,1) = SP_ice_thickness
      End If

    End If

  End Do
End Do

If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
  write (6,*) 'Ice Thickness derived from Ice Fraction '
End If

!--------------------------------------------------------------------
! Deallocate space for ice fraction
!--------------------------------------------------------------------

Call Rcf_Dealloc_Field( ice_fraction )

Return
End Subroutine Rcf_Derv_Ice_Thick
End Module Rcf_Derv_Ice_Thick_Mod
