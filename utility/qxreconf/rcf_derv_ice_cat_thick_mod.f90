! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reinitialises sea ice category fields for thickness and concentration

Module Rcf_Derv_Ice_Cat_Thick_Mod

! Description: Populate the category fields for ice thickness and ice
!              fraction based on the Grid Box Mean (GBM) of both fields.
!              Done by generating bounds for the thickness of each
!              category (pseudo level) and assinging values for fraction
!              and thickness by mapping where the GBM thicknesses fall
!              within the bounds for that category (pseudo level).
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Derv_Ice_Cat_Thick( stash_item, fields_out,            &
                                   field_count_out, hdr_out,          &
                                   ice_cat_field )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_icefrac,         stashcode_icethick,         &
    stashcode_ice_conc_cat,    stashcode_ice_thick_cat,    &
    stashcode_prog_sec

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

USE UM_ParVars, Only : &
    mype

USE decomp_params, ONLY : &
    decomp_rcf_output

USE nlsizes_namelist_mod, ONLY: &
    Nice

Implicit None

! Arguments
Type( field_type ), Pointer       :: fields_out(:)
Type( um_header_type), Intent(In) :: hdr_out
Type( field_type ), Intent(InOut), Target :: ice_cat_field
Integer, Intent(In)               :: STASH_Item
Integer, Intent(In)               :: field_count_out

! Internal variables
Type( field_type ), Pointer       ::  ice_thick
Type( field_type ), Pointer       ::  ice_frac
Type( field_type ), Pointer       ::  ice_cat_thick
Type( field_type ), Pointer       ::  ice_cat_frac

Real                              ::  rNice ! real of 1 / Nice
Real                              ::  HIN_Max(0:Nice)
Real                              ::  CC1,CC2,CC3,X1

Integer                           ::  pos   ! position in array
Integer                           ::  i,j,k ! loop index



!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
    If ( STASH_Item == stashcode_ice_conc_cat ) Then
      Write (6,*) 'Reinitialising Ice Concentrations on categories'
    Else If ( STASH_Item == stashcode_ice_thick_cat ) Then
      Write (6,*) 'Reinitialising Ice Thickness on categories'
    End If
  End If

!----------------------------------------------------------------------
! Find required fields in output dump and read them in where available
!----------------------------------------------------------------------
! GBM ice thickness; will abort if icethick not found

  Call Rcf_Locate( stashcode_prog_sec, stashcode_icethick,           &
                   fields_out,field_count_out,pos)

  ice_thick => fields_out(pos)
  Call Rcf_Alloc_Field( ice_thick )
  Call Rcf_Read_Field( ice_thick, hdr_out, decomp_rcf_output )

! GBM ice fraction; will abort if icefrac not found

  Call Rcf_Locate( stashcode_prog_sec, stashcode_icefrac,            &
                   fields_out,field_count_out,pos)

  ice_frac => fields_out(pos)
  Call Rcf_Alloc_Field( ice_frac )
  Call Rcf_Read_Field( ice_frac, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Allocate space, and pointers.
!----------------------------------------------------------------------
! Set up temporary space for one of ice_conc_cat or ice_thick_cat and
! point the other to the space passed in as ice_cat_field.

  If ( STASH_item == stashcode_ice_conc_cat ) Then
    ! Allocate Category ice thickness; will abort if icethick not found
    Call Rcf_Locate(stashcode_prog_sec, stashcode_ice_thick_cat,     &
                    fields_out, field_count_out, pos)

    ice_cat_thick => fields_out(pos)
    Call Rcf_Alloc_Field( ice_cat_thick )

    ice_cat_frac => ice_cat_field
  Else
    ! Allocate Category ice fraction; will abort if icefrac not found
    Call Rcf_Locate(stashcode_prog_sec, stashcode_ice_conc_cat,      &
                    fields_out, field_count_out, pos)

    ice_cat_frac => fields_out(pos)
    Call Rcf_Alloc_Field( ice_cat_frac )

    ice_cat_thick => ice_cat_field
  End If


!----------------------------------------------------------------------
! Calculate some of the constants
!----------------------------------------------------------------------

  rNice = 1.0 / Real(Nice)
! calcluate the thickness category offset scalar
  CC1 = 3.0 * rNice
! calcluate the thickness category adjustment scalar
  CC2 = 15.0 * CC1
! CC3 is the thickness category phase adjustment scaler
  CC3 = 3.0

  ! use different minimum thickness in special case Nice=1
  If ( Nice == 1 ) Then
    HIN_Max(0) = 0.10
  Else
    HIN_Max(0) = 0.00
  End If

!----------------------------------------------------------------------
! Loop through ice_cat_field pseudo levels
!----------------------------------------------------------------------

  Do i = 1,Nice

    ! Calculate thickness category domain.
    X1 = Real ( i - 1 ) * rNice
    HIN_Max(i) = HIN_Max(i-1) + CC1 +                                 &
                 CC2 * (1.0 + TanH( CC3 * ( X1 - 1.0 ) ) )

    !Ensure upper limit for HIN_Max is an arbitrary LARGE value
    HIN_Max(Nice) = Huge(HIN_Max(0))

    ! Map ice thickness and fraction to thier psuedolevel equivalents.
    Where ( HIN_Max(i-1) <= ice_thick % data(:,1)  .AND. &
            ice_thick % data(:,1) < HIN_Max(i) )
      ice_cat_thick % data(:,i) = ice_thick % data(:,1)
      ice_cat_frac % data(:,i) = ice_frac % data(:,1)
    Elsewhere
      ice_cat_thick % data(:,i) = 0.00
      ice_cat_frac % data(:,i) = 0.00
    End Where

  EndDo

! Desired field has now been 'populated' as the relevant pointer
! was used to point to ice_cat_field earlier, and ice_cat_field is the
! field returned to the calling routine via the argument list.

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------

  Call Rcf_Dealloc_Field( ice_thick )
  Call Rcf_Dealloc_Field( ice_frac )

  If ( STASH_item == stashcode_ice_conc_cat ) Then
    Call Rcf_Dealloc_Field( ice_cat_thick )
  Else
    Call Rcf_Dealloc_Field( ice_cat_frac )
  End If

Return

End Subroutine Rcf_Derv_Ice_Cat_Thick

End Module Rcf_Derv_Ice_Cat_Thick_mod
