! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reinitialises frozen/unfrozen soil moisture

Module Rcf_Freeze_Soil_Mod

! Description:
!     This subroutine is the interface to FREEZE_SOIL which
!     reinitialises frozen/unfrozen soil moisture fields after
!     soil moisture/temperature has been updated from ancillary.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2   20/02/01   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Freeze_Soil( fields_out, field_count_out, hdr_out,  &
                            grid_out )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_soil_moist,      &
    stashcode_soil_temp,       &
    stashcode_vol_smc_sat,     &
    stashcode_soil_suction,    &
    stashcode_clapp_hb,        &
    stashcode_unfrozen_soil,   &
    stashcode_frozen_soil,     &
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

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Implicit None

! Arguments
Type( field_type ), Pointer       :: fields_out(:)
Type( um_header_type), Intent(In) :: hdr_out
Type( grid_type ), Intent(In)     :: grid_out
Integer, Intent(In)               :: field_count_out

! Internal variables
Type( field_type ), Pointer       ::  soil_moist
Type( field_type ), Pointer       ::  soil_temp
Type( field_type ), Pointer       ::  vol_smc_sat
Type( field_type ), Pointer       ::  soil_suction
Type( field_type ), Pointer       ::  clapp_hb
Type( field_type ), Pointer       ::  unfrozen_soil
Type( field_type ), Pointer       ::  frozen_soil

Integer                           ::  pos   ! position in array

! Externals
External Freeze_Soil


!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Reinitialising frozen/unfrozen soil moisture'
End If

!----------------------------------------------------------------------
! Find required fields in output dump and read them in
!----------------------------------------------------------------------
! Soil moisture
Call Rcf_Locate( stashcode_prog_sec, stashcode_soil_moist,           &
                 fields_out, field_count_out, pos)
soil_moist => fields_out(pos)
Call Rcf_Alloc_Field( soil_moist )
Call Rcf_Read_Field( soil_moist, hdr_out, decomp_rcf_output )

! Soil temperature
Call Rcf_Locate( stashcode_prog_sec, stashcode_soil_temp,            &
                 fields_out, field_count_out, pos)
soil_temp => fields_out(pos)
Call Rcf_Alloc_Field( soil_temp )
Call Rcf_Read_Field( soil_temp, hdr_out, decomp_rcf_output )

! Volume SMC at Saturation
Call Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,          &
                 fields_out, field_count_out,pos)
vol_smc_sat => fields_out(pos)
Call Rcf_Alloc_Field( vol_smc_sat )
Call Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )

! Saturated Soil Water Suction
Call Rcf_Locate( stashcode_prog_sec, stashcode_soil_suction,         &
                 fields_out, field_count_out, pos)
soil_suction => fields_out(pos)
Call Rcf_Alloc_Field( soil_suction )
Call Rcf_Read_Field( soil_suction , hdr_out, decomp_rcf_output )

! Clapp-Hornberger "B" Coefficient
Call Rcf_Locate( stashcode_prog_sec, stashcode_clapp_hb,             &
                 fields_out, field_count_out, pos)
clapp_hb => fields_out(pos)
Call Rcf_Alloc_Field( clapp_hb )
Call Rcf_Read_Field( clapp_hb , hdr_out, decomp_rcf_output )

! Unfrozen soil moisture
Call Rcf_Locate( stashcode_prog_sec, stashcode_unfrozen_soil,        &
                 fields_out, field_count_out, pos)
unfrozen_soil => fields_out(pos)
Call Rcf_Alloc_Field( unfrozen_soil )
Call Rcf_Read_Field( unfrozen_soil, hdr_out, decomp_rcf_output )

! Frozen soil moisture
Call Rcf_Locate( stashcode_prog_sec, stashcode_frozen_soil,          &
                 fields_out, field_count_out, pos)
frozen_soil => fields_out(pos)
Call Rcf_Alloc_Field( frozen_soil )
Call Rcf_Read_Field( frozen_soil, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Call Freeze_Soil to do the calculation
!----------------------------------------------------------------------
! DEPENDS ON: freeze_soil
Call Freeze_Soil( soil_temp % level_size, grid_out % sm_levels,    &
                  clapp_hb % Data,        grid_out % soil_depths,  &
                  soil_suction % Data,    soil_moist % Data,       &
                  soil_temp % Data,       vol_smc_sat % Data,      &
                  unfrozen_soil % Data,   frozen_soil % Data   )

!----------------------------------------------------------------------
! Write out corrected soil frozen/unfrozen fractions
!----------------------------------------------------------------------
Call Rcf_Write_Field( unfrozen_soil, hdr_out, decomp_rcf_output )
Call Rcf_Write_Field( frozen_soil, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------
Call Rcf_Dealloc_Field( frozen_soil )
Call Rcf_Dealloc_Field( unfrozen_soil )
Call Rcf_Dealloc_Field( clapp_hb )
Call Rcf_Dealloc_Field( soil_suction )
Call Rcf_Dealloc_Field( vol_smc_sat )
Call Rcf_Dealloc_Field( soil_temp )
Call Rcf_Dealloc_Field( soil_moist )


Return
End Subroutine Rcf_Freeze_Soil
End Module Rcf_Freeze_Soil_Mod
