! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialisations for FLake lake model prognostics.

Module Rcf_Init_Flake_Mod

IMPLICIT NONE

! Subroutine Rcf_Init_Flake
!
! Description:
!   Initialises FLake model prognostics
!
! Method:
!   Based on land conditions of surface and soil temperatures,
!   and presence of snow.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Reconfiguration

Contains

Subroutine Rcf_Init_Flake( fields_out,  field_count_out, hdr_out, &
                           flake_field, field_stashcode )


Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

USE decomp_params, ONLY : &
    decomp_rcf_output

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_prog_sec,        &
    stashcode_tstar_tile,      &
    stashcode_snow_tile,       &
    stashcode_soil_temp,       &
    stashcode_flake_depth,     &
    stashcode_flake_t_mean,    &
    stashcode_flake_t_mxl,     &
    stashcode_flake_t_ice,     &
    stashcode_flake_h_mxl,     &
    stashcode_flake_h_ice

Use PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use UM_Parvars, Only : &
    mype

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_Lsm_Mod, Only : &
    local_lsm_out

Use c_0_dg_c, Only : &
    tm

Use nstypes, Only : &
    lake

  USE lake_mod, Only :     &
         h_snow_min_flk    &
        ,h_ice_min_flk     &
        ,h_ice_max         &
        ,lake_h_mxl_0

Implicit None

! Arguments
! Arguments
Type( field_type ), Pointer          :: fields_out(:)
Type( field_type ), Intent( InOut )  :: flake_field
Type( um_header_type ), Intent( In ) :: hdr_out

Integer, Intent( In )                :: field_count_out
Integer, Intent( In )                :: field_stashcode

! Local variables
Type( field_type ), Pointer          :: tstar_tile
Type( field_type ), Pointer          :: snow_tile
Type( field_type ), Pointer          :: soil_temp
Type( field_type ), Pointer          :: lake_depth

Real                                 :: lake_t_mean &
                                       (flake_field % level_size)
Real                                 :: lake_t_mxl  &
                                       (flake_field % level_size)
Real                                 :: lake_t_ice  &
                                       (flake_field % level_size)
Real                                 :: lake_h_mxl  &
                                       (flake_field % level_size)
Real                                 :: lake_h_ice  &
                                       (flake_field % level_size)

Real                                 :: soil_mean_temp
Integer                              :: size
Integer                              :: pos  ! field position
Integer                              :: i    ! looper
Integer                              :: j    ! looper
Integer                              :: l    ! looper

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Initialising FLake variable ',field_stashcode
End If

!----------------------------------------------------------------------
! Find tile T* temperature in output fields and read it in
!----------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_tstar_tile,         &
                 fields_out, field_count_out, pos)
tstar_tile => fields_out(pos)
Call Rcf_Alloc_Field( tstar_tile )
Call Rcf_Read_Field( tstar_tile, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Find tile snow in output fields and read it in
!----------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_snow_tile,         &
                 fields_out, field_count_out, pos)
snow_tile => fields_out(pos)
Call Rcf_Alloc_Field( snow_tile )
Call Rcf_Read_Field( snow_tile, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Find soil temperature in output fields and read it in
!----------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_soil_temp,         &
                 fields_out, field_count_out, pos)
soil_temp => fields_out(pos)
Call Rcf_Alloc_Field( soil_temp )
Call Rcf_Read_Field( soil_temp, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Find lake depth in output fields and read it in
!----------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_flake_depth,        &
                 fields_out, field_count_out, pos)
lake_depth => fields_out(pos)
Call Rcf_Alloc_Field( lake_depth )
Call Rcf_Read_Field( lake_depth, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Perform FLake init calculations
!----------------------------------------------------------------------
Do l = 1, flake_field % level_size

! calculate mean soil temperature
!
  soil_mean_temp = 0.0
  Do i = 1, Output_Grid % sm_levels
    soil_mean_temp = soil_mean_temp + soil_temp % Data(l,i) * Output_Grid % soil_depths(i)
  End Do
  soil_mean_temp= soil_mean_temp / SUM( Output_Grid % soil_depths )

!----------------------------------------------------------------------
! Decide which calculations are necessary
!----------------------------------------------------------------------
Select Case( field_stashcode )

 Case (stashcode_flake_t_mxl,stashcode_flake_t_mean)

! set mixed-layer T based on the temperature of the 1st soil level
        lake_t_mxl(l) = MAX(soil_temp % Data(l,1), tm + 0.10)

! set mean T to the mean temperature over the soil column,
! BUT constrain to between the mixed-layer temperature and freezing.
        lake_t_mean(l) = MIN(soil_mean_temp,  &
                                lake_t_mxl(l) - 0.05)
        lake_t_mean(l) = MAX(lake_t_mean(l), tm + 0.05)


 Case (stashcode_flake_h_ice,stashcode_flake_h_mxl)

! initialise the ice thickness
        IF(soil_mean_temp                    <= tm) THEN
          lake_h_ice(l) =   SUM( Output_Grid % soil_depths ) &
              * MIN(1.0,(tm - soil_mean_temp         + 0.05)/25.0)
        ELSE IF (soil_temp % Data(l, 1)      <= tm) THEN
          lake_h_ice(l) = Output_Grid % soil_depths(1) &
              * MIN(1.0,(tm - soil_temp % Data(l, 1) + 0.05)/25.0)
        ELSE IF (tstar_tile % Data(l, lake ) <= tm) THEN
          lake_h_ice(l) = 2.0 * h_ice_min_flk
        ELSE
          lake_h_ice(l) = 0.0
        END IF
! ...and bound by the Mironov & Ritter (2004) maximum
!    to avoid SQRT(-ve) in FLake_driver
        lake_h_ice(l) = MIN( lake_h_ice(l),h_ice_max )

! bound the mixed-layer depth by the available unfrozen depth
        lake_h_mxl(l) = MIN(lake_h_mxl_0,  &
                               lake_depth % Data(l,1)-lake_h_ice(l))
        lake_h_mxl(l) = MAX(lake_h_mxl(l),0.0)


 Case (stashcode_flake_t_ice)

! set the ice upper-boundary temperature depending on the
! presence of snow
        IF(snow_tile % Data(l, lake ) <= h_snow_min_flk)THEN
          lake_t_ice(l) = tstar_tile % Data(l, lake)
        ELSE
          lake_t_ice(l) = soil_temp % Data(l, 1)
        END IF

! put an upper bound on the ice T
        lake_t_ice(l) = MIN(lake_t_ice(l),tm - 0.05)

End Select


!----------------------------------------------------------------------
! write the correct data to the output field
!----------------------------------------------------------------------
Select Case( field_stashcode )

 Case (stashcode_flake_t_mxl)

   flake_field % Data(l,1) = lake_t_mxl(l)

 Case (stashcode_flake_t_mean)

   flake_field % Data(l,1) = lake_t_mean(l)

 Case (stashcode_flake_h_ice)

   flake_field % Data(l,1) = lake_h_ice(l)

 Case (stashcode_flake_h_mxl)

   flake_field % Data(l,1) = lake_h_mxl(l)

 Case (stashcode_flake_t_ice)

   flake_field % Data(l,1) = lake_t_ice(l)

End Select

    END DO


End Subroutine Rcf_Init_Flake

End Module Rcf_Init_Flake_Mod
