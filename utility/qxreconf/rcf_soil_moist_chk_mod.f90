! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Checks that soil moisture is appropriately set for vegetation

Module Rcf_Soil_Moist_Chk_Mod

!  Subroutine Rcf_Soil_Moist_Chk
!
! Description:
!   Ensures that soil moisture is consistent with vegetation fields.
!
! Method:
!    Set Soil Moisture Content (SMC) = 0 for land-ice points
!    Ensure SMC is 0.1 * volumetric wilt as minimum at non-land-ice points
!    (0.1 is fairly arbitary but matches the nudging scheme - see
!     SURF documentation)
!    Ensure SMC is saturation as maximum at non-land-ice points
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Soil_Moist_Chk( fields, field_count, grid, decomp, hdr, &
                                l_changed )

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
    stashcode_vol_smc_wilt,    &
    stashcode_vol_smc_sat,     &
    stashcode_soil_moist,      &
    stashcode_prog_sec

USE UM_ParVars, Only : &
    nproc,                  &
    mype

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

USE water_constants_mod, ONLY: rho_water

Implicit None

! Arguments
Type( field_type ), Pointer        :: fields(:)
Type( grid_type ), Intent(In)      :: grid
Type( um_header_type ), Intent(In) :: hdr
Integer, Intent(In)                :: field_count
Integer, Intent(In)                :: decomp
Logical, Intent(  Out)             :: l_changed

! Local variables
Integer                            :: pos_sat
Integer                            :: pos_wilt
Integer                            :: pos_smc
Integer                            :: i
Integer                            :: j
Integer                            :: istat
Integer                            :: count

Real                               :: smc_min
Real                               :: smc_max

Type( field_type ), Pointer        :: soil_moist
Type( field_type ), Pointer        :: wilt
Type( field_type ), Pointer        :: sat


!-------------------------------------------------------------------
! Locate required paramaters in output dump
!-------------------------------------------------------------------

Call Rcf_Locate ( stashcode_prog_sec, stashcode_vol_smc_sat,  &
                  fields, field_count, pos_sat,  .TRUE. )
Call Rcf_Locate ( stashcode_prog_sec, stashcode_vol_smc_wilt, &
                  fields, field_count, pos_wilt, .TRUE. )
Call Rcf_Locate ( stashcode_prog_sec, stashcode_soil_moist,   &
                  fields, field_count, pos_smc,  .TRUE. )

!--------------------------------------------------------------------
! Reset smc if appropriate
!--------------------------------------------------------------------

If (pos_sat  /= 0 .AND. &
    pos_wilt /= 0 .AND. &
    pos_smc  /= 0 ) Then

  soil_moist => fields( pos_smc  )
  wilt       => fields( pos_wilt )
  sat        => fields( pos_sat  )

  Call Rcf_Alloc_Field( soil_moist )
  Call Rcf_Alloc_Field( wilt       )
  Call Rcf_Alloc_Field( sat        )

  Call Rcf_Read_Field( soil_moist, hdr, decomp )
  Call Rcf_Read_Field( wilt,       hdr, decomp )
  Call Rcf_Read_Field( sat,        hdr, decomp )

  count = 0
  Do i = 1, soil_moist % level_size

    If ( wilt % data(i,1) > 0.0 ) Then
      ! Non-ice land point
      ! Make sure soil_moist lies between appropriate bounds

      Do j = 1, soil_moist % levels
        ! min is 0.1 * wilt in absolute (not volumetric) terms
        smc_min = 0.1 * wilt % data(i,1) * grid % soil_depths(j) &
                      * rho_water

        ! max is saturation in absolute terms
        smc_max = sat % data(i,1) * grid % soil_depths(j) * rho_water

        If ( soil_moist % data(i,j) < smc_min) Then
          ! Too low - set to the minimum value
          soil_moist % data(i,j) = smc_min
          count = count + 1

        Else If ( soil_moist % data(i,j) > smc_max) Then
          ! Too high - set to the maximum value
          soil_moist % data(i,j) = smc_max
          count = count + 1

        End If
     End Do

    Else If (wilt % data(i,1) == 0.0 ) Then
      ! land-ice point
      ! soil_moist should be set to 0.0 if it isn't
      Do j = 1, soil_moist % levels
        If ( soil_moist % data(i,j) /= 0.0 ) Then
          soil_moist % data(i,j) = 0.0
          count = count + 1
        End If
      End Do
    End If

  End Do
      
  Call gc_isum (1, nproc, istat, count)
  If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
    write (6,*) 'Soil Moisture: No of values reset ', &
                 count
  End If


  ! Set the changed flag to true if count is positive
  If (count > 0) Then
    l_changed = .TRUE.
  Else
    l_changed = .FALSE.
  End If

!---------------------------------------------------------------------
! Write out changed field
!---------------------------------------------------------------------

  If (l_changed) Then
    Call Rcf_Write_Field( soil_moist, hdr, decomp )
  End If

  Call Rcf_Dealloc_Field( soil_moist )
  Call Rcf_Dealloc_Field( wilt       )
  Call Rcf_Dealloc_Field( sat        )
End If

Return
End Subroutine Rcf_Soil_Moist_Chk
End Module Rcf_Soil_Moist_Chk_Mod
