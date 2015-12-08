! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Covert soil stress to moisture if soil stress was used in the 
!  interpolation in the place of soil moisture

Module Rcf_Soilstress_to_soilmoist_Mod

!  Subroutine Rcf_Soilstress_to_soilmoist
!
! Description:
!   Convert soil stress back to soil moisture after interpolation.
!
! Method:
!    1. Read soil stress, wilt, sat and cri from output dump
!    2. convert soil stress to Vol soil:
!         soil_vol=soilstress*(soil_cri-soil_wilt)+soil_wilt
!    3. Ensure SMC_vol is less or equal to soil_sat:
!         soil_vol=min(soil_vol, soil_sat)
!    4. Convert soil_vol to  soilmoist:
!         soilmoist=soil_vol*rhow*delta_z  ! delta_z is layer depth
!    5. Check if soil moisture is not between 0 and sat.
!         if it is less than 0, the set it to soil_min and
!         if it is greater than sat, then set to soil_max
!         (see soil_min, soil_max in the codes)
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Soilstress_to_soilmoist( fields, field_count, &
           grid, decomp, hdr)

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
    stashcode_vol_smc_cri,     &
    stashcode_vol_smc_sat,     &
    stashcode_soil_moist,      &
    stashcode_prog_sec

USE UM_ParVars, Only : &
    mype

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal


USE water_constants_mod, ONLY: rho_water

Implicit None

! Arguments
Type( field_type ), Pointer   :: fields(:)
Type( grid_type )             :: grid
Type( um_header_type )        :: hdr
Integer                       :: field_count
Integer                       :: decomp

! Local variables
Integer                            :: pos_sat
Integer                            :: pos_wilt
Integer                            :: pos_smc
Integer                            :: pos_cri
Integer                            :: i
Integer                            :: j

Type( field_type ), Pointer        :: soil_moist
Type( field_type ), Pointer        :: wilt
Type( field_type ), Pointer        :: sat
Type( field_type ), Pointer        :: cri

Real                               :: smc_min
Real                               :: smc_max

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
  Write (6,*) 'Converting soil stress to soil moisture...'
End If

!-------------------------------------------------------------------
! Locate required peramaters in output dump
!-------------------------------------------------------------------

Call Rcf_Locate &
    ( stashcode_prog_sec,stashcode_vol_smc_sat,  &
     fields, field_count, pos_sat,  .TRUE. )
Call Rcf_Locate &
    ( stashcode_prog_sec,stashcode_vol_smc_wilt, &
      fields, field_count, pos_wilt, .TRUE. )
Call Rcf_Locate &
    ( stashcode_prog_sec,stashcode_soil_moist,   &
      fields, field_count, pos_smc,  .TRUE. )
Call Rcf_Locate &
    ( stashcode_prog_sec,stashcode_vol_smc_cri,  &
      fields, field_count, pos_cri,  .TRUE. )

!--------------------------------------------------------------------
! Reset smc if appropriate
!--------------------------------------------------------------------

If (pos_sat  /= 0  .AND. &
    pos_wilt /= 0 .AND. &
    pos_cri  /= 0 .AND. &
    pos_smc  /= 0 ) Then

  soil_moist => fields( pos_smc )
  wilt       => fields( pos_wilt)
  sat        => fields( pos_sat )
  cri        => fields( pos_cri )

 
  Call Rcf_Alloc_Field( soil_moist )
  Call Rcf_Alloc_Field( wilt       )
  Call Rcf_Alloc_Field( sat        )
  Call Rcf_Alloc_Field( cri        )

  Call Rcf_Read_Field( soil_moist, hdr, decomp )
  Call Rcf_Read_Field( wilt ,      hdr, decomp )
  Call Rcf_Read_Field( sat  ,      hdr, decomp )
  Call Rcf_Read_Field( cri  ,      hdr, decomp )


  ! convert smc_stress to VOL SMC
  Do j = 1, soil_moist % levels
     soil_moist % data(:,j) = soil_moist % data(:,j)* &
     (cri % data(:,1)-wilt % data(:,1))+wilt % data(:,1)
         
  End do

  ! check saturation 
  Do j = 1, soil_moist % levels
     Do i=1,soil_moist % level_size
        If (soil_moist % data(i,j) .gt. sat % data(i,1)) Then
          soil_moist % data(i,j) = sat % data (i,1)
        End If
     End do   
  End do

  ! calculate smc
  Do j = 1, soil_moist % levels
     soil_moist % data(:,j) = soil_moist % data(:,j)* &
      Grid % soil_depths(j)*RHO_WATER               
  End do



  If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
    write (6,*)'Soil Moisture has been converted from soil stress.'
                 
  End If
 

  !---------------------------------------------------------------------
  ! Write out converted field
  !---------------------------------------------------------------------

  Call Rcf_Write_Field( soil_moist, hdr, decomp )

  Call Rcf_Dealloc_Field( soil_moist )
  Call Rcf_Dealloc_Field( wilt       )
  Call Rcf_Dealloc_Field( sat        )
  Call Rcf_Dealloc_Field( cri        )

End If

Return
End Subroutine Rcf_Soilstress_to_soilmoist
End Module Rcf_Soilstress_to_soilmoist_Mod

