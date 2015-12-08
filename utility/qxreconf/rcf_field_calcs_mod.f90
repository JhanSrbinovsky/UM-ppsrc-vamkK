! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Performs Source=8 field initialisation calculations

Module Rcf_Field_Calcs_Mod

!  Subroutine Rcf_Field_Calcs - field initialisation calculations.
!
! Description:
!   Some fields may have hard-coded Source=8 initialisation (certain
!   slab fields may be set by UMUI in this way also). These fields
!   are initialised by the calculations in this routine.
!
! Method:
!   Choice of method applied determined by stashcode.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Field_Calcs( fields_in, fields_out, field_count_in, &
                            field_count_out, data_source, hdr_in,  &
                            hdr_out )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use rcf_data_source_Mod, Only : &
    data_source_type,           &
    Field_Calcs

Use Submodel_Mod, Only : &
    Submodel_Ident,      &
    Atmos_IM

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_UMhead_Mod, Only : &
    um_header_type

USE decomp_params, ONLY : &
    decomp_rcf_output

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Derv_Ice_Temp_Mod, Only : &
    Rcf_Derv_Ice_Temp

Use Rcf_Derv_2D_CCA_Mod, Only : &
    Rcf_Derv_2D_CCA

Use Rcf_Derv_Ice_Thick_Mod, Only : &
    Rcf_Derv_Ice_Thick

Use Rcf_Derv_Sea_Ice_Temp_Mod, Only :&
    Rcf_Derv_Sea_Ice_Temp

USE PrintStatus_mod, Only : &
    LTimer

Use Rcf_Init_Field_On_Tiles_Mod, Only : &
    Rcf_Init_Field_On_Tiles

Use Rcf_Init_Canopy_Water_Mod, Only : &
    Rcf_Init_Canopy_Water

Use Rcf_Init_Tile_Snow_Mod, Only : &
    Rcf_Init_Tile_Snow

Use Rcf_Init_Tile_T_Mod, Only : &
    Rcf_Init_Tile_T

USE rcf_init_ml_snow_mod, ONLY : &
    rcf_init_ml_snow

Use Rcf_Init_Flake_Mod, Only : &
    Rcf_Init_Flake

Use Rcf_Calc_Gamtot_Mod, Only : &
    Rcf_Calc_Gamtot

Use Rcf_Est_Zw_Mod, Only :      &
    Rcf_Est_Zw

Use Rcf_Est_Sthzw_Mod, Only :      &
    Rcf_Est_Sthzw

Use Rcf_Fit_Fsat_Mod, Only : &
    Rcf_Fit_Fsat

Use Rcf_Calc_Fsat_Mod, Only :    &
    Rcf_Calc_Fsat

Use Rcf_Derv_Ice_Cat_Thick_Mod, Only :    &
    Rcf_Derv_Ice_Cat_Thick

Use Rcf_Derv_Adv_Winds_Mod, Only : &
    Rcf_Derv_Adv_Winds

Use Rcf_Change_Dust_Bins_Mod, Only : &
    Rcf_Change_Dust_Bins

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_prog_sec,                                   &
    stashcode_ice_conc_cat,    stashcode_ice_thick_cat,   &
    stashcode_ice_temp_cat,    stashcode_ice_snow_depth,  &
    stashcode_icethick,                                   &
    stashcode_cca,             stashcode_can_water_tile,  &
    stashcode_sea_ice_temp,                               &
    stashcode_catch_snow,      stashcode_snow_grnd,       &
    stashcode_rgrain,                                     &
    stashcode_snow_tile,       stashcode_tstar_tile,      &
    stashcode_gamtot,          stashcode_zw,              &
    stashcode_sthzw,           stashcode_fsat,            &
    stashcode_fwetl,                                      &
    stashcode_a_fsat,          stashcode_c_fsat,          &
    stashcode_a_fwet,          stashcode_c_fwet,          &
    stashcode_u_adv,           stashcode_v_adv,           &
    stashcode_snow_T_tile,     stashcode_snowdep_grd_tile, &
    stashcode_snow_ice_tile,   stashcode_snow_liq_tile,   &
    stashcode_snow_laydns_tiles,                          &
    stashcode_snowpack_bk_dens,                           &
    stashcode_nsnow_layrs_tiles,                          &
    stashcode_snow_laythk_tiles,                          &
    stashcode_dust1_mmr, stashcode_dust2_mmr,             &
    stashcode_dust3_mmr, stashcode_dust4_mmr,             &
    stashcode_dust5_mmr, stashcode_dust6_mmr,             &
    stashcode_flake_t_mean,   stashcode_flake_t_mxl,      &
    stashcode_flake_t_ice,    stashcode_flake_h_mxl,      &
    stashcode_flake_h_ice

IMPLICIT NONE

! Arguments
Integer,            Intent(In)       :: field_count_in
Integer,            Intent(In)       :: field_count_out
Type( field_type ), Pointer          :: fields_in( : )
Type( field_type ), Pointer          :: fields_out( : )
Type( data_source_type ), Pointer    :: data_source( : )
Type( um_header_type ), Intent(In)   :: hdr_in
Type( um_header_type ), Intent(In)   :: hdr_out

! Local variables
Integer                           :: i
Integer                           :: ErrorStatus
Character (Len=*), Parameter      :: RoutineName='Rcf_Field_Calcs'
Character (Len=80)                :: Cmessage

External Timer

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

!-----------------------------------------------------------------
! Loop around output fields/sources until a source=field_calcs
! is found
!-----------------------------------------------------------------
Do i = 1, field_count_out
  If ( data_source( i ) % Source == Field_Calcs ) Then

    Call Rcf_Alloc_Field( fields_out( i ) )

    Select Case( fields_out( i ) % stashmaster % model )

    Case ( Atmos_IM )

      ! -----------
      ! Atmos items
      ! -----------

      Select Case( fields_out( i ) % stashmaster % section )

      ! For the moment only section zero fields can be used here.
      Case ( stashcode_prog_sec )

        Select Case( fields_out(i) % stashmaster % item )

        Case( stashcode_cca )
          Call Rcf_Derv_2D_CCA( fields_in, fields_out, field_count_in, &
                                field_count_out, hdr_in, hdr_out ,     &
                                fields_out( i ) )

        Case( stashcode_icethick )
          Call Rcf_Derv_Ice_Thick( fields_out, field_count_out,        &
                                   hdr_out, fields_out(i) )

        Case( stashcode_sea_ice_temp )
          Call Rcf_Derv_Sea_Ice_Temp( fields_out, field_count_out,  &
                                      fields_out( i ), hdr_out )

        Case( stashcode_catch_snow       &
             ,stashcode_snow_grnd        &
             ,stashcode_rgrain             )
          Call Rcf_Init_Field_On_Tiles( fields_in, field_count_in,       &
                                      hdr_in, fields_out( i ),           &
                                      fields_out(i) % stashmaster % item )

        CASE( stashcode_snowdep_grd_tile  &
             ,stashcode_snowpack_bk_dens  &
             ,stashcode_nsnow_layrs_tiles &
             ,stashcode_snow_laythk_tiles &
             ,stashcode_snow_ice_tile     &
             ,stashcode_snow_liq_tile     &
             ,stashcode_snow_t_tile       &
             ,stashcode_snow_laydns_tiles   )
          CALL rcf_init_ml_snow( fields_out, field_count_out,       &
                                 hdr_out, fields_out( i ),          &
                                 fields_out(i) % stashmaster % item )

        Case( stashcode_flake_t_mean     &
             ,stashcode_flake_t_mxl      &
             ,stashcode_flake_t_ice      &
             ,stashcode_flake_h_mxl      &
             ,stashcode_flake_h_ice        )
          Call Rcf_Init_Flake( fields_out, field_count_out,       &
                               hdr_out, fields_out( i ),          &
                               fields_out(i) % stashmaster % item )

        Case( stashcode_can_water_tile )
          Call Rcf_Init_Canopy_Water( fields_in, field_count_in,       &
                                      hdr_in, fields_out( i ) )

        Case( stashcode_snow_tile )
          Call Rcf_Init_Tile_Snow( fields_out, field_count_out,        &
                                   hdr_out, fields_out( i ) )

        Case( stashcode_ice_temp_cat )
          Call Rcf_Derv_Ice_Temp( fields_out,field_count_out, hdr_out, &
                                 fields_out( i) )

        Case( stashcode_tstar_tile )
          Call Rcf_Init_Tile_T( fields_out, field_count_out, hdr_out,  &
                                fields_out( i ) )

! NB The fields for LSH have quite a few inter-dependencies.  It seems
! that due to the dump being in STASHcode order the dependencies should be okay
! but please check code that if another field require calculations that it isnt
! used to calculate another field unless you are really sure.
        Case( stashcode_Gamtot )
          Call Rcf_Calc_Gamtot( fields_out, field_count_out,  &
                                    i, hdr_out )
! If calculation required for the LSH fitting variables,
! must first allocate space for them.

        Case( stashcode_a_fsat, stashcode_c_fsat, &
               stashcode_a_fwet, stashcode_c_fwet)

          Call Rcf_Fit_Fsat ( fields_out, field_count_out,      &
                                   i, hdr_out )

        Case( stashcode_Sthzw )
          Call Rcf_Est_Sthzw( fields_out, field_count_out,  &
                                   fields_out( i ), hdr_out )

        Case( stashcode_Zw )
          Call Rcf_Est_Zw( fields_out, field_count_out,  &
                                   fields_out( i ), hdr_out )

        Case( stashcode_Fsat, stashcode_fwetl )
          
          Call Rcf_Calc_Fsat( fields_out, field_count_out,  &
                              i, hdr_out )

        Case( stashcode_ice_conc_cat, stashcode_ice_thick_cat )
          Call Rcf_Derv_Ice_Cat_Thick (                               &
                                fields_out( i ) % stashmaster % item, &
                                fields_out, field_count_out,          &
                                hdr_out, fields_out( i ) )

        Case( stashcode_u_adv, stashcode_v_adv )
          Call Rcf_Derv_Adv_Winds(                                   &
                               fields_out( i ) % stashmaster % item, &
                               fields_out, field_count_out,          &
                               hdr_out, fields_out( i ) )

        Case( stashcode_dust1_mmr, stashcode_dust2_mmr,                        &
              stashcode_dust3_mmr, stashcode_dust4_mmr,                        &
              stashcode_dust5_mmr, stashcode_dust6_mmr )
          Call Rcf_Change_Dust_Bins( fields_in, field_count_in,                &
                                     fields_out, field_count_out,              &
                                     fields_out( i ) % stashmaster % item,     &
                                     hdr_in, hdr_out, fields_out( i ) )

        Case Default
          ErrorStatus = 30
          Write (Cmessage, '(A, I3, A, I5)')                         &
                'No Field Calculations specified for section ',      &
                fields_out( i ) % stashmaster % section, ' item ',   &
                fields_out( i ) % stashmaster % item
          Call Ereport( RoutineName, ErrorStatus, Cmessage )

        End Select

      Case Default
        ErrorStatus = 30
        Write (Cmessage, '(A, I3)')                                  &
              'No Field Calculations specified for section ',        &
              fields_out( i ) % stashmaster % section
        Call Ereport( RoutineName, ErrorStatus, Cmessage )

      End Select

    Case Default                          ! Selection on Internal Model
      ErrorStatus = 40
      Write (Cmessage, '(A, I2)')                                    &
            'No Field Calculations specified for Submodel ',         &
            fields_out( i ) % stashmaster % model
      Call Ereport( RoutineName, ErrorStatus, Cmessage )

    End Select                            ! Selection on Internal Model
!-----------------------------------------------------------------
! Clean up and write out
!-----------------------------------------------------------------
    Call Rcf_Write_Field( fields_out( i ), Hdr_out, decomp_rcf_output )
    Call Rcf_DeAlloc_Field( fields_out( i ) )
  End If

End Do

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

Return
End Subroutine Rcf_Field_Calcs
End Module Rcf_Field_Calcs_Mod
