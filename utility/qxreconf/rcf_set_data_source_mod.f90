! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets the data_source array corresponding to the fields array.

MODULE Rcf_Set_Data_Source_Mod

IMPLICIT NONE

!  Subroutine Rcf_Set_Data_Source - sets the data_source array.
!
! Description:
!   Allocates the data_source array and sets the elements according
!   to how the data in the the associated fields array should be
!   initialised.
!
! Method:
!   Uses data from the ITEMS namelist.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

  SUBROUTINE Rcf_Set_data_source( data_source, fields_in, fields_out,          &
      field_count_in, field_count_out,                                         &
      hdr_in, hdr_out )

    USE Ereport_Mod, ONLY :                                                    &
        Ereport

    USE Rcf_Items_Mod             ! All of it

    USE Rcf_Data_Source_Mod   ! Most of It

    USE Rcf_Field_Type_Mod, ONLY :                                             &
        Field_type

    USE Rcf_Interp_Weights_Mod, ONLY :                                         &
        h_int_method,                                                          &
        nearest_neighbour,                                                     &
        h_int_active

    USE UM_ParVars, ONLY :                                                     &
        mype

    USE PrintStatus_mod, ONLY :                                                &
        PrintStatus,                                                           &
        PrStatus_Normal,                                                       &
        PrStatus_Oper

    USE Submodel_Mod, ONLY :                                                   &
        Atmos_IM

    USE Rcf_Locate_mod, ONLY :                                                 &
        Rcf_Locate

    USE Rcf_UMhead_Mod, ONLY :                                                 &
        UM_Header_type

    USE Rcf_Recon_Mod, ONLY :                                                  &
        grib,                                                                  &
        grib2ff,                                                               &
        l_interp_input_only

    USE nlsizes_namelist_mod, ONLY :                                           &
        nice

    USE Rcf_Stashcodes_Mod   ! So many of these are used the whole lot
                         ! have been included. Anything starting
                         ! with stashcode_ came from here.

    USE rad_param, ONLY :                                                      &
        r0                   !  Grain size for fresh snow (microns)

    USE snow_param, ONLY :                                                     &
        rho_snow_const       ! density of lying snow (kg per m**3)

    USE lake_mod, ONLY :                                                       &
        lake_depth_0,                                                          &
        lake_fetch_0,                                                          &
        lake_shape_0,                                                          &
        g_dt_0

    USE rcf_lsh_field_checks_mod, ONLY :                                       &
        rcf_lsh_field_checks

    USE Lookup_addresses

    USE cppxref_mod, ONLY :                                                    &
        ppx_atm_lbc_theta,                                                     &
        ppx_atm_lbc_u,                                                         &
        ppx_atm_lbc_v

    IMPLICIT NONE

! Arguments
    TYPE( data_source_type ), POINTER  :: data_source(:)
    TYPE( field_type), POINTER         :: fields_in(:)
    TYPE( field_type ) , POINTER       :: fields_out(:)
    TYPE( um_header_type ), INTENT(IN) :: hdr_in
    TYPE( um_header_type ), INTENT(IN) :: hdr_out
    INTEGER                            :: field_count_in
    INTEGER                            :: field_count_out

! Local variables
    INTEGER                            :: item_index
    INTEGER                            :: i
    INTEGER                            :: j
    INTEGER                            :: pos
    INTEGER                            :: pos_out
    INTEGER                            :: ErrorStatus
    CHARACTER (LEN=*), PARAMETER       :: RoutineName='Rcf_Set_Data_Source'
    CHARACTER (LEN=80)                 :: Cmessage

    CHARACTER (LEN=*), PARAMETER    :: form                                    &
        ="(2(i5,' '),4(i4,' '),e10.4,' ',50a)"
    CHARACTER (LEN=*), PARAMETER    :: form_c=                                 &
        "(2(a5,' '),4(a4,' '),a10,  ' ',50a)"


!---------------------------------------------------------------
! check some of the arrays needed have been set up
!---------------------------------------------------------------

    IF ( ASSOCIATED( data_source ) ) THEN
      ErrorStatus = -10
      Cmessage = 'data_source is already set up!'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

    IF (.NOT. ASSOCIATED( fields_out )  )THEN
      ErrorSTatus = 20
      Cmessage = 'Fields have not yet been set up'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

!---------------------------------------------------------------
! Allocate the space
!---------------------------------------------------------------
    ALLOCATE( data_source( field_count_out ) )

!---------------------------------------------------------------
! Find any approriate item from namelist and set values
!----------------------------------------------------------------
    DO i = 1, field_count_out
      item_index = 0
      DO j = 1, num_items
        IF ( fields_out( i ) % stashmaster % item == item_array( j ) .AND.     &
            fields_out( i ) % stashmaster % section == sctn_array( j ) ) THEN
          item_index = j
          EXIT
        END IF
      END DO

      IF ( item_index > 0 .AND. h_int_method /= nearest_neighbour ) THEN
    ! Have found the item - fill in gaps
        data_source( i ) % source      = source_array( item_index )
        data_source( i ) % Domain      = area_array  ( item_index )
        data_source( i ) % Ancil_SctnC = upas_array  ( item_index )
        data_source( i ) % Ancil_ItemC = upaa_array  ( item_index )
        data_source( i ) % RConst      = uprc_array  ( item_index )
        data_source( i ) % Ancil_File  = upaf_array  ( item_index )

      ELSE                             ! No item found - use defaults

        data_source( i ) % source      = Input_Dump
        data_source( i ) % Domain      = Whole_Grid
        data_source( i ) % Ancil_SctnC = 0
        data_source( i ) % Ancil_ItemC = 0
        data_source( i ) % RConst      = 0.0
        data_source( i ) % Ancil_File  = ' '

      END IF

 ! Setting Source for items required from input dump but not
  ! available (and known about with relevant code written)
      CALL Rcf_Locate( fields_out( i ) % stashmaster % section,                &
          fields_out( i ) % stashmaster % item,                                &
          fields_in, field_count_in, pos, .TRUE.)

      IF (l_interp_input_only) THEN
        IF (pos == 0) THEN
          WRITE(6,'(A,I4,I3)') "Error processing STASHcode: ", &
            fields_out( i ) % stashmaster % section,     &
            fields_out( i ) % stashmaster % item
          ErrorSTatus = 20
          Cmessage = 'Should only be interpolating input fields but not found.'
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )
        ELSE
          IF (data_source ( i ) % source /= Input_Dump ) THEN
            data_source( i ) % source = Input_Dump
          END IF
        END IF
      END IF

      IF ( pos == 0 .AND. data_source( i ) % source == Input_Dump ) THEN

  ! ----------------------------------
  ! This item is not in the input dump
  ! ----------------------------------

        SELECT CASE( fields_out( i ) % stashmaster % model )

        CASE ( Atmos_IM )                 ! Atmosphere Sub-Model

     ! ----------------
     ! Atmosphere items
     ! ----------------

          SELECT CASE( fields_out( i ) % stashmaster % section )

          CASE( stashcode_prog_sec )

            SELECT CASE( fields_out( i ) % stashmaster % item )

            CASE( stashcode_3d_cca )           ! 3D Convective cloud

              data_source( i ) % source = set_to_zero

            CASE( stashcode_cca )              ! 2D Convective cloud

              data_source( i ) % source = Field_Calcs

            CASE( stashcode_npp_pft_acc,     & ! Carbon accumulation fields
                stashcode_g_lf_pft_acc,                                        &
                stashcode_g_ph_lf_pft_acc,                                     &
                stashcode_rsp_w_pft_acc,                                       &
                stashcode_rsp_s_acc,                                           &
                stashcode_catch_snow,      & ! NLT canopy snow capacity
                stashcode_catch_tile,      & ! Tiled canopy capacity
                stashcode_z0_tile,         & ! Tiled roughness length
                stashcode_z0h_tile         & ! Tiled thermal roughness length
                )
              data_source( i ) % source = set_to_const
              data_source( i ) % rconst = -1.0

            CASE( stashcode_can_water_tile,  & ! Vegetation fields
                stashcode_tstar_tile,                                          &
                stashcode_ice_temp_cat,                                        &
                stashcode_ice_conc_cat,    & ! ice conc (fraction) cats
                stashcode_ice_thick_cat,   & ! ice thickness categories
                stashcode_Gamtot,          & !\ large-scale hydrology
                stashcode_Zw,              & !/ fields
                stashcode_Sthzw,           & !/
                stashcode_Fsat,            & !\.
                stashcode_Fwetl,           & !/.
                stashcode_a_fsat,         & !\.
                stashcode_c_fsat,         & !/.
                stashcode_a_fwet,         & !\.
                stashcode_c_fwet,         & !/.
                stashcode_snow_tile,                                           &
                stashcode_snowdep_grd_tile,  & ! JULES snowdepth on tiles
                stashcode_snowpack_bk_dens,  & ! JULES snowpack bulk density
                stashcode_nsnow_layrs_tiles, & ! JULES number of snow layers
                stashcode_snow_laythk_tiles, & ! JULES snow layer thicknesses
                stashcode_snow_ice_tile,     & ! JULES snow layer solid mass
                stashcode_snow_liq_tile,     & ! JULES snow layer liquid mass
                stashcode_snow_T_tile,       & ! JULES snow layer temps.
                stashcode_snow_laydns_tiles, & ! JULES snow layer density
                stashcode_flake_t_mean,     & ! FLake prognostics
                stashcode_flake_t_mxl,      & ! "
                stashcode_flake_t_ice,      & ! "
                stashcode_flake_h_mxl,      & ! "
                stashcode_flake_h_ice       & ! "
                )

              data_source( i ) % source = Field_Calcs

            CASE( stashcode_flake_depth )  ! FLake lake depth

              data_source( i ) % source = set_to_const
              data_source( i ) % rconst = lake_depth_0

            CASE( stashcode_flake_fetch )  ! FLake lake fetch

              data_source( i ) % source = set_to_const
              data_source( i ) % rconst = lake_fetch_0

            CASE( stashcode_flake_shape )  ! FLake thermocline shape factor

              data_source( i ) % source = set_to_const
              data_source( i ) % rconst = lake_shape_0

           CASE( stashcode_flake_g_over_dt ) ! (ht flux/delta Tsurf from FLake)

              data_source( i ) % source = set_to_const
              data_source( i ) % rconst = g_dt_0

            CASE( stashcode_fexp          & ! LSH conductivity decay scale
                )
              data_source( i ) % source = set_to_const
              data_source( i ) % rconst = 1.0

            CASE( stashcode_infil_max_tile,                                    &
                stashcode_snow_grnd,       & ! Snow beneath NLT canopy
                stashcode_snow_on_ice,                                         &
                stashcode_ice_snow_depth,                                      &
                stashcode_sw_down_tile,                                        &
                stashcode_sw_down,                                             &
                stashcode_lw_up_diff,                                          &
                stashcode_mmr_smoke_fr,                                        &
                stashcode_mmr_smoke_ag,                                        &
                stashcode_mmr_smoke_cl,                                        &
                stashcode_biom_surf_em,                                        &
                stashcode_biom_elev_em,                                        &
                stashcode_dust1_mmr,                                           &
                stashcode_dust2_mmr,                                           &
                stashcode_dust3_mmr,                                           &
                stashcode_dust4_mmr,                                           &
                stashcode_dust5_mmr,                                           &
                stashcode_dust6_mmr,                                           &
                stashcode_ddmfx,                                               &
                stashcode_albedo_sice,                                         &
                stashcode_albedo_land,                                         &
                stashcode_soilcarb_dpm,    & ! RothC soil C prognostic
                stashcode_soilcarb_rpm,    & ! RothC soil C prognostic
                stashcode_soilcarb_bio,    & ! RothC soil C prognostic
                stashcode_soilcarb_hum,    & ! RothC soil C prognostic
                stashcode_qcf2,            & ! second cloud ice prognostic
                stashcode_qrain,           & ! prognostic rain
                stashcode_qgraup,          & ! prognostic graupel
                stashcode_lcbase,          & ! lowest convective cloud base
                stashcode_3d_ccw,          & ! CCW profile sent to radiation
                stashcode_deep_conv_flag,  & ! deep convection flag
                stashcode_past_conv_precip,& ! past convective precipitation
                stashcode_past_conv_depth, & ! past convective depth
                stashcode_urbhgt  ,        & ! URBAN building height
                stashcode_urbhwr  ,        & ! URBAN height to width ratio
                stashcode_urbwrr  ,        & ! URBAN width ratio
                stashcode_urbdisp ,        & ! URBAN displacement height
                stashcode_urbztm  ,        & ! URBAN effective roughness length
                stashcode_urbalbwl,        & ! URBAN wall albedo
                stashcode_urbalbrd,        & ! URBAN road albedo
                stashcode_urbemisw,        & ! URBAN wall emissivity
                stashcode_urbemisr         & ! URBAN road emissivity
                )

              data_source( i ) % source = set_to_zero

            CASE( stashcode_rgrain,                                            &
                stashcode_snow_grnsiz_tiles & ! JULES snow layer grain size
                )

              data_source( i ) % source = set_to_const
              data_source( i ) % rconst = r0

            CASE( stashcode_tstar_land,                                        &
                stashcode_tstar_sea,                                           &
                stashcode_tstar_anom,                                          &
                stashcode_tstar_sice )

              data_source( i ) % source = Already_Processed

            CASE( stashcode_dctemp_tile,                                       &
                stashcode_dctemp_ssi,                                          &
                stashcode_tm_trans,                                            &
                stashcode_e_trb,                                               &
                stashcode_tsq_trb,                                             &
                stashcode_qsq_trb,                                             &
                stashcode_cov_trb,                                             &
                stashcode_zhpar_shcu )

              data_source( i ) % source = set_to_mdi
            CASE( stashcode_etadot,                                            &
                  stashcode_thetavd,                                           &
                  stashcode_dry_rho,                                           &
                  stashcode_exner_surf,                                        &
                  stashcode_psiw_surf,                                         &
                  stashcode_psiw_lid,                                          &
                  stashcode_mv,                                                &
                  stashcode_mcl,                                               &
                  stashcode_mcf,                                               &
                  stashcode_mr,                                                &
                  stashcode_mgr,                                               &
                  stashcode_mcf2 )
              ! Set Endgame prognostics to missing data.
              data_source( i ) % source = set_to_mdi

            END SELECT                       ! based on STASH item code
          END SELECT                         ! based on STASH section code

      !For files not found when reading from GRIB data
          IF ( grib .OR. grib2ff ) THEN

            SELECT CASE( fields_out( i ) % stashmaster % section )

            CASE( stashcode_prog_sec )

              SELECT CASE( fields_out( i ) % stashmaster % item )

              CASE( stashcode_sea_ice_temp,                                    &
                  stashcode_icethick,                                          &
                  stashcode_u_adv,                                             &
                  stashcode_v_adv,                                             &
                  stashcode_bulk_cf,                                           &
                  stashcode_liquid_cf,                                         &
                  stashcode_frozen_cf)

                data_source( i ) % source = Field_Calcs

              CASE( stashcode_w_adv,                                           &
                  stashcode_qcf,                                               &
                  stashcode_cc_lwp,                                            &
                  stashcode_unfrozen_soil,                                     &
                  stashcode_frozen_soil,                                       &
                  stashcode_qcl,                                               &
                  stashcode_n_turb_mixlvs,                                     &
                  stashcode_lvl_bse_dp_sc,                                     &
                  stashcode_lvl_top_dp_sc,                                     &
                  stashcode_bl_conv_flag,                                      &
                  stashcode_turb_temp,                                         &
                  stashcode_turb_humid,                                        &
                  stashcode_area_cf,                                           &
                  stashcode_sfc_zonal_cur,                                     &
                  stashcode_sfc_merid_cur,                                     &
                  stashcode_3d_cca,           & ! has to overide source=8
                  stashcode_can_water_tile,                                    &
                  stashcode_rho,              & ! rho calc'd after interp
                  stashcode_exner,                                             &
                  stashcode_can_conduct,                                       &
                  stashcode_ccb,                                               &
                  stashcode_cct,                                               &
                  stashcode_mean_canopyw,                                      &
                  stashcode_surf_z_curr,                                       &
                  stashcode_surf_m_curr,                                       &
                  stashcode_w)

                data_source( i ) % source = Set_To_Zero

              CASE( stashcode_bl_depth )

                data_source( i ) % source = set_to_const
                data_source( i ) % rconst = 500.000

              CASE( stashcode_z0 )

                data_source( i ) % source = set_to_const
                data_source( i ) % rconst = 0.500
              CASE( stashcode_theta )
                data_source( i ) % source = other_field
              END SELECT                 ! select by item code

            END SELECT                   ! select by section code

          END IF                       ! If GRIB section

        CASE Default                 ! Couldn't find a proper Internal Model

      ! The Model type didn't match one of the three defined types.
          WRITE (Cmessage, '(2A, I2, A, I3, A, I5)')                           &
              " Couldn't find a Sub-Model ID type for : ",                     &
              " Model ", fields_out( i ) % stashmaster % model,                &
              " Section ", fields_out( i ) % stashmaster % section,            &
              " Item ",  fields_out( i ) % stashmaster % item
          ErrorStatus = 25
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )

        END SELECT                   ! select by Internal model

    ! Check that Source is now set correctly otherwise, fail
        IF ( data_source( i ) % source == Input_Dump ) THEN
          WRITE ( Cmessage, '(A, I3, A, I5, A)')                               &
              'Section ',                                                      &
              fields_out( i ) % stashmaster % section,                         &
              ' Item ',                                                        &
              fields_out( i ) % stashmaster % item ,                           &
              ' : Required field is not in input dump!'
          ErrorStatus = 30
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )
        END IF

      END IF  ! If (pos == 0 ..)

  ! ---------------------------------------------------------------
  ! Some fields, even when found in the input dump, still need
  ! extra work. Mostly for fields that may have differing numbers
  ! of pseudo levels. But also fields needing to be recalculated to
  ! ensure consistancy with others.
  ! ---------------------------------------------------------------

      IF ( pos /= 0 .AND. data_source( i ) % source == Input_Dump .AND.        &
          h_int_method /= nearest_neighbour) THEN

  ! ----------------------------------
  ! This item is in the input dump
  ! ----------------------------------

        SELECT CASE( fields_out( i ) % stashmaster % model )

    ! -----------
    ! Atmos items
    ! -----------

        CASE ( Atmos_IM )

          SELECT CASE( fields_out( i ) % stashmaster % section )

          CASE( stashcode_prog_sec )

            SELECT CASE( fields_out( i ) % stashmaster % item )

            CASE( stashcode_icethick)

              IF (h_int_active) THEN
                data_source( i ) % source = Field_Calcs
              END IF

            CASE( stashcode_ice_temp_cat,                                      &
                stashcode_ice_conc_cat,                                        &
                stashcode_ice_thick_cat )

          !-------------------------------------------------------------
          ! Check number of pseudo levels in input dump is same as
          ! number of sea ice categories (nice) in output dump.
          ! If not -recalculate
          !-------------------------------------------------------------

              IF (fields_in (pos) % levels /= nice) THEN
                data_source( i ) % source = Field_Calcs
                IF ( mype == 0 .AND.                                           &
                    Printstatus >= PrStatus_Normal ) THEN
                  WRITE (6,'(A,I7,A)') "Setting source of stashcode item ",    &
                      fields_out( i ) % stashmaster % item,                    &
                      " to calculations due to difference in NICE"
                END IF
              END IF

            CASE( stashcode_ice_snow_depth )

          !-------------------------------------------------------------
          ! Check number of pseudo levels in input dump is same as
          ! number of sea ice categories (nice) in output dump.
          ! If not -set to zero
          !-------------------------------------------------------------

              IF (fields_in (pos) % levels /= nice) THEN
                data_source( i ) % source = set_to_zero
                IF ( mype == 0 .AND.                                           &
                    Printstatus >= PrStatus_Normal ) THEN
                  WRITE (6,'(A,I7,A)') "Setting source of stashcode item ",    &
                      fields_out( i ) % stashmaster % item,                    &
                      " to -set-to-zero- due to difference in NICE"
                END IF
              END IF

            CASE( stashcode_infil_max_tile,                                    &
                stashcode_sw_down_tile )

          !-------------------------------------------------------------
          ! Check number of tiles in input dump is same as
          ! number of tiles in output dump.
          ! If not -set to zero
          !-------------------------------------------------------------

              IF (fields_in (pos) % levels /=                                  &
                  fields_out( i )    % levels ) THEN
                data_source( i ) % source = set_to_zero
                IF ( mype == 0 .AND.                                           &
                    Printstatus >= PrStatus_Normal ) THEN
                  WRITE (6,'(A,I7,A)') "Setting source of stashcode item ",    &
                      fields_out( i ) % stashmaster % item,                    &
                      " to -set-to-zero- due to difference in no. of tiles"
                END IF
              END IF

            CASE( stashcode_catch_tile,                                        &
                stashcode_z0_tile,                                             &
                stashcode_z0h_tile )

          !-------------------------------------------------------------
          ! Check number of tiles in input dump is same as
          ! number of tiles in output dump.
          ! If not -set to user constant -1
          !-------------------------------------------------------------

              IF (fields_in (pos) % levels /=                                  &
                  fields_out( i )    % levels ) THEN
                data_source( i ) % source = set_to_const
                data_source( i ) % rconst = -1.000
                IF ( mype == 0 .AND.                                           &
                    Printstatus >= PrStatus_Normal ) THEN
                  WRITE (6,'(A,I7,A)') "Setting source of stashcode item ",    &
                      fields_out( i ) % stashmaster % item,                    &
                      " to -1 due to difference in no. of tiles"
                END IF
              END IF

            CASE( stashcode_can_water_tile,                                    &
                stashcode_tstar_tile,                                          &
                stashcode_snow_tile,                                           &
                stashcode_catch_snow,                                          &
                stashcode_snow_grnd,                                           &
                stashcode_rgrain,                                              &
                stashcode_snowdep_grd_tile,                                    &
                stashcode_snowpack_bk_dens,                                    &
                stashcode_nsnow_layrs_tiles,                                   &
                stashcode_snow_laythk_tiles,                                   &
                stashcode_snow_ice_tile,                                       &
                stashcode_snow_liq_tile,                                       &
                stashcode_snow_t_tile,                                         &
                stashcode_snow_laydns_tiles )

          !-------------------------------------------------------------
          ! Check number of pseudo-levels in input dump is same as
          ! number of tiles in output dump.
          ! If not -set to Field Calcs
          !-------------------------------------------------------------

              IF (fields_in (pos) % levels /=                                  &
                  fields_out( i ) % levels ) THEN
                data_source( i ) % source = Field_Calcs
                IF ( mype == 0 .AND.                                           &
                    Printstatus >= PrStatus_Normal ) THEN
                  WRITE (6,'(A,I7,A)') "Setting source of stashcode item ",    &
                      fields_out( i ) % stashmaster % item,                    &
                   " to -FieldCalcs- due to difference in no. of pseudo-levels"
                END IF
              END IF
            CASE (stashcode_dctemp_tile,                                       &
                stashcode_dctemp_ssi,                                          &
                stashcode_tm_trans)
              IF (h_int_active .OR.                                            &
                  fields_in (pos) % levels /=                                  &
                  fields_out( i ) % levels ) THEN
                data_source( i ) % source = set_to_mdi
                IF ( mype == 0 .AND.                                           &
                    Printstatus >= PrStatus_Normal ) THEN
                  WRITE (6,'(A,I7,A)') "Setting source of stashcode item ",    &
                      fields_out( i ) % stashmaster % item,                    &
                      " to -MDI- due to difference in no. of "//               &
                      " tiles or resolution"
                END IF
              END IF
            CASE( stashcode_etadot,                                            &
                  stashcode_thetavd,                                           &
                  stashcode_dry_rho,                                           &
                  stashcode_exner_surf,                                        &
                  stashcode_psiw_surf,                                         &
                  stashcode_psiw_lid,                                          &
                  stashcode_mv,                                                &
                  stashcode_mcl,                                               &
                  stashcode_mcf,                                               &
                  stashcode_mr,                                                &
                  stashcode_mgr,                                               &
                  stashcode_mcf2 )
              ! Set Endgame prognostics to missing data.
              data_source( i ) % source = set_to_mdi

            END SELECT                      ! Select on item number
          END SELECT                          ! Select on section number

        CASE Default                 ! Couldn't find a proper Internal Model

      ! The Model type didn't match one of the three defined types.
          WRITE (Cmessage, '(2A, I2, A, I3, A, I5)')                           &
              " Couldn't find a Sub-Model ID type for : ",                     &
              " Model ", fields_out( i ) % stashmaster % model,                &
              " Section ", fields_out( i ) % stashmaster % section,            &
              " Item ",  fields_out( i ) % stashmaster % item
          ErrorStatus = 70
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )

        END SELECT                            ! Select on Internal model

      END IF                                  ! Item was in Dump
    END DO

  ! ---------------------------------------------------------------
  ! Some fields still need extra work regardless of whether they're
  ! in the input dump or not.  Need to restart loop since we need to make sure
  ! all fields are set with some default values above.
  ! ---------------------------------------------------------------
    DO i = 1, field_count_out
      SELECT CASE( fields_out( i ) % stashmaster % model )

  ! -----------
  ! Atmos items
  ! -----------

      CASE ( Atmos_IM )

        SELECT CASE( fields_out( i ) % stashmaster % section )

        CASE( stashcode_prog_sec )

          SELECT CASE( fields_out( i ) % stashmaster % item )

        ! Field in output dump is one of the dust bins.
          CASE ( stashcode_dust1_mmr,                                          &
              stashcode_dust2_mmr,                                             &
              stashcode_dust3_mmr,                                             &
              stashcode_dust4_mmr,                                             &
              stashcode_dust5_mmr,                                             &
              stashcode_dust6_mmr)

          !-------------------------------------------------------------
          ! Check for presence of a dust bin in both input dump
          ! and output dump. pos and pos_out will both be non-zero.
          ! This means both dumps are using a dust scheme.
          !-------------------------------------------------------------
            CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust1_mmr,          &
                fields_in, field_count_in, pos, .TRUE.)
            CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust1_mmr,          &
                fields_out, field_count_out, pos_out, .TRUE.)
            IF ( pos /= 0 .AND. pos_out /= 0 ) THEN
          !-------------------------------------------------------------
          ! Check for presence of 3rd dust bin in both input and output
          ! dumps.
          ! If dust3 is in input dump (pos != 0) but not in output (pos_out=0)
          ! OR dust3 not in input dump (pos=0) but is in output (pos_out != 0)
          ! Then changing scheme means dust need fieldcalcs
          !-------------------------------------------------------------

              CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust3_mmr,        &
                  fields_in, field_count_in, pos, .TRUE.)
              CALL Rcf_Locate( stashcode_prog_sec, stashcode_dust3_mmr,        &
                  fields_out, field_count_out, pos_out, .TRUE.)

              IF ( ( pos /= 0 .AND. pos_out == 0 ) .OR.                        &
                  ( pos == 0 .AND. pos_out /= 0 ) ) THEN
                data_source( i ) % source = Field_Calcs
                IF ( mype == 0 .AND.                                           &
                    Printstatus >= PrStatus_Normal ) THEN
                  WRITE (6,'(A,I7,A)') "Setting source of stashcode item ",    &
                      fields_out( i ) % stashmaster % item,                    &
                      " to -FieldCalcs- due to difference in no. of dust bins"
                END IF
              END IF ! If dust bin 3 only in one of input or output dumps.

            END IF ! If dust bin 1 present in both input and output dumps.
          CASE (  stashcode_zw,                                                &
              stashcode_fexp,                                                  &
              stashcode_sthzw,                                                 &
              stashcode_fsat,                                                  &
              stashcode_fwetl,                                                 &
              stashcode_gamtot,                                                &
              stashcode_a_fsat,                                                &
              stashcode_c_fsat,                                                &
              stashcode_a_fwet,                                                &
              stashcode_c_fwet )
          ! NB. Some of these may have had data_source set above.
          ! For LSH we can use subroutine to check consistency of data sources.
            CALL rcf_lsh_field_checks ( data_source, fields_out,               &
                field_count_out, i )

          ! LSM is set independently of everything else since we
          ! need LSM for land-packed fields.  This was done before this routine 
          ! was called.
          CASE(stashcode_lsm)
            ! LSM are performed already.
            IF (data_source( i ) % source == input_dump) THEN
              data_source( i ) % source = already_processed
            END IF
          END SELECT                      ! Select on item number
        END SELECT                          ! Select on section number

      CASE default                 ! Couldn't find a proper Internal Model

    ! The Model type didn't match one of the three defined types.
        WRITE (Cmessage, '(2A, I2, A, I3, A, I5)')                             &
            " Couldn't find a Sub-Model ID type for : ",                       &
            " Model ", fields_out( i ) % stashmaster % model,                  &
            " Section ", fields_out( i ) % stashmaster % section,              &
            " Item ",  fields_out( i ) % stashmaster % item
        ErrorStatus = 80
        CALL Ereport( RoutineName, ErrorStatus, Cmessage )

      END SELECT                            ! Select on Internal model

  ! Check on Rimwidtha if an atmos LBC file for copying
      IF ( data_source( i ) % source == Input_Dump .AND.                       &
          (fields_out( i ) % stashmaster % grid_type ==ppx_atm_lbc_theta.OR.   &
          fields_out( i ) % stashmaster % grid_type == ppx_atm_lbc_u   .OR.    &
          fields_out( i ) % stashmaster % grid_type == ppx_atm_lbc_v) ) THEN

        IF ( hdr_in  % Lookup( lbrow, fields_out( i ) % dump_pos ) /=          &
            hdr_out % Lookup( lbrow, fields_out( i ) % dump_pos) ) THEN
          Cmessage = 'Rimwidth needs to be equal for input ' //                &
              'and output grids if lbcs are copied'
          ErrorStatus = 90
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )
        END IF
      END IF

    END DO

!------------------------------------------------------------------
! Print out some diagnostics if required....
!------------------------------------------------------------------
    IF (mype == 0 .AND. PrintStatus >= PrStatus_Oper ) THEN
! form   ="(2(i5,' '),4(i4,' '),e10.4,' ',50a)"
! form_c ="(2(a5,' '),4(a4,' '),a10,  ' ',50a)"

      WRITE (6,'(A)') ''
      WRITE (6,'(A)') 'Data source :'
      WRITE (6,form_c) 'Sect', 'Item', 'src', 'Dom', 'AncS', 'AncI',           &
          'Real Const', 'Ancil File'

      DO i = 1, field_count_out
        WRITE (6,form) fields_out( i ) % stashmaster % section,                &
            fields_out( i ) % stashmaster % item,                              &
            data_source( i ) % source,                                         &
            data_source( i ) % Domain,                                         &
            data_source( i ) % Ancil_SctnC,                                    &
            data_source( i ) % Ancil_ItemC,                                    &
            data_source( i ) % RConst,                                         &
            TRIM(data_source(i) % Ancil_File)
      END DO
      WRITE (6,'(A)') ''

    END IF


    RETURN
  END SUBROUTINE Rcf_Set_Data_Source
END MODULE Rcf_Set_Data_Source_Mod
