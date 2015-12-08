! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : INJULES

MODULE s_injules

  USE scm_cntl_mod, ONLY: scm_nml

  USE s_maxdim, ONLY:                                                &
    mx_rw_lng, mx_rw, mx_nlnd, mx_sm_lv, mx_ntile, dim_cs1           &
  , mx_npft, mx_ntype

  USE scm_utils, ONLY:                                               &
    rmdi, rw_lng, rw, nlnd, sm_lv, ntype, ntile, npft                &
  , zhook_in, zhook_out, jprb, lhook, dr_hook

  USE s_main_force, ONLY:                    &
    smi_opt                                  &
  , scm_gs              => gs                &
  , scm_fsmc            => fsmc              &
  , scm_frac_typ        => frac_typ          &
  , scm_smcli           => smcli             &
  , scm_sth             => sth               &
  , scm_canht           => canht             &
  , scm_lai             => lai               &
  , scm_z0_tile         => z0_tile           &
  , scm_z0h_tile        => z0h_tile          &
  , scm_rgrain          => rgrain            &
  , scm_infil_tile      => infil_tile        &
  , scm_snow_tile       => snow_tile         &
  , scm_tstar_tile      => tstar_tile        &
  , scm_catch           => catch             &
  , scm_canopy          => canopy            &
  , scm_frac_disturb    => frac_disturb      &
  , scm_npp_ft_acc      => npp_ft_acc        &
  , scm_resp_w_ft_acc   => resp_w_ft_acc     &
  , scm_resp_s_acc      => resp_s_acc        &
  , scm_cs              => cs                &
  , scm_g_leaf_phen_acc => g_leaf_phen_acc   &
  , scm_g_leaf_acc      => g_leaf_acc        &
  , scm_lw_down         => lw_down

  USE s_main_force, ONLY:                    &
    scm_lake_depth      => lake_depth        &
  , scm_lake_fetch      => lake_fetch        &
  , scm_lake_t_mean     => lake_t_mean       &
  , scm_lake_t_mxl      => lake_t_mxl        &
  , scm_lake_t_ice      => lake_t_ice        &
  , scm_lake_h_mxl      => lake_h_mxl        &
  , scm_lake_h_ice      => lake_h_ice        &
  , scm_lake_shape      => lake_shape        &
  , scm_lake_g_dt       => lake_g_dt         &
  , scm_clapp_levs      => clapp_levs        &
  , scm_sathh_levs      => sathh_levs        &
  , scm_hcap_levs       => hcap_levs         &
  , scm_hcon_levs       => hcon_levs         &
  , scm_satcon_levs     => satcon_levs       &
  , scm_smvccl_levs     => smvccl_levs       &
  , scm_smvcwt_levs     => smvcwt_levs       &
  , scm_smvcst_levs     => smvcst_levs

  IMPLICIT NONE

  PRIVATE


!=============================================================================
!
! Description:
!   Allows SCM to read in INJULES namelist from forcing file, scm_nml.
!   INJULES contains data for initialising soil moisture for the JULES code
!
! Method:
!   Namelist INJULES is defined in this module and read in by contained
!   subroutine read_injules.  Scalar variables are read directly to
!   s_main_force.  Array variables are declared as private single rank arrays
!   of fixed length (arrays cannot be varying length in namelists) and
!   reshaped before being transferred to arrays of the correct size/shape in
!   s_main_force.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!=============================================================================
!
! Declarations:

  !--------------------------------------------------
  ! Fixed arrays to read in namelist before reshaping
  !--------------------------------------------------
  REAL ::                              &
    gs       (mx_nlnd)          = rmdi &! Stomatal conductance
  , fsmc     (mx_nlnd)          = rmdi &! Soil Moisture Stress
  , frac_typ (mx_nlnd*mx_ntype) = rmdi &! Fractions of surface types
  , smcli    (mx_nlnd*mx_sm_lv) = rmdi &! Initial SMCL profile (kg/m2)
  , sth      (mx_nlnd*mx_sm_lv) = rmdi  ! Total soil moisture in layers
                                        ! as fraction of saturation

  REAL ::                              &
    canht    (mx_nlnd*mx_npft)  = rmdi &! Canopy height (m)
  , lai      (mx_nlnd*mx_npft)  = rmdi  ! Leaf area index


  REAL ::                                 &
    z0_tile    (mx_nlnd*mx_ntile) = rmdi  &! Tile roughness lengths (m).
  , z0h_tile   (mx_nlnd*mx_ntile) = rmdi  &! Tile thermal roughness lengths (m).
                                           ! These are values without
                                           ! snow cover (i.e. bare roughnesses)
  , rgrain     (mx_nlnd*mx_ntile) = rmdi  &! Snow grain size (microns)
  , infil_tile (mx_nlnd*mx_ntile) = rmdi  &! Max. surface infiltration
  , snow_tile  (mx_nlnd*mx_ntile) = rmdi  &! Lying snow on tiles (MOSES 2)
  , tstar_tile (mx_nlnd*mx_ntile) = rmdi  &! Surface tile temperature
  , catch      (mx_nlnd*mx_ntile) = rmdi  &! Surf/canopy water capacity
                                           ! (snow-free land tiles) (kg/m2)
  , canopy     (mx_nlnd*mx_ntile) = rmdi   ! Surf/canopy water
                                           ! (snow-free land tiles) (kg/m2)


  ! 2A Triffid variables
  !=====================
  REAL ::                                    &
    frac_disturb    (mx_nlnd)         = rmdi &! Fraction of gridbox in which
                                              ! vegetation is disturbed.
  , npp_ft_acc      (mx_nlnd)         = rmdi &! Acc. npp_ft
  , resp_w_ft_acc   (mx_nlnd*mx_npft) = rmdi &! Acc. resp_w_ft
  , resp_s_acc      (mx_nlnd*dim_cs1) = rmdi &! Acc. resp_s
  , cs              (mx_nlnd*dim_cs1) = rmdi &! Soil carbon (kg C/m2)
  , g_leaf_phen_acc (mx_nlnd*mx_npft) = rmdi &! Acc. leaf turnover rate
                                              ! with phenology
  , g_leaf_acc      (mx_nlnd*mx_npft) = rmdi  ! Acc. g_leaf


  REAL ::                                    &
     lw_down (mx_rw_lng*mx_rw) = rmdi     ! Surface downward LW

  ! FLake lake scheme prognostics
  !===============================
  REAL ::                        &
    lake_depth  (mx_nlnd) = rmdi &! lake depth (m)
  , lake_fetch  (mx_nlnd) = rmdi &! lake fetch (m)
  , lake_t_mean (mx_nlnd) = rmdi &! lake mean temperature (K)
  , lake_t_mxl  (mx_nlnd) = rmdi &! lake mixed-layer temperature (K)
  , lake_t_ice  (mx_nlnd) = rmdi &! lake ice upper boundary temperature (K)
  , lake_h_mxl  (mx_nlnd) = rmdi &! lake mixed-layer depth (m)
  , lake_h_ice  (mx_nlnd) = rmdi &! lake ice thickness (m)
  , lake_shape  (mx_nlnd) = rmdi &! thermocline shape factor
  , lake_g_dt   (mx_nlnd) = rmdi  ! lake ht.flx / dT (W m-2 K-1)


  ! Jules Soil Profiles
  !===============================
  REAL ::                                     &
    clapp_levs  (mx_nlnd*mx_sm_lv)     = rmdi &
  , sathh_levs  (mx_nlnd*mx_sm_lv)     = rmdi &
  , hcap_levs   (mx_nlnd*mx_sm_lv)     = rmdi &
  , hcon_levs   (mx_nlnd*(mx_sm_lv+1)) = rmdi &
  , satcon_levs (mx_nlnd*(mx_sm_lv+1)) = rmdi &
  , smvccl_levs (mx_nlnd*mx_sm_lv)     = rmdi &
  , smvcwt_levs (mx_nlnd*mx_sm_lv)     = rmdi &
  , smvcst_levs (mx_nlnd*mx_sm_lv)     = rmdi


  !---------------------------------------------------------------------------
  ! Define namelist
  !---------------------------------------------------------------------------
  NAMELIST/injules/                                                           &
    smi_opt, smcli, fsmc, sth, canht, catch, snow_tile, lai, z0_tile, z0h_tile&
  , tstar_tile, canopy, frac_typ, frac_disturb, infil_tile, rgrain, cs, gs    &
  , g_leaf_acc, g_leaf_phen_acc, npp_ft_acc, resp_w_ft_acc, resp_s_acc        &
  , lw_down, lake_depth, lake_fetch, lake_t_mean, lake_t_mxl, lake_t_ice      &
  , lake_h_mxl, lake_h_ice, lake_shape, lake_g_dt                             &
  , clapp_levs, sathh_levs, hcap_levs, hcon_levs, satcon_levs                 &
  , smvccl_levs, smvcwt_levs, smvcst_levs


  PUBLIC read_injules

!=============================================================================
CONTAINS

  SUBROUTINE read_injules

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    INTEGER :: istatus
    CHARACTER(LEN=11), PARAMETER :: routinename='read_injules'

    IF (lhook) CALL dr_hook('READ_INJULES',zhook_in,zhook_handle)

    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      WRITE(6,'(A)') ' Error opening file on unit 10 from '//routinename
      WRITE(6,'(A)') ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,'(A)') ' Icode    = 500'
    END IF

    READ  (10, injules)
    CLOSE (10)

    IF (nlnd > 0) THEN

      IF (gs(1)   /= rmdi) scm_gs   = RESHAPE( gs,   (/nlnd/) )
      IF (fsmc(1) /= rmdi) scm_fsmc = RESHAPE( fsmc, (/nlnd/) )

      IF (smcli(1) /= rmdi) scm_smcli = RESHAPE( smcli, (/nlnd,sm_lv/) )
      IF (sth(1)   /= rmdi) scm_sth   = RESHAPE( sth,   (/nlnd,sm_lv/) )
      IF (canht(1) /= rmdi) scm_canht = RESHAPE( canht, (/nlnd,npft/) )
      IF (lai(1)   /= rmdi) scm_lai   = RESHAPE( lai,   (/nlnd,npft/) )

      IF (rgrain(1) /= rmdi) scm_rgrain = RESHAPE( rgrain, (/nlnd,ntile/) )
      IF (catch(1)  /= rmdi) scm_catch  = RESHAPE( catch,  (/nlnd,ntile/) )
      IF (canopy(1) /= rmdi) scm_canopy = RESHAPE( canopy, (/nlnd,ntile/) )

      IF (frac_typ(1) /= rmdi)                                                &
        scm_frac_typ   = RESHAPE( frac_typ, (/nlnd,ntype/) )

      IF (z0_tile(1)    /= rmdi)                                              &
        scm_z0_tile      = RESHAPE( z0_tile,    (/nlnd,ntile/) )

      IF (z0h_tile(1)   /= rmdi)                                              &
        scm_z0h_tile     = RESHAPE( z0h_tile,   (/nlnd,ntile/) )

      IF (infil_tile(1) /= rmdi)                                              &
        scm_infil_tile   = RESHAPE( infil_tile, (/nlnd,ntile/) )

      IF (snow_tile(1)  /= rmdi)                                              &
        scm_snow_tile    = RESHAPE( snow_tile,  (/nlnd,ntile/) )

      IF (tstar_tile(1) /= rmdi)                                              &
        scm_tstar_tile   = RESHAPE( tstar_tile, (/nlnd,ntile/) )

      ! Triffid
      IF (frac_disturb(1)    /= rmdi)                                         &
        scm_frac_disturb      = RESHAPE( frac_disturb,    (/nlnd/) )
      IF (npp_ft_acc(1)      /= rmdi)                                         &
        scm_npp_ft_acc        = RESHAPE( npp_ft_acc,      (/nlnd/) )
      IF (resp_w_ft_acc(1)   /= rmdi)                                         &
        scm_resp_w_ft_acc     = RESHAPE( resp_w_ft_acc,   (/nlnd,npft/) )
      IF (resp_s_acc(1)      /= rmdi)                                         &
        scm_resp_s_acc        = RESHAPE( resp_s_acc,      (/nlnd,dim_cs1/) )
      IF (cs(1)              /= rmdi)                                         &
        scm_cs                = RESHAPE( cs,              (/nlnd,dim_cs1/) )
      IF (g_leaf_phen_acc(1) /= rmdi)                                         &
        scm_g_leaf_phen_acc   = RESHAPE( g_leaf_phen_acc, (/nlnd,npft/) )
      IF (g_leaf_acc(1)      /= rmdi)                                         &
        scm_g_leaf_acc        = RESHAPE( g_leaf_acc,      (/nlnd,npft/) )

      IF (lake_depth(1)  /= rmdi)                                             &
        scm_lake_depth    = RESHAPE( lake_depth,  (/nlnd/) )
      IF (lake_fetch(1)  /= rmdi)                                             &
        scm_lake_fetch    = RESHAPE( lake_fetch,  (/nlnd/) )
      IF (lake_t_mean(1) /= rmdi)                                             &
        scm_lake_t_mean   = RESHAPE( lake_t_mean, (/nlnd/) )
      IF (lake_t_mxl(1)  /= rmdi)                                             &
        scm_lake_t_mxl    = RESHAPE( lake_t_mxl,  (/nlnd/) )
      IF (lake_t_ice(1)  /= rmdi)                                             &
        scm_lake_t_ice    = RESHAPE( lake_t_ice,  (/nlnd/) )
      IF (lake_h_mxl(1)  /= rmdi)                                             &
        scm_lake_h_mxl    = RESHAPE( lake_h_mxl,  (/nlnd/) )
      IF (lake_h_ice(1)  /= rmdi)                                             &
        scm_lake_h_ice    = RESHAPE( lake_h_ice,  (/nlnd/) )
      IF (lake_shape(1)  /= rmdi)                                             &
        scm_lake_shape    = RESHAPE( lake_shape,  (/nlnd/) )
      IF (lake_g_dt(1)   /= rmdi)                                             &
        scm_lake_g_dt     = RESHAPE( lake_g_dt,   (/nlnd/) )


      IF (clapp_levs(1)  /= rmdi)                                             &
        scm_clapp_levs    = RESHAPE( clapp_levs,  (/nlnd,sm_lv/) )
      IF (sathh_levs(1)  /= rmdi)                                             &
        scm_sathh_levs    = RESHAPE( sathh_levs,  (/nlnd,sm_lv/) )
      IF (hcap_levs(1)   /= rmdi)                                             &
        scm_hcap_levs     = RESHAPE( hcap_levs,   (/nlnd,sm_lv/) )
      IF (hcon_levs(1)   /= rmdi)                                             &
        scm_hcon_levs     = RESHAPE( hcon_levs,   (/nlnd,sm_lv+1/) )
      IF (satcon_levs(1) /= rmdi)                                             &
        scm_satcon_levs   = RESHAPE( satcon_levs, (/nlnd,sm_lv+1/) )
      IF (smvccl_levs(1) /= rmdi)                                             &
        scm_smvccl_levs   = RESHAPE( smvccl_levs, (/nlnd,sm_lv/) )
      IF (smvcwt_levs(1) /= rmdi)                                             &
        scm_smvcwt_levs   = RESHAPE( smvcwt_levs, (/nlnd,sm_lv/) )
      IF (smvcst_levs(1) /= rmdi)                                             &
        scm_smvcst_levs   = RESHAPE( smvcst_levs, (/nlnd,sm_lv/) )

    END IF

    IF (lw_down(1) /= rmdi) scm_lw_down = RESHAPE( lw_down, (/rw_lng,rw/) )

    IF (lhook) CALL dr_hook('READ_INJULES',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_injules

!=============================================================================
END MODULE s_injules

