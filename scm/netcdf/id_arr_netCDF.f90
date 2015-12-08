!
! *********************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *********************************COPYRIGHT*********************************
!
! Declaration statments for possible ids/allocatble arrays in netCDF files.

MODULE id_arr_netCDF

  USE global_SCMop, ONLY: incdf

  IMPLICIT NONE

!-----------------------------------------------------------------------------
!
! Description:
!   Declares ids/arrays used when reading in netcdf files into SCM
!
! Method:
!   Not all ids/arrays may be used and no units are specified here as these
!   are dependent on the content of the netcdf file to be read.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!-----------------------------------------------------------------------------


  INTEGER(incdf) ::  &
    status = 0        ! Error status 

  ! NetCDF file dimension ids for ...
  INTEGER(incdf) ::  &
    ncid             &!... file 
  , dimid_time       &!... time
  , dimid_nlev       &!... number of level
  , dimid_nlevp1     &!... number of half levels
  , dimid_slev        !... number of soil levels

  ! Netcdf dimensions
  INTEGER(incdf) ::  &
    xtype            &
  , attnum           &
  , ns_nc            &! Number of soil levels in file
  , nt_nc            &! Number of time levels in file
  , nk_nc            &! Number of levels in file
  , nk_half_nc       &! Number of half levels in file
  , coor_lena        &! 
  , coor_lenb         !

  
  ! Ids for possible variables in netCDF files
  INTEGER(incdf) ::  &
    id_time    &! Time
  , id_date    &! Date
  , id_mon     &! Months
  , id_day     &! Days
  , id_hour    &! Hours
  , id_min     &! minutes
  , id_sec     &! Seconds
  , id_lat     &! latitude
  , id_lon     &! longitude
  , id_p       &! Pressure
  , id_t       &! Temperture
  , id_q       &! H2O mixing ratio (vapour)
  , id_ql      &! H2O mixing ratio (liquid)
  , id_qi      &! H2O mixing ratio (ice)
  , id_u       &! Zonal wind
  , id_v       &! Meridional wind
  , id_omega   &! Vertical pressure velocity
  , id_hdivt   &! Hor.  adv. of s (dry static temperature)
  , id_hdivq   &! Hor.  adv. of q
  , id_vdivt   &! Vert. adv. of s (dry static temperature)
  , id_vdivq   &! Vert. adv. of q
  , id_uadv    &! Advective u tendency
  , id_vadv    &! Advective v tendency
  , id_tadv    &! Advective T tendency
  , id_qadv    &! Advective q tendency
  , id_alb     &! Albedo
  , id_alb_snw &! Snow albedo
  , id_flx_H   &! Srf. sensible heat flux
  , id_flx_E   &! Srf. latent heat flux
  , id_osst    &! Open sst
  , id_orog    &! Orography
  , id_st      &! Soil temperature
  , id_sm      &! Soil moisture
  , id_hvc     &! High veg cover
  , id_lvc     &! Low veg cover
  , id_s_thk   &! Soil layer thickness
  , id_z0m     &! Srf roughness length (Momentum)
  , id_z0h     &! Srf roughness length (Heat)
  , id_tsskin  &! Skin temperature
  , id_qskin   &! Skin reservoir content
  , id_ts_air  &! Srf. air temperature
  , id_ps      &! Srf. pressure (Pa)
  , id_omegas  &! Srf. rate of change of pressure (dp/dt)
  , id_qs_srf  &! Srf. saturation mixing ratio
  , id_cf      &! Cloud fraction 
  , id_ug      &! Geostrophic u
  , id_vg      &! Geostrophic v
  , id_tsnow   &! Snow temperature
  , id_rho_snw &! Snow density
  , id_lsm     &! Land sea mask
  , id_sice_cf &! Sea ice fraction
  , id_sice_t  &! Sea ice temperature
  , id_sndep    ! Snow depth

  !---------------------------------------------------------------------------
  ! Arrays to receive forcing data/info from netCDF files
  !---------------------------------------------------------------------------

  INTEGER, Allocatable :: &
    file_date(:) &
  , file_secs(:) &
  , file_time(:)  ! Time levels(s)

!  REAL ::           &
!    file_lat     (:)   &! Latitude
!  , file_long    (:)   &! Longitude
!  , file_albsoil (:)   ! Albedo
!  , file_tstari 

  REAL, ALLOCATABLE ::  &
    file_t        (:,:) &! Temperature
  , file_q        (:,:) &! H2O mixing ratio (vapour)
  , file_ql       (:,:) &! H2O mixing ratio (liquid)
  , file_qi       (:,:) &! H2O mixing ratio (ice)
  , file_cf       (:,:) &! Cloud fration
  , file_u        (:,:) &! Zonal wind
  , file_v        (:,:) &! Meridional wind
  , file_uadv     (:,:) &! LS adv. of u
  , file_vadv     (:,:) &! LS adv. of v
  , file_tadv     (:,:) &! LS adv. of T
  , file_qadv     (:,:) &! LS adv. of q
  , file_ug       (:,:) &! Geostrophic zonal wind
  , file_vg       (:,:) &! Geostrophic meridional wind
  , file_st       (:,:) &! Soil layer temperature
  , file_sm       (:,:) &! Soil layer moisture
  , file_omega    (:,:) &! Rate of change of pressure (dp/dt)
  , file_divt     (:,:) &! Hor.  adv. of s (dry static temperature)
  , file_divq     (:,:) &! Hor.  adv. of q
  , file_vertdivt (:,:) &! Vert. adv. of s (dry static temperature)
  , file_vertdivq (:,:) &! Vert. adv. of q
  , file_w        (:,:) &! Vertical wind         (calculated)
  , file_th       (:,:)  ! Potential temperature (calculated)


  REAL, ALLOCATABLE, DIMENSION(:,:) :: &
    file_p_full                        &! Pressure on full levels
  , file_p_half                        &! Pressure on half levels
  , file_exner_full                     ! Exner    on full levels

  REAL, ALLOCATABLE, DIMENSION(:)   :: &
    file_coor_par_a                    &
  , file_coor_par_b                    &
  , file_dzsoil                         ! Soil layer thickness

  REAL, ALLOCATABLE, DIMENSION(:)   :: &
    file_year          &! Year
  , file_month         &! Month
  , file_day           &! Day
  , file_hour          &! Hour
  , file_min           &! Minute
  , file_z_prf         &! Height of obs levels (assumed fixed at init runtime)
  , file_z_half_prf    &! Height of obs levels (assumed fixed at init runtime)
  , file_z_full_prf    &! Height of obs levels (assumed fixed at init runtime)
  , file_layer_thk     &! Layer thickness of obs level (Using mean T of layer)
  , file_flx_H         &! Srf. sensible heat flux
  , file_flx_E         &! Srf. latent   heat flux
  , file_tsskin        &! Srf. skin temperature
  , file_ts_air        &! Srf. air temperature
  , file_p_srf         &! Srf. mean pressure
  , file_omega_srf     &! Srf. rate of change of pressure (dp/dt)
  , file_qs_srf        &! Srf. saturation mixing ratio
  , file_u_srf         &! Srf. zonal wind
  , file_v_srf         &! Srf. meridional wind
  , file_rhs_air       &! Srf. air relative humidity
  , file_w_srf         &! Srf. vertical velocity (calculated)
  , file_th_srf        &! Srf. potential temperature (calculated)
  , file_q_srf         &! Srf. mixing ratio (calculated)
  , file_p             &! Pressure levels
  , file_z0h           &! Srf. roughness length (heat)
  , file_z0m           &! Srf. roughness length (momentum)
  , file_t_skin        &! Srf. skin temperature
  , file_q_skin        &! Srf. skin moisture
  , file_t_snow        &! Snow temperature
  , file_alb_snow      &! Snow albedo
  , file_rho_snow      &! Snow density
  , file_lsm           &! Land mask/fraction
  , file_sice_frc      &! Sea-ice fraction
  , file_sice_t        &! Sea-ice temperature
  , file_osst          &! Open-sea surface temperature
  , file_orog          &! Orographic height
  , file_z_snow        &! Snow depth
  , file_alb           &! Soil albedo  
  , file_lat           &! latitude
  , file_long          &! longitude
  , file_hvc           &! High vegetation cover
  , file_lvc            ! Low vegetation cover

END MODULE id_arr_netCDF
