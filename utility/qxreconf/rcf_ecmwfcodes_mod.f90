! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ List of ECMWF magic numbers

Module Rcf_ECMWFcodes_Mod

! Description:
!    Magic numbers for ECMWF [dia|pro]gnostic codes, as defined by
!    the ECMWF version of 'Table 2' used in the GRIB description.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: Fortran 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None
Public

! ECMWF codes used in the reconfiguration

! Sea-ice Fraction
Integer, Parameter      :: ECMWFcode_icefrac             = 31

! Snow Density
Integer, Parameter      :: ECMWFcode_snow_density        = 33

! Volumetric Soil Water Layer 1
Integer, Parameter      :: ECMWFcode_soil_moist_1        = 39

! Volumetric Soil Water Layer 2
Integer, Parameter      :: ECMWFcode_soil_moist_2        = 40

! Volumetric Soil Water Layer 3
Integer, Parameter      :: ECMWFcode_soil_moist_3        = 41

! Volumetric Soil Water Layer 4
Integer, Parameter      :: ECMWFcode_soil_moist_4        = 42

! Atmospheric Tide
Integer, Parameter      :: ECMWFcode_atm_tide            = 127

! Budget Values
Integer, Parameter      :: ECMWFcode_budget              = 128

! Geopotential
Integer, Parameter      :: ECMWFcode_geopot              = 129

! Temperature
Integer, Parameter      :: ECMWFcode_T                   = 130

! U-component of Wind
Integer, Parameter      :: ECMWFcode_u                   = 131

! V-component of Wind
Integer, Parameter      :: ECMWFcode_v                   = 132

! Specific Humidity
Integer, Parameter      :: ECMWFcode_spec_hum            = 133

! Surface Pressure
Integer, Parameter      :: ECMWFcode_pstar               = 134

! Vertical Velocity
Integer, Parameter      :: ECMWFcode_w                   = 135

! Precipitable Water Content
Integer, Parameter      :: ECMWFcode_precip_water        = 137

! Vorticity
Integer, Parameter      :: ECMWFcode_vort                = 138

! Soil Temperature level 1
Integer, Parameter      :: ECMWFcode_soil_temp_1         = 139

! Soil Wetness level 1
Integer, Parameter      :: ECMWFcode_soil_wet_1          = 140

! Snow Depth
Integer, Parameter      :: ECMWFcode_snow_depth          = 141

! Large Scale Precipitation
Integer, Parameter      :: ECMWFcode_lge_precip          = 142

! Convective Precipitation
Integer, Parameter      :: ECMWFcode_conv_precip         = 143

! Snow Fall
Integer, Parameter      :: ECMWFcode_snowfall            = 144

! Boundary Layer Dissipation
Integer, Parameter      :: ECMWFcode_b_layer_diss        = 145

! Surface Flux of Sensible Heat
Integer, Parameter      :: ECMWFcode_s_flux_s_heat       = 146

! Surface Flux of Latent Heat
Integer, Parameter      :: ECMWFcode_s_flux_l_heat       = 147

! Mean Sea Level Pressure
Integer, Parameter      :: ECMWFcode_msl_press           = 148

! Log Surface Pressure
Integer, Parameter      :: ECMWFcode_log_p               = 152

! Divergence
Integer, Parameter      :: ECMWFcode_div                 = 153

! Height (Geopotential)
Integer, Parameter      :: ECMWFcode_h_geopot            = 156

! Relative Humidity
Integer, Parameter      :: ECMWFcode_rel_hum             = 157

! Tendency of Surface Pressure
Integer, Parameter      :: ECMWFcode_tend_p              = 158

! Total Cloud Cover
Integer, Parameter      :: ECMWFcode_t_cloud_cover       = 159

! U-wind at 10m
Integer, Parameter      :: ECMWFcode_u_10                = 160

! V-wind at 10m
Integer, Parameter      :: ECMWFcode_v_10                = 161

! Temperature at 2m
Integer, Parameter      :: ECMWFcode_temp_2              = 167

! Dewpoint at 2m
Integer, Parameter      :: ECMWFcode_dew_p_2             = 168

! Soil Temperature level 2
Integer, Parameter      :: ECMWFcode_soil_temp_2         = 170

! Soil Wetness level 2
Integer, Parameter      :: ECMWFcode_soil_wet_2          = 171

! Land-Sea Mask
Integer, Parameter      :: ECMWFcode_lsm                 = 172

! Surface Roughness
Integer, Parameter      :: ECMWFcode_surface_rough       = 173

! Albedo
Integer, Parameter      :: ECMWFcode_albedo              = 174

! Downwards Shortwave Radiation (surface)
Integer, Parameter      :: ECMWFcode_down_swave_surface  = 175

! Net Shortwave Radiation (surface)
Integer, Parameter      :: ECMWFcode_net_swave_surface   = 176

! Net Longwave Radiation (surface)
Integer, Parameter      :: ECMWFcode_net_lwave_surface   = 177

! Net Shortwave Radiation (top of atmosphere)
Integer, Parameter      :: ECMWFcode_net_swave_top       = 178

! Net Longwave Radiation (top of atmosphere)
Integer, Parameter      :: ECMWFcode_net_lwave_top       = 179

! U-component of Surface Wind Stress
Integer, Parameter      :: ECMWFcode_u_stress            = 180

! V-component of Surface Wind Stress
Integer, Parameter      :: ECMWFcode_v_stress            = 181

!Evaporation
Integer, Parameter      :: ECMWFcode_evap                = 182

! Soil Temperature level 3
Integer, Parameter      :: ECMWFcode_soil_temp_3         = 183

! Soil Wetness level 3
Integer, Parameter      :: ECMWFcode_soil_wet_3          = 184

! Convective Cloud Cover
Integer, Parameter      :: ECMWFcode_conv_cloud_cover    = 185

! Low Cloud Cover
Integer, Parameter      :: ECMWFcode_low_cloud_cover     = 186

! Medium Cloud Cover
Integer, Parameter      :: ECMWFcode_med_cloud_cover     = 187

! High Cloud Cover
Integer, Parameter      :: ECMWFcode_high_cloud_cover    = 188

! Sunshine Duration
Integer, Parameter      :: ECMWFcode_sun_duration        = 189

! Ozone mass mixing ratio
Integer, Parameter      :: ECMWFcode_ozone               = 203

! Skin Temperature
Integer, Parameter      :: ECMWFcode_skin_temp           = 235

! Soil Temperature level 4
Integer, Parameter      :: ECMWFcode_soil_temp_4         = 236

! Soil Wetness level 4
Integer, Parameter      :: ECMWFcode_soil_wet_4          = 237

! Cloud Liquid Water Content
Integer, Parameter      :: ECMWFcode_qcl                 = 246

! Cloud Ice Water Content
Integer, Parameter      :: ECMWFcode_qcf                 = 247

! Cloud Cover
Integer, Parameter      :: ECMWFcode_cc                  = 248

! GEMS Hydrophobic Organic Matter Aerosol (GEMS Table 210)
Integer, Parameter      :: ECMWFcode_OMFRSH              = 7

! GEMS Hydrophilic Organic Matter Aerosol (GEMS Table 210)
Integer, Parameter      :: ECMWFcode_OMAGD               = 8

! GEMS Hydrophobic Black Carbon Aerosol (GEMS Table 210)
Integer, Parameter      :: ECMWFcode_BCFRSH              = 9

! GEMS Hydrophilic Black Carbon Aerosol (GEMS Table 210)
Integer, Parameter      :: ECMWFcode_BCAGD               = 10

! Methane (GEMS Table 210)
Integer, Parameter      :: ECMWFcode_CH4                 = 62

! Nitrogen Oxides (GEMS Table 210)
Integer, Parameter      :: ECMWFcode_NOX                 = 129

! GEMS NO2 
Integer, Parameter      :: ECMWFcode_NO2                 = 121

! Carbon monoxide (GEMS Table 210)
Integer, Parameter      :: ECMWFcode_CO                  = 123

! Formaldehyde (GEMS Table 210)
Integer, Parameter      :: ECMWFcode_HCHO                = 124

! GEMS Ozone (GEMS Table 210)
Integer, Parameter      :: ECMWFcode_GO3                 = 203

End Module Rcf_ECMWFcodes_Mod
