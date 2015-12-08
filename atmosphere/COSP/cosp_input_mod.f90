! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE cosp_input_mod

! Description:
!       Input/namelist control of COSP simulator.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: COSP

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

USE missing_data_mod, ONLY: rmdi, imdi 

IMPLICIT NONE

!- Control variables set in the UMUI via the RUN_COSP namelist
LOGICAL :: l_cosp            = .FALSE.
LOGICAL :: cosp_cloudsat_sim = .FALSE.
LOGICAL :: cosp_lidar_sim    = .FALSE.
LOGICAL :: cosp_isccp_sim    = .FALSE.
LOGICAL :: cosp_misr_sim     = .FALSE.
LOGICAL :: cosp_modis_sim    = .FALSE.
LOGICAL :: cosp_rttov_sim    = .FALSE.
INTEGER :: cosp_npoints_it   = IMDI
INTEGER :: cosp_ncolumns     = IMDI
LOGICAL :: cosp_use_vgrid    = .FALSE.

!- Defaults for vertical grid if cosp_use_vgrid is set to true
INTEGER :: cosp_nlr        = 40
LOGICAL :: cosp_csat_vgrid = .TRUE.

!- Other control variables
INTEGER,PARAMETER :: COSP_NCOLUMNS_MAX = 100


!- Inputs related to radar simulations
REAL    :: cosp_radar_freq     = 94.0
INTEGER :: cosp_surface_radar  = 0
INTEGER :: cosp_use_mie_tables = 0
INTEGER :: cosp_use_gas_abs    = 1
INTEGER :: cosp_do_ray         = 0
INTEGER :: cosp_melt_lay       = 0
REAL    :: cosp_k2             = -1
LOGICAL :: cosp_use_reff       = .TRUE.
LOGICAL :: cosp_use_precipitation_fluxes = .TRUE.
!- Inputs related to lidar simulations
INTEGER :: cosp_nprmts_max_hydro = 12
INTEGER :: cosp_naero            = 1
INTEGER :: cosp_nprmts_max_aero  = 1
INTEGER :: cosp_lidar_ice_type   = 0
!- Inputs related to ISCCP simulator
INTEGER :: cosp_overlap                   = 3
INTEGER :: cosp_isccp_topheight           = 1
INTEGER :: cosp_isccp_topheight_direction = 2
REAL    :: cosp_emsfc_lw                  = 0.99
!- Inputs related to RTTOV
INTEGER :: cosp_platform   = 1
INTEGER :: cosp_satellite  = 15
INTEGER :: cosp_instrument = 0
INTEGER,PARAMETER :: COSP_NCHANNELS = 8
INTEGER :: cosp_channels(cosp_nchannels) = (/1,3,5,6,8,10,11,13/)
REAL    :: cosp_surfem(cosp_nchannels) = (/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
REAL    :: cosp_zenang = 50.0
REAL    :: cosp_co2    = 5.241E-04
REAL    :: cosp_ch4    = 9.139E-07
REAL    :: cosp_n2o    = 4.665E-07
REAL    :: cosp_co     = 2.098E-07


NAMELIST/RUN_COSP/ cosp_cloudsat_sim, cosp_lidar_sim, cosp_isccp_sim, &
                   cosp_misr_sim, cosp_modis_sim, cosp_rttov_sim,     &
                   cosp_npoints_it, cosp_ncolumns, cosp_use_vgrid,    &
                   l_cosp

END MODULE cosp_input_mod

