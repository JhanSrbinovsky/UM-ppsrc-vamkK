! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Sizes for many of the UM's main, dynamic data arrays
!
MODULE nlsizes_namelist_mod

! Description:
!  Module-based interface to the NLSIZES namelist and associated declarations.
!  Contains the sizes needed for the dynamic allocation of the main data arrays
!  within the model.
!
! Method:
!  Sizes read in via the NLSIZES namelist are SAVEd here.
!  Many of the declarations here are currently duplicated from typsize.h.
!  As an interim measure, the reconfiguration exclusively uses this module 
!  while the UM uses the older include file.
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Grids
!
! JULES variable
USE ancil_info, ONLY:                                             &
  nsmax

USE lbc_mod, ONLY:                                                &
  rimwidtha, nrim_timesa

IMPLICIT NONE

SAVE

! Main sizes of fields for each submodel
! Grid-related sizes for ATMOSPHERE submodel
INTEGER :: global_row_length    ! Points per global row
INTEGER :: global_rows          ! No of global (theta) rows
INTEGER :: model_levels         ! No of model levels
INTEGER :: land_field           ! No of land points in field
INTEGER :: ntiles               ! No of land surface tiles
INTEGER :: nice                 ! No. of sea ice thickness categories
INTEGER :: nice_use             ! No. of sea ice categories used fully
                                !  in surface exchange and radiation
                                !  (If nice>1 & nice_use=1, categories only
                                !  partially used in surface exchange)

! Physics-related sizes for ATMOSPHERE submodel
INTEGER :: wet_levels          ! No of moist-levels
INTEGER :: cloud_levels        ! No of cloud-levels
INTEGER :: st_levels           ! No of soil temperature levels
INTEGER :: sm_levels           ! No of soil moisture levels
INTEGER :: bl_levels           ! No of boundary-layer-levels
INTEGER :: ozone_levels        ! No of ozone-levels
INTEGER :: tpps_ozone_levels   ! No of tropopause-ozone-levels
INTEGER :: river_rows          ! No of rows for river routing
INTEGER :: river_row_length    ! Row length for river routing

! Dynamics-related sizes for ATMOSPHERE submodel
INTEGER :: tr_levels            ! No of tracer-levels
INTEGER :: tr_vars              ! No of passive tracers
INTEGER :: tr_lbc_vars          ! No of tracers in lbcs
INTEGER :: tr_ukca              ! No of UKCA tracers
INTEGER :: tr_lbc_ukca          ! No of UKCA tracer lbcs

! Grid related sizes for data structure
! Data structure sizes for ATMOSPHERE submodel
INTEGER :: a_len_inthd       ! Length of INTEGER header
INTEGER :: a_len_realhd      ! Length of REAL header
INTEGER :: a_len2_levdepc    ! No of LEVEL-dependent arrays
INTEGER :: a_len2_rowdepc    ! No of ROW-dependent arrays
INTEGER :: a_len2_coldepc    ! No of COLUMN-dependent arrays
INTEGER :: a_len2_flddepc    ! No of FIELD arrays
INTEGER :: a_len_extcnst     ! No of EXTRA scalar constants
INTEGER :: a_len_cfi1        ! Length of compressed fld index 1
INTEGER :: a_len_cfi2        ! Length of compressed fld index 2
INTEGER :: a_len_cfi3        ! Length of compressed fld index 3

! Data structure sizes for ATMOSPHERE ANCILLARY file control
! routines
INTEGER :: nancil_lookupsa  ! Max no of fields to be read

! Data structure sizes for ATMOSPHERE INTERFACE file control
! routines
INTEGER :: n_intf_a          ! No of atmosphere interface areas

!  Data structure sizes for ATMOSPHERE BOUNDARY file control
! routines:
! Sizes applicable to all configurations (DUMPS/FIELDSFILE)
INTEGER :: pp_len_inthd   ! Length of PP file integer header
INTEGER :: pp_len_realhd  ! Length of PP file real    header

NAMELIST/NLSIZES/                                                 &
  global_row_length, global_rows, land_field,                     &
  model_levels, wet_levels,                                       &
  cloud_levels, tr_levels, st_levels,                             &
  bl_levels, ozone_levels,                                        &
  pp_len_inthd, pp_len_realhd

NAMELIST/ATM_SIZES/                                               &
  nice, nsmax, nancil_lookupsa, nrim_timesa

END MODULE nlsizes_namelist_mod
