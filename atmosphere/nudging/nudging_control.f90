! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Module containing most important switches for the nudged model
!  All minor changes to the nudging should be made here

!  During the development of the nudging  the analyses were obtained
!  from two sources; directly from the BADC website and stored on
!  BDAN. Each had slightly different set-ups (eg different names).
!  Although BDAN data alone is expected to be used the switches remain
!  for flexibility

!  Part of the Nudged model (see nudging_main.F90)

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
MODULE nudging_control

USE Control_Max_Sizes

USE nudging_input_mod

IMPLICIT NONE

!**********************************************************

INTEGER            :: glob_debug = 0           ! Global debug flag
INTEGER, PARAMETER :: no_debug   = 0           ! No debug flag
INTEGER, PARAMETER :: OUT        = 6           ! Output unit
CHARACTER (LEN=256) :: nmessage                ! String for errors

! flag to limit large extrapolation.
! If set to true cuts off extrapolation if > 2 or < -1
! Not important in free atmosphere and recommended to be true as
! prevents strange numbers being produced near surface
LOGICAL, PARAMETER        :: l_extrapcutoff     = .TRUE.

!*********************************************************
! Define resolution of analysis data (standard N48)
INTEGER, PARAMETER        :: ana_rows           = 73
INTEGER, PARAMETER        :: ana_row_length     = 96
REAL, PARAMETER           :: ana_base_lambda    = 0.000
REAL, PARAMETER           :: ana_base_lambda_u  = 1.875
REAL, PARAMETER           :: ana_delta_lambda   = 3.750
REAL, PARAMETER           :: ana_base_phi       = -90.000
REAL, PARAMETER           :: ana_base_phi_v     = -88.750
REAL, PARAMETER           :: ana_delta_phi      =  2.500

!*********************************************************
! This section  relates to aspects of the netcdf file
! Included mainly for historical reasons. Data is now
! required to be in data format 0
! Where is the data from
!     (0=BDAN/Reading, 1 = BADC) This affects a couple of things
!  i) Data starts at top of the atmosphere in BDAN data
!     so have to flip the data to get in desired (BADC) format
!  ii) surface pressure is 2d+1 rather than just 2D so have to
!     use slightly different routines
!  iii) naming conventions differ (see below)
INTEGER, PARAMETER:: data_source           = 0

! Total no. of dimensions read from netcdf file
INTEGER, PARAMETER:: totdims               = 6

! string array containing dimension names in the netcdf file
! Note that these are the BDAN naming conventions
! The BADC conventions are slightly different
INTEGER, PARAMETER        :: dimname_length        = 11

CHARACTER(LEN=dimname_length), PARAMETER, DIMENSION(totdims)::     &
 dim_names = (/                                                    &
  'longitude  ', 'latitude   ', 'hybrid     ',                     &
  'hybrid_1   ', 'longitude_1', 'latitude_1 '                      &
/)

! Data pathname and filenames: Currently Hard Wired
! eg `ecm-e40_1deg-model-levs_yyyymmddhh_t.nc' for temperature
! Note use this convention for both BDAN and BADC
CHARACTER(LEN=*), PARAMETER :: file_model_basis =                    &
 'ecm-e40_1deg-model-levs_'                   ! filename starting
CHARACTER(LEN=*), PARAMETER :: file_pres_basis  =                    &
 'ecm-e40_1deg-model-levs_'                   ! pressure levels are
                                              ! not supported at present
CHARACTER(LEN=*), PARAMETER :: file_end   =  '_all.nc'
                                              ! filename ending

! UM character strings
CHARACTER(LEN=*), PARAMETER   :: file_um_basis = 'HP'
CHARACTER(LEN=*), PARAMETER   :: file_um_end   = '00.GLOUM6.nc'

! JRA character strings
CHARACTER(LEN=*), PARAMETER   :: file_jra_basis = 'JRA'
CHARACTER(LEN=*), PARAMETER   :: file_jra_end   = '_all.nc'

! These allow conversion from array to individual variables
! Defined here for consistency with array names
INTEGER, PARAMETER:: dimindex_lon  = 1
INTEGER, PARAMETER:: dimindex_lat  = 2
INTEGER, PARAMETER:: dimindex_lev  = 3
INTEGER, PARAMETER:: dimindex_tim  = 4
INTEGER, PARAMETER:: dimindex_ulon = 5
INTEGER, PARAMETER:: dimindex_vlat = 6

! Names of variables that we load in from the netcdf files.
! Note these are the BDAN (rather than BADC) default names
CHARACTER(LEN=*), PARAMETER :: u_name      = 'U'
CHARACTER(LEN=*), PARAMETER :: v_name      = 'V'
CHARACTER(LEN=*), PARAMETER :: temp_name   = 'T'
CHARACTER(LEN=*), PARAMETER :: surfp_name  = 'LNSP'
CHARACTER(LEN=*), PARAMETER :: usurfp_name = 'LNSP_1'
CHARACTER(LEN=*), PARAMETER :: vsurfp_name = 'LNSP_2'

!***********************************************************
! Misecllanea

! When loading non prognostic variables choose a tag value
! (eg in UKCA code variables are tagged witha  98)
! At present taking a job that uses a stochem tag (99)
INTEGER, SAVE             :: nudging_tagvalue = 99
INTEGER                   :: dbg_x,dbg_y,dbg_z
                             ! Location for debug prints

END MODULE nudging_control
