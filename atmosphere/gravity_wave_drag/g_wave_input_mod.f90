! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE g_wave_input_mod

! Description:
!       Input/namelist control of GWD scheme.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: GWD

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

USE missing_data_mod, ONLY: rmdi,imdi

IMPLICIT NONE

! Control logicals
LOGICAL :: l_gwd      = .FALSE.  ! Use orographic drag scheme
LOGICAL :: l_use_ussp = .FALSE.  ! Use spectral GWD scheme


LOGICAL :: l_gwd_40km = .TRUE.   
                              ! Turn off orographic GWD above 40km

LOGICAL :: l_nonhydro = .FALSE. 
                              ! Turn on nonhydro scheme
LOGICAL :: l_dynbeta  = .FALSE.  
                              ! Turn on dynamically adjusting angle 
                              ! of group vel in nonhydro scheme
                              
LOGICAL :: l_taus_scale = .FALSE.
                              ! Froude No. Dependent surface stress  
                              
LOGICAL :: l_fix_gwsatn = .TRUE. 
                              ! Bug fixes to gwsatn                                                         
INTEGER :: sat_scheme = 0 
                              ! Switch to determine whether to use
                              ! stress based    (sat_scheme=0) or
                              ! amplitude based (sat_scheme=1)
                              ! saturation test

! GWD variable inputs

INTEGER :: i_gwd_vn = imdi 
                              ! Switch to determine version of GWD scheme
                              ! 4 => 4A scheme
                              ! 5 => 5A scheme
INTEGER, PARAMETER :: i_gwd_vn_4a = 4
INTEGER, PARAMETER :: i_gwd_vn_5a = 5

REAL :: kay_gwave   =  rmdi   ! Surface stress constant for GW drag
REAL :: gwd_frc     =  rmdi   ! Critical Froude Number for 4A Scheme

! Maximum and minimum values for the STPH_RP scheme
! Gravity Wave drag parameters moved to stochastic_physics_run_mod.F90
! Max and min values for the critical froud number
! Max and min values for the gravity wave constant

REAL :: gwd_fsat  = rmdi      ! Critical Froude number for wave breaking
                              ! used in the amplitude based saturation
                              ! test


!New parameters for 5A version of orographic drag scheme
REAL    :: gsharp  = rmdi     ! Function of mtn sharpness (used to tune 
                              ! amplitude of orographic gwd) 
REAL    :: fbcd = rmdi        ! Flow blocking drag coefficient

LOGICAL :: l_smooth     = .FALSE.
                              ! Turn on smoothing of acceleration over  
                              ! a vertical wavelength

LOGICAL :: l_gw_heating(3) = .FALSE.
                              ! Turn on for heating due to
                              ! (1) flow blocking drag
                              ! (2) mountain-wave dissipation
                              ! (3) non-orographic gravity-wave dissipation

! Non-orographic gravity wave (USSP scheme) variable inputs
!
REAL    :: ussp_launch_factor = rmdi
                              ! Factor enhancement for wave launch 
                              ! amplitude

LOGICAL :: l_ussp_opaque = .FALSE.
                              ! OPAQUE LID FOR USSP SCHEME

      NAMELIST/run_gwd/i_gwd_vn                                         &
                      ,kay_gwave, gwd_frc, gwd_fsat                     &
                      ,gsharp, fbcd, l_smooth, l_gw_heating             &
                      ,l_ussp_opaque, ussp_launch_factor                &
                      ,l_gwd, l_use_ussp

END MODULE g_wave_input_mod

