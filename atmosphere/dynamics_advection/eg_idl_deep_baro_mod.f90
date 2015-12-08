! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Parameters for deep atmosphere baroclinic wave test

MODULE eg_idl_deep_baro_mod

USE earth_constants_mod,   ONLY: Earth_radius
USE conversions_mod,       ONLY: pi

IMPLICIT NONE

!
! Description:
!
!   Parameters for Staniforth deep atmosphere
!   baroclinic wave test
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

REAL, PARAMETER :: u0 = 35.0        ! max wind (m/s)
REAL, PARAMETER :: eg_lapse = 0.005 ! tropo lapse rate (K/m)
REAL, PARAMETER :: rp = 0.1*Earth_radius           
                                    ! e-fold radius of perturbation (m)
REAL, PARAMETER :: xc = pi/9.0      ! lon of perturbation centre (rad)
REAL, PARAMETER :: yc = 2*pi/9.0    ! lat of perturbation centre (rad)
REAL, PARAMETER :: up = 1.0         ! Magnitude of perturbation (m/s)
REAL, PARAMETER :: taper_top = 15000.0 ! Top of u perturbation

REAL, SAVE :: T0
REAL, SAVE :: H, B, C
REAL, SAVE, ALLOCATABLE :: tau1(:), tau2(:),                          &
                           tau1_int(:), tau2_int(:)

INTEGER, SAVE :: k_const_save, b_const_save

LOGICAL, SAVE :: L_shallow_save

CONTAINS

SUBROUTINE eg_init_idl_deep_baro()

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod, ONLY: tdims

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('eg_init_idl_deep_baro',zhook_in,zhook_handle)

ALLOCATE ( tau1(tdims%k_start:tdims%k_end) )
ALLOCATE ( tau2(tdims%k_start:tdims%k_end) )
ALLOCATE ( tau1_int(tdims%k_start:tdims%k_end) )
ALLOCATE ( tau2_int(tdims%k_start:tdims%k_end) )

IF (lhook) CALL dr_hook('eg_init_idl_deep_baro',zhook_out,zhook_handle)

RETURN

END SUBROUTINE eg_init_idl_deep_baro

END MODULE eg_idl_deep_baro_mod
