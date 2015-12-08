! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Options for high order interpolation in the semi-Lagrangian scheme.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE highos_mod

IMPLICIT NONE

! Tri-cubic Lagrange interpolation
INTEGER, PARAMETER :: cubicLagrange         = 1

! Tri-quintic Lagrange interpolation
INTEGER, PARAMETER :: quinticLagrange       = 2

! ECMWF quasi-cubic (linear interpolation on outer-most parts of
!                    a tri-cubic stencil)
INTEGER, PARAMETER :: ECMWF_quasiCubic      = 3

! ECMWF quasi-cubic with Bermejo & Staniforth monotonicity
INTEGER, PARAMETER :: ECMWF_mono_quasiCubic = 4

! Bi-cubic Lagrange interpolation in the horizontal; linear in
! the vertical
INTEGER, PARAMETER :: hCubic_vLin           = 5

! ECMWF quasi-cubic in horizontal; quintic Lagrange in the vertical
INTEGER, PARAMETER :: hQuasiCubic_vQuintic  = 6

! Bi-cubic Lagrange in the horizontal; quintic Lagrange in the vertical
INTEGER, PARAMETER :: hCubic_vQuintic       = 7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following are used by ENDGame only.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Bi-cubic Lagrange in horizontal; Hermite cubic, with quadratic
! derivative estimates, in the vertical:
INTEGER, PARAMETER :: hLag3_vHerm3_d2       = 8

! Bi-cubic Lagrange in horizontal; Hermite cubic, with quartic
! derivative estimates, in the vertical:
INTEGER, PARAMETER :: hLag3_vHerm3_d4       = 9

END MODULE highos_mod
