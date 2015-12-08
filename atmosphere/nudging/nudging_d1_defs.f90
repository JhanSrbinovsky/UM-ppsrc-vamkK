! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Module used to store information about variables transferred
!  to/fro D1 array.

!  Part of the Nudged model (see nudging_main.F90)

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
MODULE nudging_d1_defs

IMPLICIT NONE

! ************************************************************
! Section & Item information for prognostic and diagnostic variables
INTEGER, PARAMETER :: sect_prog  = 0           ! Section No. for prognostic

INTEGER, PARAMETER :: sect_tropht= 30          ! Section No. for climate diags

INTEGER, PARAMETER :: sect_nudge = 39          ! Section No. for nudging

! Unique indices in the code type array
! Each one refers to the location of a variable in it
INTEGER, PARAMETER :: u_index       = 2
INTEGER, PARAMETER :: v_index       = 3
INTEGER, PARAMETER :: theta_index   = 4
INTEGER, PARAMETER :: exner_index   = 406
INTEGER, PARAMETER :: p_rho_index   = 407
INTEGER, PARAMETER :: p_theta_index = 408
INTEGER, PARAMETER :: trop_pres_index = 451

INTEGER, PARAMETER :: tdiag_data_index   = 1
INTEGER, PARAMETER :: tdiag_model_index  = 2
INTEGER, PARAMETER :: tdiag_ntend_index  = 3
INTEGER, PARAMETER :: tdiag_mtend_index  = 4
INTEGER, PARAMETER :: tdiag_relax_index  = 5

INTEGER, PARAMETER :: udiag_data_index   = 6
INTEGER, PARAMETER :: udiag_model_index  = 7
INTEGER, PARAMETER :: udiag_ntend_index  = 8
INTEGER, PARAMETER :: udiag_mtend_index  = 9
INTEGER, PARAMETER :: udiag_relax_index  = 10

INTEGER, PARAMETER :: vdiag_data_index   = 11
INTEGER, PARAMETER :: vdiag_model_index  = 12
INTEGER, PARAMETER :: vdiag_ntend_index  = 13
INTEGER, PARAMETER :: vdiag_mtend_index  = 14
INTEGER, PARAMETER :: vdiag_relax_index  = 15

END MODULE nudging_d1_defs
