! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to hold SCM fixed parameters for dimensioning input variables
! read in via namelists.

MODULE s_maxdim

  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Holds fixed dimensions for variables read in via SCM namelists. Fotran 90
!   does not allow namelist array members to have varying dimensions, so they
!   must be dimensioned using fixed parameters.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!=============================================================================

  ! Parameters

  INTEGER, PARAMETER :: &
    mx_rw      = 1      &
  , mx_rw_lng  = 1      &
  , mx_nlnd    = 1      &
  , mx_mod_lv  = 300    &
  , mx_ntile   = 20     &
  , mx_npft    = 10     &
  , mx_ntype   = 20     &
  , mx_st_lv   = 10     &
  , mx_sm_lv   = 10     &
  , mx_wet_lv  = 300    &
  , mx_tr_lv   = 300    &
  , mx_tr_vars = 300    &
  , mx_o3_lv   = 300    &
  , mx_nobs    = 1500   &
  , mx_nsprog  = 15     &
  , mx_nsmax   = 3      &
  , dim_cs1    = 1      &
  , dim_cs2    = 1

!=============================================================================
END MODULE s_maxdim
