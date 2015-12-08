! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Constants used by convective diagnosis routines.
!
MODULE cv_diag_param_mod

!-----------------------------------------------------------------------
! Description:
!   Module containing parameters used by the convective diagnosis routines
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3  v8 programming standards.
!
!-----------------------------------------------------------------------

IMPLICIT NONE

! Parameters used in determining the thetav perturbation given to the parcel
! for the undilute ascent.
 
REAL, PARAMETER ::      &
  a_plume = 0.2         & ! Minimum initial parcel dthetav
 ,b_plume = 3.26        & ! Used in initial parcel perturbation cal.
 ,max_t_grad = 1.0E-3     ! Used in initial parcel perturbation cal.


! Parameter used in cumulus testing to see if a cloud layer present

REAL, PARAMETER ::      &
  sc_cftol   = 0.1        ! Cloud fraction required for a cloud layer 
                          ! to be diagnosed.

! Parameters for the calculation of T at the lifting condensation level.
! See Bolton 1980 Mon Wea Rev 108, P1046-1053 for details of empirical 
! relationship

REAL, PARAMETER ::      &
  a_bolton = 55.0       &
, b_bolton = 2840.0     &
, c_bolton = 3.5        &
, d_bolton = 4.805

END MODULE cv_diag_param_mod
