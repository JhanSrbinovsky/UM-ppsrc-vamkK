! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Global data module for variables concerned with solar incidence.

MODULE solinc_data

  IMPLICIT NONE
  SAVE

! Description:
!   Global data necessary for calculating the angle of solar incidence
!   on sloping terrain.
!
! Method:
!   Provides global data.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):

  REAL, ALLOCATABLE, DIMENSION(:,:)   :: slope_aspect, slope_angle
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: sol_bearing, f_orog, orog_corr
  REAL, ALLOCATABLE, DIMENSION(:)     :: lg_orog_corr, lg_f_orog
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: horiz_ang, horiz_aspect
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: sky
  LOGICAL :: l_orog = .FALSE.
  LOGICAL :: l_skyview = .FALSE.
  INTEGER :: horiz_limit = 30
  INTEGER, PARAMETER :: n_horiz_ang = 16

!Variables defined as global in this module need to be made private
!to individual threads: lg_orog_corr and lg_f_orog. 
!$OMP THREADPRIVATE (lg_orog_corr, lg_f_orog)

! slope_aspect: The direction faced by the mean slope - i.e. the
!               bearing of the slope normal projected on the surface
!               measured in radians clockwise from true north.
!
! slope_angle:  Angle of the mean slope measured in radians from
!               the horizontal.
!
! sol_bearing:  Mean local bearing of the sun over the timestep
!               expressed in radians clockwise from grid north.
!
! orog_corr:    correction factor for the direct solar flux
!               reaching the surface for sloping terrain.
!
! lg_orog_corr: orog_corr gathered at lit points
!
! lg_f_orog:    The extra direct solar flux at the surface due to
!               the orography correction. This is used in the
!               correction to the sw_incs calculation.
!
! f_orog:       lg_f_orog ungathered onto the full grid. This is
!               used to correct net_atm_flux.
!
! l_orog:       model switch for orography scheme
!
! horiz_ang:    Angle in radians measured from the local zenith to
!               the obscuring terrain.
!
! horiz_aspect: The local bearing for each horizon angle measured
!               in radians clockwise from true north.
!
! sky:          Sky-view correction factor for net surface LW.
!
! l_skyview:    Model switch for skyview scheme.
!
! horiz_limit:  Number of grid-lengths to furthest point considered
!               in horizon calculation.
!
! n_horiz_ang:  Number of horizon angles calculated.
!
!- End of header

END MODULE solinc_data
