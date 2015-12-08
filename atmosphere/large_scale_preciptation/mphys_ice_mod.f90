! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE mphys_ice_mod

! Description:
! Holds Ice constants required by the large-scale
! precipitation scheme
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

  IMPLICIT NONE

      !-----------------------------------------------------------------
      ! Nucleation of ice
      !-----------------------------------------------------------------

      ! Note that the assimilation scheme uses temperature thresholds
      ! in its calculation of qsat.

      ! Nucleation mass
      REAL, PARAMETER :: m0    =  1.0E-12

      ! Maximum temperature for homogenous nucleation
      REAL, PARAMETER :: thomo = -40.0

      !  1.0/Scaling quantity for ice in crystals
      REAL, PARAMETER :: qcf0  = 1.0E4    ! This is an inverse quantity

      ! Minimum allowed QCF after microphysics
      REAL,PARAMETER:: qcfmin  = 1.0E-8

      ! 1/scaling temperature in aggregate fraction calculation
      REAL, PARAMETER :: t_scaling = 0.0384

      !  Minimum temperature limit in calculation  of N0 for ice (deg C)
      REAL, PARAMETER :: t_agg_min = -45.0

      !-----------------------------------------------------------------
      ! Hallett Mossop process
      !-----------------------------------------------------------------

      ! Switch off Hallett Mossop in this version but allow
      ! functionality

      ! Min temp for production of Hallett Mossop splinters (deg C)
      REAL, PARAMETER :: hm_t_min = -8.0

      ! Max temp for production of Hallett Mossop splinters (deg C)
      REAL, PARAMETER :: hm_t_max = -3.0

      !  Residence distance for Hallett Mossop splinters (1/deg C)
      REAL, PARAMETER :: hm_decay = 1.0 / 7.0

      ! Reciprocal of scaling liquid water content for HM process
      REAL, PARAMETER :: hm_rqcl = 1.0  / 0.1E-3


END MODULE mphys_ice_mod
