! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Parameters for deep turbulence convection scheme
! 
MODULE dts_fitpars_mod

!------------------------------------------------------------------------------
! Description:
!   Module containing parameters used by the deep turbulence convection code.
!
! Method:
!   Default values have been declared where appropriate.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP 3 programming standards version vn8.2.
!------------------------------------------------------------------------------

IMPLICIT NONE
SAVE

!------------------------------------------------------------------------------
! Parameters for calculating mb, the mass flux at cloud base
! mb = a1*wstar + a2
!------------------------------------------------------------------------------

  REAL, PARAMETER  :: a1      = 0.16
  REAL, PARAMETER  :: a2      = -0.1

! mb = a3*wstar

  REAL, PARAMETER :: a3 = 0.02 ! closures for small wstar

!------------------------------------------------------------------------------
! Parameters for calculating mfr, the mass flux at the freezing level
! mfr = b1+b2*wcld
!------------------------------------------------------------------------------

  REAL, PARAMETER  :: b1 = -0.17
  REAL, PARAMETER  :: b2 = 0.043

! linear behaviour taking mfr down to zero for lower values of wall

  REAL, PARAMETER  :: b3 = 2.3e-3

!------------------------------------------------------------------------------
! Parameters for calculating the conversion between w^3 and w^2 -- 
! linking skewness to variance from eqn 13 velocity skewness to variance
!------------------------------------------------------------------------------

  REAL, PARAMETER :: c1 = 0.189
  REAL, PARAMETER :: c2 = -0.045

!------------------------------------------------------------------------------
! Parameters for calculating qicemax-relating maximum ice amount to qv/qsat
! above freezing level
!------------------------------------------------------------------------------

  REAL, PARAMETER :: d1 = 7.1e-4
  REAL, PARAMETER :: d2 = -5.4e-4

!------------------------------------------------------------------------------
! Parameters for calculating qrainmax = e1+e2*sigma 
!(may want this to be land/sea dependent)
!------------------------------------------------------------------------------

! These are the sea parameters
  REAL, PARAMETER :: e1 = -6.9e-6
  REAL, PARAMETER :: e2 = 0.004

! These are the land parameters
  REAL, PARAMETER :: land1 = -6.0e-5
  REAL, PARAMETER :: land2 = 8.57e-4

!------------------------------------------------------------------------------
! Parameters for calculating sigma, the fractional cloud area at cloud base
!------------------------------------------------------------------------------

  REAL, PARAMETER :: f1 = 0.04
  REAL, PARAMETER :: f2 = 2.1
     
!------------------------------------------------------------------------------
! used in dts_wthv,dts_qflux
!------------------------------------------------------------------------------

  REAL, PARAMETER :: h1 = 50.0  !50  s used in wthv,qflux
  REAL, PARAMETER :: h2 = 4.0   ! 3.5 !27/2/093.5 ! obsolete

  REAL, PARAMETER :: h4 = 0.0018   ! used in dts_wthv
  REAL, PARAMETER :: j1 = 0.8      ! 0.8 should be 0.5 ! used in dts_qflux

!------------------------------------------------------------------------------
! used in dts_w_variance:
!------------------------------------------------------------------------------

  REAL, PARAMETER :: wamp = 3.0 ! scale factor for ww

!------------------------------------------------------------------------------
! Temperature at which the plume water content is entirely ice
!------------------------------------------------------------------------------

  REAL, PARAMETER :: tcritplume = -10.0 ! degrees C

!------------------------------------------------------------------------------
! Factor used in the parametrisation of revp (rate of evaporation of rain)
! The 0.005 is the proportionality factor if normalising by total rain amount,
! but I need to normalise by rain production rate, so a timescale of 600s is 
! needed
! The 0.005 factor appears to be quite a good fit for several different LEM 
! simulations,  the 600s is a compromise between different runs, and not such
! a good fit...
!------------------------------------------------------------------------------
! The following constants are used in dts_rainprod
!------------------------------------------------------------------------------

  REAL, PARAMETER :: alfac = 3.5    !3. ! 0.005*600

! NEXT 2 values currently unused
  REAL, PARAMETER :: bratio = 0.9   ! Must be a value less than or equal to 1.
                                    ! 1-bratio is the proportion of relative
                                    !  humidity dependence above the lcl for 
                                    ! calculating rain evaporation
  REAL, PARAMETER :: rhconst = 0.95 ! Also must be a value less than or equal
                                    ! to 1.0 assume 90% humidity near downdrafts
                                    ! for use in rain evaporation calculation 


!------------------------------------------------------------------------------
! Used in melt routine
!------------------------------------------------------------------------------

  REAL, PARAMETER :: tsig = 7 ! K  
  REAL, PARAMETER :: taumelt = 400.0 ! s but not important, as renormalise


!------------------------------------------------------------------------------

END MODULE  dts_fitpars_mod
