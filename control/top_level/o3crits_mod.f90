MODULE o3crits_mod

IMPLICIT NONE
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! ----------------------- Header file O3CRITS  ------------------------
! Description: Parameters of ozone tropopause criteria.
!----------------------------------------------------------------------
!     ! Criteria for establishing the level of the
!     ! ozone tropopause to enable mapping of ozone concentrations
!     ! in ozone ancillary file to model levels using the
!     ! thermal tropopause level in the model for scaling.

!     ! The gradient of O3 concentration has to exceed 60 ppbv per km
!     ! i.e. 99.423E-09 (= 60ppbv) / 1000 to convert to conc per meter
      Real, Parameter :: O3_grad_crit    = 99.423E-12

!     ! The concentration of O3 has to exceed 80 ppbv
      Real, Parameter :: O3_conc_crit    = 132.56E-09

!     ! The conc. of O3 has to exceed 110 ppbv in the stratosphere
      Real, Parameter :: O3_strat_crit   = 182.27E-09
! O3CRITS end

END MODULE o3crits_mod
