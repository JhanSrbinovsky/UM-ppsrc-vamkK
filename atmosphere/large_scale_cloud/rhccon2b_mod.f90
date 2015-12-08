! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
MODULE rhccon2b_mod

IMPLICIT NONE
! RHCCON2B start
!     Comdeck for use with the RHcrit parametrization of the large-scale
!     cloud scheme, A09_2B.
!     Four constants are used to specify the variable (a function of
!     pressure) which relates the variability of the saturation variable
!     in one box to the variability over 9 climate grid boxes. The
!     variability of the saturation variable in one box is required to
!     specify RHcrit.
!     Note that the constants
!       RHC_CON1=0.522, RHC_CON2=0.122, RHC_CON3=2.5E3, RHC_CON4=1.75E4
!     are only suitable for use in the 2.5*3.75 degrees climate model:
!     these constants depend upon the size of a grid-box.
!
!     The fit is of the form:
!          A=RHC_CON1+RHC_CON2*(p-RHC_CON4)/(RHC_CON3+abs(p-RHC_CON4))
!       where p is the pressure at the layer midpoint.
!     Then,
!           sigma(s) = A * sigma(s,9)
!        where sigma(s) is the std dev of the saturation variable s in
!     one grid-box, and sigma(s,9) is the std dev over 9 boxes.
!
!     RHC_MIN and RHC_MAX are user defined limits on the values of RHc.
!
!     S. Cusack   02-09-98
!
      REAL, PARAMETER :: RHC_CON1 = 0.522
      REAL, PARAMETER :: RHC_CON2 = 0.122
      REAL, PARAMETER :: RHC_CON3 = 2.5E3
      REAL, PARAMETER :: RHC_CON4 = 1.75E4
      REAL, PARAMETER :: RHC_MIN = 0.3
      REAL, PARAMETER :: RHC_MAX = 0.98

! RHCCON2B end

END MODULE rhccon2b_mod
