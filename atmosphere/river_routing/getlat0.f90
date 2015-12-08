! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  returns latituds at (iy) in (nla)
!
!     NOTICE: this is for GPCP/Xie and Arkin data set
!     which data is located at the center of 2.5 grid box
!
!     iy from SOUTH to NORTH.
!
!     from 23.Feb.1996, by Taikan OKI
!
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: River Routing

      REAL FUNCTION getlat0(iy, nla)
      IMPLICIT NONE
!
      integer iy, nla
!
      getlat0 = 180.0 * (real(iy)-0.5) / real(nla) - 90.0
!

      END FUNCTION getlat0

