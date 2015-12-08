! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     returns longitude at (ix) in (nlo)
!
!     NOTICE: this is for GPCP/Xie and Arkin data set
!     which data is located at the center of 2.5 grid box
!
!     ix from west to east.

      REAL FUNCTION getlon0(ix, nlo)
      IMPLICIT NONE
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: River Routing
      integer ix, nlo
!
      getlon0 = 360.0 * (real(ix)-0.5) / real(nlo)
!

      END FUNCTION getlon0

