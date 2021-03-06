! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


      real function cosd(degree)
!
!     obtain sine for degree
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: River Routing
      IMPLICIT NONE
      real degree, pi, cos
      data pi/3.141593/
!
      cosd = cos (pi*degree/180.0)
!
      END FUNCTION cosd
