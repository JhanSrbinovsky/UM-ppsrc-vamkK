! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


      real function sind(degree)
!
!     obtain sine for degree
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: River Routing
      IMPLICIT NONE
      real degree, pi, sin
      data pi/3.141593/
!
      sind = sin (pi*degree/180.0)
!
      END FUNCTION sind

