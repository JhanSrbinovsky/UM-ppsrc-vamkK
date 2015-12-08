! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine: LEVEL_RLEVEL ------------------------------------------
!
!  Purpose: To return a real value even though the routine is called
!  with integer arguments.
!
!  Tested under compiler:   cft77
!  Tested under OS version: UNICOS 5.1
!
!  Model            Modification history from model version 3.0:
!
!  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
!  Logical components covered: ...
!
!  Project task: ...
!
!  External documentation:
!
!  -------------------------------------------------------------------
!  Interface and arguments: ------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      SUBROUTINE LEVEL_RLEVEL(INT_LEVEL,REAL_LEVEL,REAL_LEVEL_OUT)
      IMPLICIT NONE
      INTEGER                                                           &
     &     INT_LEVEL              !    first dimension of the lookup
      REAL                                                              &
     &     REAL_LEVEL,                                                  &
                                  !    secnd dimension of the lookup
     &     REAL_LEVEL_OUT         !    secnd dimension of the lookup
!*
      REAL_LEVEL_OUT=REAL_LEVEL


      RETURN
      END SUBROUTINE LEVEL_RLEVEL
