! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Processor parallel properties for 2d comms

      MODULE comms_2d_mod
      IMPLICIT NONE

! Description:
!
! Method:
!         The required parameters are copied from the parallel constants
!         in the routine init_comms_2d which is called from
!                                  control/top_level/setcona
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

      LOGICAL :: l_2dcomm      ! T = doing 2D comms on demand

      INTEGER :: size_2dcomm   ! no of procs involved in comms on demand
      INTEGER :: group_2dcomm  ! group ID of procs involved in comms on demand
      INTEGER :: max_comm_size ! error check size comms on demand
      INTEGER :: comu_size     ! error check size comms on demand u
      INTEGER :: comv_size     ! error check size comms on demand v
      INTEGER :: comw_size     ! error check size comms on demand vec_w
      INTEGER :: comm_size     ! error check size comms on demand thmono

      END MODULE comms_2d_mod
