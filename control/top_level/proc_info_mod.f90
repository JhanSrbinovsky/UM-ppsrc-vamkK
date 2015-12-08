! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Processor parallel properties and pointers.

      MODULE proc_info_mod

! Description:
!          Parameters needed by processors to use
!          parallel gcom routines in l2 norm calculations
!          and other calculations (max, mins, sums etc)
!
! Method:
!         The required parameters are copied from the parallel constants
!         in the routine init_proc_info which is called from
!                                  control/top_level/setcona
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

      IMPLICIT NONE

       LOGICAL ::  at_extremity(4)  ! Processor position is at north,
                                    ! south, east or west of domain

      INTEGER  ::  l_datastart(3)
      INTEGER  ::  global_row_length
      INTEGER  ::  global_rows
      INTEGER  ::  gc_proc_col_group
      INTEGER  ::  gc_proc_row_group
      INTEGER  ::  model_domain
      INTEGER  ::  me            ! processor identifier
      INTEGER  ::  n_proc
      INTEGER  ::  n_procx
      INTEGER  ::  n_procy

      END MODULE proc_info_mod
