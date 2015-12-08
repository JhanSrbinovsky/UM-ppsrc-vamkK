! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Init_proc_info
!
      SUBROUTINE Init_proc_info(                                        &
                                g_row_length, g_rows, domain,           &
                                mype, nproc, nprocx, nprocy,            &
                                proc_col_group, proc_row_group,         &
                                datastart, extremity )

       USE yomhook, ONLY: lhook, dr_hook
       USE parkind1, ONLY: jprb, jpim

       USE proc_Info_Mod

! Purpose:
!          Initialises parameters needed by processors to use
!          parallel gcom routines in l2 norm calculations
!          and other calculations (max, mins, sums etc)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      IMPLICIT NONE

! Input variables.

      LOGICAL :: extremity(4)   ! Processor position is at north,
                                ! south, east or west of domain
      INTEGER  ::  datastart(3)
      INTEGER  ::  g_row_length
      INTEGER  ::  g_rows
      INTEGER  ::  proc_col_group
      INTEGER  ::  proc_row_group
      INTEGER  ::  domain
      INTEGER  ::  mype            ! processor identifier
      INTEGER  ::  nproc
      INTEGER  ::  nprocx
      INTEGER  ::  nprocy

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('INIT_PROC_INFO',zhook_in,zhook_handle)

!  Just copy to rename for module

      at_extremity = extremity

      l_datastart = datastart
      global_row_length = g_row_length
      global_rows = g_rows
      gc_proc_col_group = proc_col_group
      gc_proc_row_group = proc_row_group
      model_domain = domain
      me = mype
      n_proc = nproc
      n_procx = nprocx
      n_procy = nprocy

      IF (lhook) CALL dr_hook('INIT_PROC_INFO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Init_proc_info

