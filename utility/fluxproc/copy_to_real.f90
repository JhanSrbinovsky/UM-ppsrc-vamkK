! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines:copy_to_real, copy_to_integer
!
! Purpose: Flux processing routines.
!          Copy_to_real: copies an integer to a real number.
!          Copy_to_integer: copies a real to an integer number.
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine copy_to_real(ValueIn, ValueOut)

      implicit none

! declaration of argument list
      REAL ValueIn    ! value to be copied (usually an integer is
                      ! passed to this routine)
      REAL ValueOut   ! output value

      ValueOut = ValueIn

      return
      END SUBROUTINE copy_to_real

!----------------------------------------------------------------------
