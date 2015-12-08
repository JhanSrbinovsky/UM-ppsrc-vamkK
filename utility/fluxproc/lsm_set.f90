! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: lsm_set
!
! Purpose: Flux processing routine.
!          Sets missing data values in a field according to the lsm
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine lsm_set (ncols, nrows, lsm, ivalue_mask,               &
     &                    rmdi, ldebug, field )

      implicit none

! declaration of argument list

      integer ncols  ! IN number of columns in field
      integer nrows  ! IN number of rows in field
      integer lsm(ncols, nrows) ! IN land / sea mask
      integer ivalue_mask ! IN integer value in lsm for which to set a
                          !    missing value in array "field"
      real rmdi           ! IN missing data value
      logical ldebug      ! IN T => debug output is produced
      real field(ncols, nrows) ! IN/OUT field


! no parameters, globals or local arrays

! declaration of local scalars
      integer i, j   ! do loop indices

! externals
      external output_debug
!----------------------------------------------------------------------

! 1. Set values

      do j = 1, nrows
        do i = 1, ncols
          if ( lsm(i, j)  ==  ivalue_mask ) then
            field(i,j) = rmdi
          end if
        end do
      end do

! 2. output selected values for debugging

      if (ldebug) then
! DEPENDS ON: output_debug
        call output_debug('values after land-sea mask imposed',         &
     &                    nrows, ncols, Field)
      end if

      return
      END SUBROUTINE lsm_set
!----------------------------------------------------------------------
