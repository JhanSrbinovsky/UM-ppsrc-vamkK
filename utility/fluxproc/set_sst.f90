! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: set_sst
!
! Purpose: Flux processing routine.
!          Convert a surface temperature field to a sea surface
!          temperature field by setting all values less than the
!          freezing point of seawater to the freezing point (-1.8degC)
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
       SUBROUTINE set_sst                                               &
     &           (ncols, nrows,                                         &
     &            fieldSST, rmdi,                                       &
     &            out_field)




      IMPLICIT NONE


!     Input:
!     ------

      integer nrows         ! IN number of rows of array
      integer ncols         ! IN number of columns in array

      real fieldSST(ncols,nrows)   ! IN array of input values (NWPfield)
      real rmdi                    ! IN missing data indicator

!     Output:
!     -------


      real out_field(ncols,nrows)  ! OUT composite array of NWP and
                                   ! Climatology values


!    Local variables
!    ---------------


      integer i        ! Loop counter over columns
      integer j        ! Loop counter over rows
      real sst_under_ice   ! SST under ice = -1.8

      parameter ( sst_under_ice = -1.8 )


! ------------------------------------------------------------------

! 1.Loop over each element in field
!     and check if SST element is less than sst_under_ice
      do j = 1,nrows            ! Loop over rows
        do i = 1,ncols         ! Loop over columns
          if ( fieldSST(i,j)  /=  rmdi ) then

            if ( fieldSST(i,j)  <   sst_under_ice ) then

! 1.2 If test is true, set SST for that element
!    else use read in field element
              out_field(i,j) = sst_under_ice
            else
              out_field(i,j) = fieldSST(i,j)
            endif

          else
            out_field(i,j) = rmdi
          endif
        enddo      ! i
      enddo        ! j
      return
      END SUBROUTINE set_sst
!----------------------------------------------------------------------
