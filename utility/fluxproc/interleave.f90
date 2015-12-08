! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: interleave
!
! Purpose: Flux processing routine.
!          To perform a check on each element of an input field.
!          If the test fails, that element shall be replaced by
!          the climatological value.
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
       SUBROUTINE interleave                                            &
     &           (ncols, nrows,                                         &
     &            fieldNWP, fieldClim,                                  &
     &            icefrac, rmdi,                                        &
     &            l_leads,out_field)




      IMPLICIT NONE


!     Input:
!     ------

      integer ncols         ! IN number of columns of array
      integer nrows         ! IN number of rows in array

      real minicefrac            ! minimum ice fraction
      real minleadsfrac          ! minimum leads fraction

      parameter ( minicefrac = 0.005 )
      parameter ( minleadsfrac = 0.005 )

      real fieldNWP(ncols,nrows)   ! IN array of input values (NWPfield)
      real fieldClim(ncols,nrows)  ! IN array of input values
                                   ! (Climatology)
      real icefrac(ncols,nrows)    ! IN array of input values (Icefrac)
      real rmdi                    ! IN missing data indicator

      logical l_leads        ! IN T => using minleadsfrac
                             !    F => using minicefrac

!     Output:
!     -------


      real out_field(ncols,nrows)  ! OUT composite array of NWP and
                                   ! Climatology values


!    Local variables
!    ---------------


      integer i        ! Loop counter over columns
      integer j        ! Loop counter over rows


! ------------------------------------------------------------------

! 1. Use l_leads to test whether using leads or ice frac
      if ( l_leads ) then
! 1.1 Loop over each element in field
!     and check if icefrac element is missing data
        do j = 1,nrows          ! Loop over rows
          do i = 1,ncols         ! Loop over columns
            if ( icefrac (i,j)  ==  rmdi ) then
              out_field (i,j) = fieldClim(i,j)
            else
              if ( (1 - icefrac(i,j))  <   minleadsfrac ) then
! 1.3 If test is true, use climatology for that element
!    else manipulate NWP field
                out_field(i,j) = fieldClim(i,j)
              else
                out_field(i,j) =                                        &
     &              fieldNWP(i,j) / ( 1 - icefrac(i,j))
              endif
            endif
          enddo      ! i
        enddo        ! j
        else
! 2.1 Loop over each element in field
!     and check if icefrac element is missing data
        do j = 1,nrows           ! Loop over rows
          do i = 1,ncols         ! Loop over columns
            if ( icefrac (i,j)  ==  rmdi ) then
              out_field (i,j) = fieldClim(i,j)
            else
              if ( icefrac(i,j)  <   minicefrac) then
! 2.3 If test is true, use climatology for that element
!     else manipulate NWP field
                out_field(i,j) = fieldClim(i,j)
              else
                out_field(i,j) = fieldNWP(i,j) / icefrac(i,j)
              endif
            endif
          enddo      ! i
        enddo
      endif
      return
      END SUBROUTINE interleave
!----------------------------------------------------------------------
