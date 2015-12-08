! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: check_header
!
! Purpose: Flux processing routine.
!          Checks to see if NWP file exists
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine check_header (StCode,Len1_Lookup,                      &
     &                         Len2_Lookup,Lookup,                      &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &                         l_present_output)


      USE lookup_addresses

      IMPLICIT NONE

! declaration of argument list
      integer StCode               ! IN Field code in lookup table
      integer Len1_Lookup,Len2_Lookup     ! IN true lengths of tables
      integer Lookup(Len1_Lookup, Len2_Lookup)  ! IN lookup tables
      logical l_present_output     ! OUT logical to check if item exists

! declaration of global parameters used
!----------------------------------------------------------------------
! comdeck: CVALTIM
! Purpose: declares local variables storing a time of validity
!          This deck is linked to AVALTIM.
!----------------------------------------------------------------------
! declarations:
! local variables defining a time of validity
      integer ValidYear, ValidMonth, ValidDay, ValidHour,               &
     &        ValidMin, ValidSec
!----------------------------------------------------------------------

! declaration of local variables
      integer i                    ! loop counter
      logical l_climate_field      ! F => not a climate field
      logical l_data_time          ! T => search using data time
                                   ! F => search using validity time

      external time_to_use
!-----------------------------------------------------------------------

! 1. Decide whether to use 1st or second header in lookup table

      l_climate_field = .false.
! DEPENDS ON: time_to_use
      call time_to_use ( StCode, l_climate_field, l_data_time)

! 2. Set l_present_output to false
      l_present_output = .false.

! 3. Loop through NWP lookup table to find match with data time
!    If match is found, set l_present_output to true and exit

      if ( l_data_time ) then

        do i = 1,Len2_lookup

          if ( Lookup(LBYRD,i)   ==  ValidYear   .and.                  &
     &      Lookup(LBMOND,i)  ==  ValidMonth  .and.                     &
     &      Lookup(LBDATD,i)  ==  ValidDay    .and.                     &
     &      Lookup(LBHRD,i)   ==  ValidHour   .and.                     &
     &      Lookup(LBMIND,i)  ==  ValidMin    .and.                     &
     &      Lookup(ITEM_CODE,i)   ==  StCode ) then
            l_present_output = .true.
            goto 9999

          endif

        enddo

      else  ! .not. l_data_time

! 4. Loop through NWP lookup table to find match with validity time
!    If match is found, set output to true and exit

        do i = 1,Len2_lookup

          if ( Lookup(LBYR,i)   ==  ValidYear   .and.                   &
     &      Lookup(LBMON,i)  ==  ValidMonth  .and.                      &
     &      Lookup(LBDAT,i)  ==  ValidDay    .and.                      &
     &      Lookup(LBHR,i)   ==  ValidHour   .and.                      &
     &      Lookup(LBMIN,i)  ==  ValidMin    .and.                      &
     &      Lookup(ITEM_CODE,i)   ==  StCode ) then
            l_present_output = .true.
            goto 9999
          endif

        enddo

      endif  !  l_data_time

9999  continue
      return
      END SUBROUTINE check_header
!----------------------------------------------------------------------
