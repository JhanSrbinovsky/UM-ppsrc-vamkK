! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: interp_time
!
! Purpose: Flux processing routine.
!          Sets fieldout = field1 * weight1 + field2 * weight2
!          and changes date in lookup table to that of
!          validity time input
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine interp_time(Int_Head, ncols, nrows, rmdi,              &
!----------------------------------------------------------------------
! comdeck: AVALTIM
! Purpose: argument list for variables storing a time of validity.
!          This deck is linked to CVALTIM.
!----------------------------------------------------------------------
     & ValidYear, ValidMonth, ValidDay, ValidHour, ValidMin, ValidSec,  &
!----------------------------------------------------------------------
     &         weight1, weight2, Field1, Field2, FieldOut)

      USE plookups_mod, ONLY: len_fixhd, len1_lookup, len_inthd,        &
                              len_realhd, itgrid, iugrid,               &
                              max_num_fc_times, max_num_clim_times,     &
                              max_num_in_flux, len2_lookuppreferred,    &
                              len2_lookupprevious, len2_lookupclimate
      USE lookup_addresses

      IMPLICIT NONE

! declaration of parameters

! declaration of argument list
      integer Int_Head(Len_IntHd) ! IN/OUT   date is changed
      integer ncols               ! IN  number of columns
      integer nrows               ! IN  number of rows
      real    rmdi                ! IN  real missing data indicator
! validity time to insert in Lookup table: intent IN
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
      real weight1   ! IN weight to give to 1st climate field
      real weight2   ! IN weight to give to 2nd climate field
      real Field1(ncols,nrows)   ! IN  1st field
      real Field2(ncols,nrows)   ! IN  2nd field
      real FieldOut(ncols,nrows) ! OUT interpolated field

! declaration of local scalars
      integer jrow  ! row number
      integer icol  ! column number
!----------------------------------------------------------------------

! 1. do time interpolation

      do jrow = 1, nrows
        do icol = 1, ncols
          if ( Field1 (icol, jrow)  /=  rmdi .and.                      &
     &         Field2 (icol, jrow)  /=  rmdi       ) then

            FieldOut (icol, jrow) = weight1 * Field1 (icol, jrow)       &
     &                            + weight2 * Field2 (icol, jrow)

          else
            FieldOut (icol, jrow) = rmdi

          end if
        end do    ! icol
      end do      ! jrow

! 2. set validity time in integer lookup table

      Int_Head(LBYR)  = ValidYear
      Int_Head(LBMON) = ValidMonth
      Int_Head(LBDAT) = ValidDay
      Int_Head(LBHR)  = ValidHour
      Int_Head(LBMIN) = ValidMin

      return
      END SUBROUTINE interp_time
!----------------------------------------------------------------------
