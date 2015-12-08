! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: amend_lookup_flux
!                    amend_lookup_grid
!                    amend_lookup
!
! Purpose: Flux processing Routine.
!          Amends grid information in Int_Head and Real_Head
!          Also amends packing and level information
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine amend_lookup_flux ( Int_Head, Real_Head,               &
     &                          StCode, FFCode, PPCode, IVTOffHr )

! Purpose:  Amends field code, packing and level information in
!           the real and integer headers for a field

      USE plookups_mod, ONLY: len_fixhd, len1_lookup, len_inthd,        &
                             len_realhd, itgrid, iugrid,                &
                             max_num_fc_times, max_num_clim_times,      &
                             max_num_in_flux, len2_lookuppreferred,     &
                             len2_lookupprevious, len2_lookupclimate
      USE lookup_addresses

      IMPLICIT NONE

! declaration of parameters

! declaration of argument list
      integer Int_Head(Len_IntHd) ! IN/OUT int part of lookup table
      real Real_Head(Len_RealHd)  ! IN/OUT real part of lookup table

! field codes etc to insert in integer header that is output
      integer StCode   ! IN stash code
      integer FFCode   ! IN Met O 8 field code
      integer PPCode   ! IN PP package code
      integer IVTOffHr ! IN offset of validity time from reference

! no other variables  used

! no externals

!----------------------------------------------------------------------

! 1. Set the field codes
      Int_Head(ITEM_CODE) = StCode
      Int_Head(LBFC)      = PPCode
      Int_Head(LBTYP)     = FFCode

! 2. Set packing code and levels codes
      Int_Head(LBPACK) = 0
      Int_Head(LBLEV) = 8888
      Real_Head(BLEV - Len_IntHd) = 0.0

! 3. Set the forcast time
      Int_Head(LBFT) = IVTOffHr

      return
      END SUBROUTINE amend_lookup_flux
!----------------------------------------------------------------------
