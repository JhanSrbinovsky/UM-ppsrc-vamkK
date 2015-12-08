! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL -----------------------------------------------------------
!LL Stash routine
!LL purpose: Generate extra data for the timeseries.
!LL This extra data provides information about what processing was done
!LL to produce the timeseries. This information will hopefully be of som
!LL use to users doing further processing of the timeseries data.
!LL This deck contains two subroutines
!LL (1) EXTRA_TS_INFO : which generates the codes and sets up the space
!LL                   : for the extra data.
!LL
!LL (2) EXTRA_MAKE_VECTOR: which computes the long/latt ht domain info
!LL                   : and puts that into the correct place in the
!LL                   : extra data
!LL Routines are  called by stmulspa1.
!LL
!LL To some extent this routine has much in common with the
!LL multi_spatial routine but as it has a different function
!LL viz generate info on timeseries rather than generating a single time
!LL for the timeseries it is coded separately.
!LL when modifying multi_spatial be sure also to modify this routine and
!LL vice versa
!LL
!LL Programming Standard: UM DOC Paper3, Verion 4 (05/02/92)
!LL
!LL System Component Covered: D711
!LL
!LL System Task:C4
!LL

!*L Interface and arguments ------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: STASH


!*L Interface and arguments: -----------------------------


!LL    Subroutine: STUFF_INT-----------------------------------------
!LL
!LL    Purpose: To put the binary representation of an integer into a
!LL    real variable through hidden equivalencing via argument passing
!LL

      SUBROUTINE STUFF_INT(array_out,data_in)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      implicit none
      real array_out    ! OUT: Real array element in calling routine
      real data_in      ! IN: Integer data in calling routine

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('STUFF_INT',zhook_in,zhook_handle)
      array_out=data_in
      IF (lhook) CALL dr_hook('STUFF_INT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE STUFF_INT
