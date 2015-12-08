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
!LL 16/3/92 Written by Simon Tett
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
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
      SUBROUTINE EXTRA_TS_INFO(extra_data,extra_data_len,no_records)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      implicit none

      integer no_records ! IN how many timeseries records are there ?
      integer extra_data_len ! IN  size of extra data required
      real extra_data(extra_data_len) ! OUT the extra data array
!*L LOCAL PARAMETERS
!LL --------------------------------------------------------------------
      integer no_extra_blocks ! how many blocks of extra data we got ?
      parameter(no_extra_blocks=6) ! 6 words to describe
!*L Subroutines called
       EXTERNAL stuff_int ! put an integer into a real
!*L Local variables
!LL -------------------------------------------------
      integer record_len ! size of block for extra data
      integer hdr(no_extra_blocks) ! the headers for each block
! order is lat, long, 2nd lat, 2nd long, first level, 2nd level
      data hdr/3,4,5,6,7,8/ ! codes for above
      integer addr ! address in array for writting/reading data
      integer i ! loop count

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!L-------------------------------------------------------------------
      IF (lhook) CALL dr_hook('EXTRA_TS_INFO',zhook_in,zhook_handle)
      record_len=no_records+1 ! how much info in a block
      addr=1
      DO i=1,no_extra_blocks ! put the headers into the extra data

! We cannot do extra_data(addr)=1000*no_records+hdr(i) as this will
! put in a floating point conversion of the integer. We actually want
! to save the binary representation of the integer which is done
! by STUFF_INT
! DEPENDS ON: stuff_int
        CALL STUFF_INT(extra_data(addr),                                &
     &    1000*no_records+hdr(i))
        addr=addr+record_len
      ENDDO
      IF (lhook) CALL dr_hook('EXTRA_TS_INFO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE EXTRA_TS_INFO


!*L Interface and arguments: -----------------------------


!LL    Subroutine: STUFF_INT-----------------------------------------
!LL
!LL    Purpose: To put the binary representation of an integer into a
!LL    real variable through hidden equivalencing via argument passing
!LL
!LL  Model              Modification history from model version 3.0:
!LL version  date
!LL   5.3    14/05/01   Changes data type of data_in to integer.
!LL                     E.Leung
!LL   5.4    10/04/02   Reverse the above change. S.D.Mullerworth

