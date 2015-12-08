! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Reduce storage requirements by removing land points from
!          fields within ocean fieldsfiles. This is done using run
!          length encoding to compress the sequences of missing data
!          values that represent the land points. Using this method,
!          a sequence of missing data values are mapped into a single
!          missing data value followed by a value indicating the
!          length of the run of missing data values.
!
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

      subroutine runlen_encode(unpacked,unpacked_size,packed,           &
     &                         packed_size,packed_end,bmdi,             &
     &                         icode,cmessage)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      integer i,j
      integer nmdi
      integer tmdi
      integer other

      integer unpacked_size    ! IN Size of input field used to hold
                               !    unpacked data.
      integer packed_size      ! IN Size of field used to hold packed
                               !    data.
      integer packed_end       ! OUT Contains the size of the encoded
                               !     field on return.
      integer icode, ocode     ! IN/OUT Error codes >0 indicates error
                               !                    =0 success.

      real unpacked(unpacked_size) ! IN Input field representing
                               !    unpacked data.
      real packed(packed_size) ! OUT output field representing run
                               !     length decoded fields.
      real bmdi                ! IN Real missing data values

      CHARACTER(LEN=80) cmessage  ! OUT Returned error message if icode >0.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('RUNLEN_ENCODE',zhook_in,zhook_handle)
      j       = 1
      packed_end = 1
      nmdi    = 0
      tmdi    = 0
      other   = 0

!    Check if an error has already been encountered, and get out
!    if it has.
      ocode = 0
      if (icode  >   0) then
         goto 999
      else if (icode  <   0) then
         ocode = icode
         icode = 0
      end if

      do i=1,unpacked_size

        if (unpacked(i)  ==  bmdi) then
          nmdi = nmdi + 1
          tmdi = tmdi + 1
        else
          if (nmdi  >   0) then
            if ((tmdi + other)  <=  unpacked_size) then
              packed(j) = bmdi
              packed(j+1) = nmdi
              j = j + 2
              packed_end = packed_end + 2
              nmdi = 0
            else
              write (6,*)                                               &
     &        'RLENCODE : ERROR : Ocean Run length encoding failed'
              cmessage = 'Ocean Run length encoding failed'
              icode = 1
              goto 999
            end if
          end if
          packed(j) = unpacked(i)
          packed_end    = packed_end + 1
          other      = other + 1
          j          = j + 1
        end if
      end do
      if (nmdi  >   0) then
        if ((tmdi + other)  <=  unpacked_size) then
          packed(j) = bmdi
          packed(j+1) = nmdi
          j = j + 2
          packed_end = packed_end + 2
        else
          write (6,*)                                                   &
     &      'RLENCODE : ERROR : Ocean Run length encoding failed'
          cmessage = 'Ocean Run length encoding failed'
          icode = 1
          goto 999
        end if
      end if

 999  continue
!     If we have found an error, leave it in icode.  If no error
!     occurred then check if the original input value of icode was]
!     non-zero (a previous untrapped error/warning), and copy this
!     back into ICODE before eaving the routine.
      IF (icode  ==  0 .and. ocode  /=  0) then
         icode = ocode
      END IF

      packed_end = packed_end - 1

      IF (lhook) CALL dr_hook('RUNLEN_ENCODE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE runlen_encode



!----------------------------------------------------------------------
