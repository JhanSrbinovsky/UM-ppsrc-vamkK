! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   SUBROUTINE COEX,COEX2,CMPS,XPND,INSTIN,EXTRIN -----------------
!LL
!LL   PURPOSE:   PACK TO AND UNPACK FROM WGDOS FORMAT
!LL
!LL   (Note that for optimal performance the routines INSTIN and
!LL    EXTRIN are inline expanded (enable option 8 on fpp)
!LL
!LL
!LL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
!LL   STANDARD B, VERSION 2, DATED 18/01/90
!LL
!LL  Logical component number: S72
!LL
!LL   SYSTEM TASK: P7
!LLEND-------------------------------------------------------------

!     This is a portable replacement for ishft, which behaves as the UM
!     expects it too when shifting more than 64 bits.
!
!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: STASH


      function ishift(num, shift)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      implicit none
        integer num, shift
        integer ishift

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle
        IF (lhook) CALL dr_hook('ISHIFT',zhook_in,zhook_handle)
        if(abs(shift) >= kind(num)*8) then
          ishift = 0
        else
          ishift = ishft(num, shift)
        end if
        IF (lhook) CALL dr_hook('ISHIFT',zhook_out,zhook_handle)
        RETURN
      end function ishift
