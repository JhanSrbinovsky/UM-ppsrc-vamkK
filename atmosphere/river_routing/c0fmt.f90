! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      CHARACTER(LEN=32) FUNCTION c0fmt(indat, ncol)

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook

      IMPLICIT NONE
      
!
!     convert integer into ncol characters with 0 head.
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: River Routing
      character cdum*32, cd(0:9), cdig(32)*1
      integer indat, idat, ncol, i, idig(32)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      data cd/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/

      IF (lhook) CALL dr_hook('C0FMT',zhook_in,zhook_handle)


!
      idat = indat
!      write (6, *) 'idat, ncol=', idat, ncol
      cdum = ''
      if (ncol <= 31) then
        idig(ncol+1) = 0
        do i = ncol, 1, -1
          idig(i) = int(aint(real(idat)/(10**(i-1))))
          cdig(i) = cd(int(aint(real(idat)/(10**(i-1)))))
!          write(6, *) cdig(i)
          idat = idat - idig(i)*(10**(i-1))
        end do
      end if
!
      write (cdum,'(32a1)') (cdig(i), i= ncol, 1, -1)
!      write(6, *) "cdum = ", cdum
      c0fmt = cdum
!      write(6, *) "c0fmt = ", c0fmt

      IF (lhook) CALL dr_hook('C0FMT',zhook_out,zhook_handle)

      END FUNCTION c0fmt
