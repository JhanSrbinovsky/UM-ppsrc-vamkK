! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE IOERROR----------------------------------------
!
!    Purpose: Prints out a message after using buffer in/out when
!             either a return code < 0.0 is encountered
!             by UNIT function or value returned by LENGTH
!             differs from length of I/O request.
!
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Dump I/O

SUBROUTINE ioerror(string,error,len_io,len_io_req)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER  :: len_io      ! Number of 64-bit words transferred as
                        !  registered  by LENGTH function
INTEGER  :: len_io_req  ! Number of 64-bit words requested for
                        !  transfer via BUFFER IN/OUT
CHARACTER(len=*) :: string

REAL     :: error       ! Error code returned by UNIT function

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! -------------------------------------------------------------

IF (lhook) CALL dr_hook('IOERROR',zhook_in,zhook_handle)
WRITE(6,'('' **FATAL ERROR WHEN READING/WRITING MODEL DUMP**'')')
WRITE(6,'('' '',A)') TRIM(string)
WRITE(6,'('' Error code = '',F6.2)') error
WRITE(6,'('' Length requested            = '',I9)') len_io_req
WRITE(6,'('' Length actually transferred = '',I9)') len_io
WRITE(6,'('' Fatal error codes are as follows:'')')
WRITE(6,'(A)')                                                    &
' -1.0 Mismatch between actual and requested data length'
WRITE(6,'(''  0.0 End-of-file was read'')')
WRITE(6,'(''  1.0 Error occurred during read'')')
WRITE(6,'(''  2.0 Other disk malfunction'')')
WRITE(6,'(''  3.0 File does not exist'')')
WRITE(6,'('' ***********************************************'')')

WRITE(0,'(//)')
WRITE(0,'('' **FATAL ERROR WHEN READING/WRITING MODEL DUMP**'')')
WRITE(0,'('' '',A)') TRIM(string)
WRITE(0,'('' Error code = '',F6.2)') error
WRITE(0,'('' Length requested            = '',I9)') len_io_req
WRITE(0,'('' Length actually transferred = '',I9)') len_io
WRITE(0,'('' Fatal Error codes are as follows:'')')
WRITE(0,'(A)')                                                    &
' -1.0 Mismatch between actual and requested data length'
WRITE(0,'(''  0.0 End-of-file was read'')')
WRITE(0,'(''  1.0 Error occurred during read'')')
WRITE(0,'(''  2.0 Other disk malfunction'')')
WRITE(0,'(''  3.0 File does not exist'')')
WRITE(0,'('' ***********************************************'')')

IF (lhook) CALL dr_hook('IOERROR',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ioerror
