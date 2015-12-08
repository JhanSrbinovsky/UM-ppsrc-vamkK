! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE READ_FLH --------------------------------------
!
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Purpose:
!             Reads in the fixed length header from file attached to
!             unit NFTIN.
!
!    Documentation:
!             Unified Model Documentation Paper No F3
!
!             Code Owner: See Unified Model Code Owners HTML page
!             This file belongs in section: Dump I/O

SUBROUTINE read_flh (nftin,fixhd,len_fixhd,                       &
           icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IO
IMPLICIT NONE

INTEGER                                                           &
 nftin                                                            &
                  ! IN  Unit no of dump
,len_fixhd                                                        &
                  ! IN  Length of fixed length header
,fixhd(len_fixhd) ! OUT Fixed length header

INTEGER  icode    ! OUT Return code; successful=0, error > 0
CHARACTER(LEN=80)                                                    &
 cmessage         ! OUT Error message if ICODE > 0

! Local variables:---------------------------------------------
INTEGER len_io
REAL a

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! -------------------------------------------------------------

IF (lhook) CALL dr_hook('READ_FLH',zhook_in,zhook_handle)
icode=0
cmessage=' '

!  1. Buffer in fixed length header record


CALL buffin (nftin,fixhd,len_fixhd,len_io,a)

!  2. Check for I/O errors
IF (a /= -1.0 .OR. len_io /= len_fixhd) THEN
! DEPENDS ON: ioerror
  CALL ioerror('buffer in of fixed length header',a,len_io        &
               ,len_fixhd)
  cmessage='READ_FLH: I/O error'
  icode=1
END IF

IF (lhook) CALL dr_hook('READ_FLH',zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_flh
