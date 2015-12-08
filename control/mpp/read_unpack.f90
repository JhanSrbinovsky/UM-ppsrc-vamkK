! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP
SUBROUTINE read_unpack(buf,isize,unpack_size,lookup,fixhd12,      &




 icode,cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE PrintStatus_mod
USE lookup_addresses

IMPLICIT NONE

! Description: Unpacking codes from READFL1A & RDMULT1A is
!              combined into this new deck.  This enables
!              readin data to be unpacked for both UM and SX


! Subroutine Arguments:

INTEGER :: pack_code
INTEGER, INTENT(IN) :: fixhd12      !   IN : 12th element of fixed length header
INTEGER, INTENT(IN) :: isize        !   IN :
INTEGER, INTENT(IN) :: unpack_size  !   IN : Size of data when unpacked

INTEGER, INTENT(OUT) :: icode       !   OUT: return code






REAL, INTENT(INOUT) ::  buf(unpack_size) ! INOUT: Holds field to be unpacked
INTEGER, INTENT(INOUT) :: lookup(64)     ! INOUT: LOOKUP header from dump

CHARACTER(LEN=80) ::    cmessage

!local variables
INTEGER ::  i 
INTEGER ::  num_unpack_values   ! used to detect error
INTEGER ::  len_full_word 
INTEGER ::  process_size






REAL ::  tempbuf(unpack_size) ! Holds field while it's unpacked from buf
REAL ::  probuf(isize),extbuf(isize)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------


! ----------------------------------------------------------------

IF (lhook) CALL dr_hook('READ_UNPACK',zhook_in,zhook_handle)
len_full_word=64
pack_code=MOD((lookup(lbpack)),10)

IF (pack_code == 2) THEN ! 32 Bit CRAY packing
  IF (lookup(data_type)  ==  1) THEN ! if it's REAL
! DEPENDS ON: expand32b
    CALL expand32b( lookup(lblrec) , buf, fixhd12 )
  ELSE ! not a real field, so not allowed to be compressed
    WRITE(6,*) 'READ_UNPACK : Error trying to uncompress field'
    WRITE(6,*) 'This field cannot be compressed/uncompressed ',   &
               'at it is type ',lookup(data_type)
    icode=200
    cmessage='Failure uncompressing field'
    GO TO 9999
  END IF  ! if it's REAL


ELSE IF (pack_code /= 0) THEN
  icode=6
  cmessage='READ_UNPACK: packing type not supported'
  IF (lhook) CALL dr_hook('READ_UNPACK',zhook_out,zhook_handle)
  RETURN
END IF




 9999 CONTINUE

IF (lhook) CALL dr_hook('READ_UNPACK',zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_unpack
