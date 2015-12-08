! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Reads in a number of fields from UM format file

! Subroutine Interface:

SUBROUTINE readflds_serial (nftin, number_of_fields,              &
                            first_field, lookup, len1_lookup_arg, &
                            d1, fixhd,                            &
                            expand, icode, cmessage)

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE IO
  USE ereport_mod, ONLY : ereport
  USE lookup_addresses

  IMPLICIT NONE

! Description:

!  Reads in NUMBER_OF_FIELDS fields from file on unit NFTIN,
!  starting at field number FIRST_FIELD. The data is returned
!  in the D1 array.

!  ==============================================================
!  ========    WARNING    WARNING   WARNING    ==================
!  ==============================================================
!  This routine is set up to read any number of fields but it only
!  works if number of fields is set to 1. Needs investigating why
!  it fails with more than one level.
!  ===============================================================
!  ===============================================================

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dump I/O

! Subroutine Arguments:

  INTEGER :: nftin               ! IN: unit number to read data from
  INTEGER :: number_of_fields    ! IN: number of fields to read in
  INTEGER :: first_field         ! IN: first field to read in
  INTEGER :: len1_lookup_arg     ! IN: first dimension of LOOKUP table
  INTEGER :: lookup(len1_lookup_arg,*)  ! IN: lookup table starting
                                 !     at field 1
  INTEGER :: fixhd(*)            ! IN: fixed length header

  INTEGER :: expand        ! IN: (=1 if WGDOS or RLE packed data
                           !      is to be expanded)
  INTEGER :: icode         ! OUT: return code

  REAL :: d1(*)            ! OUT: array to return the data in

  CHARACTER(LEN=80) :: cmessage ! OUT: Error message if ICODE <> 0

! COMMON blocks and PARAMETERs

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

! Local variables

  INTEGER :: k                ! loop over fields to read in
  INTEGER :: pack_code        ! packing code for field
  INTEGER :: field_start      ! location of field on disk
  INTEGER :: data_size        ! number of words of data on disk
                              ! (including padding for WFIO)
  INTEGER :: data_read_size   ! number of words to read from disk
  INTEGER :: data_full_size   ! number of words after any unpacking
  INTEGER :: len_io           ! number of words read from disk
  INTEGER :: field_item       ! Item number of field
  INTEGER :: field_sect       ! Section number of field
  INTEGER :: field_model      ! Model ID of field
  INTEGER :: i                ! loop index

  REAL    :: a_io             ! Return code from BUFFIN

  CHARACTER(LEN=*), PARAMETER :: routinename = 'READFLDS_SERIAL'

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------

  IF (lhook) CALL dr_hook('READFLDS_SERIAL',zhook_in,zhook_handle)
IF (fixhd(12) < 403) THEN
  WRITE (6,*) 'READFLDS_SERIAL: file created by UM version ', fixhd(12)
  icode = 10
  cmessage = 'READFLDS_SERIAL: Cannot read files from before vn4.3'

  CALL ereport ( routinename, icode, cmessage )
END IF

DO k = first_field, first_field+number_of_fields-1

  field_item  = MOD(lookup(42,k),1000)
  field_sect  = (lookup(42,k)-field_item)/1000
  field_model = lookup(45,k)
  pack_code   = MOD((lookup(lbpack,k)),10)

!   ---------------------------------------
!   Determine location of the field on disk
!   ---------------------------------------

  field_start = lookup(lbegin,k) ! position of field in file
  IF (field_start <= 0) THEN
    WRITE (6,*) 'READFLDS_SERIAL: start address =',field_start
    icode = 20
    cmessage = 'READFLDS_SERIAL: start address of field not given'
    CALL ereport ( routinename, icode, cmessage )
  END IF

!   --------------------
!   Determine data sizes
!   --------------------

! DATA_SIZE : contains the number of words of data used to store the
! field on disk (needs to be halved if 32 bit packing has been used)

  IF (pack_code == 2) THEN
    data_size = (lookup(lblrec,k)+1)/2
  ELSE
    data_size = lookup(lblrec,k)
  END IF

! DATA_FULL_SIZE : is the number of words required to store the field
! in memory after any unpacking is done.

! This is to give buf the correct size in RDUNPCK, as
! buf will be the final expanded size of whole field
! including extra data
! WARNING LBEXT - may be -32768 MISSING VALUE !

  IF ((pack_code == 4) .AND. (lookup(lbext, k) > 0)) THEN
    data_full_size = MAX(lookup(lbrow, k)*lookup(lbnpt, k)  &
         +lookup(lbext, k) ,lookup(lblrec,k))
  ELSE
    data_full_size = MAX(lookup(lbrow, k)*lookup(lbnpt, k)  &
         ,lookup(lblrec,k))
  END IF

  IF ( (lookup(lbrow,k) < 0) .OR. (lookup(lbnpt,k) < 0) ) THEN
    data_full_size = lookup(lblrec,k)
  END IF

! DATA_READ_SIZE : contains the number of words to data that need to
! be read in for a field. Each field has extra words of dummy data
! added at the end to ensure each field starts on a disk sector
! boundary. The last field on a dump does not have these extra words
! added

  IF (k /= (first_field+number_of_fields-1)) THEN
    data_read_size = lookup(lbnrec,k)
  ELSE
    data_read_size = data_size
  END IF

  IF (data_read_size < 0) THEN
    WRITE (6,*) 'READFLDS_SERIAL: number of words to read =',   &
           data_read_size
    icode = 30
    cmessage = 'READFLDS_SERIAL: number of words to read not given'
    CALL ereport ( routinename, icode, cmessage )
  END IF

! data_full_size needs to be at least as big as data_read_size since
! it is used to dimension the BUF array (in READ_MULTI ?)

  data_full_size = MAX(data_full_size, data_read_size)

!   -------------------------------------------
!   Move file pointer to the start of the field
!   -------------------------------------------


  CALL setpos (nftin, field_start, icode)
  IF (icode /= 0) THEN
    WRITE (6,*) 'READFLDS_SERIAL - SETPOS failed to move to ',  &
                 field_start,' on unit ',nftin
    WRITE (6,*) 'SETPOS returned error code ',icode
    icode = 40
    cmessage = 'SETPOS failed in READFLDS_SERIAL. See output.'
    CALL ereport ( routinename, icode, cmessage )
  END IF

!   ------------------------
!   Read the field from disk
!   ------------------------

  IF (pack_code == 2) THEN

!     Data is packed using CRAY 32 bit method - note
!     that we need to read in 2*data_read_size 32 bit
!     words using BUFFIN32

! DEPENDS ON : buffin32_f77
    CALL buffin32_f77 (nftin,d1,2*data_read_size,len_io,a_io)

!     And then halve LEN_IO to satisfy tests against
!     data_read_size

    len_io = len_io/2

  ELSE

!     Data is not 32-bit packed

    CALL buffin (nftin,d1(1:data_read_size),data_read_size,len_io,a_io)

  END IF

!   --------------------------------
!   Check that data been read in OK?
!   --------------------------------

  IF ((a_io /= -1.0) .OR. (len_io /= data_read_size)) THEN
    WRITE (6,*) 'READFLDS_SERIAL : Error in call to BUFFIN'
    WRITE (6,*) 'LEN_IO : ',len_io
    WRITE (6,*) 'A_IO : ',a_io
    WRITE (6,*) 'Attempting to read field (Model,Section,Item) ',  &
                 field_model, field_sect, field_item
    WRITE (6,'(10(E10.5,1X))') d1(1:100)
    icode = 50
    cmessage = 'READFLDS_SERIAL: Failure reading in field'
    CALL ereport ( routinename, icode, cmessage )
  END IF

!   ----------------------------------------------
!   If it's a compressed REAL field, expand it out
!   ----------------------------------------------

! DEPENDS ON: READ_UNPACK
  CALL read_unpack (d1, data_read_size, data_full_size,  &
                    lookup(1,k), fixhd(12),              &
                    expand,                              &
                    icode, cmessage)

  IF (icode /= 0) THEN
    WRITE(6,*)'READFLDS_SERIAL: Failure unpacking field'
    CALL ereport ( routinename, icode, cmessage )
  END IF

END DO ! K : loop over fields

IF (lhook) CALL dr_hook('READFLDS_SERIAL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE readflds_serial
