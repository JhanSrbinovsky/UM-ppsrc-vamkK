! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Purpose: To set the LBEGIN and LBNREC fields in the LOOKUP Headers
!             for VN 16 Type Dumpfiles - addressed by location and
!             length which are rounded up the 'io_field_padding'. The
!             data start is also rounded up to the 'io_data_alignment'
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Dump I/O

!    Interface and arguments: ------------------------------------------
SUBROUTINE set_dumpfile_address(fixhd, len_fixhd,                 &
                                lookup, len1_lookup, len2_lookup, &
                                number_of_data_words_in_memory,   &
                                number_of_data_words_on_disk,     &
                                disk_address)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE io_configuration_mod, ONLY : &

    io_field_padding,        &
    io_data_alignment

USE PrintStatus_mod
USE UM_ParVars
USE lookup_addresses

IMPLICIT NONE

INTEGER                                                           &
 len_fixhd                                                        &
                                 ! IN  Length of fixed length
                                 !     header
,len1_lookup                                                      &
                                 ! IN  1st dim of lookup
,len2_lookup                                                      &
                                 ! IN  2nd dim of lookup
,number_of_data_words_in_memory                                   &
                                 ! OUT Number of Data Words
                                 !     in memory
,number_of_data_words_on_disk                                     &
                                 ! OUT Number of data words
                                 !     on disk
,disk_address                    ! OUT Current rounded disk
                                 !     address and final data
                                 !     length

INTEGER                                                           &
 fixhd(len_fixhd)                                                 &
                                 !IN Fixed length header
,lookup(len1_lookup,len2_lookup) !IN/OUT PP lookup tables


INTEGER                                                           &
 disk_length                                                      &
                                 ! current data length on disk
,i                                                                &
                                 ! Loop Index
,old_fixhd_160                   ! Original value of fixhd(160)
                                 ! checking as the new addresses
                                 ! are computed

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Including cprintst so that we can use the PrintStatus variable

IF (lhook) CALL dr_hook('SET_DUMPFILE_ADDRESS',zhook_in,zhook_handle)

IF (fixhd(160) <  0) THEN 
  IF (lhook) CALL dr_hook('SET_DUMPFILE_ADDRESS',zhook_out,zhook_handle)
  RETURN
END IF

!--check that the initial data address has been rounded up
!  to a sector boundary - REMEMBER all the code removes
!  one from this address because addresses start at zero.
IF ((fixhd(160)-1)  /=                                            &
  (((fixhd(160)-1)+io_data_alignment-1)/io_data_alignment)*       &
     io_data_alignment) THEN
!--save the current initial data address
  old_fixhd_160=fixhd(160)
!--round up the initial disk address
  fixhd(160)=(((fixhd(160)-1)+io_data_alignment-1)/               &
   io_data_alignment)*io_data_alignment+1
  IF (mype == 0) THEN
    IF (printstatus >= prstatus_diag) THEN
      WRITE(6,'(/,A,i10,A,i10)')                                  &
 'SET_DUMPFILE_ADDRESS: Start of Data Address on disk reset from '&
      ,old_fixhd_160-1,' to ', fixhd(160)-1
    END IF
  END IF
END IF

!--adjust the Dumpfile version Number
!      if(fixhd(1) <  16) fixhd(1)=16

!--count the number of words on disk and in memory
number_of_data_words_on_disk=0
number_of_data_words_in_memory=0

!--find the initial data location on disk
disk_address=fixhd(160)-1

!--loop over all the entries and alter the addresses and lengths
DO i=1, len2_lookup
!--check for a PP type file with an incomplete lookup table
  IF (lookup(1, i) == -99) EXIT
!--check for packing to 32-bits
    IF (lookup(lbpack,i)-                                         &
      ((lookup(lbpack,i)/10)*10) == 2) THEN
      disk_length=(lookup(lblrec,i)+1)/2
    ELSE
      disk_length=lookup(lblrec,i)
    END IF
!--count the number of words
    number_of_data_words_on_disk=                                 &
    number_of_data_words_on_disk+disk_length
    number_of_data_words_in_memory=                               &
    number_of_data_words_in_memory+lookup(lblrec,i)
!--round up the length to a number of sectors
    disk_length=((disk_length+io_field_padding-1)/                &
    io_field_padding)*io_field_padding
!--set the disk address
    lookup(lbegin,i)=disk_address
!--set the disk length
    lookup(lbnrec,i)=disk_length
!--increment the disk address
    disk_address=disk_address+lookup(lbnrec,i)
END DO

IF (mype == 0) THEN
  IF (printstatus >= prstatus_diag) THEN
!--find the number of bytes in a word
    CALL word_length(i)
!--print the diagnostic message
    WRITE(6,'(/''SET_DUMPFILE_ADDRESS: Dumpfile LOOKUP Address'', '      // &
        ' '' and Lengths Rewritten:''//                                 '// &
        ' i10,'' Words Stored as Data Length in FIXHD(161)''/           '// &
        ' i10,'' Words Used in Memory for Data''/                       '// &
        ' i10,'' Words Used on Disk for Data''/                         '// &
        ' i10,'' Words Used on Disk for Data after Rounding''/          '// &
        ' i10,'' Words Used on Disk in Total for the File'',''  ('',i11,'// &
        ' '' Bytes)''/)') fixhd(161), number_of_data_words_in_memory,       &
        number_of_data_words_on_disk, disk_address-fixhd(160),          &
        disk_address, disk_address*i
  END IF
END IF
IF (lhook) CALL dr_hook('SET_DUMPFILE_ADDRESS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE set_dumpfile_address
