
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine Interface
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dump I/O

!    Purpose: Reads in model obs dump on unit NFTIN and checks model
!             and dump dimensions for consistency.
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Documentation: Unified Model Documentation Paper No F3

SUBROUTINE readacobs(nftin,fixhd,len_fixhd                        &
 ,inthd,len_inthd                                                 &
 ,realhd,len_realhd                                               &
 ,levdepc,len1_levdepc,len2_levdepc                               &
 ,rowdepc,len1_rowdepc,len2_rowdepc                               &
 ,coldepc,len1_coldepc,len2_coldepc                               &
 ,flddepc,len1_flddepc,len2_flddepc                               &
 ,extcnst,len_extcnst                                             &
 ,dumphist,len_dumphist                                           &
 ,cfi1,len_cfi1                                                   &
 ,cfi2,len_cfi2                                                   &
 ,cfi3,len_cfi3                                                   &
 ,lookup,len1_lookup,len2_lookup                                  &
 ,len_data,d1,                                                    &
  icode,cmessage                                                  &
 ,ipt                                                             &
                )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IO
USE UM_ParVars
USE lookup_addresses
USE Submodel_Mod

IMPLICIT NONE

INTEGER ipt
INTEGER                                                           &
 nftin                                                            &
               !IN Unit no of dump
,len_fixhd                                                        &
               !IN Length of fixed length header
,len_inthd                                                        &
               !IN Length of integer header
,len_realhd                                                       &
               !IN Length of real header
,len1_levdepc                                                     &
               !IN 1st dim of level dep consts
,len2_levdepc                                                     &
               !IN 2nd dim of level dep consts
,len1_rowdepc                                                     &
               !IN 1st dim of row dep consts
,len2_rowdepc                                                     &
               !IN 2nd dim of row dep consts
,len1_coldepc                                                     &
               !IN 1st dim of column dep consts
,len2_coldepc                                                     &
               !IN 2nd dim of column dep consts
,len1_flddepc                                                     &
               !IN 1st dim of field dep consts
,len2_flddepc                                                     &
               !IN 2nd dim of field dep consts
,len_extcnst                                                      &
               !IN Length of extra constants
,len_dumphist                                                     &
               !IN Length of history block
,len_cfi1                                                         &
               !IN Length of comp field index 1
,len_cfi2                                                         &
               !IN Length of comp field index 2
,len_cfi3                                                         &
               !IN Length of comp field index 3
,len1_lookup                                                      &
               !IN 1st dim of lookup
,len2_lookup   !IN 2nd dim of lookup

INTEGER                                                           &
 len_data                                                         &
                !IN Length of model data
,icode          !OUT Return code; successful=0
                !                 error > 0

CHARACTER(LEN=80)                                                    &
 cmessage       !OUT Error message if ICODE > 0

INTEGER                                                           &
 fixhd(len_fixhd)                                                 &
                  !IN Fixed length header
,inthd(len_inthd)                                                 &
                  !IN Integer header
,lookup(len1_lookup,len2_lookup)                                  &
                  !IN PP lookup tables

,cfi1(len_cfi1+1)                                                 &
                  !IN Compressed field index no 1
,cfi2(len_cfi2+1)                                                 &
                  !IN Compressed field index no 2
,cfi3(len_cfi3+1) !IN Compressed field index no 3

REAL                                                              &
 realhd(len_realhd)                                               &
                                      !IN Real header
,levdepc(1+len1_levdepc*len2_levdepc)                             &
                                      !IN Lev dep consts
,rowdepc(1+len1_rowdepc*len2_rowdepc)                             &
                                      !IN Row dep consts
,coldepc(1+len1_coldepc*len2_coldepc)                             &
                                      !IN Col dep consts
,flddepc(1+len1_flddepc*len2_flddepc)                             &
                                      !IN Field dep consts
,extcnst(len_extcnst+1)                                           &
                                      !IN Extra constants
,dumphist(len_dumphist+1)                                         &
                                      !IN History block
,d1(len_data)     !IN Real equivalence of data block

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

! Local variables:---------------------------------------------
INTEGER start_block                                               &
                     ! Pointer to current position in file
,len_io                                                           &
                     ! No of 64-bit words buffered in
,k,i                                                              &
                     ! Loop counts
,ipts                ! No of 64-bit words requested to be
                     ! buffered in
REAL a               ! Error code returned by UNIT

INTEGER real_start_block                                          &
                         ! Real disk address
 , l                                                              &
                         ! loop counter
 , word_address                                                   &
                         ! word address on disk of the record
 , um_sector_ipts                                                 &
                         ! number fo words to read, rounded up
                         ! to a sector size
 , l_ipts                ! local value of ipts for address calc.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER :: ipts_to_read
INTEGER :: word_address_after_record
!--------------------------------------------------------------

IF (lhook) CALL dr_hook('READACOBS',zhook_in,zhook_handle)

IF (mype  ==  0) THEN
  WRITE(6,'(/,'' READING ACOBS FILE ON UNIT'',I3)')nftin
  WRITE(6,'('' #####################################'',/)')
END IF

icode=0
cmessage=' '

!  1. Read in all header records and check for consistency.
!     START_BLOCK points to position of model data block
!     on return

! DEPENDS ON: readhead
CALL readhead(nftin,fixhd,len_fixhd,                              &
              inthd,len_inthd,                                    &
              realhd,len_realhd,                                  &
              levdepc,len1_levdepc,len2_levdepc,                  &
              rowdepc,len1_rowdepc,len2_rowdepc,                  &
              coldepc,len1_coldepc,len2_coldepc,                  &
              flddepc,len1_flddepc,len2_flddepc,                  &
              extcnst,len_extcnst,                                &
              dumphist,len_dumphist,                              &
              cfi1,len_cfi1,                                      &
              cfi2,len_cfi2,                                      &
              cfi3,len_cfi3,                                      &
              lookup,len1_lookup,len2_lookup,                     &
              len_data,                                           &
              start_block,icode,cmessage)

IF (icode >  0) GO TO 9999


!  2. Buffer in model data one field at a time for
!     conversion from 32-bit to 64-bit numbers

IF (fixhd(160) >  0) THEN

! Check for error in file pointers
  real_start_block=start_block
  IF (start_block /= fixhd(160)) THEN
! If new format Dumpfile, we must reset the start address
    IF ((lookup(lbnrec,1) == 0 .AND. lookup(lblrec,1) >  0) .OR.  &
! Ocean ACOBS Files (?)
     ((lookup(lbnrec,1) == imdi) .OR. (lookup(lbegin,1) == imdi)) &
     .OR.                                                         &
! Prog lookups in dump before vn3.2:
     ((lookup(lbnrec,1) == imdi) .AND. (fixhd(12) <= 301))) THEN
      cmessage='READACOBS: Addressing conflict'
      icode=1
! DEPENDS ON: poserror
      CALL poserror('model data',                                 &
                    start_block,160,fixhd(160))
      GO TO 9999
    ELSE
      real_start_block=fixhd(160)
    END IF
  END IF

!      Move to start of data.

  CALL setpos (nftin,fixhd(160)-1,icode)

! Loop over number of fields in data block
  DO k=1,fixhd(152)

    IF (lookup(lblrec,k) >  0) THEN   !  Any data for this field?

! Test whether data stored as 32-bit on disk
      IF (MOD((lookup(lbpack,k)),10) == 2) THEN
        ipts=(lookup(lblrec,k)+1)/2
      ELSE
        ipts=lookup(lblrec,k)
      END IF

! Compute word address in file from which to begin I/O

! Old Format dumpfiles
      IF ((lookup(lbnrec,k) == 0) .OR.                            &
! Ocean ACOBS Files (?)
    ((lookup(lbnrec,k) == imdi) .OR. (lookup(lbegin,k) == imdi))  &
    .OR.                                                          &
! Prog lookups in dump before vn3.2:
    ((lookup(lbnrec,k) == imdi) .AND. (fixhd(12) <= 301))) THEN
! Dump and ancillary files
        word_address=1
        IF (k >  1) THEN
          DO l=2,k
            IF (MOD(lookup(lbpack,l-1),10) == 2) THEN
              l_ipts=(lookup(lblrec,l-1)+1)/2
            ELSE
              l_ipts=(lookup(lblrec,l-1))
            END IF
            word_address=word_address+l_ipts
          END DO
        END IF
        word_address=fixhd(160)+word_address-2
        um_sector_ipts=ipts
        ipts_to_read=um_sector_ipts
      ELSE

! PP type files and new format Dumpfiles (vn4.4 onwards)
        word_address=lookup(lbegin,k)
! Use the stored round-up value
        um_sector_ipts=lookup(lbnrec,k)
! We only want to read the data length, excluding any sector padding
        ipts_to_read=lookup(lblrec,k)
      END IF

! Position file pointer

      CALL setpos(nftin,word_address,icode)

! Read data into final position
! Check that data_type is valid no: 1 to 3 or -1 to -3
  IF ((lookup(data_type,k) >= 1 .AND. lookup(data_type,k) <= 3)   &
  .OR.(lookup(data_type,k) <= -1 .AND. lookup(data_type,k) >= -3))&
      THEN
        ipts=um_sector_ipts
        word_address_after_record=um_sector_ipts+word_address
        CALL buffin(nftin,d1(lookup(naddr,k):),ipts_to_read,len_io,a)
        IF ((a /= -1.0) .OR. (len_io /= ipts_to_read)) THEN
          WRITE(6,*)'ERROR READING DUMP ON UNIT ',nftin
          icode=2
          cmessage='READACOBS: BAD BUFFIN OF DATA'
! DEPENDS ON: ioerror
          CALL ioerror('BUFFER IN FROM READACOBS',a,len_io,ipts)
          GO TO 9999
        END IF
! Move the file pointer to the file address after any sector 
! size based padding that was in the record
        CALL setpos(nftin,word_address_after_record,icode)
! Error in lookup(data_type,k)
      ELSE
        icode=3
        WRITE(cmessage,*)'READACOBS:  Invalid code (',lookup(data_type,k),&
            ') in LOOKUP(DATA_TYPE,K)'
      END IF

    END IF  !  Skip to here if no data for this field

    start_block=start_block+lookup(lblrec,k)
    real_start_block=real_start_block+um_sector_ipts

  END DO

  IF (mype  ==  0) THEN
    WRITE(6,'('' '')')
    IF (fixhd(5) >= 6 .AND. fixhd(5) <= 8) THEN ! AC/Var Obs/ Cx file
      WRITE(6,'('' OBSERVATION DATA'')')
    ELSE
      WRITE(6,'('' MODEL DATA'')')
    END IF
    WRITE(6,'('' '',i8,'' words long'')')fixhd(161)
  END IF ! mype  ==  0

END IF

IF (mype  ==  0) THEN
  WRITE(6,'('' '')')
  WRITE(6,'(A,i9,A,I3)')' INITIAL DATA SUCCESSFULLY READ -',      &
  start_block,' WORDS FROM UNIT',nftin
  IF (real_start_block /= start_block) THEN
    WRITE(6,'(/'' Number of Words Read from Disk was '',i9)')     &
          real_start_block
  END IF
END IF ! mype  ==  0

9999 CONTINUE

IF (lhook) CALL dr_hook('READACOBS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE readacobs
