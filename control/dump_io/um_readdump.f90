
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE UM_READDUMP---------------------------------------
!
!    Purpose: Reads in model dump on unit NFTIN and checks model
!             and dump dimensions for consistency.
!
!             Code Owner: See Unified Model Code Owners HTML page
!             This file belongs in section: Dump I/O

! Subroutine Interface

SUBROUTINE um_readdump(nftin,fixhd,len_fixhd                      &
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
 ,mpp_lookup,mpp_len1_lookup                                      &
 ,submodel_id,n_objs_d1,d1_addr                                   &
 ,len_data,d1                                                     &
 ,read_header                                                     &
  )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IO
USE ereport_mod, ONLY : ereport
USE PrintStatus_mod
USE UM_ParVars
USE Decomp_DB
USE lookup_addresses
USE Submodel_Mod

IMPLICIT NONE

 INTEGER                                                          &
  nftin                                                           &
                  !IN Unit number of dump
, len_fixhd                                                       &
                  !IN Length of fixed length header
, len_inthd                                                       &
                  !IN Length of integer header
, len_realhd                                                      &
                  !IN Length of real header
, len1_levdepc                                                    &
                  !IN 1st dim of level dep consts
, len2_levdepc                                                    &
                  !IN 2nd dim of level dep consts
, len1_rowdepc                                                    &
                  !IN 1st dim of row dep consts
, len2_rowdepc                                                    &
                  !IN 2nd dim of row dep consts
, len1_coldepc                                                    &
                  !IN 1st dim of column dep consts
, len2_coldepc                                                    &
                  !IN 2nd dim of column dep consts
, len1_flddepc                                                    &
                  !IN 1st dim of field dep consts
, len2_flddepc                                                    &
                  !IN 2nd dim of field dep consts
, len_extcnst                                                     &
                  !IN Length of extra constants
, len_dumphist                                                    &
                  !IN Length of history block
, len_cfi1                                                        &
                  !IN Length of comp field index 1
, len_cfi2                                                        &
                  !IN Length of comp field index 2
, len_cfi3                                                        &
                  !IN Length of comp field index 3
, len1_lookup                                                     &
                  !IN 1st dim of lookup
, len2_lookup                                                     &
                  !IN 2nd dim of lookup
, mpp_len1_lookup                                                 &
                  !IN 1st dim of MPP lookup
, submodel_id                                                     &
                  !IN submodel of dump
, n_objs_d1                                                       &
                  !IN number of objects (3D fields) in D1
, len_data        !IN length of model data

INTEGER                                                           &
  fixhd(len_fixhd)                                                &
                     !IN Fixed length header
, inthd(len_inthd)                                                &
                     !IN Integer header
, lookup(len1_lookup,len2_lookup)                                 &
                     !IN PP lookup tables
, cfi1(len_cfi1+1)                                                &
                     !IN Compressed field index no 1
, cfi2(len_cfi2+1)                                                &
                     !IN Compressed field index no 2
, cfi3(len_cfi3+1)                                                &
                     !IN Compressed field index no 3

, mpp_lookup(mpp_len1_lookup,len2_lookup)
                     !OUT Local processor lookup

REAL                                                              &
  realhd(len_realhd)                                              &
                     !IN Real header
, levdepc(1+len1_levdepc*len2_levdepc)                            &
                                       !IN Lev dep consts
, rowdepc(1+len1_rowdepc*len2_rowdepc)                            &
                                       !IN Row dep consts
, coldepc(1+len1_coldepc*len2_coldepc)                            &
                                       !IN Col dep consts
, flddepc(1+len1_flddepc*len2_flddepc)                            &
                                       !IN Field dep consts
, extcnst(len_extcnst+1)                                          &
                           !IN Extra constants
, dumphist(len_dumphist+1)                                        &
                           !IN History block

, d1(len_data)       !OUT Local subdomain of dump

LOGICAL                                                           &
 read_header         !IN  True if header is to be read in

! Parameters required for dimensioning the D1_ADDR array
! D1_ADDR start
      ! Information for accessing D1 addressing array
      ! Number of items of info needed for each object and maximum
      ! number of objects in D1 -

      ! Number of items of information in D1 addressing array
      INTEGER,PARAMETER:: D1_LIST_LEN=17

! Names of items in D1 addressing array. Update D1_LIST_LEN above if
! items added

      ! Prognostic, Diagnostic, Secondary or other
      INTEGER,PARAMETER:: d1_object_type    = 1 ! Internal model id
      INTEGER,PARAMETER:: d1_imodl          = 2  ! Internal model id
      INTEGER,PARAMETER:: d1_section        = 3  ! Section
      INTEGER,PARAMETER:: d1_item           = 4  ! Item
      INTEGER,PARAMETER:: d1_address        = 5  ! Address in D1
      INTEGER,PARAMETER:: d1_length         = 6  ! Record length
      INTEGER,PARAMETER:: d1_grid_type      = 7  ! Grid type
      INTEGER,PARAMETER:: d1_no_levels      = 8  ! Number of levels

      ! Stash list number for diags. -1 for progs
      INTEGER,PARAMETER:: d1_stlist_no      = 9

      ! Pointer to dump header lookup table
      INTEGER,PARAMETER:: d1_lookup_ptr     = 10

      INTEGER,PARAMETER:: d1_north_code     = 11 ! Northern row
      INTEGER,PARAMETER:: d1_south_code     = 12 ! Southern row
      INTEGER,PARAMETER:: d1_east_code      = 13 ! Eastern row
      INTEGER,PARAMETER:: d1_west_code      = 14 ! Western row
      INTEGER,PARAMETER:: d1_gridpoint_code = 15 ! gridpoint info
      INTEGER,PARAMETER:: d1_proc_no_code   = 16 ! Processing Code
      INTEGER,PARAMETER:: d1_halo_type      = 17 ! Halo width type

      ! Types of items for d1_type

      INTEGER,PARAMETER:: prognostic = 0
      INTEGER,PARAMETER:: diagnostic = 1
      INTEGER,PARAMETER:: secondary  = 2
      INTEGER,PARAMETER:: other      = 3

! D1_ADDR end

INTEGER                                                           &
  d1_addr(d1_list_len,n_objs_d1)
                     ! IN D1 addressing info.

! Include files/common blocks/parameters
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

INTEGER                                                           &
  start_block                                                     &
                ! first word of field data in dump (unused)
, object_index                                                    &
                ! pointer to entry in D1_ADDR
, level                                                           &
                ! level number of multi-level field
, d1_item_code                                                    &
                ! sec/item in d1_addr converted into single code
, number_of_fields                                                &
                ! total number of fields to read in
, field_start                                                     &
                ! start address of a field in the file
, data_size                                                       &
                ! number of words of data on disk for a field
, data_read_size                                                  &
                ! total number of words to read for a field
, data_full_size                                                  &
                ! total number of words after any unpacking
, len_io                                                          &
                ! number of words of data successfully read
, k                                                               &
                ! loop counter over fields
, orig_decomp                                                     &
                ! current decomposition on entry
, address                                                         &
                ! start address of field in D1
, local_len     ! number of words of data put into D1 on this
                ! processor


LOGICAL                                                           &
  packed_field  ! TRUE if a field has been packed to 32 bits

! Error reporting
INTEGER       icode       ! =0 normal exit; >0 error exit
CHARACTER(LEN=256) cmessage    ! Error message
CHARACTER(LEN=*) routinename

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
PARAMETER (   routinename='UM_READDUMP')



!--------------------------------------------------------------

IF (lhook) CALL dr_hook('UM_READDUMP',zhook_in,zhook_handle)
icode=0
cmessage=''

IF (mype  ==  0) THEN
  WRITE(6,'(/,'' READING UNIFIED MODEL DUMP ON UNIT'',I3)')nftin
  WRITE(6,'('' #####################################'',/)')
END IF

! Change to the relevant decomposition type for this dump

orig_decomp=current_decomp_type

IF (submodel_id  ==  a_im) THEN
  IF (current_decomp_type  /=  decomp_standard_atmos)             &
  CALL change_decomposition(decomp_standard_atmos,icode)

ELSE  ! unsupported decomposition type
  WRITE(6,*)                                                      &
    'UM_READDUMP : Could not change to decomposition required ',  &
    'for submodel type ',submodel_id
  icode=1
  cmessage='Unsupported submodel for MPP code'

  CALL ereport(routinename,icode,cmessage)
END IF

IF (icode  /=  0) THEN
  icode=2
  WRITE(6,*) 'UM_READDUMP : Error - Could not set decomposition ',&
             'for selected submodel.'
  cmessage='Unsupported decomposition selected for MPP code'

  CALL ereport(routinename,icode,cmessage)
END IF

! Read in the header records, and do a consistency check

IF (read_header) THEN

! DEPENDS ON: readhead
  CALL readhead(nftin,fixhd,len_fixhd,                            &
                inthd,len_inthd,                                  &
                realhd,len_realhd,                                &
                levdepc,len1_levdepc,len2_levdepc,                &
                rowdepc,len1_rowdepc,len2_rowdepc,                &
                coldepc,len1_coldepc,len2_coldepc,                &
                flddepc,len1_flddepc,len2_flddepc,                &
                extcnst,len_extcnst,                              &
                dumphist,len_dumphist,                            &
                cfi1,len_cfi1,                                    &
                cfi2,len_cfi2,                                    &
                cfi3,len_cfi3,                                    &
                lookup,len1_lookup,len2_lookup,                   &
                len_data,                                         &
                start_block,icode,cmessage)

  IF (icode  /=  0) THEN
    WRITE(6,*) 'UM_READDUMP : Error reading dump header ',        &
               'on unit ',nftin
    WRITE(6,*) 'Return code from READHEAD was ',icode,            &
               ' and error message was ',cmessage
    icode=3
    cmessage='Error reading dump header'

    CALL ereport(routinename,icode,cmessage)
  END IF

END IF  ! IF (READ_HEADER)


IF (fixhd(160)  >   0) THEN ! If there is data to read

! Loop over fields and read into D1

  number_of_fields=fixhd(152)
  address=1
  object_index=1
  level=1
  local_len=0

  DO k=1,number_of_fields  ! loop over fields to read in
     
    ! Use previous value of local_len to increment address
    address = address+local_len
    ! Reset local_len
    local_len = 0
    mpp_lookup(p_lblrec,k)=0
    mpp_lookup(p_naddr,k)=address

    IF (lookup(lblrec,k)  >   0) THEN ! If there's data in
                                      ! the field

! Check that DATA_TYPE is valid no: +/-1 to +/-3
      IF (( ABS(lookup(data_type,k))  >=  1) .AND.                &
          ( ABS(lookup(data_type,k))  <=  3)) THEN

! Set "packed_field" to .TRUE. if 32bit packing has been used
        packed_field=(MOD((lookup(lbpack,k)),10) == 2)

! Check that the diagnostic in the dump matches that expected
! from D1_ADDR

        IF (d1_addr(d1_object_type,object_index)  ==  diagnostic) THEN

          d1_item_code= (d1_addr(d1_section,object_index)*1000) + &
                         d1_addr(d1_item,object_index)
          IF (lookup(item_code,k)  /=  d1_item_code) THEN
            WRITE(6,*)                                            &
              'UM_READDUMP : Dump field ',k,                      &
              ' does not match STASH request for item ',          &
              d1_addr(d1_item,object_index),                      &
              ' section ',d1_addr(d1_section,object_index)
            WRITE(6,*) 'Expected code ',lookup(item_code,k)
            cmessage='UM_READDUMP Dump does not match STASH list'
            icode=4

            CALL ereport(routinename,icode,cmessage)
          END IF ! IF (LOOKUP(ITEM_CODE,K)  /=  d1_item_code)
        END IF ! IF (D1_ADDR(d1_object_type,object_index)  ==
              !     diagnostic)

! Set up the location of the field on disk and how much data
! needs to be read in
        field_start=lookup(lbegin,k) ! position of field in file

! data_size contains the number of words to data used to store
! the field on disk
        IF (packed_field) THEN
          data_size=(lookup(lblrec,k)+1)/2
        ELSE
          data_size=lookup(lblrec,k)
        END IF

! data_read_size contains the number of words to data that need to
! be read in for a field. Each field has extra words of dummy data
! added at the end to ensure each field starts on a disk sector
! boundary. The last field on a dump does not have these extra words
! added
        IF (k  /=  number_of_fields) THEN
          data_read_size=lookup(lbnrec,k)
        ELSE
          data_read_size=data_size
        END IF

! This is the max of number of words required to store the field in
! memory after any unpacking is done, and number of words required
! to read in the data.

        data_full_size=MAX(lookup(lblrec,k),data_read_size)

! Move file pointer to the start of the field

        CALL setpos(nftin,field_start,icode)
        IF (icode  /=  0) THEN
          WRITE(6,*)                                              &
          'UM_READDUMP - SETPOS failed to move file pointer to ', &
          field_start,' on unit ',nftin
          WRITE(6,*) 'SETPOS returned error code ',icode
          icode=5
          cmessage='SETPOS failed while reading dump. See output.'

          CALL ereport(routinename,icode,cmessage)
        END IF

        IF (address > len_data) THEN
          IF (d1_addr(d1_length,object_index) > 0) THEN
            icode=5
            WRITE(6,'(A,I5,A,I9,A,I9,A,I9)') &
                       "Error for field ", k, " with calculated address ",    &
                       address, " and maximum D1 address of ", len_data,      &
                       " and D1 length of ", d1_addr(d1_length,object_index)
            cmessage='Error calculating address in D1'
            CALL ereport(routinename,icode,cmessage)
          ELSE
            ! A zero sized field (land packed maybe?) set to maximum.
            address = len_data
          END IF
        END IF
         
! DEPENDS ON: um_read_multi
        CALL um_read_multi(nftin,d1(address),data_read_size,         &
                           data_full_size,                           &
                           len_io,local_len,lookup(1,k),fixhd(12),   &
                           d1_addr(1,object_index),1,                &
                           icode,cmessage)

        mpp_lookup(p_lblrec,k)=local_len

        IF (icode  /=  0) THEN
          WRITE(6,*)                                              &
            'UM_READDUMP - Error while attempting to read field ',&
            k,' of ',number_of_fields,' from unit ',              &
            nftin
          WRITE(6,*) 'Return code from UM_READ_MULTI was ',icode, &
            'and error message was ',cmessage
          WRITE(6,*) 'Field Information: '
          WRITE(6,*) 'Section ',d1_addr(d1_section,object_index), &
                     ' Item ',d1_addr(d1_item,object_index)
          WRITE(6,*) 'Disk address : ',field_start
          WRITE(6,*) 'D1 address : ',address
          WRITE(6,*) 'Number of words requested : ',data_read_size
          WRITE(6,*) 'Number of words returned : ',len_io

          icode=6
          cmessage='Error reading field from dump'

          CALL ereport(routinename,icode,cmessage)

        END IF ! If an error was detected reading the field

      ELSE ! Error in LOOKUP(DATA_TYPE,K)

        IF (( fixhd(5)  <   6) .OR.                               &
                                     ! Not AC, Var or Cx
            ( fixhd(5)  >   8)) THEN
! DEPENDS ON: pr_look
          CALL pr_look(                                           &
                  lookup,lookup,len1_lookup,k)
        END IF

        WRITE(6,*) 'um_readdump : failure for field ',k,' of ',   &
                   number_of_fields
        WRITE(6,*) 'LOOKUP(DATA_TYPE,K)= ',lookup(data_type,k)
        icode=7
        cmessage='Invalid data type ( LOOKUP(DATA_TYPE,K) )'

        CALL ereport(routinename,icode,cmessage)

      END IF ! Check the LOOKUP(DATA_TYPE,K) is valid

    END IF ! If there was data in the field

    level=level+1
    IF (level  >   d1_addr(d1_no_levels,object_index)) THEN
      level=1
      object_index=object_index+1
    END IF

  END DO ! K : loop over fields to read in

  IF (printstatus  >=  prstatus_normal) THEN
    IF (mype  ==  0) THEN
      WRITE(6,*) 'Data successfully read'
      WRITE(6,*) fixhd(161),' words read from unit ',nftin
      IF ((fixhd(5)  >=  6) .AND.                                 &
                                   ! AC/Var
          (fixhd(5)  <=  8)) THEN ! Obs/Cx
        WRITE(6,*) '(Observational data)'
      ELSE
        WRITE(6,*) '(Model data)'
      END IF
    END IF ! IF (mype  ==  0)
  END IF ! IF (PrintStatus  >=  PrStatus_Normal)

END IF ! IF (FIXHD(160)  >   0)

! Reset to original decomposition type
CALL change_decomposition(orig_decomp,icode)

IF (lhook) CALL dr_hook('UM_READDUMP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE um_readdump
