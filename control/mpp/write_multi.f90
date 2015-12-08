! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Parallel UM interface to BUFFOUT

! Subroutine Interface:
SUBROUTINE write_multi(nft,d1,isize,len_io,local_len,             &
                       lookup,fixhd12,compbuf,                    &
                       addr_info,icode,cmessage)
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE IO
  USE UM_ParVars
  USE Decomp_DB
  USE lookup_addresses

  IMPLICIT NONE
  
! Description:
!  This routine provides an interface to BUFFOUT for the parallel
!  Unified Model. It is used where each process must write out a
!  local section of a global field.

! Method:
!  Each processor sends its local part of the global field to PE 0
!  which assembles all the parts, and then writes them to disk.
!  Fields which are compressed to land points are expanded before
!  sending to PE 0, PE 0 then compresses the global field before
!  writing it to disk.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine Arguments:

INTEGER, INTENT(IN) :: nft             !   IN : FORTRAN unit number
INTEGER, INTENT(IN) :: isize           !   IN : no. of words to write out
INTEGER, INTENT(OUT) :: len_io         !  OUT : no. of words written out
INTEGER, INTENT(OUT) :: local_len      !  OUT : size of the local field written out
INTEGER, INTENT(IN) :: lookup(64)      !   IN : LOOKUP header from dump
INTEGER, INTENT(IN) :: fixhd12         !   IN : 12th. element of fixed length header
                                       !        required for packing fields

! Required for dimensioning ADDR_INFO
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

INTEGER, INTENT(IN) ::  addr_info(d1_list_len) ! IN addressing info about field
REAL             :: d1(*)                      !   IN : Array to write out
REAL             :: compbuf(*)                 !   IN : Workspace for compressing field

INTEGER, INTENT(OUT) ::  icode                 !  OUT : Return code
CHARACTER(LEN=80) ::   cmessage                !  OUT : Error message

! Parameters and Common blocks

! CSMID start
!
! Description:
!    Hold parameters defining internal model identifiers and submodel
!    data partition (ie main D1 data array and consequent dump), both
!    short and long form.
!
!
! Declarations:
!
      ! Internal models
      INTEGER,PARAMETER:: A_IM      = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: ATMOS_IM  = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: O_IM      = 2 ! Ocean internal model
      INTEGER,PARAMETER:: OCEAN_IM  = 2 ! Ocean internalmodel
      INTEGER,PARAMETER:: S_IM      = 3 ! Slab internal model
      INTEGER,PARAMETER:: SLAB_IM   = 3 ! Slab internal model
      INTEGER,PARAMETER:: W_IM      = 4 ! Wave internal model
      INTEGER,PARAMETER:: WAVE_IM   = 4 ! Wave internal model
!      INTEGER,PARAMETER:: I_IM      = 5 ! Sea=ice internal model
!      INTEGER,PARAMETER:: SEAICE_IM = 5 ! Sea=ice internal model
      ! New dynamics (Charney-Phillips grid)
!      INTEGER,PARAMETER:: N_IM      = 6 ! ND internal model
!      INTEGER,PARAMETER:: NATMOS_IM = 6 ! ND internal model
      ! Small Executables
      INTEGER,PARAMETER:: X_IM      = 7 ! SX indicator

      ! Submodels
      INTEGER,PARAMETER:: A_SM      = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: ATMOS_SM  = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: O_SM      = 2 ! Ocean submodel
      INTEGER,PARAMETER:: OCEAN_SM  = 2 ! Ocean submodel
      INTEGER,PARAMETER:: W_SM      = 4 ! Wave submodel
      INTEGER,PARAMETER:: WAVE_SM   = 4 ! Wave submodel
      ! New dynamics (Charney-Phillips grid)
!      INTEGER,PARAMETER:: N_SM      = 6 ! ND submodel
!      INTEGER,PARAMETER:: NATMOS_SM = 6 ! ND submodel
      ! Small Executables
      INTEGER,PARAMETER:: X_SM      = 7 ! SX indicator

! CSMID end

REAL ::  buf(isize*2)      ! Buffer for holding data to be written.
                           ! Factor of two is incase the data is
                           ! packed (ISIZE is the size on disk)
!DIR$ CACHE_ALIGN buf
INTEGER  ::  i             ! loop counter

REAL  ::  buf_icode           ! return code from buffout

INTEGER ::  orig_decomp    ! decomposition on entry
INTEGER ::  new_decomp     ! decomposition to change to

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ------------------------------------------------------------------
IF (lhook) CALL dr_hook('WRITE_MULTI',zhook_in,zhook_handle)
len_io=isize
local_len=0

orig_decomp=current_decomp_type
new_decomp=orig_decomp

IF ((addr_info(d1_imodl)  ==  atmos_im) .AND.                     &
    (orig_decomp  /=  decomp_standard_atmos)) THEN

  new_decomp=decomp_standard_atmos
END IF

IF (new_decomp  /=  orig_decomp) THEN

  CALL change_decomposition(new_decomp,icode)

  IF (icode  /=  0) THEN
    WRITE(6,*) 'WTMULT : Error changing to decomposition ',       &
      new_decomp
    WRITE(6,*)                                                    &
      'Attempting to write field (Model,Section,Item) ',          &
      addr_info(d1_imodl),                                        &
      addr_info(d1_section),                                      &
      addr_info(d1_item)
    icode=1
    cmessage='Failure changing decomposition'
    GO TO 9999
  END IF

END IF

! Gather the field from the local D1 array to buf


! DEPENDS ON: general_gather_field
CALL general_gather_field(                                        &
  d1,buf,local_len,lookup(lblrec),1,                              &
  addr_info,0,                                                    &
  -1,icode,cmessage)

  IF (icode  /=  0) THEN
    WRITE(6,*) 'WTMULT : Call to GENERAL_GATHER_FIELD failed'
    WRITE(6,*) 'Return code was ',icode
    WRITE(6,*) 'Error message was ',cmessage
    WRITE(6,*) 'Field number ',lookup(item_code)
    WRITE(6,*) 'dimensions ',lookup(lbnpt),' x ',lookup(lbrow)
    WRITE(6,*) 'Grid type ',addr_info(d1_grid_type)
    WRITE(6,*) 'Field was not written out'

    icode=300
    cmessage='Failure to gather field'
    GO TO 9999
  END IF

! ------------------------------------------------------------------
! And finally the code to write the global field in array buf
! out to disk.

IF (mype  ==  0) THEN
!       Does this field need to be compressed?
  IF(MOD((lookup(lbpack)),10)  ==  2) THEN
    IF(lookup(data_type)  ==  1) THEN
! DEPENDS ON: pack21
      CALL pack21(lookup(lblrec),buf,                             &
                  compbuf)
    END IF
  ELSE ! no compression required - just do a copy
    DO i=1,lookup(lblrec)
      compbuf(i)=buf(i)
    END DO
  END IF

! Now write out the global field

  IF(MOD((lookup(lbpack)),10)  ==  2) THEN
    IF(lookup(data_type)  ==  1) THEN
! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*ISIZE 32 bit words using BUFFOUT32_F77 (as the array is 64 bit)
! DEPENDS ON : buffout32_f77
      CALL buffout32_f77(nft,compbuf(1:2*isize),2*isize,len_io,buf_icode)
! And then halve LEN_IO to satisfy tests against ISIZE
      len_io = len_io/2
    END IF
  ELSE
! For non-packed data
    CALL buffout(nft,compbuf(1:isize),isize,len_io,buf_icode)
  END IF
  IF ((buf_icode  /=  -1.0) .OR. (len_io  /=  isize)) THEN
    WRITE(6,*) 'WTMULT : Error in call to BUFFOUT'
    WRITE(6,*) 'LEN_IO : ',len_io
    WRITE(6,*) 'IOSTAT : ',buf_icode
    WRITE(6,*)                                                    &
      'Attempting to read field (Model,Section,Item) ',           &
      addr_info(d1_imodl),                                        &
      addr_info(d1_section),                                      &
      addr_info(d1_item)
    icode=400
    cmessage='Failure writing out field'
    GO TO 9999
  END IF

END IF ! am I PE 0 ?


! If the field was compressed for writing on disk, we need to compress
! and expand the field in memory. This ensures the same field exists in
! memory that would exist if this dump was read back in.

IF(MOD((lookup(lbpack)),10)  ==  2) THEN
  IF(lookup(data_type)  ==  1) THEN
! DEPENDS ON: pack21
    CALL pack21(local_len,d1,compbuf)
! DEPENDS ON: expand21
    CALL expand21(local_len,compbuf,d1)
  END IF
END IF

IF (new_decomp  /=  orig_decomp) THEN  ! change back

  CALL change_decomposition(orig_decomp,icode)

END IF
 9999 CONTINUE  ! point to jump to if there is a failure

IF (lhook) CALL dr_hook('WRITE_MULTI',zhook_out,zhook_handle)
RETURN
END SUBROUTINE write_multi
