
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Parallel UM interface to BUFFIN - primarily for small execs

! Subroutine Interface:
SUBROUTINE read_multi(nft,d1,isize,unpack_size,len_io,local_len,  &
                      lookup,fixhd12,                             &
                      addr_info,n_levels,                         &




                      icode,cmessage)
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IO
USE UM_ParVars
USE Decomp_DB
USE lookup_addresses
USE cppxref_mod, ONLY: ppx_atm_ozone, ppx_atm_tzonal

IMPLICIT NONE

! Description:
!  This routine provides an interface to BUFFIN for the parallel
!  Unified Model. It is used where each process must read in a
!  local section of a global field.

! Method:
!  PE 0 reads in the global field, and then distributes the
!  relevant parts of it to each processor.
!  Fields compressed to land points are expanded by PE 0, and
!  recompressed after being received by the relevant processor.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine Arguments:

INTEGER, INTENT(IN) ::  nft        ! FORTRAN unit number
INTEGER, INTENT(IN) ::  isize      ! no. of words to be read in (global field)
INTEGER, INTENT(IN) ::  unpack_size ! no. of words after any unpacking (global)
INTEGER, INTENT(OUT) :: len_io     ! no. of words read in (global field)
INTEGER, INTENT(OUT) :: local_len  !  no. of words in local field






INTEGER, INTENT(IN) :: fixhd12     ! 12th element of fixed length header
INTEGER, INTENT(IN) :: n_levels    ! Number of levels
INTEGER, INTENT(OUT) :: icode      ! Return code

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

INTEGER, INTENT(INOUT) :: lookup(64)  ! LOOKUP header from dump
INTEGER, INTENT(INOUT) :: addr_info(d1_list_len)   
                                   ! addressing info about field

CHARACTER(LEN=80) ::  cmessage     ! Error message
REAL, INTENT(OUT) ::  d1(*)        ! Array to read data in to

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

! Local variables
INTEGER :: info          ! GCOM return code
INTEGER :: i             ! loop counter
INTEGER :: map(n_levels) ! Which processor holds which level

INTEGER :: orig_decomp  ! decomposition on entry
INTEGER :: new_decomp   ! decomposition to change to

INTEGER :: grid_code    ! grid code
REAL    :: buf_icode    ! return code from BUFFIN

REAL    :: buf(unpack_size)       ! buffer for reading the field into
REAL    :: local_buf(unpack_size) ! and distributed - this is too big

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ------------------------------------------------------------------

IF (lhook) CALL dr_hook('READ_MULTI',zhook_in,zhook_handle)
len_io=isize
local_len=0

grid_code = addr_info(d1_grid_type)

orig_decomp=current_decomp_type
new_decomp=orig_decomp

IF ((addr_info(d1_imodl)  ==  atmos_im) .AND.                     &
    (orig_decomp  /=  decomp_standard_atmos)) THEN

  new_decomp=decomp_standard_atmos

ELSE IF ((addr_info(d1_imodl)  ==  ocean_im) .AND.                &
        (addr_info(d1_object_type)  ==  prognostic) .AND.         &
        (orig_decomp  /=  decomp_standard_ocean)) THEN

  new_decomp=decomp_standard_ocean

ELSE IF ((addr_info(d1_imodl)  ==  ocean_im) .AND.                &
        (addr_info(d1_object_type)  /=  prognostic) .AND.         &
        (orig_decomp  /=  decomp_nowrap_ocean)) THEN

  new_decomp=decomp_nowrap_ocean

END IF

IF (new_decomp  /=  orig_decomp) THEN

  CALL change_decomposition(new_decomp,icode)

  IF (icode  /=  0) THEN
    WRITE(6,*) 'READ_MULTI : Error changing to decomposition ',   &
      new_decomp
    WRITE(6,*)                                                    &
      'Attempting to read field (Model,Section,Item) ',           &
      addr_info(d1_imodl),                                        &
      addr_info(d1_section),                                      &
      addr_info(d1_item)
    icode=1
    cmessage='Failure changing decomposition'
    GO TO 9999
  END IF

END IF

! First thing to do is to read the field in to PE 0

IF (mype  ==  0) THEN
  CALL set_unit_bcast_flag(nft)! Switch off IO broadcasting
  IF (MOD((lookup(lbpack)),10)  ==  2) THEN
! Data is packed using CRAY 32 bit method - note that we need to read
! in 2*ISIZE 32 bit words using BUFFIN32_F77 (as the array is 64 bit)
! DEPENDS ON : buffin32_f77
    CALL buffin32_f77(nft,buf,2*isize,len_io,buf_icode)
! And then halve LEN_IO to satisfy tests against ISIZE

    len_io = len_io/2
  ELSE
! For non-packed data
    CALL buffin(nft,buf,isize,len_io,buf_icode)
  END IF
  CALL clear_unit_bcast_flag(nft)! Restore broadcast behaviour

!       Has the data been read in OK?
  IF ((buf_icode  /=  -1.0) .OR. (len_io  /=  isize)) THEN
    WRITE(6,*) 'READ_MULTI : Error in call to BUFFIN'
    WRITE(6,*) 'LEN_IO : ',len_io
    WRITE(6,*) 'IOSTAT : ',buf_icode
    WRITE(6,*)                                                    &
      'Attempting to read field (Model,Section,Item) ',           &
      addr_info(d1_imodl),                                        &
      addr_info(d1_section),                                      &
      addr_info(d1_item)
    icode=100
    cmessage='Failure reading in field'
    GO TO 9999
  END IF

! If it's a compressed REAL field, expand it out (except for fields
! we want to distribute in 32-bit)
! DEPENDS ON: read_unpack
  CALL read_unpack(buf,isize,unpack_size,lookup,fixhd12,            &
                   icode,cmessage)

  IF (icode /= 0) THEN
    cmessage='Failure unpacking field'
    GO TO 9999
  END IF
END IF ! IF (mype  ==  0)

! For atmosphere zonal ozone fields - set to zonal grid type
IF (( grid_code  ==  ppx_atm_ozone) .AND.           &
    ( lookup(lbnpt)  ==  1)) THEN

  addr_info(d1_grid_type)=ppx_atm_tzonal

END IF
! Create map showing which processor contains which level.
! At present, all levels are stored on PE 0
DO i=1,n_levels
  map(i)=0
END DO


! Now decompose the field in buf to the local D1 arrays
! DEPENDS ON: general_scatter_field
CALL general_scatter_field(d1, buf, local_len, lookup(lblrec),        &
                           n_levels, addr_info, map, icode, cmessage)

IF (icode  ==  1) THEN
  WRITE(6,*) 'READ_MULTI : Call to GENERAL_SCATTER_FIELD failed'
  WRITE(6,*) 'Return code was ',icode
  WRITE(6,*) 'Error message was ',cmessage
  WRITE(6,*) 'Field number ',lookup(item_code)
  WRITE(6,*) 'Grid type ', grid_code
  WRITE(6,*) 'Dimensioned : ',                                    &
               lookup(lbnpt),' x ',lookup(lbrow)
  icode=300
  cmessage='Failure decomposing field'
  GO TO 9999
END IF

! A sync can help avoid MPI congestion here.
CALL gc_gsync(nproc, icode)


IF (new_decomp  /=  orig_decomp) THEN  ! change back
  CALL change_decomposition(orig_decomp,icode)
END IF

9999  CONTINUE

IF (lhook) CALL dr_hook('READ_MULTI',zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_multi
