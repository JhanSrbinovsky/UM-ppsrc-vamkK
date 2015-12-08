
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Write out fields to a UM format file.

SUBROUTINE writflds ( unit,                                       &
                                    ! in
                      numfields,                                  &
                                    ! in
                      position,                                   &
                                    ! in
                      lookup,                                     &
                                    ! in
                      len1lookup,                                 &
                                    ! in
                      d1,                                         &
                                    ! in
                      buflen,                                     &
                                    ! in
                      fixhd,                                      &
                                    ! in
                      icode,                                      &
                                    ! out
                      cmessage )    ! out

! Description:

!   Buffers out NumFields fields from D1 to a UM format file on unit
!   Unit, starting at field number Position.


! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dump I/O

! Declarations:

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IO
USE UM_ParVars
USE io_configuration_mod, ONLY: io_field_padding
USE lookup_addresses
USE cppxref_mod, ONLY: &
    ppx_atm_rim,       &
    ppx_ocn_rim,       &
    ppx_wam_rim,       &
    ppx_halo_type,     &    
    ppx_grid_type
USE Submodel_Mod

IMPLICIT NONE

! Include files:
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

! Subroutine arguments:

INTEGER,       INTENT(IN) ::                                      &
  len1lookup            ! Lookup dimension.

INTEGER,       INTENT(IN) ::                                      &

  unit,                                                           &
                        ! Unit number of UM file.
  numfields,                                                      &
                        ! Number of fields to write out.
  position,                                                       &
                        ! Field number from which to begin I/O.
  lookup(len1lookup,*)  ! Lookup table for UM file.

REAL,          INTENT(IN) ::                                      &

  d1(*)                 ! Array containing fields to be written.

INTEGER,       INTENT(IN) ::                                      &

  buflen,                                                         &
                        ! Length of I/O buffer.
  fixhd(*)              ! Fixed-length header for UM file.

INTEGER,       INTENT(OUT) ::                                     &

  icode                 ! Return code. >0 => error.

CHARACTER(LEN=80), INTENT(OUT) ::                                     &

  cmessage              ! Error message.

! Local variables:

INTEGER :: i, j, k,                                               &
           lenio,                                                 &
           ipts,                                                  &
                           ! No. of values to be written to disk.
           wordaddress,                                           &
                           ! Address from which to begin I/O.
           l_ipts,                                                &
                           ! Record length during index search.
           um_sector_ipts,                                        &
                           ! No. of words to write, rounded up.
           ipts_write      ! No. of words written to disk.

REAL :: buf( ( (buflen+io_field_padding-1)/io_field_padding ) *       &
             io_field_padding )

!dir$ cache_align buf

INTEGER :: address,                                               &
                           ! Start address of D1 field.
           locallen,                                              &
           item,                                                  &
           section,                                               &
           model,                                                 &
           jcode,                                                 &
                           ! Return code.
           fld_type,                                              &
           grid_type,                                             &
           halo_type,                                             &
           fake_d1_addr(d1_list_len) ! Fake D1_addr record to
                                     ! be fed to write_multi.

! Integer functions called:

INTEGER  exppxi, get_fld_type

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header ------------------------------------------------------

!----------------------------------------------------------------------
! [1]: Initialize.
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('WRITFLDS',zhook_in,zhook_handle)
address = 1

icode    = 0
j        = 0
cmessage = ' '

!----------------------------------------------------------------------
! [2]: Write fields.
!----------------------------------------------------------------------

DO k = POSITION, POSITION + numfields - 1

 ! See whether data is stored in 32-bit format on disk.
  IF (MOD(lookup(lbpack, k), 10) == 2) THEN
    ipts = ( lookup(lblrec, k) + 1 ) / 2
  ELSE
    ipts = ( lookup(lblrec, k) + 1 )
  END IF

  IF ( lookup(lbnrec, k) == 0    .OR.                             &
                                      ! Old format dumps.
       lookup(lbnrec, k) == imdi .OR.                             &
                                      ! \ Ocean ACOBS
       lookup(lbegin, k) == imdi .OR.                             &
                                      ! / files?
       ( lookup(lbnrec, k) == imdi                                &
                                      ! \ Prognostic lookup
         .AND. fixhd(12) < 301 )                                  &
                                      ! / in pre-vn3.2 dump.
     ) THEN

    wordaddress = 1
    DO i = 2, k
      IF (MOD(lookup(lbpack, i-1), 10) == 2) THEN
        l_ipts = ( lookup(lblrec, i-1) + 1 ) / 2
      ELSE
        l_ipts =   lookup(lblrec, i-1)
      END IF
      wordaddress = wordaddress + l_ipts
    END DO
    wordaddress = wordaddress + fixhd(160) - 2
    um_sector_ipts = ipts

  ELSE ! PP type files and post vn4.3 dumps.

    wordaddress = lookup(lbegin, k)

   ! Use the stored rounded-up value:
    um_sector_ipts = lookup(lbnrec, k)

  END IF

  ipts_write = um_sector_ipts

! Position file pointer:

  CALL setpos (unit, wordaddress, icode)

!----------------------------------------------------------------------
! [2.2]: MPP write.
!----------------------------------------------------------------------

! Set up fake D1_addr record:

  DO i = 1, d1_list_len
    fake_d1_addr(i) = 0
  END DO

  item    = MOD(lookup(item_code, k), 1000)
  section = ( lookup(item_code, k) - item ) / 1000
  model   = lookup(model_code, k)

! DEPENDS ON: exppxi
  halo_type = exppxi ( model, section, item, ppx_halo_type,       &
                       jcode, cmessage )

! DEPENDS ON: exppxi
  grid_type = exppxi ( model, section, item, ppx_grid_type,       &
                       jcode, cmessage )

  IF (jcode /= 0) THEN
    WRITE (6,*) ''
    WRITE (6,*) 'WRITFLDS: Failed to get PPXREF info.'
    WRITE (6,*) ''
    WRITE (6,*) '  Model ID:      ', model
    WRITE (6,*) '  Section:       ', section
    WRITE (6,*) '  Item:          ', item
    WRITE (6,*) '  Error code:    ', jcode
    WRITE (6,*) '  Error message: ', cmessage
    WRITE (6,*) ''
    icode    = 1
    cmessage = 'WRITFLDS: Failed to get PPXREF info.'
    IF (lhook) CALL dr_hook('WRITFLDS',zhook_out,zhook_handle)
    RETURN
  END IF

  fake_d1_addr(d1_object_type) = prognostic
  fake_d1_addr(d1_imodl)       = model
  fake_d1_addr(d1_section)     = section
  fake_d1_addr(d1_item)        = item
  fake_d1_addr(d1_halo_type)   = halo_type

 ! Grid type: for LBCs we need some special logic...
  IF (lookup(lbhem, k) == 99) THEN

    IF ( lookup(model_code, k) == atmos_im ) THEN
      fake_d1_addr(d1_grid_type) = ppx_atm_rim
    ELSE
      icode = 2
      WRITE (6,*) ''
      WRITE (6,*) 'WRITFLDS: Cannot process LBC for model type',  &
                             lookup(model_code, k)
      WRITE (6,*) ''
      cmessage = 'Cannot write LBCs for this model type.'
      IF (lhook) CALL dr_hook('WRITFLDS',zhook_out,zhook_handle)
      RETURN
    END IF

  ELSE ! Not an LBC.

    fake_d1_addr(d1_grid_type) = grid_type

  END IF

! DEPENDS ON: get_fld_type
  fld_type = get_fld_type(grid_type)

  fake_d1_addr(d1_length)    = lasize(1, fld_type, halo_type) *   &
                               lasize(2, fld_type, halo_type)
  fake_d1_addr(d1_no_levels) = 1

! Write field:

! DEPENDS ON: write_multi
  CALL write_multi (                                              &
    unit,         d1(address), um_sector_ipts,                    &
                                                 ! in
    lenio,        locallen,                                       &
                                                 ! out
    lookup(1,k),  fixhd(12),   buf,                               &
                                                 ! in
    fake_d1_addr,                                                 &
                                                 ! in
    jcode,        cmessage )                     ! out

  address = address + locallen
  IF(locallen == 0)THEN
    address=address+lookup(lblrec,k)
  END IF

! Check for errors:

  IF (lenio /= um_sector_ipts) THEN

    WRITE (6,*) ''
    WRITE (6,*) 'WRITFLDS: Error writing field number ', k,       &
                         ' on unit ', unit
    WRITE (6,*) ''
    WRITE (6,*) '  write_multi error code:    ', jcode
    WRITE (6,*) '  write_multi error message: ', cmessage
    WRITE (6,*) ''
    icode    = jcode + 1
    cmessage = 'WRITFLDS: Error from write_multi.'
    IF (lhook) CALL dr_hook('WRITFLDS',zhook_out,zhook_handle)
    RETURN

  END IF

END DO ! k


IF (lhook) CALL dr_hook('WRITFLDS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE writflds
