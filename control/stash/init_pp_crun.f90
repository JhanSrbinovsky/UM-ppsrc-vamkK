! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Initialises buffer space for lookups and headers for PP files
!   on a Crun.
!
! Method:
!    Simply allocates space, reads in the data and then attaches
!    data to the pointers in the PP buffer structures.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

SUBROUTINE Init_PP_Crun( ftn_unit, filename, &
    len1_lookup, pp_len2_lookup,              &
    filetype)
  
  USE Model_file, ONLY :     &
      InitHeader,            &
      attachHeader,          &
      LoadHeader,            &
      setHeader,             &
      FixedHeader,           &
      InitLookups,           &
      AttachLookups,         &
      model_file_open
  USE IOS, ONLY :            &
      IOS_MergeLookup
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE filenamelength_mod, ONLY : filenamelength
  USE IO
  USE PrintStatus_mod
  USE ereport_mod, ONLY : ereport
  USE UM_ParVars
  USE nlstcall_mod, ONLY : pp_len2_look, &
                           model_status

  USE chsunits_mod, ONLY : nunits

  IMPLICIT NONE

! Subroutine arguments
  INTEGER :: ftn_unit
  INTEGER :: len1_lookup
  INTEGER :: pp_len2_lookup
  INTEGER :: len_fixhd
  
  CHARACTER(LEN=filenamelength) :: filename
  
  CHARACTER(LEN=1) :: filetype
  
! Local variables
  INTEGER :: icode
  INTEGER :: len_io
  INTEGER :: isize
  INTEGER :: address_ipplook
  INTEGER :: dummy
  INTEGER :: step
  INTEGER :: i
  INTEGER :: wordsRemaining
  
  INTEGER, PARAMETER :: current_io_pe=0
  
  REAL :: ierr

! Pointer for the fixed length header
  INTEGER, POINTER :: pp_fixhd(:)

! Pointer for the lookup table
  INTEGER, POINTER :: ipplook(:,:)

  CHARACTER(LEN=*), PARAMETER :: RoutineName='Init_PP_Crun'
  CHARACTER(LEN=80)           :: Cmessage

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  

  IF (lhook) CALL dr_hook('INIT_PP_CRUN',zhook_in,zhook_handle)

! first we need to open the preattached file

  CALL Model_File_Open(ftn_unit, filename, LEN_TRIM(filename), 1, 1, icode)
  IF (icode /= 0) THEN
    cmessage='Error opening file'    
    CALL Ereport(RoutineName, icode, cmessage)
  END IF
      
! Is the file a pp file ?
  IF ((filetype == 'p') .OR. (filetype == 'c')) THEN

! allocate the arrays
    ALLOCATE( ipplook(len1_lookup,pp_len2_lookup) )
    
!  attach the fixed length header
    CALL initHeader(ftn_unit, FixedHeader)
    CALL loadHeader(ftn_unit, FixedHeader, icode=icode)
    pp_fixhd=>attachHeader(ftn_unit, FixedHeader)

    IF (icode /= 0) THEN
      icode = -10
      CALL Ereport(RoutineName, icode, 'Error: Failed to load fixed header.')
!     In a non-operational run, skip to the end of this routine.
      IF (MODEL_STATUS /= 'Operational') GO TO 9999
    END IF

    ! If we are using an async stash server we need to propogate 
    ! the fixed length header to the server
    CALL setHeader(ftn_unit,FixedHeader)

! init the lookup table
    step=1
    CALL initLookups(ipplook, ftn_unit, len1_lookup,               &
        pp_len2_lookup,address_ipplook,step)
    
! the position of the lookup table is pp_fixhd(150)-1
    IF (mype == 0) THEN
      address_ipplook = pp_fixhd(150)-1
    END IF
    
    step = 2
    CALL initLookups(ipplook, ftn_unit, dummy, dummy,              &
        address_ipplook, step)

! attach the lookup table
    IF (mype == current_io_pe) THEN
      ipplook=>attachLookups(ftn_unit)
    END IF

! read the lookup table

    CALL setpos(ftn_unit, address_ipplook, icode )
    IF (icode /= 0) THEN
      Cmessage = 'Error in setpos'
      CALL Ereport(RoutineName, icode, cmessage)
    END IF
    
! Only current_io_pe should do the read as only this PE has 
! memory allocated for this operation!
    IF (mype == current_io_pe) THEN
      CALL set_unit_bcast_flag(ftn_unit)
      wordsRemaining=readableWords(ftn_unit)
      IF (printstatus>=prstatus_oper)             &
          WRITE(6,'(A,I3,A,I10,A)')               &
          'init_pp_crun: unit ',ftn_unit,' has ', &
          wordsRemaining,' readable words for lookup read'
      
      IF (wordsRemaining>= len1_lookup*pp_len2_lookup) THEN
        CALL buffin(ftn_unit, ipplook,                         &
            len1_lookup*pp_len2_lookup, len_io, ierr )

        CALL IOS_MergeLookup(ftn_unit, ipplook, &
            len1_lookup*pp_len2_lookup)
      ELSE
        icode = -20
        Cmessage = 'Error reading lookup table'
        
        CALL Ereport(RoutineName, icode, cmessage)
!     In a non-operational run, skip to the end of this routine.
        IF (model_status /= 'Operational') THEN
          Go To 9999
        END IF
      END IF
      CALL clear_unit_bcast_flag(ftn_unit)            
    END IF ! current_io_pe

! Nullify ipplook and pp_fixhd for safety
    NULLIFY(ipplook)
    NULLIFY(pp_fixhd)
    
  END IF ! is_ppfile
      
9999 CONTINUE

  IF (lhook) CALL dr_hook('INIT_PP_CRUN',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE init_pp_crun
