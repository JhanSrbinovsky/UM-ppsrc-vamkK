! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

MODULE Rcf_Read_Multi_Mod

  USE lookup_addresses

  IMPLICIT NONE

CONTAINS

! Subroutine Rcf_Read_Multi - parallel interface to buffin
!
! Description:
!  This routine provides an interface to BUFFIN for the parallel
!  reconfiguration. It is used where each process must read in a
!  local section of a global field.
!
! Method:
!  PE 0 reads in the global field, and then distributes the
!  relevant parts of it to each processor.
!  Fields compressed to land points are expanded by PE 0, and
!  recompressed after being received by the relevant processor.

  SUBROUTINE Rcf_Read_Multi(  NFT, D1, ISIZE, EXPAND_SIZE, LEN_IO,  &
      LOCAL_LEN, IOSTAT, LOOKUP, FIXHD12,                           &
      Stash_Record )
    
    USE IO, ONLY :            &
        Buffin,               &
        set_unit_bcast_flag,  &
        clear_unit_bcast_flag

    USE UM_ParVars, ONLY : &
        mype,           &
        nproc,          &
        glsizep

    USE Ereport_Mod, ONLY : &
        Ereport

    USE PrintStatus_Mod, ONLY :     &
        PrintStatus,                &
        PrStatus_Normal

    USE Rcf_Average_Polar_Mod, ONLY : &
        Rcf_Average_Polar

    USE Rcf_General_Scatter_Field_Mod, ONLY : &
        Rcf_General_Scatter_Field

    USE Rcf_Ppx_Info_Mod, ONLY : &
        STM_Record_Type

    USE lookup_addresses

    IMPLICIT NONE

    ! Subroutine Arguments:
    INTEGER, INTENT(IN)   :: nft        ! Fortain unit number
    INTEGER, INTENT(IN)   :: isize      ! no. of words to read (glob)
    INTEGER, INTENT(IN)   :: expand_size! size of expanded field
    INTEGER, INTENT(IN)   :: Fixhd12    ! 12th element of fixed header
    INTEGER, INTENT(IN)   :: Lookup(64) ! lookup table
    INTEGER, INTENT(OUT)  :: Len_IO     ! no. of words read in (glob)
    INTEGER, INTENT(OUT)  :: Local_Len  ! no. of local field words

    REAL,    INTENT(OUT)  :: IOstat     ! Return Code
    REAL,    INTENT(OUT)  :: D1(*)      ! Array to read data into

    TYPE( STM_Record_Type ), INTENT(IN) :: Stash_Record

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
    INTEGER      :: info            ! GCOM return code
    INTEGER      :: dimx            ! x dimension returned from COEX
    INTEGER      :: dimy            ! y dimension returned from COEX
    INTEGER      :: idum            ! dummy for COEX
    INTEGER      :: ErrorStatus     ! Error code returned from COEX
    LOGICAL      :: averaged        ! Have polar rows been averaged?
    REAL         :: buf( expand_size ) ! Buffer for reading data into
    REAL         :: buf2( expand_size) ! Temp. buffer for unpacking.
    CHARACTER (Len=80) :: Cmessage

    INTEGER, PARAMETER           :: len_full_word = 64 ! Output word size
                                                       ! for COEX
    CHARACTER (LEN=*), PARAMETER :: RoutineName = 'Rcf_Read_Multi'

! ------------------------------------------------------------------

    iostat      = -1.0
    len_io      = isize
    local_len   = 0
    errorstatus = 0

    ! First thing to do is to read the field in to PE 0

    CALL set_unit_bcast_flag(nft) ! Switch off IO broadcasting

    IF (mype == 0) THEN

      ! Check for 32-bit packed data and call the correct buffin version
      ! accordingly.  This should also handle 32-bit packed lateral boundary
      ! data.
      IF (MOD((lookup(lbpack)),10) == 2) THEN
        ! Data is packed using CRAY 32 bit method

        IF (lookup(lbhem) == 99) THEN
          ! For packed lateral boundary data ISIZE is the 
          ! number of 32-bit words to be read
          len_io = isize
! DEPENDS ON : buffin32_f77
          CALL buffin32_f77(nft,buf,isize,len_io,iostat)
        ELSE
          ! For other packed fields we need to read in 2*ISIZE 32 bit words
          ! using BUFFIN32
! DEPENDS ON : buffin32_f77
          CALL buffin32_f77(nft,buf,2*isize,len_io,iostat)
          ! And then halve len_io to satisfy tests against ISIZE
          len_io = len_io/2
        ENDIF
      ELSE
        ! For non-packed data
        CALL buffin(nft,buf,isize,len_io,iostat)
      ENDIF
      
      ! We must check to see if it is a 32 bit field on disk, and if
      ! so, expand it before doing anything with it.
      IF (MOD((lookup(lbpack)),10) == 2) THEN
        IF (lookup(data_type) == 1) THEN
          ! For special case of lateral boundary data, the length
          ! is given by ISIZE.
          IF (lookup(lbhem) == 99) THEN
! DEPENDS ON: expand32b
            CALL expand32b( isize , buf, fixhd12 )
          ELSE
! DEPENDS ON: expand32b
            CALL expand32b( lookup(lblrec) , buf, fixhd12 )
          ENDIF
        ELSE
          iostat=100
        ENDIF

        ! What about WGDOS packing?
      ELSE IF (MOD((lookup(lbpack)),10) == 1) THEN

        ! temporary copy
        buf2(:) = buf(:)
! DEPENDS ON: coex
        CALL coex( buf, expand_size, buf2, isize, dimx, dimy,  &
            idum, idum, .false., rmdi, len_full_word,   &
            errorstatus, cmessage )

        IF ( ErrorStatus /= 0 ) Then
          CALL Ereport( RoutineName, ErrorStatus, Cmessage )
        ENDIF
      ENDIF

      ! Lets do a polar row check/averaging (on p-grids) for read data only
      IF (lookup(lbhem) == 0 .AND. lookup(lbegin) == 1 .AND.  &
          Stash_Record % grid_type <= 3 ) THEN
        CALL Rcf_Average_Polar( buf, glsizep(2), glsizep(1), .TRUE., averaged)

        !  Print a warning if we've done averaging
        IF (averaged .AND. PrintStatus >= PrStatus_Normal) THEN
          WRITE(6,*) 'Field had had polar rows averaged in Rcf_Read_Multi'
        END IF
      END IF

    END IF  ! (mype == 0)

    CALL clear_unit_bcast_flag(nft)! Restore broadcast behaviour

    ! Now we can distribute it out to the other processes
    CALL Rcf_General_Scatter_Field( D1, buf, &
        local_len, lookup(lblrec), stash_record, 0 )

    RETURN
  END SUBROUTINE Rcf_Read_Multi

END MODULE Rcf_Read_Multi_Mod
