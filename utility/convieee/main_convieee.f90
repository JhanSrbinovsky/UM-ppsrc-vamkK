! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! PROGRAM MAIN_CONVIEEE ------------------
!
!
!
!    Purpose: Converts a dump, ancillary or fieldsfile
!             (atmosphere or ocean) between different
!             into 32-bit or 64-bit IEEE format or vice-versa.
!             The following conversions are supported:-
!               On a IEEE machine
!                 64-bit IEEE to 32-bit IEEE
!                 64-bit IEEE to 64-bit IEEE
!             In either case, WGDOS data will be unpacked if
!             requested.
!
!             MAIN_CONVIEEE reads in fixed length and integer
!             headers of UM file to be converted, extracts dimensions
!             of file and then passes these values to
!             subroutine CONVIEEE.
!
!    Documentation: UM Doc Paper F5
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Small execs
      PROGRAM MAIN_CONVIEEE
      USE IO
        USE ereport_mod, ONLY : ereport, ereport_finalise
        USE UM_Config, ONLY : &
            appInit, &  
            exe_convieee
        IMPLICIT NONE

      INTEGER                                                           &
     & FIXHD(256)                                                       &
                      !Space for fixed length header
     &,INTHD(100)     !Space for integer header

      INTEGER                                                           &
     & LEN_FIXHD                                                        &
                      !Length of fixed length header on input file
     &,LEN_INTHD                                                        &
                      !Length of integer header on input file
     &,LEN_REALHD                                                       &
                      !Length of real header on input file
     &,LEN1_LEVDEPC                                                     &
                      !1st dim of lev dependent consts on input file
     &,LEN2_LEVDEPC                                                     &
                      !2nd dim of lev dependent consts on input file
     &,LEN1_ROWDEPC                                                     &
                      !1st dim of row dependent consts on input file
     &,LEN2_ROWDEPC                                                     &
                      !2nd dim of row dependent consts on input file
     &,LEN1_COLDEPC                                                     &
                      !1st dim of col dependent consts on input file
     &,LEN2_COLDEPC                                                     &
                      !2nd dim of col dependent consts on input file
     &,LEN1_FLDDEPC                                                     &
                      !1st dim of field dependent consts on input file
     &,LEN2_FLDDEPC                                                     &
                      !2nd dim of field dependent consts on input file
     &,LEN_EXTCNST                                                      &
                      !Length of extra consts on input file
     &,LEN_DUMPHIST                                                     &
                      !Length of history header on input file
     &,LEN_CFI1                                                         &
                      !Length of index1 on input file
     &,LEN_CFI2                                                         &
                      !Length of index2 on input file
     &,LEN_CFI3                                                         &
                      !Length of index3 on input file
     &,LEN1_LOOKUP                                                      &
                      !1st dim of LOOKUP on input file
     &,LEN2_LOOKUP                                                      &
                      !2nd dim of LOOKUP on input file
     &,LEN_DATA                                                         &
                      !Length of data on input file
     &,ROW_LENGTH                                                       &
                      !No of points E-W on input file
     &,P_ROWS                                                           &
                      !No of p-rows on input file
     &,P_FIELD                                                          &
                      !No of p-points per level on input file
     &,MAX_FIELD_SIZE !Maximum field size on file

      INTEGER                                                           &
     & LEN_IO                                                           &
                !Length of I/O returned by BUFFER IN
     &,I                                                                &
                !Loop index
     &,NFTIN                                                            &
                !Unit number of input UM dump
     &,NFTOUT                                                           &
                !Unit number of output IEEE dump
     &,ERR                                                              &
                !Return code from OPEN
     &,IEEE_TYPE                                                        &
                 ! Output file precision
     &,ICODE                                                            &
                !Return code from setpos
     &,LEN                                                              &
                !Length of string returned by PXFGETARG
     &,IERR     !Return code from PXFGETARG

      CHARACTER(LEN=5) INPREC
      CHARACTER(LEN=80)                                                     &
     & STRING    ! Character string holding command line arg

      REAL A    !Return code from BUFFIN; -1.0 = O.K.

      integer wgdos_expand

      CHARACTER(LEN=100) DUMMY_ENV   ! replaces deprecated string that would 
                                     ! hold the executable path
      INTEGER ME_GC,NPROC_GC


      DUMMY_ENV='dummy path'

      CALL GC_INIT(DUMMY_ENV,ME_GC,NPROC_GC)
! Initialise print status for standard output
! DEPENDS ON: initprintstatus
      CALL InitPrintStatus
      CALL appInit(exe_convieee)
      CALL ioInit()

!--select no WGDOS expansion
      wgdos_expand=0

!L 0. Read in precision of output file
      CALL FORT_GET_ENV("PRECISION",9,INPREC,5,IERR)
      READ(INPREC,'(a5)') string
      LEN = LEN_TRIM(string)
      IF(LEN /= 2.OR.IERR /= 0)THEN
        IEEE_TYPE=32
        if(len == 3 .and. ierr == 0) then
          if(string == '32e' .or. string == '32E') then
            wgdos_expand=1
          else if(string == '64e' .or. string == '64E') then
            ieee_type=64
            wgdos_expand=1
          else
            WRITE(6,*)'Unsupported word length ',STRING

            CALL EREPORT('MAIN_CONVIEEE', 1000,                         &
     &       'Unsupported word length')

          endif
        endif
      ELSE
        IF(STRING == '32')THEN
          IEEE_TYPE=32
        ELSEIF(STRING == '64')THEN
          IEEE_TYPE=64
        ELSE
          WRITE(6,*)'Unsupported word length ',STRING

          CALL EREPORT('MAIN_CONVIEEE', 1002,                           &
     &     'Unsupported word length')

        ENDIF
      ENDIF
!
      IF(WGDOS_EXPAND == 0) THEN
        WRITE(6,'(/''Conversion to IEEE  '',i2,''-bit Format'',         &
     &    '' with no expansion of WGDOS Fields''/)') ieee_type
      ELSE
        WRITE(6,'(/''Conversion to IEEE  '',i2,''-bit Format'',         &
     &    '' with expansion of WGDOS Fields''/)') ieee_type
      END IF


!L 1. Assign unit numbers

      NFTIN=20
      NFTOUT=21

      WRITE(6,'(20x,''FILE STATUS'')')
      WRITE(6,'(20x,''==========='')')
!     CALL OPEN(1,'PPXREF',6,0,0,ERR)
      CALL FILE_OPEN(NFTIN,'FILE1',5,0,0,ERR)
      CALL FILE_OPEN(NFTOUT,'FILE2',5,1,0,ERR)


!L 2. Buffer in fixed length header record

      CALL BUFFIN(NFTIN,FIXHD,256,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= 256)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header of input dump',  &
     &  A,LEN_IO,256)

      CALL EREPORT('MAIN_CONVIEEE', 1003,                               &
     & 'Buffer in of fixed length header of input dump wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,256
        IF(FIXHD(I) <  0)FIXHD(I)=0
      ENDDO

! Input file dimensions (ensure sizes are >= 0 for allocating)
      LEN_FIXHD    = 256
      LEN_INTHD    = MAX(FIXHD(101),0)
      LEN_REALHD   = MAX(FIXHD(106),0)
      LEN1_LEVDEPC = MAX(FIXHD(111),0)
      LEN2_LEVDEPC = MAX(FIXHD(112),0)
      LEN1_ROWDEPC = MAX(FIXHD(116),0)
      LEN2_ROWDEPC = MAX(FIXHD(117),0)
      LEN1_COLDEPC = MAX(FIXHD(121),0)
      LEN2_COLDEPC = MAX(FIXHD(122),0)
      LEN1_FLDDEPC = MAX(FIXHD(126),0)
      LEN2_FLDDEPC = MAX(FIXHD(127),0)
      LEN_EXTCNST  = MAX(FIXHD(131),0)
      LEN_DUMPHIST = MAX(FIXHD(136),0)
      LEN_CFI1     = MAX(FIXHD(141),0)
      LEN_CFI2     = MAX(FIXHD(143),0)
      LEN_CFI3     = MAX(FIXHD(145),0)
      LEN1_LOOKUP  = MAX(FIXHD(151),0)
      LEN2_LOOKUP  = MAX(FIXHD(152),0)
      LEN_DATA     = MAX(FIXHD(161),0)


!L 3. Buffer in integer constants from dump

       CALL BUFFIN(NFTIN,INTHD,LEN_INTHD,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= LEN_INTHD)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of integer constants in input dump',    &
     &  A,LEN_IO,LEN_INTHD)

      CALL EREPORT('MAIN_CONVIEEE', 1004,                               &
     & 'Buffer in of integer constants in input dump wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,LEN_INTHD
        IF(INTHD(I) <  0)INTHD(I)=0
      ENDDO

       ROW_LENGTH=INTHD(6)
       P_ROWS=INTHD(7)
       P_FIELD=ROW_LENGTH*P_ROWS

!L Extract maximum field size from LOOKUP header
! DEPENDS ON: find_max_field_size
      CALL FIND_MAX_FIELD_SIZE(                                         &
     &      NFTIN,LEN1_LOOKUP,LEN2_LOOKUP,FIXHD,MAX_FIELD_SIZE,         &
     &      wgdos_expand)
! Rewind file
      CALL SETPOS(NFTIN,0,ICODE)

!L 4. Call CONVIEEE

! DEPENDS ON: convieee
      CALL CONVIEEE(LEN_FIXHD,LEN_INTHD,LEN_REALHD,                     &
     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  P_ROWS,ROW_LENGTH,                                              &
     &  NFTIN,NFTOUT,IEEE_TYPE,                                         &
     &  MAX_FIELD_SIZE, WGDOS_EXPAND)

      CALL file_close(NFTIN,'FILE1',5,0,0,ERR)
      CALL file_close(NFTOUT,'FILE2',5,0,0,ERR)

      CALL ioShutdown()

      CALL ereport_finalise( )

      CALL gc_exit()

      END PROGRAM MAIN_CONVIEEE
