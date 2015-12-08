! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  PROGRAM MAIN_CONVPP --------------------------------------------
!LL
!LL  Purpose: Converts a UM file into PP format.
!LL
!LL  Programming standards:
!LL
!LL  Logical components covered:
!LL
!LL  System Tasks: F3,F4,F6
!LL
!LL  Documentation: UM Doc Paper F5
!LL
!LL  -----------------------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Small execs
      PROGRAM MAIN_CONVPP

      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE IO
      USE ereport_mod, ONLY : ereport, ereport_finalise
      USE UM_Config, ONLY : &
          appInit, &
          exe_convpp
      IMPLICIT NONE

      INTEGER, PARAMETER :: NmelistInput = 5

      CHARACTER(LEN=filenamelength) :: arg1,arg2  ! Filenames


      INTEGER                                                           &
     & FIXHD(256)                                                       &
                         !Space for fixed length header
     &,INTHD(100)        !Space for integer header

      INTEGER                                                           &
     & LEN_FIXHD                                                        &
                      !Length of fixed length header on input file
     &,LEN_INTHD                                                        &
                      !Length of integer header on input file
     &,JOC_NO_SEAPTS                                                    &
                      !Number of points in compressed array
     &,LEN_OCFLD                                                        &
                      !Length of uncompressed ocean field
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
                !Unit number of input UM dump 1
     &,ERR                                                              &
                !Return code from OPEN

     &,ICODE    !Return code from setpos
      REAL A    !BUFFER IN UNIT function

      CHARACTER(LEN=100) DUMMY_ENV  ! replaces deprecated string that would
                                    ! hold the executable path
      INTEGER ME_GC,NPROC_GC

      DUMMY_ENV='dummy path'

      CALL GC_INIT(DUMMY_ENV,ME_GC,NPROC_GC)
! Initialise print status for standard output
! DEPENDS ON: initprintstatus
      CALL InitPrintStatus
      CALL appInit(exe_convpp)
      CALL ioInit()


! DEPENDS ON: ini_ppheader
      CALL INI_PPHEADER(NmelistInput)

!L 1. Assign unit numbers

      NFTIN=20

      WRITE(6,'(20x,''FILE STATUS'')')
      WRITE(6,'(20x,''==========='')')

      CALL FILE_OPEN(20,'FILE1',5,0,0,ERR)
      CALL GET_FILE(10,ARG2,filenamelength,ICODE)
      OPEN(10,FILE=ARG2,FORM='UNFORMATTED')

!L 2. Buffer in fixed length header record

      CALL BUFFIN(NFTIN,FIXHD,256,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= 256)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header of input dump',  &
     &  A,LEN_IO,256)

      CALL EREPORT('MAIN_CONVPP', 1000,                                 &
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

      CALL EREPORT('MAIN_CONVPP', 1001,                                 &
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
       CALL FIND_MAX_FIELD_SIZE(                                        &
     &      NFTIN,LEN1_LOOKUP,LEN2_LOOKUP,FIXHD,MAX_FIELD_SIZE)

! Calculate sizes of compressed and uncompressed ocean fields
      JOC_NO_SEAPTS=INTHD(11)
      IF(FIXHD(2) == 2)THEN
        LEN_OCFLD=INTHD(6)*INTHD(7)*INTHD(8)
      ELSE
        LEN_OCFLD=0
      ENDIF
! Rewind file
      CALL SETPOS(NFTIN,0,ICODE)

      IF(FIXHD(2) == 1)THEN

! Atmospheric dump
! DEPENDS ON: atmos_convpp
      CALL ATMOS_CONVPP (LEN_FIXHD,LEN_INTHD,LEN_REALHD,                &

     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  NFTIN,MAX_FIELD_SIZE)

      ELSEIF (FIXHD(2) == 2)THEN

! Oceanic dump
! DEPENDS ON: ocean_convpp
      CALL OCEAN_CONVPP (LEN_FIXHD,LEN_INTHD,LEN_REALHD,                &
     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  NFTIN,MAX_FIELD_SIZE,JOC_NO_SEAPTS,LEN_OCFLD)
      ENDIF

      CALL ioShutdown()

      CALL ereport_finalise( )

      CALL gc_exit()

      END PROGRAM MAIN_CONVPP
