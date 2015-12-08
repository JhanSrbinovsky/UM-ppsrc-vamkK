! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  MAIN PROGRAM FOR MERGE UTILITY---------------------------------
!
!  Purpose:
!          This program was primarily written to merge two boundary
!          datasets for use with the mesoscale model model in .
!          test mode. It has been extended to cope with the merging
!          of any two SEQUENTIAL datasets in unified model format,
!          provided they are of the same type and resolution.
!          A namelist allows the user to specify the point of merging
!          and in the case of time series to merge the datasets at the
!          point the times overlap.
!
!          MAIN_MERGE reads in fixed length ,integer and lookup
!          headers of UM files to be merged, extracts dimensions
!          of each file, sets the dimensions of the merged file
!          to that of the first input file and then passes these
!          values to subroutine MERGE.
!
!  Documentation:
!          UM Doc Paper 3
!
!  -----------------------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Small execs
      PROGRAM MAIN_MERGE

      USE filenamelength_mod, ONLY:                                     &
          filenamelength
      USE IO
      USE ereport_mod, ONLY : ereport, ereport_finalise
      USE UM_Config, ONLY : &
          appInit, &
          exe_merge
      IMPLICIT NONE

      INTEGER                                                           &
     & FIXHD1(256)                                                      &
                         !Space for fixed length header file 1
     &,INTHD1(100)                                                      &
                         !Space for integer header file 1
     &,LOOKUP1(128)      !Space for lookup record file 1

      INTEGER                                                           &
     & FIXHD2(256)                                                      &
                         !Space for fixed length header file 2
     &,INTHD2(100)                                                      &
                         !Space for integer header file 2
     &,LOOKUP2(128)      !Space for lookup record file 2

      INTEGER                                                           &
     & LEN_FIXHD1                                                       &
                      !Length of fixed length header on file 1
     &,LEN_INTHD1                                                       &
                      !Length of integer header on file 1
     &,LEN_REALHD1                                                      &
                      !Length of real header on file 1
     &,LEN1_LEVDEPC1                                                    &
                      !1st dim of lev dependent consts on file 1
     &,LEN2_LEVDEPC1                                                    &
                      !2nd dim of lev dependent consts on file 1
     &,LEN1_ROWDEPC1                                                    &
                      !1st dim of row dependent consts on file 1
     &,LEN2_ROWDEPC1                                                    &
                      !2nd dim of row dependent consts on file 1
     &,LEN1_COLDEPC1                                                    &
                      !1st dim of col dependent consts on file 1
     &,LEN2_COLDEPC1                                                    &
                      !2nd dim of col dependent consts on file 1
     &,LEN1_FLDDEPC1                                                    &
                      !1st dim of field dependent consts on file 1
     &,LEN2_FLDDEPC1                                                    &
                      !2nd dim of field dependent consts on file 1
     &,LEN_EXTCNST1                                                     &
                      !Length of extra consts on file 1
     &,LEN_DUMPHIST1                                                    &
                      !Length of history header on file 1
     &,LEN_CFI11                                                        &
                      !Length of index1 on file 1
     &,LEN_CFI21                                                        &
                      !Length of index2 on file 1
     &,LEN_CFI31                                                        &
                      !Length of index3 on file 1
     &,LEN1_LOOKUP1                                                     &
                      !1st dim of LOOKUP on file 1
     &,LEN2_LOOKUP1                                                     &
                      !2nd dim of LOOKUP on file 1
     &,LEN_DATA1                                                        &
                      !Length of data on file 1
     &,ROW_LENGTH1                                                      &
                      !No of points E-W on file 1
     &,P_ROWS1                                                          &
                      !No of p-rows on file 1
     &,P_FIELD1       !No of p-points per level on file 1

      INTEGER                                                           &
     & LEN_FIXHD2                                                       &
                      !Length of fixed length header on file 2
     &,LEN_INTHD2                                                       &
                      !Length of integer header on file 2
     &,LEN_REALHD2                                                      &
                      !Length of real header on file 2
     &,LEN1_LEVDEPC2                                                    &
                      !1st dim of lev dependent consts on file 2
     &,LEN2_LEVDEPC2                                                    &
                      !2nd dim of lev dependent consts on file 2
     &,LEN1_ROWDEPC2                                                    &
                      !1st dim of row dependent consts on file 2
     &,LEN2_ROWDEPC2                                                    &
                      !2nd dim of row dependent consts on file 2
     &,LEN1_COLDEPC2                                                    &
                      !1st dim of col dependent consts on file 2
     &,LEN2_COLDEPC2                                                    &
                      !2nd dim of col dependent consts on file 2
     &,LEN1_FLDDEPC2                                                    &
                      !1st dim of field dependent consts on file 2
     &,LEN2_FLDDEPC2                                                    &
                      !2nd dim of field dependent consts on file 2
     &,LEN_EXTCNST2                                                     &
                      !Length of extra consts on file 2
     &,LEN_DUMPHIST2                                                    &
                      !Length of history header on file 2
     &,LEN_CFI12                                                        &
                      !Length of index1 on file 2
     &,LEN_CFI22                                                        &
                      !Length of index2 on file 2
     &,LEN_CFI32                                                        &
                      !Length of index3 on file 2
     &,LEN1_LOOKUP2                                                     &
                      !1st dim of LOOKUP on file 2
     &,LEN2_LOOKUP2                                                     &
                      !2nd dim of LOOKUP on file 2
     &,LEN_DATA2                                                        &
                      !Length of data on file 2
     &,ROW_LENGTH2                                                      &
                      !No of points E-W on file 2
     &,P_ROWS2                                                          &
                      !No of p-rows on file 2
     &,P_FIELD2       !No of p-points per level on file 2

      INTEGER                                                           &
     & LEN_FIXHD3                                                       &
                      !Length of fixed length header on file 3
     &,LEN_INTHD3                                                       &
                      !Length of integer header on file 3
     &,LEN_REALHD3                                                      &
                      !Length of real header on file 3
     &,LEN1_LEVDEPC3                                                    &
                      !1st dim of lev dependent consts on file 3
     &,LEN2_LEVDEPC3                                                    &
                      !2nd dim of lev dependent consts on file 3
     &,LEN1_ROWDEPC3                                                    &
                      !1st dim of row dependent consts on file 3
     &,LEN2_ROWDEPC3                                                    &
                      !2nd dim of row dependent consts on file 3
     &,LEN1_COLDEPC3                                                    &
                      !1st dim of col dependent consts on file 3
     &,LEN2_COLDEPC3                                                    &
                      !2nd dim of col dependent consts on file 3
     &,LEN1_FLDDEPC3                                                    &
                      !1st dim of field dependent consts on file 3
     &,LEN2_FLDDEPC3                                                    &
                      !2nd dim of field dependent consts on file 3
     &,LEN_EXTCNST3                                                     &
                      !Length of extra consts on file 3
     &,LEN_DUMPHIST3                                                    &
                      !Length of history header on file 3
     &,LEN_CFI13                                                        &
                      !Length of index1 on file 3
     &,LEN_CFI23                                                        &
                      !Length of index2 on file 3
     &,LEN_CFI33                                                        &
                      !Length of index3 on file 3
     &,LEN1_LOOKUP3                                                     &
                      !1st dim of LOOKUP on file 3
     &,LEN2_LOOKUP3                                                     &
                      !2nd dim of LOOKUP on file 3
     &,LEN_DATA3                                                        &
                      !Length of data on file 3
     &,ROW_LENGTH3                                                      &
                      !No of points E-W on file 3
     &,P_ROWS3                                                          &
                      !No of p-rows on file 3
     &,P_FIELD3       !No of p-points per level on file 3


      INTEGER      ERR          !Return code from open
      INTEGER      I            !Loop index
      INTEGER      LEN_IO       !Length of I/O returned by BUFFER IN
      INTEGER      MAX_LEN1     !Length of longest data record in file 1
      INTEGER      MAX_LEN2     !Length of longest data record in file 2
      INTEGER      NFTIN1       !Unit number of input UM file 1
      INTEGER      NFTIN2       !Unit number of input UM file 2

      INTEGER      ErrorStatus  !Error code returned from FILE_OPEN
      INTEGER      OpenStatus   !Error code returned from GET_FILE
      INTEGER      ICODE        !Error code returned from SETPOS

      REAL         A            !BUFFER IN UNIT function

      CHARACTER(LEN=filenamelength) :: filename     
                                !Name of user preSTASHmaster file


      CHARACTER(LEN=100) DUMMY_ENV   ! replaces deprecated string that would
                                     ! hold the executable path
      INTEGER ME_GC,NPROC_GC

      DUMMY_ENV='dummy path'

      CALL GC_INIT(DUMMY_ENV,ME_GC,NPROC_GC)
! Initialise print status for standard output
! DEPENDS ON: initprintstatus
      CALL InitPrintStatus
      CALL appInit(exe_merge)
      CALL ioInit()
!L 1. Assign unit numbers

      NFTIN1=20
      NFTIN2=21

      WRITE(6,*) " "
      WRITE(6,*)' MERGE UTILITY'
      WRITE(6,*)' -------------'
      WRITE(6,*)' '

      WRITE(6,'(20x,''FILE STATUS'')')
      WRITE(6,'(20x,''==========='')')
      CALL FILE_OPEN(NFTIN1,'FILE1',5,0,0,ERR)
      CALL FILE_OPEN(NFTIN2,'FILE2',5,0,0,ERR)

!L 2. Buffer in fixed length header record from file 1

      CALL BUFFIN(NFTIN1,FIXHD1,256,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= 256)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header of input file',  &
     &  A,LEN_IO,256)

      CALL EREPORT('MAIN_MERGE', 1000,                                  &
     & 'buffer in of fixed length header of input file wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,256
        IF(FIXHD1(I) <  0)FIXHD1(I)=0
      ENDDO

! Input file dimensions
      LEN_FIXHD1=256
      LEN_INTHD1=FIXHD1(101)
      LEN_REALHD1=FIXHD1(106)
      LEN1_LEVDEPC1=FIXHD1(111)
      LEN2_LEVDEPC1=FIXHD1(112)
      LEN1_ROWDEPC1=FIXHD1(116)
      LEN2_ROWDEPC1=FIXHD1(117)
      LEN1_COLDEPC1=FIXHD1(121)
      LEN2_COLDEPC1=FIXHD1(122)
      LEN1_FLDDEPC1=FIXHD1(126)
      LEN2_FLDDEPC1=FIXHD1(127)
      LEN_EXTCNST1=FIXHD1(131)
      LEN_DUMPHIST1=FIXHD1(136)
      LEN_CFI11=FIXHD1(141)
      LEN_CFI21=FIXHD1(143)
      LEN_CFI31=FIXHD1(145)
      LEN1_LOOKUP1=FIXHD1(151)
      LEN2_LOOKUP1=FIXHD1(152)
      LEN_DATA1=FIXHD1(161)

!L 3. Buffer in fixed length header record from file 2

      CALL BUFFIN(NFTIN2,FIXHD2,256,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= 256)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header of input file',  &
     &  A,LEN_IO,256)

      CALL EREPORT('MAIN_MERGE', 1001,                                  &
     & 'buffer in of fixed length header of input file')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,256
        IF(FIXHD2(I) <  0)FIXHD2(I)=0
      ENDDO

! Input file dimensions
      LEN_FIXHD2=256
      LEN_INTHD2=FIXHD2(101)
      LEN_REALHD2=FIXHD2(106)
      LEN1_LEVDEPC2=FIXHD2(111)
      LEN2_LEVDEPC2=FIXHD2(112)
      LEN1_ROWDEPC2=FIXHD2(116)
      LEN2_ROWDEPC2=FIXHD2(117)
      LEN1_COLDEPC2=FIXHD2(121)
      LEN2_COLDEPC2=FIXHD2(122)
      LEN1_FLDDEPC2=FIXHD2(126)
      LEN2_FLDDEPC2=FIXHD2(127)
      LEN_EXTCNST2=FIXHD2(131)
      LEN_DUMPHIST2=FIXHD2(136)
      LEN_CFI12=FIXHD2(141)
      LEN_CFI22=FIXHD2(143)
      LEN_CFI32=FIXHD2(145)
      LEN1_LOOKUP2=FIXHD2(151)
      LEN2_LOOKUP2=FIXHD2(152)
      LEN_DATA2=FIXHD2(161)


!L 4. Buffer in integer constants from file 1

       CALL BUFFIN(NFTIN1,INTHD1,FIXHD1(101),LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= FIXHD1(101))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of integer constants in input file 1',  &
     &  A,LEN_IO,FIXHD1(101))

      CALL EREPORT('MAIN_MERGE', 1002,                                  &
     & 'buffer in of integer constants in input file 1 wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,FIXHD1(101)
        IF(INTHD1(I) <  0)INTHD1(I)=0
      ENDDO

       ROW_LENGTH1=INTHD1(6)
       P_ROWS1=INTHD1(7)
       P_FIELD1=ROW_LENGTH1*P_ROWS1

!L 5. Buffer in integer constants from file 2

       CALL BUFFIN(NFTIN2,INTHD2,FIXHD2(101),LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= FIXHD2(101))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of integer constants in input file 2',  &
     &  A,LEN_IO,FIXHD2(101))

      CALL EREPORT('MAIN_MERGE', 1003,                                  &
     & 'buffer in of integer constants in input file 2 wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,FIXHD2(101)
        IF(INTHD2(I) <  0)INTHD2(I)=0
      ENDDO

!L 6. Cause abort if files obviously different

      ROW_LENGTH2=INTHD2(6)
      P_ROWS2=INTHD2(7)
      P_FIELD2=ROW_LENGTH2*P_ROWS2

      IF(P_FIELD1 /= P_FIELD2)THEN
       WRITE(6,*)'COMPARE: ERROR Dumps are at different resolutions'

       CALL EREPORT('MAIN_MERGE', 1004,                                 &
     &  'Dumps are at different resolutions')

      ENDIF

!L 7. Buffer in lookup table from file 1 and find largest record

      MAX_LEN1=0
      DO I=1,FIXHD1(152)
        CALL SETPOS(NFTIN1,FIXHD1(150)-1+64*(I-1),ICODE)
        CALL BUFFIN(NFTIN1,LOOKUP1,FIXHD1(151)                          &
            ,LEN_IO,A)

! Check for I/O errors
        IF(A /= -1.0.OR.LEN_IO /= FIXHD1(151))THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer in of lookup table in input file 1',     &
              A,LEN_IO,FIXHD1(151))

          CALL EREPORT('MAIN_MERGE', 1005,                              &
              'buffer in of lookup table in input file 1 wrong size')
          
        ENDIF
          
        MAX_LEN1=MAX0(LOOKUP1(15),MAX_LEN1)
        
      ENDDO

!L 8. Buffer in lookup table from file 2 and find largest record

      MAX_LEN2=0
      DO I=1,FIXHD2(152)
        CALL SETPOS(NFTIN2,FIXHD2(150)-1+64*(I-1),ICODE)
        CALL BUFFIN(NFTIN2,LOOKUP2,FIXHD2(151)                          &
            ,LEN_IO,A)

! Check for I/O errors
          IF(A /= -1.0.OR.LEN_IO /= FIXHD2(151))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of lookup table in input file 2',       &
     &  A,LEN_IO,FIXHD2(151))

            CALL EREPORT('MAIN_MERGE', 1006,                            &
     &       'buffer in of lookup table in input file 2 wrong size')

          ENDIF

         MAX_LEN2=MAX0(LOOKUP2(15),MAX_LEN2)

       ENDDO

! Enlargen size of I/O buffer if smaller than field.
        IF(P_FIELD1 <  MAX_LEN1)                                        &
     &    P_FIELD1=MAX_LEN1
        IF(P_FIELD2 <  MAX_LEN2)                                        &
     &    P_FIELD2=MAX_LEN2

! Rewind files
      CALL SETPOS(NFTIN1,0,ICODE)
      CALL SETPOS(NFTIN2,0,ICODE)

!L 9. Output file dimensions.

! Set equal to the dimensions of file 1 initially.
! Most will not need to be changed for merged file.
      LEN_FIXHD3=256
      LEN_INTHD3=FIXHD1(101)
      LEN_REALHD3=FIXHD1(106)
      LEN1_LEVDEPC3=FIXHD1(111)
      LEN2_LEVDEPC3=FIXHD1(112)
      LEN1_ROWDEPC3=FIXHD1(116)
      LEN2_ROWDEPC3=FIXHD1(117)
      LEN1_COLDEPC3=FIXHD1(121)
      LEN2_COLDEPC3=FIXHD1(122)
      LEN1_FLDDEPC3=FIXHD1(126)
      LEN2_FLDDEPC3=FIXHD1(127)
      LEN_EXTCNST3=FIXHD1(131)
      LEN_DUMPHIST3=FIXHD1(136)
      LEN_CFI13=FIXHD1(141)
      LEN_CFI23=FIXHD1(143)
      LEN_CFI33=FIXHD1(145)
      LEN1_LOOKUP3=FIXHD1(151)
      LEN2_LOOKUP3=FIXHD1(152)+FIXHD2(152)
      LEN_DATA3=FIXHD1(161)+FIXHD2(161)
      P_FIELD3=P_FIELD1

!L 10. Call MERGE

! DEPENDS ON: merge
      CALL MERGE(LEN_FIXHD1,LEN_INTHD1,LEN_REALHD1,                     &
     &  P_ROWS1,P_ROWS2,P_ROWS3,                                        &
     &  ROW_LENGTH1,ROW_LENGTH2,ROW_LENGTH3,                            &
     &  LEN1_LEVDEPC1,LEN2_LEVDEPC1,LEN1_ROWDEPC1,                      &
     &  LEN2_ROWDEPC1,LEN1_COLDEPC1,LEN2_COLDEPC1,                      &
     &  LEN1_FLDDEPC1,LEN2_FLDDEPC1,LEN_EXTCNST1,                       &
     &  LEN_DUMPHIST1,LEN_CFI11,LEN_CFI21,LEN_CFI31,                    &
     &  LEN1_LOOKUP1,LEN2_LOOKUP1,LEN_DATA1,P_FIELD1,                   &
     &  LEN_FIXHD2,LEN_INTHD2,LEN_REALHD2,                              &
     &  LEN1_LEVDEPC2,LEN2_LEVDEPC2,LEN1_ROWDEPC2,                      &
     &  LEN2_ROWDEPC2,LEN1_COLDEPC2,LEN2_COLDEPC2,                      &
     &  LEN1_FLDDEPC2,LEN2_FLDDEPC2,LEN_EXTCNST2,                       &
     &  LEN_DUMPHIST2,LEN_CFI12,LEN_CFI22,LEN_CFI32,                    &
     &  LEN1_LOOKUP2,LEN2_LOOKUP2,LEN_DATA2,P_FIELD2,                   &
     &  LEN_FIXHD3,LEN_INTHD3,LEN_REALHD3,                              &
     &  LEN1_LEVDEPC3,LEN2_LEVDEPC3,LEN1_ROWDEPC3,                      &
     &  LEN2_ROWDEPC3,LEN1_COLDEPC3,LEN2_COLDEPC3,                      &
     &  LEN1_FLDDEPC3,LEN2_FLDDEPC3,LEN_EXTCNST3,                       &
     &  LEN_DUMPHIST3,LEN_CFI13,LEN_CFI23,LEN_CFI33,                    &
     &  LEN1_LOOKUP3,LEN2_LOOKUP3,LEN_DATA3,P_FIELD3                    &
     & ,NFTIN1,NFTIN2)

      CALL ioShutdown()

      CALL ereport_finalise( )

      END PROGRAM MAIN_MERGE
