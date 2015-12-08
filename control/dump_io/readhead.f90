! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE READHEAD---------------------------------------
!
!    Programming standard: Unified Model Documentation Paper No 3
!
!    Purpose: Reads in model dump header records on unit NFTIN and
!             checks model and dump dimensions for consistency.
!
!    Documentation: Unified Model Documentation Paper No F3
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Dump I/O

SUBROUTINE readhead(nftin,fixhd,len_fixhd,                        &
                                                 ! Intent (In)
                    inthd,len_inthd,                              &
                    realhd,len_realhd,                            &
                    levdepc,len1_levdepc,len2_levdepc,            &
                    rowdepc,len1_rowdepc,len2_rowdepc,            &
                    coldepc,len1_coldepc,len2_coldepc,            &
                    flddepc,len1_flddepc,len2_flddepc,            &
                    extcnst,len_extcnst,                          &
                    dumphist,len_dumphist,                        &
                    cfi1,len_cfi1,                                &
                    cfi2,len_cfi2,                                &
                    cfi3,len_cfi3,                                &
                    lookup,len1_lookup,len2_lookup,len_data,      &
                    start_block,icode,cmessage)   ! Intent (Out)
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE IO
USE PrintStatus_mod
USE UM_ParVars
USE lookup_addresses

IMPLICIT NONE

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
               !IN 2ndt dim of level dep consts
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
,start_block                                                      &
                !OUT Pointer to position of each block.
                !Should point to start of model data block on exit
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
,dumphist(len_dumphist+1)             !IN History block

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
INTEGER k
INTEGER len_io
INTEGER fixhd_152    !  Original value of FIXHD(152)
LOGICAL l_a_dump
LOGICAL l_ff    ! FieldsFile Logical
REAL a

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! -------------------------------------------------------------

IF (lhook) CALL dr_hook('READHEAD',zhook_in,zhook_handle)
icode=0
cmessage=' '

!  1. Buffer in fixed length header record


CALL buffin(nftin,fixhd,len_fixhd,len_io,a)


! Check for I/O errors
IF (a /= -1.0 .OR. len_io /= len_fixhd) THEN
! DEPENDS ON: ioerror
  CALL ioerror('buffer in of fixed length header',a,len_io        &
               ,len_fixhd)
  cmessage='READHEAD: I/O error'
  icode=1
  IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
  RETURN
END IF

start_block=len_fixhd+1

fixhd_152 = fixhd(152)    !  Store original value

!     Test if atmos dump read in
l_a_dump = fixhd(5) == 1 .AND. fixhd(2) == 1                      &
           .AND. len_data /= imdi

!      Test if FieldsFile read in
l_ff = fixhd(5) == 3

IF &
    (l_a_dump) &
    THEN
  IF (fixhd(152) /= len2_lookup) THEN
!XX       WRITE (6,*) 'FIXHD(152) being reset from ',FIXHD(152),' to ',
!XX  *    LEN2_LOOKUP
    fixhd(152) = len2_lookup
  END IF
END IF

! Check validity of data and print out fixed header information

IF (mype  ==  0) THEN
! DEPENDS ON: pr_fixhd
CALL pr_fixhd(fixhd,len_fixhd,len_inthd,len_realhd,len1_levdepc   &
,len2_levdepc,len1_rowdepc,len2_rowdepc,len1_coldepc,len2_coldepc &
,len1_flddepc,len2_flddepc,len_extcnst,len_dumphist,len_cfi1      &
,len_cfi2,len_cfi3,len1_lookup,len2_lookup,len_data               &
,icode,cmessage)

IF (icode >  0) THEN 
  IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
  RETURN
END IF

END IF
!  2. Buffer in integer constants

IF (fixhd(100) >  0) THEN

! Check for error in file pointers
  IF (fixhd(100) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('integer constants',start_block,100,fixhd(100))
    cmessage='READHEAD: Addressing conflict'
    icode=2
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffin(nftin,inthd,fixhd(101),len_io,a)

! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= fixhd(101)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of integer constants',a,len_io,       &
                 fixhd(101))
    cmessage='READHEAD: I/O error'
    icode=3
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(101)

END IF

!  3. Buffer in real constants

IF (fixhd(105) >  0) THEN

! Check for error in file pointers
  IF (fixhd(105) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('real constants',start_block,105,fixhd(105))
    cmessage='READHEAD: Addressing conflict'
    icode=4
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

! Check for I/O errors

  CALL buffin(nftin,realhd,fixhd(106),len_io,a)

  IF (a /= -1.0 .OR. len_io /= fixhd(106)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of real constants',a,len_io,          &
                  fixhd(106))
    cmessage='READHEAD: I/O error'
    icode=5
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(106)

END IF

!  4. Buffer in level dependent constants

IF (fixhd(110) >  0 .AND. len1_levdepc /= 0) THEN

! Check for error in file pointers
  IF (fixhd(110) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('level dependent constants',                    &
          start_block,110,fixhd(110))
    cmessage='READHEAD: Addressing conflict'
    icode=6
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffin(nftin,levdepc,fixhd(111)*fixhd(112),len_io,a)

! Check for I/O errors
  IF(a /= -1.0 .OR. len_io /= fixhd(111)*fixhd(112)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of level dependent constants',a,      &
                  len_io,fixhd(111)*fixhd(112))
    cmessage='READHEAD: I/O error'
    icode=7
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(111)*fixhd(112)

  IF (printstatus >= prstatus_oper) THEN
    IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' LEVEL DEPENDENT CONSTANTS'')')
  WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(111)*fixhd(112)
    END IF ! mype  ==  0
  END IF
END IF

!  5. Buffer in row dependent constants

IF (fixhd(115) >  0 .AND. len1_rowdepc /= 0) THEN

! Check for error in file pointers
  IF (fixhd(115) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('row dependent constants',                      &
          start_block,115,fixhd(115))
    cmessage='READHEAD: Addressing conflict'
    icode=8
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffin(nftin,rowdepc,fixhd(116)*fixhd(117),len_io,a)

! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= fixhd(116)*fixhd(117)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of row dependent constants',a,len_io, &
                  fixhd(116)*fixhd(117))
    cmessage='READHEAD: I/O error'
    icode=9
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(116)*fixhd(117)

  IF (printstatus >= prstatus_oper) THEN
    IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' ROW DEPENDENT CONSTANTS'')')
  WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(116)*fixhd(117)
    END IF ! mype  ==  0
  END IF
END IF

!  6. Buffer in column dependent constants

IF (fixhd(120) >  0 .AND. len1_coldepc /= 0) THEN

! Check for error in file pointers
  IF (fixhd(120) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('column dependent constants',                   &
          start_block,120,fixhd(120))
    cmessage='READHEAD: Addressing conflict'
    icode=10
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffin(nftin,coldepc,fixhd(121)*fixhd(122),len_io,a)

! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= fixhd(121)*fixhd(122)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of column dependent constants',a,     &
                  len_io,fixhd(121)*fixhd(122))
    cmessage='READHEAD: I/O error'
    icode=11
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(121)*fixhd(122)

  IF (printstatus >= prstatus_oper) THEN
    IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' COLUMN DEPENDENT CONSTANTS'')')
  WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(121)*fixhd(122)
    END IF ! mype  ==  0
  END IF
END IF

!  7. Buffer in constants stored as fields

IF (fixhd(125) >  0 .AND. len1_flddepc /= 0) THEN

! Check for error in file pointers
  IF (fixhd(125) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('fields of constants',                          &
    start_block,125,fixhd(125))
    cmessage='READHEAD: Addressing conflict'
    icode=12
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffin(nftin,flddepc,fixhd(126)*fixhd(127),len_io,a)

! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= fixhd(126)*fixhd(127)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of field dependent constants',a,      &
                  len_io,fixhd(126)*fixhd(127))
    cmessage='READHEAD: I/O error'
    icode=13
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(126)*fixhd(127)

  IF (printstatus >= prstatus_oper) THEN
    IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' FIELD DEPENDENT CONSTANTS'')')
  WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(126)*fixhd(127)

    END IF ! mype  ==  0
  END IF
END IF

!  8. Buffer in extra constants

IF (fixhd(130) >  0 .AND. len_extcnst /= 0) THEN

! Check for error in file pointers
  IF (fixhd(130) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('extra constants',                              &
    start_block,130,fixhd(130))
    cmessage='READHEAD: Addressing conflict'
    icode=14
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffin(nftin,extcnst,fixhd(131),len_io,a)

! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= fixhd(131)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in extra constants',a,len_io,            &
                  fixhd(131))
    cmessage='READHEAD: I/O error'
    icode=15
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(131)

  IF (printstatus >= prstatus_oper) THEN
    IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' EXTRA CONSTANTS'')')
      WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(131)
    END IF ! mype  ==  0
  END IF
END IF

!  9. Buffer in temporary history block

IF (fixhd(135) >  0 .AND. len_dumphist /= 0) THEN

! Check for error in file pointers
  IF (fixhd(135) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('history',                                      &
    start_block,136,fixhd(136))
    cmessage='READHEAD: Addressing conflict'
    icode=16
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffin(nftin,dumphist,fixhd(136),len_io,a)

! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= fixhd(136)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of history file',a,len_io,            &
                  fixhd(136))
    cmessage='READHEAD: I/O error'
    icode=17
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(136)

  IF (printstatus >= prstatus_oper) THEN
    IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' TEMPORARY HISTORY BLOCK'')')
      WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(136)

    END IF ! mype  ==  0
  END IF
END IF

!  10. Buffer in compressed field index1

IF (fixhd(140) >  0 .AND. len_cfi2 /= 0) THEN

! Check for error in file pointers

  IF (fixhd(140) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('compressed field index1',                      &
          start_block,140,fixhd(140))
    cmessage='READHEAD: Addressing conflict'
    icode=18
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffin(nftin,cfi1,fixhd(141),len_io,a)

! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= fixhd(141)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of compressed index1',a,len_io,       &
                  fixhd(141))
    cmessage='READHEAD: I/O error'
    icode=19
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(141)

  IF (printstatus >= prstatus_oper) THEN
    IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' COMPRESSED FIELD INDEX NO 1'')')
      WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(141)
    END IF ! mype  ==  0
  END IF
END IF

!  11. Buffer in compressed field index2

IF (fixhd(142) >  0 .AND. len_cfi2 /= 0) THEN

! Check for error in file pointers
  IF (fixhd(142) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('compressed field index2',                      &
          start_block,142,fixhd(142))
    cmessage='READHEAD: Addressing conflict'
    icode=20
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffin(nftin,cfi2,fixhd(143),len_io,a)

! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= fixhd(143)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of compressed index2',a,len_io,       &
                  fixhd(143))
    cmessage='READHEAD: I/O error'
    icode=21
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(143)

  IF (printstatus >= prstatus_oper) THEN
    IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' COMPRESSED FIELD INDEX NO 2'')')
      WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(143)
    END IF ! mype  ==  0
  END IF
END IF

!  12. Buffer in compressed field index3

IF (fixhd(144) >  0 .AND. len_cfi3 /= 0) THEN

! Check for error in file pointers
  IF (fixhd(144) /= start_block) THEN
! DEPENDS ON: poserror
    CALL poserror('compressed field index3',                      &
          start_block,144,fixhd(144))
    cmessage='READHEAD: Addressing conflict'
    icode=22
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF


  CALL buffin(nftin,cfi3,fixhd(145),len_io,a)

! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= fixhd(145)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of compressed index3',a,len_io,       &
                  fixhd(145))
    cmessage='READHEAD: I/O error'
    icode=23
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

  start_block=start_block+fixhd(145)

  IF (printstatus >= prstatus_oper) THEN
    IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' COMPRESSED FIELD INDEX NO 3'')')
      WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(145)
    END IF ! mype  ==  0
  END IF
END IF

!  13. Buffer in lookup table

IF (fixhd(150) >  0) THEN

! Supress checking if not full dump
  IF (len_dumphist /= 0) THEN
! Check for error in file pointers
    IF (fixhd(150) /= start_block) THEN
! DEPENDS ON: poserror
      CALL poserror('lookup table',                               &
            start_block,150,fixhd(150))
      cmessage='READHEAD: Addressing conflict'
      icode=24
      IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF

! Move to start of Look Up Table

  CALL setpos(nftin,fixhd(150)-1,icode)

! Read in fields from LOOKUP table

  CALL buffin(nftin,lookup(:,:),fixhd(151)*fixhd(152),len_io,a)

! Check for I/O errors
  IF (a /= -1.0 .OR. len_io /= fixhd(151)*fixhd(152)) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer in of lookup table',a,len_io,            &
                  fixhd(151)*fixhd(152))
    cmessage='READHEAD: I/O error'
    icode=25
    IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
    RETURN
  END IF

! Point to start of data section ( Use original FIXHD(152) )
  start_block=start_block+fixhd(151)*fixhd_152

  IF (printstatus >= prstatus_oper) THEN
    IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' LOOKUP TABLE'')')
  WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(151)*fixhd(152)

      IF (fixhd(152) <  fixhd_152) THEN
        WRITE(6,'('' '')')
      WRITE(6,'('' '',i6,'' Entries in Look Up Table.'')') fixhd_152
        WRITE(6,'('' '',i6,'' Entries read in.'')') fixhd(152)
      END IF

    END IF ! mype  ==  0
  END IF
!---------------------------------------------------------------
! Reset LOOKUP(45) if not set
!---------------------------------------------------------------

  DO k=1,len2_lookup
    IF (lookup(45,k) == 0 .OR. lookup(45,k) == imdi) THEN

      IF (printstatus >= prstatus_normal) THEN
        WRITE(6,*) 'WARNING: User defined field found - ',      &
                   'STASH code : ', lookup(42,k)
        WRITE(6,*) ' Setting internal model number to atmosphere.'
      END IF
      lookup(45,k)=1

    END IF

  END DO
!---------------------------------------------------------------
!  Reset LOOKUP headers if dump created earlier than vn2.8
!---------------------------------------------------------------

  IF (fixhd(12) <  208) THEN
! DEPENDS ON: newpack
    CALL newpack(lookup,len1_lookup,len2_lookup)
  END IF

! Check LOOKUP for consistency with PARAMETER statements
  IF (lookup(lbnrec,1) == 0 .OR.                                  &
!        Prog lookups in dump before vn3.2:
     (lookup(lbnrec,1) == imdi .AND. fixhd(12) <= 301)) THEN
    IF (len_data /= imdi) THEN
! DEPENDS ON: chk_look
      CALL chk_look(fixhd,lookup,len1_lookup,len_data,            &
                    icode,cmessage)
    END IF
  END IF

END IF

IF (lhook) CALL dr_hook('READHEAD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE readhead
