! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE WRITHEAD---------------------------------------
!
!    Purpose: Writes out model dump header records on unit NFTOUT and
!             checks model and dump dimensions for consistency.
!             32-bit IEEE output option supported
!             64-bit IEEE output option supported
!
!
!    Programming standard:
!             Unified Model Documentation Paper No 3
!
!    Documentation:
!             Unified Model Documentation Paper No F3
!
!             Code Owner: See Unified Model Code Owners HTML page
!             This file belongs in section: Dump I/O


MODULE writhead_mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE writhead(nftout,fixhd,len_fixhd,                                  &
      inthd,len_inthd,                                                         &
      realhd,len_realhd,                                                       &
      levdepc,len1_levdepc,len2_levdepc,                                       &
      rowdepc,len1_rowdepc,len2_rowdepc,                                       &
      coldepc,len1_coldepc,len2_coldepc,                                       &
      flddepc,len1_flddepc,len2_flddepc,                                       &
      extcnst,len_extcnst,                                                     &
      dumphist,len_dumphist,                                                   &
      cfi1,len_cfi1,                                                           &
      cfi2,len_cfi2,                                                           &
      cfi3,len_cfi3,                                                           &
      lookup,len1_lookup,len2_lookup,len_data,                                 &
      output_bits,                                                             &
      start_block,icode,cmessage)

    USE parkind1, ONLY: jprb, jpim
    USE yomhook, ONLY: lhook, dr_hook
    USE io
    USE UM_types
    USE ereport_mod, ONLY : ereport
    USE PrintStatus_mod
    USE UM_ParVars
    USE um_config
    IMPLICIT NONE

    INTEGER, INTENT(IN) ::                                                     &
        nftout,                 &! Unit no of dump
        len_fixhd,              &! Length of fixed length header
        len_inthd,              &! Length of integer header
        len_realhd,             &! Length of real header
        len1_levdepc,           &! 1st dim of level dep consts
        len2_levdepc,           &! 2ndt dim of level dep consts
        len1_rowdepc,           &! 1st dim of row dep consts
        len2_rowdepc,           &! 2nd dim of row dep consts
        len1_coldepc,           &! 1st dim of column dep consts
        len2_coldepc,           &! 2nd dim of column dep consts
        len1_flddepc,           &! 1st dim of field dep consts
        len2_flddepc,           &! 2nd dim of field dep consts
        len_extcnst,            &! Length of extra constants
        len_dumphist,           &! Length of history block
        len_cfi1,               &! Length of comp field index 1
        len_cfi2,               &! Length of comp field index 2
        len_cfi3,               &! Length of comp field index 3
        len1_lookup,            &! 1st dim of lookup
        len2_lookup,            &! 2nd dim of lookup
        output_bits,            &! IEEE precision (use only by ieee)
        len_data                 !Length of model data

    INTEGER, INTENT(OUT) :: start_block   !Pointer to position of each block.
                                !Should point to start of model data block
                                !on exit

    INTEGER, INTENT(OUT) :: icode      ! Return code; successful=0
                                       ! error > 0

    CHARACTER(LEN=80) :: cmessage       ! Error message if ICODE > 0

    INTEGER, INTENT(IN) ::                                                     &
        fixhd(len_fixhd),       &! Fixed length header
        inthd(len_inthd),       &! Integer header
        lookup(len1_lookup,len2_lookup), &! PP lookup tables
        cfi1(len_cfi1+1),       &! Compressed field index no 1
        cfi2(len_cfi2+1),       &! Compressed field index no 2
        cfi3(len_cfi3+1)         ! Compressed field index no 3

    REAL, INTENT(IN) ::                                                        &
        realhd(len_realhd),                   &! Real header
        levdepc(1+len1_levdepc*len2_levdepc), &! Lev dep consts
        rowdepc(1+len1_rowdepc*len2_rowdepc), &! Row dep consts
        coldepc(1+len1_coldepc*len2_coldepc), &! Col dep consts
        flddepc(1+len1_flddepc*len2_flddepc),                                  &
                                 ! Field dep consts
        extcnst(len_extcnst+1), &! Extra constants
        dumphist(len_dumphist+1) ! History block

! Local variables
    INTEGER                              :: i
    INTEGER                              :: len_io
    REAL                                 :: a
    INTEGER(KIND=integer32), ALLOCATABLE :: tmp_idata32(:)
    REAL(KIND=real32), ALLOCATABLE       :: tmp_rdata32(:)
    REAL, ALLOCATABLE                    :: tmp_rdata_default(:)
    INTEGER(KIND=integer32)              :: lookup_o32(len1_lookup,len2_lookup)
    INTEGER(KIND=jpim), PARAMETER        :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER        :: zhook_out = 1
    REAL(KIND=jprb)                      :: zhook_handle
    CHARACTER (LEN=*), PARAMETER         :: routinename = 'WRITHEAD'


    IF (lhook) CALL dr_hook('WRITHEAD',zhook_in,zhook_handle)
    icode=0
    cmessage=' '

    !  1. Buffer out fixed length header record

    IF (output_bits == 32) THEN
      ALLOCATE (tmp_idata32(len_fixhd))
      tmp_idata32(:)=fixhd(:)
      CALL buffout(nftout,tmp_idata32,len_fixhd,len_io,a)
      DEALLOCATE (tmp_idata32)
    ELSE
      CALL buffout(nftout,fixhd(:),len_fixhd,len_io,a)
    END IF


    ! Check for I/O errors
    IF (a /= -1.0 .OR. len_io /= len_fixhd) THEN
      ! DEPENDS ON: ioerror
      CALL ioerror('buffer out of fixed length header',a,len_io,               &
          len_fixhd)
      cmessage='WRITHEAD: I/O error'
      icode=1
      IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
      RETURN
    END IF

    start_block=len_fixhd+1

    ! Check validity of data and print out fixed header information
    IF (printstatus >= prstatus_oper) THEN

      ! DEPENDS ON: pr_fixhd
      CALL pr_fixhd(fixhd,len_fixhd,len_inthd,len_realhd,len1_levdepc,         &
          len2_levdepc,len1_rowdepc,len2_rowdepc,len1_coldepc,len2_coldepc,    &
          len1_flddepc,len2_flddepc,len_extcnst,len_dumphist,len_cfi1,         &
          len_cfi2,len_cfi3,len1_lookup,len2_lookup,len_data,                  &
          icode,cmessage)

      IF (icode >  0) THEN
        cmessage="WRITHEAD: Error returned from PR_FIXHD"
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF
    END IF

    !  2. Buffer out integer constants

    IF (fixhd(100) >  0) THEN

      ! Check for error in file pointers
      IF (fixhd(100) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('integer constants',start_block,100,fixhd(100))
        cmessage='WRITHEAD: Addressing conflict'
        icode=2
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      IF (output_bits == 32) THEN
        ALLOCATE (tmp_idata32(fixhd(101)))
        tmp_idata32(:)=inthd(:)
        CALL buffout(nftout,tmp_idata32,fixhd(101),len_io,a)
        DEALLOCATE (tmp_idata32)
      ELSE
        CALL buffout(nftout,inthd(:),fixhd(101),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(101)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of integer constants',a,len_io,               &
            fixhd(101))
        cmessage='WRITHEAD: I/O error'
        icode=3
        IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
        RETURN
      END IF

      start_block=start_block+fixhd(101)

    END IF

    !  3. Buffer out real constants

    IF (fixhd(105) >  0) THEN

      ! Check for error in file pointers
      IF (fixhd(105) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('real constants',start_block,105,fixhd(105))
        cmessage='WRITHEAD: Addressing conflict'
        icode=4
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      ! Check for I/O errors

      IF (output_bits == 32) THEN
        ALLOCATE (tmp_rdata32(fixhd(106)))
        tmp_rdata32(:)=realhd(:)
        CALL buffout(nftout,tmp_rdata32,fixhd(106),len_io,a)
        DEALLOCATE (tmp_rdata32)
      ELSE
        CALL buffout(nftout,realhd(:),fixhd(106),len_io,a)
      END IF

      IF (a /= -1.0 .OR. len_io /= fixhd(106)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of real constants',a,len_io,                  &
            fixhd(106))
        cmessage='WRITHEAD: I/O error'
        icode=5
        IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
        RETURN
      END IF

      start_block=start_block+fixhd(106)

    END IF

    !  4. Buffer out level dependent constants

    IF (fixhd(110) >  0) THEN

      ! Check for error in file pointers
      IF (fixhd(110) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('level dependent constants',                             &
            start_block,110,fixhd(110))
        cmessage='WRITHEAD: Addressing conflict'
        icode=6
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      IF (output_bits == 32) THEN
        ALLOCATE (tmp_rdata32(fixhd(111)*fixhd(112)))
        tmp_rdata32(:)=levdepc(:)
        CALL buffout(nftout,tmp_rdata32,fixhd(111)*fixhd(112),len_io,a)
        DEALLOCATE (tmp_rdata32)
      ELSE
        CALL buffout(nftout,levdepc(:),fixhd(111)*fixhd(112),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(111)*fixhd(112)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of level dependent constants',a,              &
            len_io,fixhd(111)*fixhd(112))
        cmessage='WRITHEAD: I/O error'
        icode=7
        IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
        RETURN
      END IF

      start_block=start_block+fixhd(111)*fixhd(112)

      IF (printstatus >= prstatus_oper) THEN
        IF (mype  ==  0) THEN
          WRITE(6,'('' '')')
          WRITE(6,'('' LEVEL DEPENDENT CONSTANTS'')')
          WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(111)*fixhd(112)
        END IF ! if mype  ==  0
      END IF

    END IF

    !  5. Buffer out row dependent constants

    IF (fixhd(115) >  0) THEN

      ! Check for error in file pointers
      IF (fixhd(115) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('row dependent constants',                               &
            start_block,115,fixhd(115))
        cmessage='WRITHEAD: Addressing conflict'
        icode=8
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      IF (output_bits == 32) THEN
        ALLOCATE (tmp_rdata32(fixhd(116)*fixhd(117)))
        tmp_rdata32(:)=rowdepc(:)
        CALL buffout(nftout,tmp_rdata32,fixhd(116)*fixhd(117),len_io,a)
        DEALLOCATE (tmp_rdata32)
      ELSE
        CALL buffout(nftout,rowdepc(:),fixhd(116)*fixhd(117),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(116)*fixhd(117)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of row dependent constants',a,len_io,         &
            fixhd(116)*fixhd(117))
        cmessage='WRITHEAD: I/O error'
        icode=9
        IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
        RETURN
      END IF


      start_block=start_block+fixhd(116)*fixhd(117)

      IF (printstatus >= prstatus_oper) THEN
        IF (mype  ==  0) THEN
          WRITE(6,'('' '')')
          WRITE(6,'('' ROW DEPENDENT CONSTANTS'')')
          WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(116)*fixhd(117)
        END IF ! if mype  ==  0
      END IF

    END IF

    !  6. Buffer out column dependent constants

    IF (fixhd(120) >  0) THEN

      ! Check for error in file pointers
      IF (fixhd(120) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('column dependent constants',                            &
            start_block,120,fixhd(120))
        cmessage='WRITHEAD: Addressing conflict'
        icode=10
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      IF(output_bits == 32)THEN
        ALLOCATE (tmp_rdata32(fixhd(121)*fixhd(122)))
        tmp_rdata32(:)=coldepc(:)
        CALL buffout(nftout,tmp_rdata32,fixhd(121)*fixhd(122),len_io,a)
        DEALLOCATE (tmp_rdata32)
      ELSE
        CALL buffout(nftout,coldepc(:),fixhd(121)*fixhd(122),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(121)*fixhd(122)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of column dependent constants',a,             &
            len_io,fixhd(121)*fixhd(122))
        cmessage='WRITHEAD: I/O error'
        icode=11
        IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
        RETURN
      END IF

      start_block=start_block+fixhd(121)*fixhd(122)

      IF (printstatus >= prstatus_oper) THEN
        IF (mype  ==  0) THEN
          WRITE(6,'('' '')')
          WRITE(6,'('' COLUMN DEPENDENT CONSTANTS'')')
          WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(121)*fixhd(122)
        END IF ! if mype  ==  0
      END IF

    END IF

    !  7. Buffer out constants stored as fields

    IF (fixhd(125) >  0) THEN

      ! Check for error in file pointers
      IF (fixhd(125) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('fields of constants',                                   &
            start_block,125,fixhd(125))
        cmessage='WRITHEAD: Addressing conflict'
        icode=12
        IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
        RETURN
      END IF

      IF (output_bits == 32) THEN
        ALLOCATE (tmp_rdata32(fixhd(126)*fixhd(127)))
        tmp_rdata32(:)=flddepc(:)
        CALL buffout(nftout,tmp_rdata32,fixhd(126)*fixhd(127),len_io,a)
        DEALLOCATE (tmp_rdata32)
      ELSE
        CALL buffout(nftout,flddepc(:),fixhd(126)*fixhd(127),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(126)*fixhd(127)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of field dependent constants',a,              &
            len_io,fixhd(126)*fixhd(127))
        cmessage='WRITHEAD: I/O error'
        icode=13
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      start_block=start_block+fixhd(126)*fixhd(127)

      IF (printstatus >= prstatus_oper) THEN
        IF (mype  ==  0) THEN
          WRITE(6,'('' '')')
          WRITE(6,'('' FIELD DEPENDENT CONSTANTS'')')
          WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(126)*fixhd(127)
        END IF ! if mype  ==  0
      END IF

    END IF

    !  8. Buffer out extra constants

    IF (fixhd(130) >  0) THEN

      ! Check for error in file pointers
      IF (fixhd(130) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('extra constants',                                       &
            start_block,130,fixhd(130))
        cmessage='WRITHEAD: Addressing conflict'
        icode=14
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      IF (output_bits == 32) THEN
        ALLOCATE (tmp_rdata32(fixhd(131)))
        tmp_rdata32(:)=extcnst(:)
        CALL buffout(nftout,tmp_rdata32,fixhd(131),len_io,a)
        DEALLOCATE (tmp_rdata32)
      ELSE
        CALL buffout(nftout,extcnst(:),fixhd(131),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(131)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out extra constants',a,len_io,                    &
            fixhd(131))
        cmessage='WRITHEAD: I/O error'
        icode=15
        IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
        RETURN
      END IF

      start_block=start_block+fixhd(131)

      IF (printstatus >= prstatus_oper) THEN
        IF (mype  ==  0) THEN
          WRITE(6,'('' '')')
          WRITE(6,'('' extra constants'')')
          WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(131)
        END IF ! if mype  ==  0
      END IF

    END IF

    !  9. Buffer out temporary history block

    IF(fixhd(135) >  0)THEN

      ! Check for error in file pointers
      IF (fixhd(135) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('history',                                               &
            start_block,136,fixhd(136))
        cmessage='WRITHEAD: Addressing conflict'
        icode=16
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      IF (output_bits == 32) THEN
        ALLOCATE (tmp_rdata32(fixhd(136)))
        tmp_rdata32(:)=dumphist(:)
        CALL buffout(nftout,tmp_rdata32,fixhd(136),len_io,a)
        DEALLOCATE (tmp_rdata32)
      ELSE
        CALL buffout(nftout,dumphist(:),fixhd(136),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(136)  )THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of history file',a,len_io,                    &
            fixhd(136))
        cmessage='WRITHEAD: I/O error'
        icode=17
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      start_block=start_block+fixhd(136)

      IF (printstatus >= prstatus_oper) THEN
        IF (mype  ==  0) THEN
          WRITE(6,'('' '')')
          WRITE(6,'('' temporary history block'')')
          WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(136)
        END IF ! if mype  ==  0
      END IF

    END IF

    !  10. Buffer out compressed field index1

    IF (fixhd(140) >  0) THEN

      ! Check for error in file pointers

      IF (fixhd(140) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('compressed field index1',                               &
            start_block,140,fixhd(140))
        cmessage='WRITHEAD: Addressing conflict'
        icode=18
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      IF (output_bits == 32) THEN
        ALLOCATE (tmp_idata32(fixhd(141)))
        tmp_idata32(:)=cfi1(:)
        CALL buffout(nftout,tmp_idata32,fixhd(141),len_io,a)
        DEALLOCATE (tmp_idata32)
      ELSE
        CALL buffout(nftout,cfi1,fixhd(141),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(141)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of compressed index1',a,len_io,               &
            fixhd(141))
        cmessage='WRITHEAD: I/O error'
        icode=19
        IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
        RETURN
      END IF

      start_block=start_block+fixhd(141)

      IF (printstatus >= prstatus_oper) THEN
        IF (mype  ==  0) THEN
          WRITE(6,'('' '')')
          WRITE(6,'('' compressed field index no 1'')')
          WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(141)
        END IF ! if mype  ==  0
      END IF

    END IF

    !  11. Buffer out compressed field index2

    IF (fixhd(142) >  0) THEN

      ! Check for error in file pointers
      IF (fixhd(142) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('compressed field index2',                               &
            start_block,142,fixhd(142))
        cmessage='WRITHEAD: Addressing conflict'
        icode=20
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      IF (output_bits == 32) THEN

        ALLOCATE (tmp_idata32(fixhd(143)))
        tmp_idata32(:)=cfi2(:)
        CALL buffout(nftout,tmp_idata32,fixhd(143),len_io,a)
        DEALLOCATE (tmp_idata32)
      ELSE
        CALL buffout(nftout,cfi2,fixhd(143),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(143)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of compressed index2',a,len_io,               &
            fixhd(143))
        cmessage='WRITHEAD: I/O error'
        icode=21
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      start_block=start_block+fixhd(143)

      IF (printstatus >= prstatus_oper) THEN
        IF (mype  ==  0) THEN
          WRITE(6,'('' '')')
          WRITE(6,'('' compressed field index no 2'')')
          WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(143)
        END IF ! if mype  ==  0
      END IF

    END IF

    !  12. Buffer out compressed field index3

    IF (fixhd(144) >  0) THEN

      ! Check for error in file pointers
      IF (fixhd(144) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('compressed field index3',                               &
            start_block,144,fixhd(144))
        cmessage='WRITHEAD: Addressing conflict'
        icode=22
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      IF (output_bits == 32) THEN
        ALLOCATE (tmp_idata32(fixhd(145)))
        tmp_idata32(:)=cfi3(:)
        CALL buffout(nftout,tmp_idata32,fixhd(145),len_io,a)
        DEALLOCATE (tmp_idata32)
      ELSE
        CALL buffout(nftout,cfi3,fixhd(145),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(145)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of compressed index3',a,len_io,               &
            fixhd(145))
        cmessage='WRITHEAD: I/O error'
        icode=23
        IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
        RETURN
      END IF

      start_block=start_block+fixhd(145)

      IF (printstatus >= prstatus_oper) THEN
        IF (mype  ==  0) THEN
          WRITE(6,'('' '')')
          WRITE(6,'('' compressed field index no 3'')')
          WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(145)
        END IF ! if mype  ==  0
      END IF

    END IF

    !  13. Buffer out lookup table

    IF (fixhd(150) >  0) THEN
      IF ( getExeType() == exe_convieee ) THEN

        IF (start_block /= fixhd(150)) THEN
          IF (start_block >  fixhd(150)) THEN
            WRITE(6,'(3A,2I12)')                                               &
                '**** ERROR - Current Disk Address is greater',                &
                ' than the Address in the Fixed Length Header for the',        &
                ' Lookup Table ***** - values are', start_block-1, fixhd(150)-1
            cmessage = 'CONVIEEE: Fixed length Header Error'
            icode = 26
            CALL ereport(routinename,icode,cmessage)
          ELSE
            WRITE(6,'(3A,I10,A,I10,A)')                                        &
                '**** WARNING - Current Disk Address does not',                &
                ' match the Address in the Fixed Length Header for the',       &
                ' Lookup Table ***** - current address altered from',          &
                start_block-1,' to ',fixhd(150)-1,                             &
                ' to match the Fixed Length Header'
            start_block=fixhd(150)
            IF (output_bits == 32) THEN
              CALL setpos32(nftout, start_block-1)
            ELSE
              CALL setpos(nftout, start_block-1)
            END IF
          END IF
        END IF
      END IF
      ! Check for error in file pointers
      IF (fixhd(150) /= start_block) THEN
        ! DEPENDS ON: poserror
        CALL poserror('lookup table',start_block,150,fixhd(150))
        cmessage='WRITHEAD: Addressing conflict'
        icode=24
        CALL ereport("WRITHEAD", icode, cmessage)
      END IF

      IF (output_bits == 32) THEN

        ALLOCATE (tmp_rdata32(46:64))
        ALLOCATE (tmp_rdata_default(46:64))

        DO i=1,fixhd(152)
          lookup_o32(1:45,i)  = lookup(1:45,i)   ! Inntegers

          ! Copy real data to real storage
          CALL um_memcpy_f(tmp_rdata_default(46),lookup(46,i),19)
          ! cast to 32 bit.
          tmp_rdata32(:) = tmp_rdata_default(:)
          ! Copy 32 bit reals into integer output storage
          CALL um_memcpy32(lookup_o32(46,i),tmp_rdata32(46),19)

          IF (fixhd(5) == 6.OR.fixhd(5) == 7.OR.                               &
                                ! 6=ACOBS 7=VAROBS
              fixhd(5) == 8.OR.fixhd(5) == 10)THEN ! 8=CX   10=OBSTORE
            lookup_o32(65:128,i)=lookup(65:128,i) ! More Integers
          END IF

        END DO

        DEALLOCATE (tmp_rdata32)
        DEALLOCATE (tmp_rdata_default)

        CALL buffout(nftout,lookup_o32,fixhd(151)*fixhd(152),len_io,a)
      ELSE
        CALL buffout(nftout,lookup,fixhd(151)*fixhd(152),len_io,a)
      END IF

      ! Check for I/O errors
      IF (a /= -1.0 .OR. len_io /= fixhd(151)*fixhd(152)) THEN
        ! DEPENDS ON: ioerror
        CALL ioerror('buffer out of lookup table',a,len_io,                    &
            fixhd(151)*fixhd(152))
        cmessage='WRITHEAD: I/O error'
        icode=25
        IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)
        RETURN
      END IF

      start_block=start_block+fixhd(151)*fixhd(152)

      IF (printstatus >= prstatus_oper) THEN
        IF (mype  ==  0) THEN
          WRITE(6,'('' '')')
          WRITE(6,'('' LOOKUP TABLE'')')
          WRITE(6,'('' '',i8,'' 64-bit words long'')')fixhd(151)*fixhd(152)
        END IF ! if mype  ==  0
      END IF

      ! No consistency checks for parallel code. The LOOKUP headers don't
      ! match the data layout in memory within a PE.
    END IF

    IF (lhook) CALL dr_hook('WRITHEAD',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE writhead

END MODULE writhead_mod
