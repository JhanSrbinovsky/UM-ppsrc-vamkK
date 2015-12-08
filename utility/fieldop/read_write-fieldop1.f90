! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
SUBROUTINE read_write_fieldop1(i_dm,                                           &
    pp_unit1,                                                                  &
    pp_unit2,                                                                  &
    len1_lookup,                                                               &
    len2_lookup,                                                               &
    len_fixhd,                                                                 &
    fixhd,                                                                     &
    lookup,                                                                    &
    fixhd2,                                                                    &
    lookup2,                                                                   &
    rookup2,                                                                   &
    len1_lookup2,                                                              &
    len2_lookup2,                                                              &
    op,                                                                        &
    rookup,                                                                    &
    entry_no,                                                                  &
    entry_no2,                                                                 &
    data_add2,                                                                 &
    data_add1,                                                                 &
    model_flag,                                                                &
    nfields,                                                                   &
    tfields,                                                                   &
    llev,                                                                      &
    ignore,                                                                    &
    pp_unit_out,                                                               &
    max_len,                                                                   &
    divisor,                                                                   &
    amdi,                                                                      &
    l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,                                    &
    l13,l14,l15,l16,l17,l18,l19,l20,                                           &
    stash1,stash2,stash3,stash4,stash5,                                        &
    stash6,stash7,stash8,stash9,stash10,                                       &
    stash11,stash12,stash13,stash14,stash15,                                   &
    stash16,stash17,stash18,stash19,stash20,                                   &
    icode,cmessage)

  USE ereport_mod, ONLY : ereport
  USE lookup_addresses

  IMPLICIT NONE

! Description:
!     Accesses the data (packed or unpacked) from one or two model
! dumps or direct access fieldsfiles and write out to a new file
! the difference, sum or product of the data values. (if a single
! datafile is read the data is divided by an integer). The output file
! is a copy of the first input file with the fields overwritten by the
! differenced/meaned etc data.

! Subroutine arguments
!   Scalar arguments with intent(in):
  INTEGER                                                                      &
      pp_unit1,                                                                &
                              ! unit no of required fieldsfile/dump
      pp_unit2,                                                                &
                              ! unit no of required fieldsfile/dump
      pp_unit_out,                                                             &
                              ! unit no of output file
      len_fixhd,                                                               &
      i_dm,                                                                    &
                              ! num_values rounded to an even no
                              ! used to dimension the output array
      data_add1,                                                               &
                              ! The word address of the data.
      data_add2,                                                               &
                              ! The word address of the data.
      len1_lookup,                                                             &
                              ! First dimension of the lookup
      len1_lookup2,                                                            &
                              ! First dimension of the lookup
      len2_lookup,                                                             &
                           ! Size of the lookup on the file
      len2_lookup2,                                                            &
                           ! Size of the lookup on the file
      max_len,                                                                 &
      comp_accry1,                                                             &
      comp_accry2,                                                             &
      entry_no,                                                                &
                              ! Lookup entry no of the Field.
      entry_no2,                                                               &
                              ! Lookup entry no of the Field.
      stash1,stash2,stash3,stash4,stash5,                                      &
      stash6,stash7,stash8,stash9,stash10,                                     &
      stash11,stash12,stash13,stash14,stash15,                                 &
      stash16,stash17,stash18,stash19,stash20,                                 &
      l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,                                  &
      l13,l14,l15,l16,l17,l18,l19,l20,                                         &
      divisor,                                                                 &
      exppxi

  REAL                                                                         &
      rookup(len1_lookup,len2_lookup),                                         &
                                            ! Real lookup
      rookup2(len1_lookup2,len2_lookup2),                                      &
                                            ! Real lookup
      amdi

  LOGICAL                                                                      &
      model_flag             !IN True => dumps, False => fieldsfile

!   Array  arguments with intent(in):
  INTEGER                                                                      &
      fixhd(len_fixhd),                                                        &
                                           ! fixed header (file1)
      fixhd2(len_fixhd),                                                       &
                                           ! fixed header (file2)
      lookup(len1_lookup,len2_lookup),                                         &
                                           ! integer lookup (file1)
      lookup2(len1_lookup2,len2_lookup2)  ! integer lookup (file2)

  CHARACTER                                                                    &
      op*(8)

  CHARACTER                                                                    &
      exppxc*(36)

!   ErrorStatus
  INTEGER                                                                      &
      icode

  CHARACTER                                                                    &
      cmessage*80

! Local parameters:
  INTEGER max_len_ilabel  ! max length of INT part of pp header
  PARAMETER(max_len_ilabel=45)

  INTEGER max_len_rlabel  ! max length of REAL part of pp header
  PARAMETER(max_len_rlabel=32)

! Local scalars:
  INTEGER                                                                      &
      i,                                                                       &
                       ! local counter
      iwa,                                                                     &
                       ! Word address in call setpos (file1).
      iwa2,                                                                    &
                       ! Word address in call setpos (file2).
      n_rows_out,                                                              &
                       ! No of rows of data in field.
      n_cols_out,                                                              &
                       ! No. of columns of data in field.
      len_ilabel,                                                              &
                       ! number of values in ilabel
      len_rlabel      ! number of values in rlabel

  LOGICAL                                                                      &
      packed,                                                                  &
                   ! indicates whether the data is packed
      cont,                                                                    &
                   ! indicates whether to operate on field
      nfields,                                                                 &
      tfields,                                                                 &
      llev,                                                                    &
      ignore,                                                                  &
      lop

! Local dynamic arrays:
  INTEGER                                                                      &
      ilabel(max_len_ilabel),                                                  &
                                    ! holds integer part of lookup
      ilabel2(max_len_ilabel)      ! holds integer part of lookup

  REAL                                                                         &
      field(max_len),                                                          &
                                   ! array holding data
      field0(max_len),                                                         &
                                   ! array holding data
      rlabel(max_len_rlabel),                                                  &
                                ! holds real part of lookup
      rlabel2(max_len_rlabel)  ! holds real part of lookup

!- End of header
  comp_accry1=0
  comp_accry2=0
      ! Fields with stashcodes are written directly to output
      ! fieldsfile or dump.
  IF ((((lookup(42,entry_no) == stash1).OR.                                    &
      (lookup(42,entry_no) == stash2)                                          &
      .OR.(lookup(42,entry_no) == stash3).OR.                                  &
      (lookup(42,entry_no) == stash4)                                          &
      .OR.(lookup(42,entry_no) == stash5).OR.                                  &
      (lookup(42,entry_no) == stash6)                                          &
      .OR.(lookup(42,entry_no) == stash7).OR.                                  &
      (lookup(42,entry_no) == stash8)                                          &
      .OR.(lookup(42,entry_no) == stash9).OR.                                  &
      (lookup(42,entry_no) == stash10)                                         &
      .OR.(lookup(42,entry_no) == stash11).OR.                                 &
      (lookup(42,entry_no) == stash12)                                         &
      .OR.(lookup(42,entry_no) == stash13).OR.                                 &
      (lookup(42,entry_no) == stash14)                                         &
      .OR.(lookup(42,entry_no) == stash15).OR.                                 &
      (lookup(42,entry_no) == stash16)                                         &
      .OR.(lookup(42,entry_no) == stash17).OR.                                 &
      (lookup(42,entry_no) == stash18)                                         &
      .OR.(lookup(42,entry_no) == stash19).OR.                                 &
      (lookup(42,entry_no) == stash20))                                        &
      .AND.(nfields))                                                          &
      .OR.                                                                     &
      (.NOT.((lookup(42,entry_no) == stash1).OR.                               &
      (lookup(42,entry_no) == stash2)                                          &
      .OR.(lookup(42,entry_no) == stash3).OR.                                  &
      (lookup(42,entry_no) == stash4)                                          &
      .OR.(lookup(42,entry_no) == stash5).OR.                                  &
      (lookup(42,entry_no) == stash6)                                          &
      .OR.(lookup(42,entry_no) == stash7).OR.                                  &
      (lookup(42,entry_no) == stash8)                                          &
      .OR.(lookup(42,entry_no) == stash9).OR.                                  &
      (lookup(42,entry_no) == stash10)                                         &
      .OR.(lookup(42,entry_no) == stash11).OR.                                 &
      (lookup(42,entry_no) == stash12)                                         &
      .OR.(lookup(42,entry_no) == stash13).OR.                                 &
      (lookup(42,entry_no) == stash14)                                         &
      .OR.(lookup(42,entry_no) == stash15).OR.                                 &
      (lookup(42,entry_no) == stash16)                                         &
      .OR.(lookup(42,entry_no) == stash17).OR.                                 &
      (lookup(42,entry_no) == stash18)                                         &
      .OR.(lookup(42,entry_no) == stash19).OR.                                 &
      (lookup(42,entry_no) == stash20))                                        &
      .AND.(tfields))                                                          &
      .OR.(ignore))                                                            &
      THEN

    WRITE(*,*)'FIELD NO.',entry_no,'DIRECTLY TRANSFERED'
    lop=.FALSE.

  ELSE IF(fixhd(5) == 3.AND.lookup(42,entry_no) == 30)THEN
    WRITE(*,*)'FIELD NO.',entry_no,                                            &
        'LAND-SEA MASK: DIRECTLY TRANSFERED'
    lop=.FALSE.
  ELSE

    IF (((llev .AND.                                                           &
        ((lookup(33,entry_no) == l1).OR.                                       &
        (lookup(33,entry_no) == l2)                                            &
        .OR.(lookup(33,entry_no) == l3).OR.                                    &
        (lookup(33,entry_no) == l4)                                            &
        .OR.(lookup(33,entry_no) == l5).OR.                                    &
        (lookup(33,entry_no) == l6)                                            &
        .OR.(lookup(33,entry_no) == l7).OR.                                    &
        (lookup(33,entry_no) == l8)                                            &
        .OR.(lookup(33,entry_no) == l9).OR.                                    &
        (lookup(33,entry_no) == l10)                                           &
        .OR.(lookup(33,entry_no) == l11).OR.                                   &
        (lookup(33,entry_no) == l12)                                           &
        .OR.(lookup(33,entry_no) == l13).OR.                                   &
        (lookup(33,entry_no) == l14)                                           &
        .OR.(lookup(33,entry_no) == l15).OR.                                   &
        (lookup(33,entry_no) == l16)                                           &
        .OR.(lookup(33,entry_no) == l17).OR.                                   &
        (lookup(33,entry_no) == l18)                                           &
        .OR.(lookup(33,entry_no) == l19).OR.                                   &
        (lookup(33,entry_no) == l20))).AND.                                    &
        (lookup(33,entry_no) /= 0)) .OR. (.NOT.llev))THEN

      WRITE(*,*)'FIELD NO.',entry_no,'OPERATED ON'
      lop=.TRUE.


    ELSE
      WRITE(*,*)'FIELD NO.',entry_no,'DIRECTLY TRANSFERED'
      lop=.FALSE.

    END IF

  END IF

      ! Fieldsfiles contain packed data which when expanded,
      ! operated upon and repacked occupy a different amount of
      ! space. It is therefore necessary to write out data in a
      ! fieldsfile to the new addresses even although only some
      ! of the fields may be changed.
  IF ((lop) .OR. (fixhd(5) == 3)) THEN

    DO i=1,i_dm           ! field is initialised.
      field(i) =0.0
    END DO

    packed=.FALSE.

        ! Access the 1st fieldsfile/dump.
! DEPENDS ON: readff_fieldop1
    CALL readff_fieldop1(pp_unit1,                                             &
                                     ! IN Unit no.
        field,                                                                 &
                                     ! OUT Data field.
        i_dm,                                                                  &
                                     ! IN Size of data field (rounded).
        entry_no,                                                              &
                                     ! IN position of field in lookup.
        ilabel,                                                                &
                                     ! OUT Integer part of lookup.
        rlabel,                                                                &
                                     ! OUT Real part of lookup.
        len2_lookup,                                                           &
                                     ! IN
        len1_lookup,                                                           &
                                     ! IN
        len_fixhd,                                                             &
        fixhd,                                                                 &
                                     ! IN Fixed header
        lookup,                                                                &
                                     ! IN Integer part of lookup.
        rookup,                                                                &
                                     ! IN Real part of lookup.
        data_add1,                                                             &
                                     ! IN Start address of data.
        model_flag,                                                            &
                                     ! IN TRUE -dump, FALSE -fieldsfile
        max_len_ilabel,                                                        &
                                     ! IN
        max_len_rlabel,                                                        &
                                     ! IN
        max_len,                                                               &
        comp_accry1,                                                           &
                                     ! OUT Packing accuracy of field
        len_ilabel,                                                            &
                                     ! IN
        len_rlabel,                                                            &
                                     ! IN
        iwa,                                                                   &
                                     ! OUT Word address of data field
                                     !     in call to setpos.
        icode,cmessage)  ! IN

    IF (lop) THEN ! Arithmetic operation performed on field.

      IF (op  /=  'idiv    ') THEN

        DO i=1,i_dm

          field0(i) =field(i) ! Write data from file1 to data0.
          field(i)  =0.0

        END DO ! i

        ! Access the 2nd fieldsfile/dump.
! DEPENDS ON: readff_fieldop1
        CALL readff_fieldop1(pp_unit2,                                         &
                                       ! IN Unit no.
            field,                                                             &
                                       ! OUT Data field corresponding
                                       !     to field accessed in file1.
            i_dm,                                                              &
                                       ! IN Size of data field (rounded)
            entry_no2,                                                         &
                                       ! IN position of field in lookup2
            ilabel2,                                                           &
                                       ! OUT Integer part of lookup.
            rlabel2,                                                           &
                                       ! OUT Real part of lookup.
            len2_lookup2,                                                      &
                                       ! IN
            len1_lookup2,                                                      &
                                       ! IN
            len_fixhd,                                                         &
            fixhd2,                                                            &
                                       ! IN Fixed header
            lookup2,                                                           &
                                       ! IN Integer part of lookup.
            rookup2,                                                           &
                                       ! IN Real part of lookup.
            data_add2,                                                         &
                                       ! IN Start address of data.
            model_flag,                                                        &
                                       ! IN TRUE dump, FALSE fieldsfile
            max_len_ilabel,                                                    &
                                       ! IN
            max_len_rlabel,                                                    &
                                       ! IN
            max_len,                                                           &
                                       ! IN Max no of points of a field
            comp_accry2,                                                       &
                                       ! OUT Packing accuracy of field
            len_ilabel,                                                        &
                                       ! IN
            len_rlabel,                                                        &
                                       ! IN
            iwa2,                                                              &
                                       ! OUT Word address of data field
                                       !     in call to setpos.
            icode,cmessage)  ! IN

          ! The data has now been read in and has 1) read in as packed
          ! and then un-packed or 2) The data was never packed at all.
          ! If packed field will have lblrec/2 values if a DUMP and
          ! LBLREC values if a pp_file. If the data is not packed field
          ! will have the no of data points length lbrow*lbnpt+lbext if
          ! a pp_file and lblrec if a dump file.

          ! For a dump lblrec will hold origonal no of data points.
          ! For a pp_file lblrec will hold the no of CRAY words needed
          ! to hold the data (if un-packed also no of data points)

          ! Difference, sum or multiply the data in fields.
        IF (op  ==  'subtract') THEN
          WRITE(*,*)'subtract',entry_no
          DO i = 1,i_dm
            IF (lookup(data_type,entry_no) == 3) THEN
              field(i) = field(i)
            ELSE
              IF (field(i) /= amdi) THEN
                field(i) = field0(i) - field(i)
              END IF
            END IF
          END DO

        ELSE IF (op  ==  'add     ') THEN
          DO i = 1,i_dm
            IF (lookup(data_type,entry_no) == 3) THEN
              field(i) = field(i)
            ELSE
              IF (field(i) /= amdi) THEN
                field(i) = field0(i) + field(i)
              END IF
            END IF
          END DO

        ELSE IF (op  ==  'multiply') THEN
          DO i = 1,i_dm
            IF (lookup(data_type,entry_no) == 3) THEN
              field(i) = field(i)
            ELSE
              IF (field(i) /= amdi) THEN
                field(i) = field0(i) * field(i)
              END IF
            END IF
          END DO

        ELSE

          WRITE(6,*)'Not a valid operation'

          CALL ereport('READ_WRITE', 1000,                                     &
              cmessage)


        END IF
      ELSE

        ! Divide data in field from a single input file by an integer.
        DO i = 1,i_dm
          IF (lookup(data_type,entry_no) == 3) THEN
            field(i) = field(i)
          ELSE
            IF (field(i) /= amdi) THEN
              field(i) = field(i) / divisor
            END IF
          END IF
        END DO

      END IF

    END IF

    IF(icode /= 0) RETURN

    n_rows_out=lookup(18,entry_no)
    n_cols_out=lookup(19,entry_no)

        ! Write data to output dump/fieldsfile. Data written to original
        ! positions in 1st file.
! DEPENDS ON: writeff
    CALL writeff(pp_unit_out,                                                  &
                                     ! IN Unit no of output file.
        field,                                                                 &
                                     ! IN Output data after arith op.
        i_dm,                                                                  &
                                     ! IN Size of data field (rounded).
        entry_no,                                                              &
                                     ! IN pos. of field in lookup table.
        data_add1,                                                             &
        lookup,                                                                &
                                     ! IN
        len_fixhd,                                                             &
        fixhd,                                                                 &
                                   ! IN
        len2_lookup,                                                           &
                                   ! IN
        len1_lookup,                                                           &
                                     ! IN
        n_rows_out,                                                            &
                                     ! IN
        n_cols_out,                                                            &
                                     ! IN
        packed,                                                                &
                                     ! IN FALSE - unpacked data
        max_len,                                                               &
                                   ! IN Max no of points of a field in f
        comp_accry1,                                                           &
                                   ! IN accuracy at which field packed
        op,                                                                    &
                                     ! IN Operation type.
        icode,cmessage) ! IN


  END IF

  RETURN
END SUBROUTINE read_write_fieldop1

