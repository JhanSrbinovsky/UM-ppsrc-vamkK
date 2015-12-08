! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine interface:
SUBROUTINE readff_fieldop1(pp_unit1,                                           &
    field,                                                                     &
    i_dm,                                                                      &
    entry_no,                                                                  &
    ilabel,                                                                    &
    rlabel,                                                                    &
    pp_len2_lookup,                                                            &
    len1_lookup,                                                               &
    len_fixhd,                                                                 &
    pp_fixhd,                                                                  &
    lookup,                                                                    &
    rookup,                                                                    &
    data_add1,                                                                 &
    model_flag,                                                                &
    max_len_ilabel,                                                            &
    max_len_rlabel,                                                            &
    max_len,                                                                   &
    pppak,                                                                     &
    len_ilabel,                                                                &
    len_rlabel,                                                                &
    iwa,                                                                       &
    icode,cmessage)
  USE lookup_addresses

  IMPLICIT NONE
!
!
! Description: To read a direct access PP file.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <appropriate code>
! System Task:              <appropriate code>
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! 1.0 Global variables (*CALLed COMDECKs etc...):
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

! Subroutine arguments
!   Scalar arguments with intent(in):
  INTEGER                                                                      &
      len1_lookup,                                                             &
                              ! first dimension of the lookup
      pp_len2_lookup,                                                          &
                              ! secnd dimension of the lookup
      pp_unit1,                                                                &
                              ! unit no of required fieldsfile
      i_dm,                                                                    &
                              ! Size of data field (rounded)
      max_len_rlabel,                                                          &
                              ! max sixe of rlabel
      max_len_ilabel,                                                          &
                              ! max sixe of ilabel
      max_len,                                                                 &
      data_add1,                                                               &
                              ! The word address of the data.
      entry_no,                                                                &
                              ! Lookup entry no of the Field.
      len_fixhd,                                                               &
      lookup(len1_lookup,pp_len2_lookup) ! integer lookup

  REAL                                                                         &
      rookup(len1_lookup,pp_len2_lookup) ! real lookup

  LOGICAL                                                                      &
      model_flag             ! True => Dump False =>Fieldsfile

!   Array  arguments with intent(in):
  INTEGER                                                                      &
      pp_fixhd(len_fixhd)    ! fixed header

!   Scalar arguments with intent(out):
  INTEGER                                                                      &
      len_rlabel,                                                              &
                              ! actual size of rlabel
      len_ilabel,                                                              &
                              ! actual size of ilabel
      pppak

!   Array  arguments with intent(out):
  INTEGER                                                                      &
      ilabel(max_len_ilabel) ! integer part of lookup

  REAL                                                                         &
      field(i_dm),                                                             &
                           ! array holding final output data.
      rlabel(max_len_rlabel) ! real part of lookup

!   ErrorStatus
  INTEGER                                                                      &
      icode                  ! error code

  CHARACTER                                                                    &
      cmessage*80          ! error message

! Local scalars:
  INTEGER                                                                      &
      i,j,                                                                     &
                              ! Local counters
      pack_type,                                                               &
                              ! packing type N1 of LBPACK
      num_cray_words,                                                          &
                              ! number of words for field
      nvals,                                                                   &
                              ! number of points in a data field
      iwa,                                                                     &
                              ! Word address in call setpos
      length_of_data,                                                          &
                              ! Length of a particular field
      addr,                                                                    &
                              ! Address of a field in the data store
      pos_rlabel,                                                              &
                              ! position of first REAL in PPhdr
      pack_type_i            ! packing type N1 of LBPACK

  REAL                                                                         &
      amdi                   ! Missing data indicator for lookup

!- End of header

  amdi=rookup(bmdi,entry_no)
  IF (amdi /= rmdi) WRITE(*,*)' NON STANDARD MISSING DATA USED'

  pack_type = MOD(lookup(lbpack,entry_no),10)

      ! Reading a model type dump
      ! A model dump has no direct addressing only relative.

  IF(model_flag) THEN

! Old Format dumpfiles
    IF((lookup(lbnrec,entry_no) == 0) .OR.                                     &
! Prog lookups in dump before vn3.2:
        ((lookup(lbnrec,entry_no) == imdi) .AND.                               &
        (pp_fixhd(12) <= 301))) THEN

      IF(pack_type == 2) THEN            ! 32 bit packing.

        num_cray_words = (lookup(lblrec,entry_no)+1)/2
      ELSE IF (pack_type >  0) THEN

        num_cray_words = lookup(lblrec,entry_no)/2
      ELSE

        num_cray_words = lookup(lblrec,entry_no)
      END IF

      nvals = lookup(lblrec,entry_no) ! No of data points
      addr=data_add1

      IF (entry_no >  1) THEN

        DO i =1,entry_no-1

          pack_type_i = MOD(lookup(lbpack,i),10)
          IF (pack_type_i  ==  2) THEN ! 32 Bit packed

            length_of_data = (lookup(lblrec,i)+1)/2
          ELSE

            length_of_data = lookup(lblrec,i)
          END IF

          addr = addr + length_of_data

        END DO ! i
      ELSE       !  If the first entry.

        addr = data_add1
        IF (pack_type  ==  2) THEN ! 32 Bit packed

          length_of_data = (lookup(lblrec,1)+1)/2
        ELSE

          length_of_data=lookup(lblrec,1)
        END IF

        WRITE(*,*)'  length_of_data  ',length_of_data

      END IF

      iwa=addr  ! Not -1 as this is already done in dump

    ELSE
! New format Dumpfiles (vn4.4 onwards)

      IF(pack_type == 2) THEN            ! 32 bit packing.
        num_cray_words=(lookup(lblrec,entry_no)+1)/2
      ELSE IF(pack_type >  0) THEN
        num_cray_words=lookup(lblrec,entry_no)/2
      ELSE
        num_cray_words=lookup(lblrec,entry_no)
      END IF
      iwa = lookup(lbegin,entry_no)
      nvals = lookup(lbrow,entry_no) * lookup(lbnpt,entry_no)
    END IF
  ELSE ! Reading a PP type file.

    num_cray_words = lookup(lblrec,entry_no) ! PP type file
    iwa = lookup(lbegin,entry_no)
    nvals = lookup(lbrow,entry_no) * lookup(lbnpt,entry_no)                    &
        + lookup(lbext,entry_no)

  END IF

  107 FORMAT(' ENTRY NO=',i5,'num_cray_words= ',i6,'nvals=',i6)

  IF (i_dm  <   num_cray_words) THEN

    icode = num_cray_words
    cmessage ='readff  Idim to small icode holds correct value'
    GO TO 9999

  END IF

  icode=0
! DEPENDS ON: read_rec_fieldop1
  CALL read_rec_fieldop1(field,                                                &
                                     ! OUT array holding data
      num_cray_words,                                                          &
                                     ! IN No of CRAY words holding data
      iwa,                                                                     &
                                  ! IN WORD address of field to be read
      pp_unit1,                                                                &
                                  ! IN unit no of the file
      max_len,                                                                 &
      icode,                                                                   &
                                  ! IN/OUT
      pack_type)

  2212 FORMAT('  FIELDS FILE NUMBER ',i2,'  ON UNIT',i2,2x,'BEING read')

  IF (icode == 0) THEN

    pos_rlabel = MOD(lookup(lbrel,entry_no),100)

          ! Treat lookup(45) (submodel identifier) as an integer.
    pos_rlabel=46


    len_rlabel=1+len1_lookup-pos_rlabel
    len_ilabel=len1_lookup-len_rlabel

    DO i=1,len_ilabel
      ilabel(i)=lookup(i,entry_no)
    END DO

!         check for valid release number
    IF (ilabel(lbrel) <  1) THEN

      WRITE(*,*)' resetting LBREL from',ilabel(lbrel),' to 3'
      ilabel(lbrel)=3

    END IF

    DO i=1,len_rlabel
      rlabel(i)=rookup(i+pos_rlabel-1,entry_no)
    END DO

  END IF

        ! At this point field holds the data either packed or un-packed
        ! Is the packing indicator set and is un-packing required?
        ! If so then the data is temp un-packed into a work array of
        ! length idim
  IF (pack_type >  0) THEN       ! Is the field packed.

! DEPENDS ON: un_pack_fieldop1
    CALL un_pack_fieldop1(pack_type,                                           &
                                       ! IN packing type N1 of LBPACK
        i_dm,                                                                  &
                                       ! IN length of unpacked pp buffer
        field,                                                                 &
                                       ! IN/OUT I/P contains packed data
                                       ! Output contains un-packed data.
        num_cray_words,                                                        &
                                         ! IN length of input field
        ilabel,                                                                &
                                       ! IN holds integer part of lookup
        len_ilabel,                                                            &
                                       ! IN length of ilabel array
        amdi,                                                                  &
                                       ! IN Missing data indicator.
        pp_fixhd,                                                              &
                                       ! IN PPfile fixed length header
        len_fixhd,                                                             &
        pppak,                                                                 &
        icode,cmessage)  ! IN/OUT

  ELSE IF(lookup(data_type,entry_no) == 2) THEN !Fld is integer

! DEPENDS ON: integer_to_real_fieldop1
    CALL integer_to_real_fieldop1(i_dm,                                        &
                                       ! IN full unpacked size of field
        field,                                                                 &
                                       ! IN contains integer data.
        field,                                                                 &
                                       ! OUT contains Real data.
        nvals,                                                                 &
                                       ! IN no of values in field
        max_len,                                                               &
        ilabel,                                                                &
                                       ! IN/OUT integer part of lookup
        icode)  ! IN/OUT error code

  END IF

  9999 CONTINUE
  100 FORMAT(//,32x,'   ARRAY        ',//,32(16f5.0/))
  101 FORMAT(//,32x,'   lookup       ',//,32(16i5/))
  103 FORMAT('   LENIN  ',i12)

  RETURN
END SUBROUTINE readff_fieldop1
