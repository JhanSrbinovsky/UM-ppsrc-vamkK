! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
SUBROUTINE writeff(pp_unit_out,                                                &
    field,                                                                     &
    i_dm,                                                                      &
    entry_no,                                                                  &
    data_add,                                                                  &
    lookup,                                                                    &
    len_fixhd,                                                                 &
    fixhd,                                                                     &
    len2_lookup,                                                               &
    len1_lookup,                                                               &
    n_rows_out,                                                                &
    n_cols_out,                                                                &
    packed,                                                                    &
    max_len,                                                                   &
    comp_accry,                                                                &
    op,                                                                        &
    icode,cmessage)
  USE io
  USE io_configuration_mod, ONLY: io_field_padding
  USE ereport_mod, ONLY : ereport
  USE Atmos_Max_Sizes
  USE UM_ParParams
  USE Decomp_DB, ONLY : decompose
  USE lookup_addresses

  IMPLICIT NONE
!
!
! Description: To ouput a field to a UM dump or fieldsfile, with the
!              data written in packed (wgdos,grib or cray 32 bits) or
!              unpacked form.
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
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
! Description: This include file contains information needed when
!              generating variable horizontal grid data in the
!              STASH extra data vector. Introduced UM 5.4 - R. Hill
!===================================================================
      LOGICAL :: X_VAR_GRID ! Whether variable grid in E-W direction
      LOGICAL :: Y_VAR_GRID ! and/or in S-N direction

      INTEGER :: VAR_GRID_TYPE ! 0 = none
                               ! 1 = T grid
                               ! 2 = U/V grid

      ! Grid boundaries for T and U,V
      REAL :: X_BOUNDARY(ROW_LENGTH_MAX+1,2)
      REAL :: Y_BOUNDARY(ROWS_MAX+1,2)

      ! Grid Points for T and U,V
      REAL :: X_GRID(ROW_LENGTH_MAX,2)
      REAL :: Y_GRID(ROWS_MAX,2)

      COMMON /OVARGRID/ X_VAR_GRID,Y_VAR_GRID                           &
     & ,X_BOUNDARY,Y_BOUNDARY,X_GRID,Y_GRID,VAR_GRID_TYPE

      ! The following parameters correspond to the extra data
      ! vector descriptors expected, for e.g., in PV-WAVE
      ! plotting routines (e.g. decode_extra.pro). There are
      ! numerous other areas of code where these integer
      ! descriptors must be handled (e.g. FIELDCOS, PPI2H, FTT)
      ! So it is not a trivial matter to introduce new code
      ! descriptors. Furthermore, ieee -32 will destroy these
      ! integers so PP data must always be processed via the
      ! long winded route: QXFIELDCOS -> PUTIBM -> FTT/PPI2H.
      ! (This thoroughly unsatisfactory state of affairs may
      ! be correctable with developments to ieee and convpp).
      INTEGER,PARAMETER :: x_coord_vector=1
                                     ! Indicates that an extra
                                     ! data vector gives LBNPT
                                     ! x-coordinate values
      INTEGER,PARAMETER :: y_coord_vector=2
                                     ! Indicates that an extra
                                     ! data vector gives LBROW
                                     ! Y-coordinate values
      INTEGER,PARAMETER :: x_lbnd_vector=12
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! x-boundary values
      INTEGER,PARAMETER :: x_ubnd_vector=13
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! x-boundary values
      INTEGER,PARAMETER :: y_lbnd_vector=14
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! y-boundary values
      INTEGER,PARAMETER :: y_ubnd_vector=15
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! y-boundary values

! Subroutine arguments
!   Scalar arguments with intent(in):
  INTEGER                                                                      &
      n_rows_out,                                                              &
      n_cols_out,                                                              &
      len_fixhd,                                                               &
      pp_unit_out,                                                             &
      len2_lookup,                                                             &
      len1_lookup,                                                             &
      entry_no,                                                                &
      grib_packing,                                                            &
      max_len,                                                                 &
      comp_accry,                                                              &
      i_dm,                                                                    &
      data_add,                                                                &
      exppxi

  CHARACTER                                                                    &
      op*(8)

  CHARACTER                                                                    &
      exppxc*(36)

  LOGICAL                                                                      &
      packed

!   Array  arguments with intent(in):
  INTEGER                                                                      &
      lookup(len1_lookup,len2_lookup),                                         &
      fixhd(len_fixhd),                                                        &
      ifield(max_len)

  REAL                                                                         &
      field(max_len)

!   ErrorStatus
  INTEGER                                                                      &
      icode

  REAL                                                                         &
      a

  CHARACTER                                                                    &
      cmessage*80

! Local scalars:
  INTEGER                                                                      &
      comp_accrcy,                                                             &
                          !
      num_words,                                                               &
                          !
      pack_type,                                                               &
                          ! Packing type N1 of LBPACK
      len_io,                                                                  &
                          !
      isize,                                                                   &
      tot_levels,                                                              &
                          ! total number of levels
      i                  !

!- End of header

  icode       = 0
  num_words   = -99
  pack_type   = MOD(lookup(lbpack,entry_no),10)
  IF (pack_type >  0) packed =.TRUE.

      ! Method of GRIB packing - use width method, with simple packing
      ! to be similar to the ECMWF MARS archive.
  grib_packing=6

  IF((lookup(44,entry_no) <  0) .OR.                                           &
      (lookup(44,entry_no) >  100)) THEN

    lookup(44,entry_no) = 0
    lookup(44,entry_no) = lookup(44,entry_no) + 1
  ELSE

    lookup(44,entry_no) = lookup(44,entry_no) + 1
  END IF

  CALL setpos(pp_unit_out,fixhd(150)+(entry_no-1)                              &
      *len1_lookup-1,icode)
  CALL buffout(pp_unit_out,lookup(1:,entry_no),fixhd(151),                     &
      len_io,a)

      ! Check for I/O errors
  IF (a  /=  -1.0 .OR. len_io  /=  fixhd(151)) THEN

! DEPENDS ON: ioerror
    CALL ioerror('buffer out of lookup table',a,len_io,                        &
        fixhd(151))
    cmessage='FIELDOP: I/O error'
    icode=25
    RETURN
  END IF

  IF (pack_type  == 1 .OR.pack_type  ==  4 )THEN
    isize=(((((i_dm+io_field_padding-1)/io_field_padding)*                     &
        io_field_padding) +1)/2) * 2
      ! Deactivate varible grid extra data processing.
    var_grid_type = 0

       ! Data packed using WGDOS method and written to O/P file.
! DEPENDS ON: fldop_pp_file
    CALL fldop_pp_file(field,                                                  &
                                    ! IN Array to store expanded data
        isize,                                                                 &
                                    ! IN length of pp buffer
        num_words,                                                             &
                                    ! IN No of 64bit words of data
        rmdi,                                                                  &
                                    ! IN Missing data
        comp_accry,                                                            &
                                   ! IN PPXREF accuracy code.
        i_dm,                                                                  &
                                   ! IN length of pp buffer
        pp_unit_out,                                                           &
                                    ! IN Unit no of O/P field.
        data_add,                                                              &
                                   ! IN Word address of data (file1)
        n_cols_out,                                                            &
                                    ! IN
        n_rows_out,                                                            &
                                   ! IN
        packed,                                                                &
                                    ! IN TRUE - packing required.
        pack_type,                                                             &
                                    ! IN WGDOS packed data.
        lookup,                                                                &
                                   ! IN lookup headers of file1.
        len1_lookup,                                                           &
                                   ! IN
        len2_lookup,                                                           &
                                   ! IN
        entry_no,                                                              &
                                    ! IN
        icode,cmessage)

    IF (icode >  0) THEN
      cmessage='FIELDOP : Error in FLDOP_PP_FILE'

      CALL ereport('WRITEFF', icode,                                           &
          cmessage)

    END IF

  ELSE IF (pack_type == 3) THEN

       ! Data compressed using the GRIB method, written to O/P file.
! DEPENDS ON: grib_file
    CALL grib_file(len1_lookup,                                                &
                                       ! IN
        len2_lookup,                                                           &
                                       ! IN
        lookup,                                                                &
                                       ! IN
        lookup,                                                                &
                                       ! IN
        entry_no,                                                              &
                                       ! IN Posn of field in lookup.
        field,                                                                 &
                                       ! IN Unpacked output data.
        max_len,                                                               &
                                       ! IN Length of pp buffer
        max_len,                                                               &
                                       ! IN
        num_words,                                                             &
                                       ! IN No of 64bit words of data
        pp_unit_out,                                                           &
                                       ! IN Unit no of O/P field.
        pp_unit_out,                                                           &
                                       ! IN Word address of record.
        grib_packing,                                                          &
                                       !
        icode,cmessage) ! IN

    IF (icode >  0) THEN
      cmessage='FIELDOP : Error in GRIB_FILE'

      CALL ereport('WRITEFF', icode,                                           &
          cmessage)

    END IF

  ELSE IF ((pack_type == 0).OR.(pack_type == 2)) THEN
        ! Update lookup header data lengths and addressing for
        ! unpacked data in fieldsfile.

    IF ((fixhd(5) == 3).AND.(pack_type == 0)) THEN
      IF (entry_no  ==  1) THEN
        lookup(29,entry_no) = data_add
      ELSE
        lookup(29,entry_no) = lookup(29,entry_no-1)                            &
            + lookup(30,entry_no-1)
      END IF
      lookup(40,entry_no) = lookup(29,entry_no)
    END IF

    CALL decompose(n_cols_out, n_rows_out,0,0,1)
    IF (lookup(data_type,entry_no)  ==  2) THEN

      DO i=1,i_dm
        ifield(i)=field(i)
      END DO

! DEPENDS ON: writflds
      CALL writflds(pp_unit_out,1,entry_no,lookup,len1_lookup,                 &
          ifield,lookup(lblrec,entry_no),fixhd,                                &
          icode,cmessage)
    ELSE

       ! Data unpacked.
! DEPENDS ON: writflds
      CALL writflds(pp_unit_out,1,entry_no,lookup,len1_lookup,                 &
          field,lookup(lblrec,entry_no),fixhd,                                 &
          icode,cmessage)
    END IF

    IF (icode >  0) THEN
      cmessage='FIELDOP : Error in MODEL DUMP'

      CALL ereport('WRITEFF', icode,                                           &
          cmessage)

    END IF

  ELSE
    cmessage='FIELDOP : Pack type not supported'

    CALL ereport('WRITEFF', 1003,                                              &
        cmessage)

  END IF

  CALL setpos(pp_unit_out,fixhd(150)+(entry_no-1)                              &
      *len1_lookup-1,icode)
  CALL buffout(pp_unit_out,lookup(1:,entry_no),fixhd(151),                     &
      len_io,a)

      ! Check for I/O errors
  IF (a  /=  -1.0 .OR. len_io  /=  fixhd(151)) THEN

! DEPENDS ON: ioerror
    CALL ioerror('buffer out of lookup table',a,len_io,                        &
        fixhd(151))
    cmessage='FIELDOP: I/O error'
    icode=25
    RETURN
  END IF

  IF (op  ==  'add     ') THEN

    fixhd(15) = 100
  ELSE IF (op  ==  'subtract') THEN

    fixhd(15) = 200
  ELSE IF (op  ==  'multiply') THEN

    fixhd(15) = 300
  ELSE IF (op  ==  'idiv') THEN

    fixhd(15) = 400
  END IF

  CALL setpos(pp_unit_out,0,icode)
  CALL buffout(pp_unit_out,fixhd(1:),len_fixhd,len_io,a)

      ! Check for I/O errors
  IF(a  /=  -1.0 .OR. len_io  /=  len_fixhd) THEN
! DEPENDS ON: ioerror
    CALL ioerror('buffer out of fixed length header',a,len_io,                 &
        len_fixhd)
    cmessage='FIELDOP: I/O error'
    icode=1
    RETURN
  END IF
  RETURN
END SUBROUTINE writeff
