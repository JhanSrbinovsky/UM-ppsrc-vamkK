! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs

SUBROUTINE un_pack_fieldop1(pack_type,                                         &
    npoints,                                                                   &
    pdata,                                                                     &
    num_cray_words,                                                            &
    ilabel,                                                                    &
    len_ilabel,                                                                &
    amdi,                                                                      &
    pp_fixhd,                                                                  &
    len_fixhd,                                                                 &
    pppak,                                                                     &
    icode,cmessage)
  USE ereport_mod, ONLY : ereport
  USE lookup_addresses

  IMPLICIT NONE
!
! Description: To unpack data from the input array pdata and return
!              the data in pdata.
!
! Subroutine arguments
!   Scalar arguments with intent(in):
  INTEGER                                                                      &
      npoints,                                                                 &
                               ! full unpacked size of a pdata
      max_len,                                                                 &
      num_cray_words,                                                          &
                            ! length of input pdata
      len_fixhd,                                                               &
      len_ilabel           ! length of ilabel array

  REAL                                                                         &
      amdi                 ! Missing data indicator.

!   Scalar arguments with intent(in):
  INTEGER                                                                      &
      pp_fixhd(len_fixhd)  ! PPfile fixed length header

!   Scalar arguments with intent(in/out):
  INTEGER                                                                      &
      pack_type            ! Type of packing used

!   Array  arguments with intent(in/out):
  INTEGER                                                                      &
      ilabel(len_ilabel)


!   ErrorStatus
  INTEGER                                                                      &
      icode

  CHARACTER                                                                    &
      cmessage*80

! Local scalars:
  INTEGER                                                                      &
      num_unpack_values,                                                       &
                              ! Number of numbers originally packed
      i,                                                                       &
                              ! loop counter
      ixx,                                                                     &
                              ! Returned X dimension from COEX
      iyy,                                                                     &
                              ! Returned Y dimension from COEX
      idum,                                                                    &
                              ! Dummy variable
      pppak                  ! Packing acc

! Local parameters:
  INTEGER len_full_word   ! The length of a FULL_WORD
  PARAMETER(len_full_word=64)

! Local arrays:
  REAL                                                                         &
      field(npoints),                                                          &
                            !WORK array used for un_packing
      pdata(npoints)       ! Input contains packed data.

! Function & Subroutine calls:

!- End of header

  IF (pack_type == 1) THEN     ! WGDOS packing

! DEPENDS ON: coex
    CALL coex(field,                                                           &
                                   ! OUT
        npoints,                                                               &
                                   ! IN
        pdata,                                                                 &
                                   ! IN
        npoints,                                                               &
                                   ! IN
        ixx,iyy,                                                               &
                                   ! OUT
        idum,                                                                  &
        pppak,                                                                 &
                                   ! OUT
        .FALSE.,                                                               &
                                   ! IN
        amdi,                                                                  &
                                   ! IN
        len_full_word,                                                         &
                                   ! IN
        icode,                                                                 &
                                   ! OUT
        cmessage)        ! OUT

    num_unpack_values = ixx * iyy
    ilabel(lblrec) = ilabel(lbrow) * ilabel(lbnpt) + ilabel(lbext)
  ELSE IF (pack_type  ==  2) THEN !  32 Bit CRAY packing


    num_cray_words = num_cray_words*2
! DEPENDS ON: expand21
    CALL expand21(num_cray_words,                                              &
                                             ! IN
        pdata,                                                                 &
                                             ! IN
        field)                 ! OUT

    num_unpack_values = num_cray_words


  ELSE IF (pack_type  ==  3) THEN !  GRIB packing

    WRITE(6,*) 'Grib unpacking only supported on NEC SX6.'

    CALL ereport('UN_PACK', 1000,                                              &
        cmessage)


  ELSE IF (pack_type  ==  4) THEN ! Run length encoded

! DEPENDS ON: runlen_decode
    CALL runlen_decode(field,                                                  &
                                   ! OUT
        npoints,                                                               &
                                   ! IN
        pdata,                                                                 &
                                   ! IN
        ilabel(lblrec),                                                        &
                                          ! IN
        amdi,                                                                  &
                                   ! IN
        icode,                                                                 &
                                   ! OUT
        cmessage)        ! OUT

    num_unpack_values =  npoints
    ilabel(lblrec) = npoints
  ELSE

    icode=6
    cmessage=' UNPACK - packing type not yet supported'
  END IF

      ! Write unpacked data back into array pdata.
  DO i =1,num_unpack_values
    pdata(i) = field(i)
  END DO

  ilabel(data_type) =1                ! data must now be real
  ilabel(lbpack)    =ilabel(lbpack)-pack_type ! data no
  pack_type         =0                        ! longer packed

  RETURN
END SUBROUTINE un_pack_fieldop1
