! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
SUBROUTINE integer_to_real_fieldop1(unused,integer_field,field,nvals,          &
    max_len,ilabel,icode)
  USE lookup_addresses

  IMPLICIT NONE
!
! Description: Converts integer data into real.
!

! Subroutine arguments
!   Scalar arguments with intent(in):
  INTEGER                                                                      &
      unused,                                                                  &
                            !IN full unpacked size of a field
      max_len,                                                                 &
      nvals                !IN no of values in an input field

!   Array  arguments with intent(in):
  INTEGER                                                                      &
      integer_field(max_len)  ! contains integer data.

!   Scalar arguments with intent(out):
  INTEGER                                                                      &
      icode                !OUT error code

!   Array arguments with intent(out):
  INTEGER                                                                      &
      ilabel(44)           !OUT integer part of lookup

  REAL                                                                         &
      field(max_len)          !OUT contains Real data.

! Local scalars:
  INTEGER                                                                      &
      i                    ! loop counter

!- End of header

  DO  i =1,nvals
    field(i) = integer_field(i)
  END DO

  ilabel(data_type) =1       ! The data type must now be real
  icode=0

  RETURN
END SUBROUTINE integer_to_real_fieldop1
