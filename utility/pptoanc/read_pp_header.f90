! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:
!
! Subroutine interface:
      subroutine read_pp_header (ftin1,pp_int,pp_real,ibm_to_cray,      &
     &                           l_bit_32)

      USE um_types
      IMPLICIT NONE
!
! Description:
! This reads the pp headers from a Fortan unformatted file.
!
! Method:
! Data is read either into a 32 bit word, or a 64 bit depending on
! the l_bit_32 or ibm_to_cray logicals. If data is in IBM format it
! is converted to IEEE before passing back. If data is IEEE 32 bit,
! it is converted to IEEE 64 bit.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
! History:
! Version   Date     Comment
! -------   ----     -------
!          16/06/94  Original code. Dave Robinson
! 4.4      14/8/97   Consolidated in UM  Ian Edmond
! 6.1      10/03/04  Port to NEC - subroutined. Paul Selwood
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!
! Declarations:


! Subroutine arguments
!   Scalar arguments with intent(in):
      integer (kind=integer64) :: ftin1






      logical :: ibm_to_cray
      logical :: l_bit_32


!   Array  arguments with intent(in):
      integer (kind=integer64) :: pp_int(45)
      real    (kind=real64)    :: pp_real(19)

! local scalars
      integer (kind=integer64) :: i
      integer (kind=integer64) ::ier

! local arrays
      integer (kind=integer64) :: pp_buffer(32)
      integer (kind=integer32) :: pp_int_32(45)
      real    (kind=real32)    :: pp_real_32(19)

! functions called
      integer (kind=integer64) ::ibm2ieee

! End of header

      if (ibm_to_cray) then

! Read in the PP header
        read(ftin1) pp_buffer

! Convert Integer part of header (Words 1-45)
        ier = ibm2ieee(2,45,pp_buffer,0,pp_int,1,64,32)

! Convert Real part of header (Words 46-64)
        ier = ibm2ieee(3,19,pp_buffer(23),32,pp_real,1,64,32)

      else

! Read in the PP header
        If (L_bit_32) Then
          read(ftin1) pp_int_32,pp_real_32
          Do i = 1, 45
            pp_int(i) = pp_int_32(i)
          End Do

          DO i = 1, 19
            pp_real(i) = pp_real_32(i)
          END DO
        Else
           ! read in straight 64 bit data
           read(ftin1) pp_int, pp_real
        End If  ! 32 bit

      End If


      return
      END SUBROUTINE read_pp_header
