! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:
      subroutine readdata(rows,columns,ftin1,ibm_to_cray,               &
     &                    len_extra,l_bit_32,l_skip,                    &
     &                    data, extra_data)

      USE um_types
  IMPLICIT NONE
!
! Description:
! Reads PP data from a Fortran unformatted file.
!
! Method:
! Data is read either into a 32 bit word, or a 64 bit depending on
! the l_bit_32 or ibm_to_cray logicals. If data is in IBM format it
! is converted to IEEE before passing back. If data is IEEE 32 bit,
! it is converted to IEEE 64 bit.
!
! If l_skip is true, the data is read, but nothing is passed back.
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
! History:
! Version   Date     Comment
! -------   ----     -------
!          16/06/94  Original code. Dave Robinson
! 4.4      14/8/97   Consolidated in UM  Ian Edmond
! 6.1      10/03/04  Port to NEC - subroutined. Paul Selwood.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!
! Declarations:
!
! Subroutine arguments
!   Scalar arguments with intent(in):


      integer (kind=integer64):: rows
      integer (kind=integer64):: columns
      integer (kind=integer64):: ftin1
      integer (kind=integer64):: len_extra







      logical :: ibm_to_cray
      logical :: l_bit_32
      logical :: l_skip


      real (kind=real64) :: data(rows*columns)
      real (kind=real64) :: extra_data(len_extra+1)

! Local scalars:
      integer (kind=integer64) ::  i
      integer (kind=integer64) :: ierr

! local arrays:
      real (kind=real64) :: field2(rows*columns)
      real (kind=real64) :: extra_data2(len_extra+1)

      real (kind=real32) :: extra_data1(len_extra)
      real (kind=real32) :: field1(rows*columns)

! externals
      integer (kind=integer64), external :: ibm2ieee
! End of header

      if (ibm_to_cray .OR. l_bit_32) then     ! 32 bit read
        if (len_extra >  0) then
          read (ftin1) field1,(extra_data1(i),i=1,len_extra)
        else
          read (ftin1) field1
        endif

        ! If not skipping, copy data to output in relevant
        ! format
        if (.NOT. l_skip) then
          if (ibm_to_cray) then
            ierr=ibm2ieee(3, rows*columns, field1, 0_8, data,           &
     &                    1_8, 64_8, 32_8)
            ierr=ibm2ieee(3, rows*columns, extra_data1, 0_8, extra_data,&
     &                   1_8, 64_8, 32_8)

          else      ! straight 32 bit data

            Do i = 1, rows*columns
              data(i) = field1(i)
            End Do

            Do i = 1, len_extra
              extra_data(i) = extra_data1(i)
            End Do

          end if ! ibm_to_cray
        end if   ! l_skip


      else  ! 64 bit read
        If (l_skip) Then   ! read into junk array
          if (len_extra >  0) then
            read (ftin1) field2,(extra_data2(i),i=1,len_extra)
          else
            read (ftin1) field2
          endif
        Else              ! read into output array
          if (len_extra >  0) then
            read (ftin1) data,(extra_data(i),i=1,len_extra)
          else
            read (ftin1) data
          endif
        End If
      end if

      return
      END SUBROUTINE readdata
!
! Subroutine interface:
