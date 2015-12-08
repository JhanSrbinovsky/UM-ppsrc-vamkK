! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine interface:
      subroutine calc_len_cfi(ftin2,cols_nowrap,len1_rowdepc,           &
     &     nlevels,len_cfi,fldsizelev,ibm_to_cray,                      &
     &     add_wrap_pts,l_bit_32,icode)

      USE um_types, ONLY: real32
      implicit none
!
! Description:
!           this subroutine calculates the dimensions of the
!           compression arrays. It is a subset of the subroutine
!           calc_cfi_and_fld
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
! Subroutine arguments
!   Scalar arguments with intent(in):

      integer ftin2           ! (in) unit number for levels dataset
      integer cols_nowrap     ! (in) number of points east-west
                              !      (without wrap points)
      integer len1_rowdepc    ! (in) number of points north-south
      integer nlevels         ! (in) number of points in vertical


      logical ibm_to_cray   ! T => input pp data is in IBM number
                            !      format and needs to be converted to
                            !      run on the Cray.
      logical add_wrap_pts  ! T => add wrap points to the output file
      logical l_bit_32

      CHARACTER(LEN=80) levels

      integer icode      ! error code

! Array arguments with intent(in):

      integer len_cfi(3) ! (out) total number of sea segments
      integer fldsizelev(nlevels) !(out) no. of points on each level
                                  ! of compressed field


! Local Scalars

      integer columns    ! no. of columns in levels dataset
      integer rows       ! no. of rows in levels dataset
      integer i,j,k      ! local loop indices
      integer ierr       ! return code from ibm2ieee
      integer count      ! local counter for points in a sea segment
      integer seg_count  ! local counter for number of sea segments
      integer last_count ! local counter for calcultaing points in a
                         ! sea segment
! Local dynamic arrays:

      integer pp_int(45)     ! integer part of levels lookup header
      real pp_real(19)       ! real part of levels lookup header

      real(KIND=real32) levels_in(cols_nowrap*len1_rowdepc)
                              ! temp array for levels dataset to be
                              ! converted to cray number format
      real levels_array(cols_nowrap,len1_rowdepc)
                              ! array of ocean levels

! Function & Subroutine calls:
      integer ibm2ieee

!- End of header

!L 1. Take required dimensions from levels dataset

!L 1.2 Obtain columns and rows by reading header

! DEPENDS ON: read_pp_header
      call read_pp_header(ftin2,pp_int,pp_real,ibm_to_cray,l_bit_32)

      rows = pp_int(18)
      columns = pp_int(19)

      print*,'rows = ',rows
      print*,'columns = ',columns

!L 1.3 Check the dimensions and read the levels_array.

      if (len1_rowdepc  /=  rows) then
         write(6,*)'wrong number of rows in SIZES namelist'
         write(6,*)'len1_rowdepc should equal rows in levels dataset'
         write(6,*)'resubmit'
         icode = 222
         go to 9999      ! Jump out
      end if

      if (cols_nowrap  /=  columns) then
         write(6,*)'wrong number of columns in SIZES namelist'
         write(6,*)'len1_coldepc should equal columns in levels dataset'
         write(6,*)'resubmit'
         icode = 223
         go to 9999      ! Jump out
      end if

!L Do number conversion if required

      if (ibm_to_cray) then
        read(ftin2) levels_in
        ierr=ibm2ieee(3,rows*columns,levels_in,0,levels_array,          &
     &                   1,64,32)
        ierr = 0
      else
        read(ftin2)levels_array
      end if

      close(ftin2)

!L 2. Calculate len_cfi and fldsizelev

!L 2.1 Loop over the points in the field to calculate the number of
!L segments

      count=0
      seg_count=0
      last_count=0

      do k=1,nlevels
       do j=1,rows

        if (k  <=  levels_array(1,j)) then
          count = count + 1
          seg_count = seg_count + 1
        end if

       do i=2,cols_nowrap

        if (k  <=  levels_array(i,j)) then
          count = count + 1
        end if

        if ((k  >   levels_array(i-1,j)) .and.                          &
     &                   (k <= levels_array(i,j))) then
          seg_count=seg_count+1
        end if

      END DO
      END DO

        fldsizelev(k) = count - last_count
        print*,'k = ',k
        print*,'fldsizelev(k) = ',fldsizelev(k)

        last_count = count

      END DO

      len_cfi(1) = seg_count
      len_cfi(2) = seg_count
      len_cfi(3) = rows * nlevels

      print*,'len_cfi(1) = ',len_cfi(1)
      print*,'len_cfi(2) = ',len_cfi(2)
      print*,'len_cfi(3) = ',len_cfi(3)


 9999 continue
      return
      END SUBROUTINE calc_len_cfi
