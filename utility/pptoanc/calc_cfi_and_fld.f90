! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
! Subroutine interface:
      subroutine calc_cfi_and_fld(ftin2,nlevels,len1_coldepc,           &
            cols_nowrap,len1_rowdepc,len1_flddepc,len2_flddepc,         &
            fields_const,fldsizelev,len_cfi,cfi1,cfi2,cfi3,compress,    &
            flddepc,ibm_to_cray,add_wrap_pts,imdi,l_bit_32,icode)

USE filenamelength_mod, ONLY :                                          & 
    filenamelength

      implicit none
!
! Description:
!     this subroutine calculates the compression arrays:
!     cfi1(len_cfi(1)), cfi2(len_cfi(2)) and cfi3(len_cfi(3))
!     using an array of numbers of ocean levels at each point:
!     levels_array(len1_coldepc,len1_rowdepc)
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
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Subroutine arguments
!   Scalar arguments with intent(in):
      integer ftin2      ! (in) unit numbr for levels dataset
      integer nlevels    ! (in) number of points in vertical

      integer len1_coldepc ! (in) 1st dimension of col_dep_consts
      integer cols_nowrap  ! (in) no. of points east-west
      integer len1_rowdepc ! (in) 1st dimension of row_dep_consts
      integer len1_flddepc ! (in) 1st dimension of fields_const
      integer len2_flddepc ! (in) 2nd dimension of fields_const
      integer imdi         ! (in) integer missing data indicator
      integer icode        ! error code

      logical compress     ! T => the dump is to be compressed
      logical flddepc      ! T => fields_const are wanted in the dump
      logical ibm_to_cray  ! T => input pp data is in IBM number
                           !      format and needs to be converted to
                           !      run on the Cray.

      logical add_wrap_pts ! T => add wrap points to the output field
      logical l_bit_32


!   Array arguments with intent(in):


      integer fldsizelev(nlevels) ! number of points on each compressed
                                  ! level

      integer len_cfi(3) ! (in) total number of sea segments

      integer cfi1(len_cfi(1))  ! (out) index array for compressed array
      integer cfi2(len_cfi(2))  ! (out) index array for expanded array
      integer cfi3(len1_rowdepc,nlevels)  ! (out)
                     ! contains number of first sea
                     ! segment in each row at each levelc
                     ! if there is a sea segment in the row
                     ! contains number of next sea segment
                     ! otherwise

      real fields_const(len1_flddepc,len2_flddepc) ! (out) array for
                           ! fields of constants

! Local scalars :

      integer columns  ! no. of points east-west
      integer rows     ! no. of points north-south
      integer i,j,k    ! local loop indices
      integer count    ! local counter for points in a sea segment
      integer seg_count! local counter for number of sea segments

      CHARACTER(LEN=filenamelength) :: levels

! Local dynamic arrays :

      integer pp_int(45)

      real    pp_real(19)
      real temp_levels_array(len1_coldepc,len1_rowdepc)
                                  ! local array of ocean levels
      real levels_array(cols_nowrap,len1_rowdepc)
                                  ! local array of ocean levels
      real dummy(1)


!- End of header

!L 1. Read the fields_const from levels dataset

!L 1.1 Read the data from levels dataset
      call get_file(ftin2,LEVELS,filenamelength,icode)
      levels=trim(levels)
      OPEN( unit=ftin2, file=levels, form="unformatted" )

! DEPENDS ON: read_pp_header
      call read_pp_header(ftin2,pp_int,pp_real,ibm_to_cray,l_bit_32)

      rows = pp_int(18)
      columns = pp_int(19)

      print*,'rows = ',rows
      print*,'columns = ',columns

!L 1.4 Read in levels_array and check the dataset is on the
!L same grid as the input pp fields.  If add_wrap_pts and flddepc
!L then add wrap points to the levels dataset.  The compression indices
!L are the same for an output dump with or without wrap points.

! DEPENDS ON: readdata
      CALL readdata( rows, columns, ftin2, ibm_to_cray, 0,              &
     &               l_bit_32, .FALSE., levels_array, dummy)
      close(ftin2)

      if (add_wrap_pts) then

        if (len1_rowdepc  /=  rows) then
           write(6,*)'wrong number of rows in SIZES namelist'
           write(6,*)'len1_rowdepc should equal rows in levels dataset'
           write(6,*)'resubmit'
           icode = 222
           go to 9999      ! Jump out
        end if

        if (len1_coldepc  /=  columns+2) then
           write(6,*)'wrong number of columns in SIZES namelist'
           write(6,*)'len1_coldepc should equal columns+2 in '//&
               'levels dataset'
           write(6,*)'resubmit'
           icode = 223
           go to 9999      ! Jump out
        end if

        if (len1_coldepc*len1_rowdepc  /=  len1_flddepc) then
           write(6,*)'len1_flddepc should equal '//&
               'len1_coldepc*len1_rowdepc'
           write(6,*)'resubmit'
           icode = 224
           go to 9999      ! Jump out
        end if

        do j = 1,rows
          do i = 1,columns
            temp_levels_array(i,j) = levels_array(i,j)
          enddo
        enddo

        do j = 1,rows
          temp_levels_array(columns+1,j)=temp_levels_array(1,j)
          temp_levels_array(columns+2,j)=temp_levels_array(2,j)
        enddo

        do j = 1,rows
         do i = 1,len1_coldepc
           fields_const(i+(j-1)*len1_coldepc,1) = temp_levels_array(i,j)
         enddo
        enddo

      else

        if (len1_rowdepc  /=  rows) then
           write(6,*)'wrong number of rows in SIZES namelist'
           write(6,*)'len1_rowdepc should equal rows in levels dataset'
           write(6,*)'resubmit'
           icode = 222
           go to 9999      ! Jump out
        end if

        if (len1_coldepc  /=  columns) then
           write(6,*)'wrong number of columns in SIZES namelist'
           write(6,*)'len1_coldepc should equal columns in '//&
               'levels dataset'
           write(6,*)'resubmit'
           icode = 223
           go to 9999      ! Jump out
        end if

        if (len1_coldepc*len1_rowdepc  /=  len1_flddepc) then
           write(6,*)'len1_flddepc should equal '//&
               'len1_coldepc*len1_rowdepc'
           write(6,*)'resubmit'
           icode = 224
           go to 9999      ! Jump out
        end if

        do j = 1,rows
         do i = 1,len1_coldepc
           fields_const(i+(j-1)*len1_coldepc,1) = levels_array(i,j)
         enddo
        enddo

      end if

!L 2.1 Initialise cfi3 array and create the compression indices.

      if (compress) then

      do k=1,nlevels
        do j=1,rows
          cfi3(j,k)=imdi
        END DO
      END DO

      count=0
      seg_count=0

      do k=1,nlevels
        do j=1,rows
!
!     if the first element in a row is sea, a new segment is starting,
!     so count and seg_count are both incremented, and cfi1 and
!     cfi2 have new entries.  Columns is used here instead of
!     len1_coldepc as the index to this array expects that in oa_pack.
!
          if (k <= levels_array(1,j)) then
            count=count+1
            seg_count=seg_count+1
            cfi1(seg_count)=count
            cfi2(seg_count)=1+(j-1)*columns+(k-1)*columns*rows
            cfi3(j,k)=seg_count
           end if

             do i=2,columns
!
!     if present point is sea, add one to count
!
               if (k <= levels_array(i,j)) then
                  count=count+1
               end if
!
!     if present point is sea and previous point is land,
!     a new segment is starting, so seg_count is incremented
!     and cfi1 and cfi2 have new entries. Columns is used here instead
!     of len1_coldepc as the index to this array expects that
!     in oa_pack.
!
               if ((k >  levels_array(i-1,j)).and.                      &
     &                          (k <= levels_array(i,j))) then
                 seg_count=seg_count+1
                 cfi1(seg_count)=count
                 cfi2(seg_count)=i+(j-1)*columns+(k-1)*columns*rows
!
!     if cfi3(j,k) has not been reset,
!     then the present segment must be the first in the row
!
                    if (cfi3(j,k) == imdi) then
                      cfi3(j,k)=seg_count
                    end if
               end if
      END DO

!
!     if there is no sea segment in the row,
!     then set cfi3 to seg_count+1
!
            if (cfi3(j,k) == imdi) then
               cfi3(j,k)=seg_count+1
            end if

      END DO
      END DO

      end if
9999  continue
      return
      END SUBROUTINE calc_cfi_and_fld
