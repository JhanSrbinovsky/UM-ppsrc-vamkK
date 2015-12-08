! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
! Subroutine interface:
!
! Subroutine interface:
      subroutine anc_head(pp_int,pp_real,rows,columns,fieldsize,nfields,&
     & field_types,n_times,nlevels,n_freq_waves,n_dir_waves,no_cmp,     &
     & len1_levdepc,len2_levdepc,len1_rowdepc,len2_rowdepc,len1_coldepc,&
     & len2_coldepc,len1_flddepc,len2_flddepc,len_extcnst,len_cfi,      &
     & tracer_grid,add_wrap_pts,periodic,single_time,ibm_to_cray,       &
     & compress,levdepc,rowdepc,coldepc,flddepc,extcnst,wave,           &
     & lfh,fixhd,len_intc,int_const,len_realc,real_const,icode)

      USE conversions_mod, ONLY: pi
      USE check_iostat_mod

      implicit none
!
! Description:
!              Creates the dump/ancillary file header.
!            (Fixed length header,integer constants and real constants)
!
!
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
      integer rows
      integer columns
      integer fieldsize
      integer nfields
      integer field_types
      integer n_times
      integer nlevels
      integer n_freq_waves
      integer n_dir_waves
      integer no_cmp      ! no. of total compressed points in compressed
                          ! array

      integer len1_levdepc
      integer len2_levdepc
      integer len1_rowdepc
      integer len2_rowdepc
      integer len1_coldepc
      integer len2_coldepc
      integer len1_flddepc
      integer len2_flddepc
      integer len_extcnst

      integer lfh
      integer len_intc
      integer len_realc

      integer icode

!L logical choices (IN)

      logical tracer_grid
      logical add_wrap_pts
      logical periodic
      logical single_time
      logical ibm_to_cray
      logical compress
      logical levdepc
      logical rowdepc
      logical coldepc
      logical flddepc
      logical extcnst

      logical wave   ! T for creating wave dump

!   Array  arguments with intent(in):
      integer pp_int(45)
      integer len_cfi(3)
      real pp_real(19)

!   Array  arguments with intent(out):

      integer fixhd(lfh)
      integer int_const(len_intc)

      real real_const(len_realc)

! Local Scalars

      integer ipos
      integer i,j
      integer errorstatus

      integer fvhh,fvdd,fvmm,fvyy   ! hour,day,month,year first validity
                                    ! time
      integer lvhh,lvdd,lvmm,lvyy   ! hour,day,month,year last validity
                                    ! time
      integer ivhh,ivdd,ivmm,ivyy   ! hour,day,month,year interval

      logical year360         ! true for 360-day calendar
      logical l_first_vt
      logical l_last_vt

      integer cdays ! century days
      integer chours ! century hours
      integer new_cdays ! new century days
      integer new_chours! new century hours
      integer ihr,idy,imn

! Function & Subroutine calls:

      integer FIND_NAMELIST

!- End of header

      namelist /first_vt/ fvhh,fvdd,fvmm,fvyy
      namelist /last_vt/  lvhh,lvdd,lvmm,lvyy
      namelist /interval/ year360,ivhh,ivdd,ivmm,ivyy


!L 1. Set fixed header

!L 1.0 Initialise to missing data
! (* set dimensions of all arrays to 1 *)

! DEPENDS ON: init_flh
      call init_flh (fixhd,lfh)

!L 1.1 First 9 elements in header DEFAULTS

      fixhd(2)=2               ! indicator for the ocean
      fixhd(3)=4               ! depth coordinates
      fixhd(4)=0               ! global grid
      fixhd(5)=4               ! ancillary fields dataset
      fixhd(8)=2               ! Calendar indicator
      fixhd(9)=1               ! Indicator for grid staggering

!L 1.2 Set dates

      if (periodic) then
        fixhd(10)=2
      else if (single_time) then
        fixhd(10)=0
      else
        fixhd(10)=1
      end if

      fixhd(12)=401                      ! UM Version number

! (* first validity time *)
! (* fixhd(21-27)       *)

      fvhh = 0
      fvdd = 0
      fvmm = 0
      fvyy = 0

      rewind(5)
! DEPENDS ON: find_namelist
      I=FIND_NAMELIST(5,"FIRST_VT")

      If(I == 0)then
        READ (UNIT=5, NML=FIRST_VT, IOSTAT=ErrorStatus)
        CALL check_iostat(errorstatus, "namelist FIRST_VT")
      Else
        write(6,*)'Cannot find namelist FIRST_VT'
      End if
      write (6,*) ' '
      write (6,*) 'FIRST_VT namelist is set up as follows:-'
      write (6,*) ' '
      write (6,first_vt)

!     Test if first VT has been provided in namelist
      l_first_vt = .not.                                                &
     &  (fvhh == 0 .and. fvdd == 0 .and. fvmm == 0 .and. fvyy == 0)

      if (l_first_vt) then  !  First VT given in namelist

        fixhd(21) = fvyy
        fixhd(22) = fvmm
        fixhd(23) = fvdd
        fixhd(24) = fvhh
        fixhd(25) = 0
        fixhd(26) = 0
        fixhd(27) = 0

      else if (single_time) then

        fixhd(21:27) = 0    !  Set First VT to zero

      else  !  Get first VT from first PP Header

        fixhd(21) = pp_int(1)
        fixhd(22) = pp_int(2)
        fixhd(23) = pp_int(3)
        fixhd(24) = pp_int(4)
        fixhd(25) = pp_int(5)
        fixhd(26) = 0
        fixhd(27) = pp_int(6)

      endif

      year360=.false.
      ivhh = 0
      ivdd = 0
      ivmm = 0
      ivyy = 0

! (* interval is read from namelist*)
! (* fixhd(35-41)       *)

      rewind(5)
! DEPENDS ON: find_namelist
      I=FIND_NAMELIST(5,"INTERVAL")

      If(I == 0)then
        READ (UNIT=5, NML=INTERVAL, IOSTAT=ErrorStatus)
        CALL check_iostat(errorstatus, "namelist INTERVAL")
      Else
        write(6,*)'Cannot find namelist INTERVAL'
      End if

      write (6,*) ' '
      write (6,*) 'INTERVAL namelist is set up as follows:-'
      write (6,*) ' '
      write (6,interval)

      fixhd(35) = ivyy
      fixhd(36) = ivmm
      fixhd(37) = ivdd
      fixhd(38) = ivhh
      fixhd(39) = 0
      fixhd(40) = 0
      fixhd(41) = 0

! (* last validity time *)
! (* fixhd(28-34)       *)

      lvhh = 0
      lvdd = 0
      lvmm = 0
      lvyy = 0

      rewind(5)
! DEPENDS ON: find_namelist
      I=FIND_NAMELIST(5,"LAST_VT")

      If(I == 0)then
        READ (UNIT=5, NML=LAST_VT, IOSTAT=ErrorStatus)
        CALL check_iostat(errorstatus, "namelist LAST_VT")
      Else
        write(6,*)'Cannot find namelist LAST_VT'
      End if

      write (6,*) ' '
      write (6,*) 'LAST_VT namelist is set up as follows:-'
      write (6,*) ' '
      write (6,last_vt)

!     Test if last VT has been provided in namelist
      l_last_vt = .not.                                                 &
     &  (lvhh == 0 .and. lvdd == 0 .and. lvmm == 0 .and. lvyy == 0)

      if (year360) then

        fixhd(8) = 2    !  360 day calendar

        if (l_last_vt) then  ! Last VT given in namelist

          fixhd(28) = lvyy
          fixhd(29) = lvmm
          fixhd(30) = lvdd
          fixhd(31) = lvhh
          fixhd(32) = 0
          fixhd(33) = 0
          fixhd(34) = 0

        else if (single_time) then

          fixhd(28:34) = 0    !  Set Last VT to zero

        else  !  calculate last VT from first VT and Interval

          fixhd(33)=fixhd(26)  ! seconds
          fixhd(32)=fixhd(25)  ! minutes

          ihr=fixhd(24)+ivhh*(n_times-1)
          fixhd(31)=mod(ihr,24)

          idy=fixhd(23)+ivdd*(n_times-1)+ihr/24
          fixhd(30)=mod(idy-1,30)+1

          imn=fixhd(22)+ivmm*(n_times-1)+(idy-1)/30
          fixhd(29)=mod(imn-1,12)+1

          fixhd(28)=fixhd(21)+ivyy*(n_times-1)+(imn-1)/12

          fixhd(34)=(fixhd(29)-1)*30+fixhd(30)

        endif

      else   !  365 calander files

        fixhd(8) = 1    !  365 day calendar

        if (l_last_vt) then  !  Last VT given in namelist

          fixhd(28) = lvyy
          fixhd(29) = lvmm
          fixhd(30) = lvdd
          fixhd(31) = lvhh
          fixhd(32) = 0
          fixhd(33) = 0
          fixhd(34) = 0

        else if (single_time) then

          fixhd(28:34) = 0    !  Set Last VT to zero

        else  !  calculate last VT from first VT and Interval

!         Check First VT and Interval first
!         First VT is OK if FIXHD(21,22,23) are all set.
!         Interval is OK if only IVHH and/or IVDD are used.
          if (fixhd(21) <= 0 .or. fixhd(22) <= 0 .or.                   &
     &        fixhd(23) <= 0 .or. ivmm >  0 .or. lvyy >  0) then
            write (6,*) ' '
            write (6,*) ' ERROR in ANC_HEAD. Last Validity Time ?'
            write (6,*) ' Last VT cant be calculated from first VT.'
            write (6,*) ' Rerun job with last VT in LAST_VT namelist.'
            write (6,*) ' '
            icode = 1
            go to 9999   !  Return
          endif

! (* calculate century day and hour for first validity time)
! DEPENDS ON: date31
          call date31(fixhd(23),fixhd(22),fixhd(21),cdays)
          chours=(cdays-1)*24+fixhd(24)

!          write(6,*)'cdays=',cdays
!          write(6,*)'chours=',chours

! (* add time interval)
          new_chours=chours+(n_times-1)*(ivhh+ivdd*24)

! (* convert to new century day)
          new_cdays=1+new_chours/24

! (* convert to actual date)
! DEPENDS ON: date13
          call date13(new_cdays,fixhd(30),fixhd(29),fixhd(28))
          fixhd(31)=new_chours-24*(new_cdays-1)

          fixhd(32) = 0
          fixhd(33) = 0
          fixhd(34) = 0

        endif

      end if

      WRITE(6,'('' '')')
      WRITE(6,'('' Validity Times (VT) in Ancillary File.'')')
      WRITE(6,'('' '')')
      WRITE(6,'(''               Year  Month Day Hour Min  Sec DayNo    &
     &'')')
      WRITE(6,'('' First VT    ='',7I5)')(FIXHD(I),I=21,27)
      WRITE(6,'('' Last  VT    ='',7I5)')(FIXHD(I),I=28,34)
      WRITE(6,'('' VT Interval ='',7I5)')(FIXHD(I),I=35,41)
      WRITE(6,'('' '')')

!L 1.3 Set pointers to starts of sections in ancillary file

      ipos = lfh + 1   ! position of start of integer consts

! (* integer constants location *)
      fixhd(100)= ipos
      fixhd(101)=len_intc
      ipos = ipos + len_intc

! (* real constants location *)
      fixhd(105)=ipos
      fixhd(106)=len_realc
      ipos = ipos + len_realc

! (* levels dependent constants*)

      if (levdepc) then
        fixhd(110) = ipos
        fixhd(111) = len1_levdepc
        fixhd(112) = len2_levdepc
        ipos = ipos + len1_levdepc*len2_levdepc
      endif

      if (rowdepc) then
        fixhd(115) = ipos
        fixhd(116) = len1_rowdepc
        fixhd(117) = len2_rowdepc
        ipos = ipos + len1_rowdepc*len2_rowdepc
      endif

      if (coldepc) then
        fixhd(120) = ipos
        fixhd(121) = len1_coldepc
        fixhd(122) = len2_coldepc
        ipos = ipos + len1_coldepc*len2_coldepc
      endif

      if (flddepc) then
        fixhd(125) = ipos
        fixhd(126) = len1_flddepc
        fixhd(127) = len2_flddepc
        ipos = ipos + len1_flddepc*len2_flddepc
      endif

      if (extcnst) then
        fixhd(130) = ipos
        fixhd(131) = len_extcnst
        ipos = ipos + len_extcnst
      endif

! (* compressed field indices *)
      if (compress .and.  .not.wave) then
        fixhd(140) = ipos
        fixhd(141) = len_cfi(1)
        ipos = ipos + len_cfi(1)

        fixhd(142) = ipos
        fixhd(143) = len_cfi(2)
        ipos = ipos + len_cfi(2)

        fixhd(144) = ipos
        fixhd(145) = len_cfi(3)
        ipos = ipos + len_cfi(3)

      end if

! (* location of lookup table *)
      fixhd(150)=ipos
      fixhd(151)=64
      fixhd(152)=nfields

! for wave dump - number of fields in dump *
      fixhd(153)=nfields

! (* location of data *)
      ipos=ipos+64*nfields
      fixhd(160)=ipos

!   fixhd(161) is set in ancfld (it is updated after each field
!   is read )

!L 2. Set Integer Constants

      do j=1,len_intc
        int_const(j)=imdi
      enddo

      int_const(3)=n_times

      if (add_wrap_pts) then
        int_const(6)=columns+2
      else
        int_const(6)=columns
      end if

! When the UM reads the ancillary files (see RPANCO1A) it checks that
! the number of rows in the model tracer grid (JMT) matches the number
! of rows declared in the integer consts; the number of rows in the
! velocity grid is one less  than that in the tracer grid.
!

      if (tracer_grid) then
        int_const(7)=rows
      else
        int_const(7)=rows+1
      end if

      int_const(8) = nlevels

!CMHWAVES
      if(wave) then
       int_const(8) = n_freq_waves
       int_const(9) = n_dir_waves
       int_const(10)= fieldsize
      endif
!CMHWAVES

      if (compress) then
        int_const(11) = no_cmp
      end if

      int_const(15)=field_types

!L 3. Set real constants

      do j=1,len_realc
        real_const(j)=rmdi
      enddo

! (* grid spacing)
      real_const(1)=pp_real(17)
      real_const(2)=abs(pp_real(15))

! (* lat of first row (3) and long of first point on row (4)
      if (tracer_grid) then
        real_const(3)=pp_real(14)+pp_real(15)
        real_const(4)=pp_real(16)+pp_real(17)
      else
        real_const(3)=pp_real(14)+0.5*pp_real(15)
        real_const(4)=pp_real(16)+0.5*pp_real(17)
      end if

! (* test value of the start longitude *)
      if (real_const(4) <  0.0) then
        real_const(4)=real_const(4)+360.0
      else if (real_const(4) >= 360.0) then
        real_const(4)=real_const(4)-360.0
      end if

! (* lat and long of pseudo north pole)
      real_const(5)=pp_real(11)
      real_const(6)=pp_real(12)

! WAVES
! direction increment
       if(wave) then
        real_const(13)=2.*pi/n_dir_waves
       endif

 9999 continue
      return
      END SUBROUTINE anc_head
!
! Subroutine interface:
!
! Subroutine interface:
!
! Subroutine interface:
!
! Subroutine interface:
! Purpose: Works out the lookup tables for the dump/ancillary    *
!           file header from the pp fields                       *
!
! Subroutine interface:

!+ Skip namelists in f90 compiled UM code removing need for
!+ assign -f 77 g:sf  in script
!
! Subroutine Interface:
