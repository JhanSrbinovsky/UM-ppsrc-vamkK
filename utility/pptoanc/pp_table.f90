! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Works out the lookup tables for the dump/ancillary    *
!           file header from the pp fields                       *
!
! Subroutine interface:
      subroutine pp_table(pp_int,pp_real,nfields,lookup,rlookup,        &
     & fieldsize,n,levn,m,runtot,number_of_codes,field_code,            &
     & stash_code,add_wrap_pts,compress,pack32,wave,len1_levdepc,       &
     & len2_levdepc,lev_dep_consts,len_realc,real_const,icode)

      USE c_model_id_mod, ONLY: model_id
      USE lookup_addresses

      IMPLICIT NONE
!
! Description:
!            Works out the lookup tables for the dump/ancillary
!            file header from the pp fields
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

      integer nfields   ! dimension for lookup tables

      integer fieldsize  ! size of field to be stored in anc file
      integer n          ! field number
      integer levn       ! level number
      integer m          ! field type
      integer number_of_codes
      integer len1_levdepc
      integer len2_levdepc
      integer len_realc

      logical add_wrap_pts
      logical compress
      logical pack32
      logical wave  ! T for wave dump creation


!   Array  arguments with intent(in):

      integer pp_int(45)
      real pp_real(19)
      integer field_code(number_of_codes)  !stash and field codes
      integer stash_code(number_of_codes)  !input by user


!   Scalar arguments with intent(in/out):
      integer runtot  ! start address for this field on input and
                      ! for next field on output
      integer icode

!   Array  arguments with intent(in/out):

      integer lookup(45,nfields)
      real rlookup(46:64,nfields)

      real lev_dep_consts(1+len1_levdepc*len2_levdepc)
      real real_const(len_realc)


! Local Scalars
      integer i

! Function & Subroutine calls:
      Integer get_um_version_id
!- End of header

      do i=1,45
        lookup(i,n) = 0
      enddo

      do i=46,64
        rlookup(i,n) = 0.0
      enddo

      lookup(lbyr,n)   = pp_int(1)   ! lbyr
      lookup(lbmon,n)  = pp_int(2)   ! lbmon
      lookup(lbdat,n)  = pp_int(3)   ! lbdat
      lookup(lbhr,n)   = pp_int(4)   ! lbhr
      lookup(lbmin,n)  = pp_int(5)   ! lbmin
      lookup(lbsec,n)  = pp_int(6)   ! lbsec
      lookup(lbyrd,n)  = pp_int(7)   ! lbyrd
      lookup(lbmond,n) = pp_int(8)   ! lbmond
      lookup(lbdatd,n) = pp_int(9)   ! lbdatd
      lookup(lbhrd,n)  = pp_int(10)  ! lbhrd
      lookup(lbmind,n) = pp_int(11)  ! lbmind
      lookup(lbsecd,n) = pp_int(12)  ! lbsecd

      lookup(lbtim,n)  = pp_int(13) ! lbtim
      lookup(lbft,n)   = pp_int(14) ! lbft
      lookup(lbcode,n) = pp_int(16) ! lbcode
      lookup(lbhem,n)  = pp_int(17) ! lbhem

!L 1.0 Obtain rows and columns depending on compress and
!L add_wrap_pts.

      if (add_wrap_pts) then
       if (compress) then
        lookup(lbrow,n) = 0  ! no rows if data is compressed
        lookup(lbnpt,n) = 0  ! no columns if data compressed
       else
        lookup(lbrow,n) = pp_int(18) ! lbrow
        lookup(lbnpt,n) = pp_int(19)+2 ! lbnpt
       end if

      else

       if (compress) then
        lookup(lbrow,n) = 0  ! no rows if data is compressed
        lookup(lbnpt,n) = 0  ! no columns if data compressed
       else
        lookup(lbrow,n) = pp_int(18) ! lbrow
        lookup(lbnpt,n) = pp_int(19) ! lbnpt
       end if

      endif
!L
      if (compress) then

       if(.not. wave) then

!C      compression using CFI
        lookup(lbpack,n) = 00110  ! compression

       else

!C      for wave dump compression using ls mask
        lookup(lbpack,n) = 00220  ! compression to sea points

       endif

      else
        lookup(lbpack,n) = 00000  ! no compression
      end if

      if (pack32) then
        lookup(lbpack,n) = lookup(lbpack,n) + 2  ! lbpack
      end if

      lookup(lblrec,n) = fieldsize  ! lblrec
      lookup(lbext,n)  = pp_int(20) ! lbext
      lookup(lbrel,n)  = 3          ! lbrel
      lookup(lbfc,n)   = pp_int(23) ! lbfc
      lookup(lbproc,n) = pp_int(25) ! lbproc
      lookup(lbvc,n)   = pp_int(26) ! lbvc
      lookup(lbegin,n) = runtot     ! lbegin

      lookup(lblev,n)  = levn       ! lblev (level number)

!MH for spectral wave energy the required value is already in pp-header
      if(pp_int(23) == 351) then
       lookup(lblev,n) = pp_int(33) ! wave model freq number
       lookup(44,n)    = pp_int(44) ! wave model dir number
      endif

      lookup(lbproj,n) = pp_int(lbproj)
      lookup(lbtyp,n)  = pp_int(lbtyp)
      lookup(lblev,n)  = pp_int(lblev)
      lookup(lbplev,n) = pp_int(lbplev)

! DEPENDS ON: get_um_version_id
      lookup(lbsrce,n) = get_um_version_id(model_id)
      lookup(data_type,n) = pp_int(data_type)  !   data type

      if (lookup(data_type,n) <  1 .or. lookup(data_type,n) >  3) then
        write (6,*) '********** WARNING ****************** '
        write (6,*) ' Data Type= ',lookup(data_type,n),' for field ',n, &
     &              ' is not recognised.'
        write (6,*) ' Either correct PP Header or set it through',      &
     &              ' the HEADER_DATA namelist.'
        write (6,*) '********** WARNING ****************** '
      endif

      lookup(naddr,n) = runtot     ! start address in data
      lookup(item_code,n) = pp_int(item_code)

      lookup(model_code,n) = pp_int(model_code) ! sub model identifier

      write (6,*) 'Field No ',n,' PP Field Code = ',lookup(lbfc,n),     &
     &                          ' Stash Code  = ',lookup(item_code,n)

      runtot=runtot+fieldsize

      rlookup(blev,n)   = pp_real(52-45) ! blev / hybrid lev 'B' value
      rlookup(brlev,n)  = pp_real(53-45) ! brlev
      rlookup(bhlev,n)  = pp_real(54-45) ! bhlev / hybrid lev 'A' value
      rlookup(bhrlev,n) = pp_real(55-45) ! bhrlev
      rlookup(bplat,n)  = pp_real(56-45) ! bplat
      rlookup(bplon,n)  = pp_real(57-45) ! bplon
      rlookup(bgor,n)   = pp_real(58-45) ! bgor
      rlookup(bzy,n)    = pp_real(59-45) ! bzy
      rlookup(bdy,n)    = pp_real(60-45) ! bdy
      rlookup(bzx,n)    = pp_real(61-45) ! bzx
      rlookup(bdx,n)    = pp_real(62-45) ! bdx

      rlookup(bmdi,n) = rmdi           ! bmdi
      rlookup(bmks,n) = pp_real(64-45) ! bmks

! for spectral wave energy the required value is set from real_const
      if(pp_int(23) == 351) then
        rlookup(blev,n) = lev_dep_consts(pp_int(33)) ! wave model freq
        rlookup(bhlev,n)= (pp_int(44)-1)*real_const(13) ! direction
      endif

 9999 continue
      return
      END SUBROUTINE pp_table
