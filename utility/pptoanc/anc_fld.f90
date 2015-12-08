! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
! Subroutine interface:
      subroutine anc_fld(ftin2,ftout,nolevsmax,number_of_codes,         &
     &  max_n_pp_files,len_cfi,fldsizelev,                              &
     &  fld_types,n_times,nlevels,n_pp_files,stash_code,field_code,     &
     &  nlevs_code,unit_no,n_freq_waves,n_dir_waves,len_intc,len_realc, &
     &  len_extra,len1_levdepc,len2_levdepc,len1_rowdepc,len2_rowdepc,  &
     &  len1_coldepc,len2_coldepc,len1_flddepc,len2_flddepc,            &
     &  len_extcnst,rmdi_input,                                         &
     &  add_wrap_pts,periodic,single_time,ibm_to_cray,compress,wave,    &
     &  levdepc,rowdepc,coldepc,flddepc,extcnst,pack32,pphead,          &
     &  grid_of_tracer,field_order,lwfio,L_bit_32,                      &
     &  len2_lookup_max,cols_nowrap,icode)

      USE filenamelength_mod, ONLY :                                    &
            filenamelength
      USE IO
      USE um_types
      USE check_iostat_mod
      USE writhead_mod
      USE UM_ParParams
      USE Atmos_Max_Sizes
      USE lookup_addresses

      IMPLICIT NONE
!
! Description:
!          Main subroutine. Creates the ancillary/dump header,
!          lookup tables and writes the header using WRITHEAD.
!          Calls dataw which writes the data out
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

      integer ftin2                  ! input unit for mask file used
                                     ! for fields consts and
                                     ! compression indices.
      integer ftout                  ! unit number for output ancillary
                                     ! file

      integer nolevsmax              ! max number of levels; dimensions
                                     !  fldsizelev array
      integer number_of_codes        ! max number of stash/field codes
      integer max_n_pp_files         ! max number of input pp files



      integer fld_types      ! number of field types in I/O files
      integer n_times          ! number of time periods in I/O files
      integer nlevels          ! number of levels (default = 1)
      integer n_pp_files       ! number of input pp files


      integer n_freq_waves  ! number of wave frequencies
      integer n_dir_waves   ! number of wave directions

      integer len_intc      !  Actual  length of integer constants array
      integer len_realc     !  Actual  length of real constants array

      integer len_extra     ! length of extra data (minimum value = 0)

      integer len1_levdepc   ! Actual 1st dimension of lev_dep_consts
      integer len2_levdepc   ! Actual 2nd dimension of lev_dep_consts

      integer len1_rowdepc   ! Actual 1st dimension of row_dep_consts
      integer len2_rowdepc   ! Actual 2nd dimension of row_dep_consts

      integer len1_coldepc   ! Actual 1st dimension of col_dep_consts
      integer len2_coldepc   ! Actual 2nd dimension of col_dep_consts

      integer len1_flddepc   ! Actual 1st dimension of fields_const
      integer len2_flddepc   ! Actual 2nd dimension of fields_const

      integer len_extcnst    ! Actual 1st dimension of fields_const

      integer len2_lookup_max ! maximum 2nd dimension of the lookup
                              ! table
      integer cols_nowrap     ! no. of columns in field without wrap
      integer icode           ! error code variable

      real    rmdi_input     ! real missing data indicator
                             ! in input pp field

      logical add_wrap_pts   ! T => adds wrapping columns
                             !      e.g. for global grid
      logical periodic       ! T => periodic in time
                             !      e.g. climate field
      logical single_time    ! T => all fields input valid at one time
      logical ibm_to_cray    ! T => input pp data is in IBM number
                             !      format and needs to be converted to
                             !      run on the Cray.
                             !      (Only use if running on Cray)
      logical compress       ! T => fields are packed into ancillary
                             !      field compressed field indices are
                             !      calculated
      logical wave           ! T => a wave dump is to be created
      logical levdepc        ! T => if level dependent constants array
                             !      required
      logical rowdepc        ! T => if row dependant constants are
                             !      required
      logical coldepc        ! T => if column dependant constants are
                             !      required
      logical flddepc        ! T => if fields of constants are
                             !      required
      logical extcnst        ! T => if fields of constants are
                             !      required
      logical pack32         ! T => use 32 bit Cray numbers
      logical pphead         ! T => print out pp headers read in


      logical field_order    ! T => input pp fields ordered by time.
                             !      i.e. different months in input
                             !         files, same fields in all files
                             ! F => inout pp fields ordered by fields.
                             !      i.e. different fields in input
                             !         files, all months in all files
      logical lwfio          ! T => set the LBEGIN and LBNREC fields
                             !      in the LOOKUP Headers for VN 16
                             !      Type Dumpfiles.
                             ! F => Old dumpfiles

      logical l_bit_32       ! T => pp file input is 32 bit
                             ! F => pp file input is 64 bit


!   Array  arguments with intent(in):

      integer len_cfi(3)          ! lengths of compressed field indices
      integer fldsizelev(nolevsmax)! size of packed field on each level
      integer stash_code(number_of_codes) ! array of stash codes
      integer field_code(number_of_codes) ! array of field codes
      integer nlevs_code(number_of_codes) ! array of levels depending
                                          ! on field code
      integer unit_no(number_of_codes)    ! array of unit numbers for
                                          ! input
      logical grid_of_tracer(number_of_codes) ! T => fields are on a
                                              ! tracer grid


! Local parameters:

      integer len_look_user     ! No. of changes to the lookup table
                                ! made by the user
      integer len1_lookup       ! 1st dimension of the lookup
      integer len1_lookup_all   ! Dimension of the whole lookup array
      integer lfh               ! length of the fixed length header

      integer max_len_intc      ! Max dimension of integer constants
      integer max_len_realc     ! Max dimension of real constants
      integer max_len1_levdepc  ! Max 1st dimension of lev_dep_consts
      integer max_len2_levdepc  ! Max 2nd dimension of lev_dep_consts
      integer max_len1_rowdepc  ! Max 1st dimension of row_dep_consts
      integer max_len2_rowdepc  ! Max 2nd dimension of row_dep_consts
      integer max_len1_coldepc  ! Max 1st dimension of col_dep_consts
      integer max_len2_coldepc  ! Max 2nd dimension of col_dep_consts
      integer max_len_extcnst   ! Max dimension of extra_const

      integer len_dumphist      ! Actual dimension of dumphist

      parameter (len_look_user = 50)
      parameter (len1_lookup = 45)
      parameter (len1_lookup_all = 64)
      parameter (lfh=256)

      parameter (max_len_intc=40)
      parameter (max_len_realc=40)
      parameter (max_len1_levdepc=100)
      parameter (max_len2_levdepc=5)
      parameter (max_len1_rowdepc=rows_max)
      parameter (max_len2_rowdepc=5)
      parameter (max_len1_coldepc=row_length_max)
      parameter (max_len2_coldepc=5)
      parameter (max_len_extcnst=1000)


      parameter (len_dumphist=1)

! Local Scalars


      integer ftin1            ! unit number for input pp fields

      integer len_data         ! length of the data record
      integer start_block      ! position of start for WRITHEAD

      integer n_sea_points     ! number of sea points for wave dump

      integer rows             ! no. of rows in input pp field
      integer columns          ! no. of columns in input pp field

      integer i,j              ! loop counters
      integer levn             ! level number
      integer m,n              ! loop counters
      integer np               ! Number of points in pp field

      integer len2_lookup      ! 2nd dimension of lookup table
                               ! total number of fields in output file
      integer len2_step        ! calculation step for len2_lookup

      integer nlevs_this_code  ! # of levels for this field code and
                               ! limit for do -loop over levels(waves)

      integer fieldn           ! present field number
      integer fieldsize        ! size of field when it is stored
                               ! in output data set
      integer runtot           ! running total of start address
                               ! in data array of present field


      integer no_cmp           ! total no. of compressed points in
                               ! compressed array
      integer ipos             ! position counter

      integer irow_number

      integer                                                           &
       disk_address                                                     &
                                       ! Current rounded disk address
      ,number_of_data_words_on_disk                                     &
                                       ! Number of data words on disk
      ,number_of_data_words_in_memory  ! Number of Data Words in memory

      logical tracer_grid   ! T => field is on a tracer grid
                            ! F =>field is on a velocity grid

      logical t_compress    ! compress argument for DATA subroutine
                            ! used for wave dump LSmask set f whatever
                            ! compress is


      CHARACTER(LEN=80) cmessage    ! error message from WRITHEAD
      CHARACTER(LEN=filenamelength) :: ancfile, levels

! Local arrays dimensioned by parameters:

      ! arrays to overwrite integer lookup tables

      integer ifld_int(len_look_user)  ! int lookup fields to change
      integer item_int(len_look_user)  ! item number to change
      integer ival_int(len_look_user)  ! integer value to use

      ! arrays to overwrite real lookup tables

      integer ifld_real(len_look_user)  ! lookup fields to change
      integer item_real(len_look_user)  ! item number to change
      real    rval_real(len_look_user)  ! real value to use

      integer fixhd(lfh)          ! fixed length header

      integer int_const(max_len_intc)  ! integer constants
      real real_const(max_len_realc)   ! real constants

      real lev_dep_consts(1+max_len1_levdepc*max_len2_levdepc)
      real row_dep_consts(1+max_len1_rowdepc*max_len2_rowdepc)
      real col_dep_consts(1+max_len1_coldepc*max_len2_coldepc)
      real extra_const(max_len_extcnst)
      real dumphist(len_dumphist)

! Local dynamic arrays:

      integer pp_int(45)
      integer lookup(45,len2_lookup_max)   ! Integer part of lookup
                                           ! table array

      real pp_real(19)
      real rlookup(46:64,len2_lookup_max)  ! Integer part of lookup
                                           ! table array

      integer lookup_all(len1_lookup_all,len2_lookup_max)
                                    ! Whole lookup table array

      integer cfi1(len_cfi(1))   ! compressed field index
      integer cfi2(len_cfi(2))   ! arrays
      integer cfi3(len_cfi(3))

      integer n_pp_flds(max_n_pp_files)  ! Number of pp fields array

      logical lsmask(len1_coldepc*len1_rowdepc) ! land sea mask
                                                ! for wave dump
      real fields_const(len1_flddepc,len2_flddepc)
                                  ! array for fields of constants
      real dummy(1)
      
      INTEGER               :: ErrorStatus  
                                 ! Return code : 0 Normal Exit : >0 Error
      
! Function & Subroutine calls:
      integer FIND_NAMELIST

!- End of header

      namelist /header_data/ fixhd,int_const,real_const,                &
                             lev_dep_consts,row_dep_consts,             &
                             col_dep_consts,extra_const,                &
                             ifld_int, item_int, ival_int,              &
                             ifld_real, item_real, rval_real

!L 0. Preliminaries

!L 0.1 Check sizes namelist dimensions do not exceed the maximum
!L dimensions and intialise arrays.


       if (len_intc >  max_len_intc) then
         write (6,*) ' len_intc in namelist is too big.'
         write (6,*) ' Max value allowed is ',max_len_intc
         write (6,*) ' Increase MAX_LEN_INTC in program'
         icode = 1
         go to 9999   !  Return
       endif

       if (len_realc >  max_len_realc) then
         write (6,*) ' len_realc in namelist is too big.'
         write (6,*) ' Max value allowed is ',max_len_realc
         write (6,*) ' Increase MAX_LEN_REALC in program'
         icode = 2
         go to 9999   !  Return
       endif

       if (len1_levdepc >  max_len1_levdepc) then
         write (6,*) ' len1_levdpec in namelist is too big.'
         write (6,*) ' Max value allowed is ',max_len1_levdepc
         write (6,*) ' Increase MAX_LEN1_LEVDEPC in program'
         icode = 3
         go to 9999   !  Return
       endif

       if (len2_levdepc >  max_len2_levdepc) then
         write (6,*) ' len2_levdpec in namelist is too big.'
         write (6,*) ' Max value allowed is ',max_len2_levdepc
         write (6,*) ' Increase MAX_LEN2_LEVDEPC in program'
         icode = 4
         go to 9999   !  Return
       endif

      if (len1_rowdepc >  max_len1_rowdepc) then
         write (6,*) ' len1_rowdpec in namelist is too big.'
         write (6,*) ' Max value allowed is ',max_len1_rowdepc
         write (6,*) ' Increase MAX_LEN1_ROWDEPC in program'
         icode = 5
      endif

      if (len2_rowdepc >  max_len2_rowdepc) then
         write (6,*) ' len2_rowdpec in namelist is too big.'
         write (6,*) ' Max value allowed is ',max_len2_rowdepc
         write (6,*) ' Increase MAX_LEN2_ROWDEPC in program'
         icode = 6
         go to 9999   !  Return
       endif

      if (len1_coldepc >  max_len1_coldepc) then
         write (6,*) ' len1_coldpec in namelist is too big.'
         write (6,*) ' Max value allowed is ',max_len1_coldepc
         write (6,*) ' Increase MAX_LEN1_COLDEPC in program'
         icode = 7
         go to 9999   !  Return
       endif

      if (len2_coldepc >  max_len2_coldepc) then
         write (6,*) ' len2_coldepc in namelist is too big.'
         write (6,*) ' Max value allowed is ',max_len2_coldepc
         write (6,*) ' Increase MAX_LEN2_COLDEPC in program'
         icode = 8
         go to 9999   !  Return
       endif

      if (len_extcnst >  max_len_extcnst) then
         write (6,*) ' len_extcnst in namelist is too big.'
         write (6,*) ' Max value allowed is ',max_len_extcnst
         write (6,*) ' Increase MAX_LEN_EXTCNST in program'
         icode = 11
         go to 9999   !  Return
      endif

      ! Initialise namelist arrays.
      do n=1,max_len2_levdepc
        do m=1,max_len1_levdepc
          lev_dep_consts(m+(n-1)*max_len1_levdepc)=0.0
        end do
      end do

      do n=1,max_len2_rowdepc
        do m=1,max_len1_rowdepc
          row_dep_consts(m+(n-1)*max_len1_rowdepc)=0.0
        end do
      end do

      do n=1,max_len2_coldepc
        do m=1,max_len1_coldepc
          col_dep_consts(m+(n-1)*max_len1_coldepc)=0.0
        end do
      end do

      do n=1,max_len_extcnst
        extra_const(n)=0.0
      end do

      do n=1,len_dumphist
        dumphist(n)=0.0
      end do

!L
!L 0.2 Read wave land/sea mask
!L
      if (wave .and. compress) then

        write(6,*) 'reading in landsea mask for waves from pp dataset'

        call get_file(ftin2,LEVELS,filenamelength,icode)
        levels=trim(levels)
        OPEN( unit=ftin2, file=levels, form="unformatted" )

! DEPENDS ON: read_pp_header
        CALL read_pp_header(ftin2,pp_int,pp_real,ibm_to_cray,           &
                            l_bit_32)

! DEPENDS ON: readdata
        CALL readdata( len1_rowdepc,len1_coldepc,ftin2,ibm_to_cray,     &
                       0, l_bit_32, .FALSE., lsmask, dummy)

        close(ftin2)

! reset so sea points are true
!
        n_sea_points=0

        do i=1,len1_coldepc*len1_rowdepc

         lsmask(i)=.not. lsmask(i)

         if(lsmask(i)) then
           n_sea_points=n_sea_points+1
         end if

        enddo

        print*,'n_sea_points set to ', n_sea_points
        fldsizelev(1)=n_sea_points

      endif     ! wave .and. compress

!L 0.3 Calculate number of fields (len2_lookup), which depends on
!L     nlevs_code.

      len2_lookup = 0
      len2_step = 0

      do i =1,fld_types
          len2_step = n_times * nlevs_code(i)
          len2_lookup = len2_lookup + len2_step
      end do

      print*,' '
      print*,'len2_lookup = ',len2_lookup

!L 0.4 Initialise arrays

      do n=1,max_n_pp_files
        n_pp_flds(n)=0
      enddo

      do n = 1, len_look_user
        ifld_int(n) = imdi
        item_int(n) = imdi
        ival_int(n) = imdi
        ifld_real(n) = imdi
        item_real(n) = imdi
        rval_real(n) = rmdi
      end do

      do i = 1,len_dumphist
        dumphist(i)= 0.0
      end do

!L 0.5 Read StashMaster files

      IROW_NUMBER=0
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_A',IROW_NUMBER,                     &
        ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_O',IROW_NUMBER,                     &
        ICODE,CMESSAGE)

!L 1. Read headers of all input data

!L 1.0 Open the UM/ancillary file

      call get_file(ftout,ANCFILE,filenamelength,icode)
      call file_open (ftout,ANCFILE,filenamelength,1,1,icode)
      if (icode >  0) then
        write (6,*) ' Problem with opening Ancillary File on Unit, '
        ICODE = 2
        go to 9999  ! Return
      endif

! note these values are used in pp_table so need to be set here
! before the do-loops
!     * set default lev_dep_consts for WAVE frequency using the
!       factor 1.1
!     * as set in real_const(15) by namelist
!     * need to set frmin as lev_dep_consts(1) in namelist input
!     * pick up the CO value from namelist input HEADER value

      if (wave) then

        rewind(5)
! DEPENDS ON: find_namelist
        I=FIND_NAMELIST(5,"HEADER_DATA")

        If(I == 0)then
          READ (UNIT=5, NML=HEADER_DATA, IOSTAT=ErrorStatus)
          CALL check_iostat(errorstatus, "namelist HEADER_DATA")
        Else
          write(6,*)'Cannot find namelist HEADER_DATA'
        End if

        print*,'real_const(15) is',real_const(15)
        print*,'lev_dep_consts(1,1)=',lev_dep_consts(1)

         do m=2,len1_levdepc
          lev_dep_consts(m) = real_const(15)*lev_dep_consts(m-1)
         enddo

      endif

!L 1.1 Read through all data sets calculating the ancillary
!L    file headers.  Loop over n_times then fld_types.

      runtot=1  ! points to start point in data array for next field
      fieldn=0  ! field number counter

      do n=1,n_times

      do m=1,fld_types

      fieldn = fieldn + 1

!L 1.2 do steps which are independent of the level first
!L Read the pp header for each field

      if (n_pp_files == 1) then
        ftin1= unit_no(1)
      else
        if (field_order) then
          ftin1= unit_no(n)
        else
          ftin1= unit_no(m)
        endif
      endif

! DEPENDS ON: read_pp_header
      call read_pp_header (ftin1,pp_int,pp_real,ibm_to_cray,l_bit_32)

      n_pp_flds(ftin1-29) = n_pp_flds(ftin1-29)+1
      write (6,*) 'Field No ',fieldn,' read in. From PP File ',ftin1-29,&
                  ' Field No ',n_pp_flds(ftin1-29)


      if (pphead) then
        write (6,*) 'pp_int for field ',fieldn
        write (6,*) (pp_int(j),j=1,45)
        write (6,*) 'pp_real for field ',fieldn
        write (6,*) (pp_real(j),j=1,19)
      endif


!L 1.3 Extract the data dimensions and determine whether field
!L      is on tracer or velocity grid and how many levels this
!L      field has.

      rows          = pp_int(lbrow)
      columns       = pp_int(lbnpt)
      np            = pp_int(lblrec)
      len_extra     = MAX(pp_int(lbext), 0)  ! extra data in pp-field
      pp_int(lbext )= 0                      ! get rid of extra data

      do i = 1, number_of_codes
        if ( pp_int(item_code)  ==  stash_code(i) ) then

!L        Get grid and number of levels for this stash code
          tracer_grid = grid_of_tracer(i)
          nlevs_this_code = nlevs_code(i)
          write (6,*) ' PP Code ',pp_int(lbfc),' tracer_grid ',         &
          tracer_grid,' nlevs_this_code ',nlevs_this_code

!L        Check field code set ; if not, set from FIELD_CODE
          if (pp_int(lbfc) /= field_code(i)) then
            write (6,*) 'Field No ',fieldn,' Field code',               &
            ' incorrect or not set. Reset from ',pp_int(lbfc),          &
            ' to ',field_code(i)
            pp_int(lbfc) =  field_code(i)
          endif

          go to 8
        end if
      end do

      write (6,*) ' WARNING from subroutine ANC_FLD'
      write (6,*) ' Stash code ', pp_int(item_code),' in PP Header ',   &
        ' was not found in STASH_CODE of CODES namelist.'

8     continue

!L 1.4 Start loop over levels

      do levn= 1, nlevs_this_code

        if(levn  /=  1) then

          fieldn = fieldn + 1

! DEPENDS ON: read_pp_header
          call read_pp_header(ftin1,pp_int,pp_real,ibm_to_cray,l_bit_32)

          n_pp_flds(ftin1-29) = n_pp_flds(ftin1-29)+1

          write (6,*) 'Field No ',fieldn,' read in. From PP File '      &
         ,ftin1-29, ' Field No ',n_pp_flds(ftin1-29)

            if (pphead) then
             write (6,*) 'pp_int for field ',fieldn
             write (6,*) (pp_int(j),j=1,45)
             write (6,*) 'pp_real for field ',fieldn
             write (6,*) (pp_real(j),j=1,19)
            endif

        end if   ! levn  /=  1

!L 1.5 Set t_compress depending on wave and nlevs_this_code.
!L     Don't compress fields of only one level.

      if (compress .and. nlevs_this_code  /=  1) then

         t_compress = .true.

      elseif (wave .and. pp_int(lbfc) == 38) then

        print*,'re-setting compress false for lsmask'
        print*,'in ppheader section of anc_fld'
        t_compress=.false.

      elseif (wave .and.  pp_int(lbfc) /= 38                            &
                         .and. nlevs_this_code  ==  1) then

        t_compress = .true.

      else

         t_compress = .false.

      endif

!L 1.6 Calculate fieldsize

      if (add_wrap_pts) then
       if (t_compress) then
        if(.not.wave) then
          fieldsize= fldsizelev(levn)
        else
          fieldsize=n_sea_points
          if(pp_int(lbfc) == 38) then  ! LS mask data field
           print*,'fieldsize set uncomp for lsmask'
           fieldsize=rows*columns
          endif
        endif
       else
         fieldsize=rows*(columns+2)
       endif
      else
        if (t_compress) then
         if (.not. wave) then
           fieldsize= fldsizelev(levn)
         else
           fieldsize=n_sea_points
           if(pp_int(lbfc) == 38) then  ! LS mask data field
            print*,'fieldsize set uncomp for lsmask'
            fieldsize=rows*columns
           endif
         endif
        else
         fieldsize=rows*columns
        endif
       endif

!L 1.7 Set the fixed header, integer and real constants

      if (fieldn == 1) then

       icode = 0

!L Calculate no_cmp

       no_cmp = 0
       do i = 1, nlevels   ! do not use levn in this loop
          no_cmp = no_cmp + fldsizelev(i)
       end do

!
! note - anc_head only uses compress if .not. wave
!
! DEPENDS ON: anc_head
      CALL anc_head(pp_int,pp_real,rows,columns,fieldsize,len2_lookup,  &
       fld_types,n_times,nlevels,n_freq_waves,n_dir_waves,no_cmp,       &
       len1_levdepc,len2_levdepc,len1_rowdepc,len2_rowdepc,len1_coldepc,&
       len2_coldepc,len1_flddepc,len2_flddepc,len_extcnst,len_cfi,      &
       tracer_grid,add_wrap_pts,periodic,single_time,ibm_to_cray,       &
       t_compress,levdepc,rowdepc,coldepc,flddepc,extcnst,wave,         &
       lfh,fixhd,len_intc,int_const,len_realc,real_const,icode)

      if (icode >  0) go to 9999  !  Error detected ; Return

! initialise header for length of data in ancillary file
! (reset to zero from mdi value set in anc_head)
      fixhd(161) = 0

      end if    ! fieldn  ==  1

! accumulate indicator of length of data in ancillary file
      fixhd(161)=fixhd(161)+fieldsize

!L 1.8 Set the lookup table for this field

!CMH note for waves - to use frequency information in lev-dep-consts
!CMH need to set before call pp_table as well as in proper place.
!CMH so do before the loops

! DEPENDS ON: pp_table
      CALL pp_table(pp_int,pp_real,len2_lookup,lookup,rlookup,          &
       fieldsize,fieldn,levn,m,runtot,number_of_codes,field_code,       &
       stash_code,add_wrap_pts,t_compress,pack32,wave,len1_levdepc,     &
       len2_levdepc,lev_dep_consts,len_realc,real_const,icode)

      if (icode >  0) go to 9999  !  Error detected ; Return

!L 1.9 Read past the data part of this pp field
! DEPENDS ON: readdata
      call readdata(rows,columns,ftin1,ibm_to_cray,len_extra,l_bit_32,  &
     &              .TRUE., dummy, dummy)

      END DO ! end of loop over levels

      END DO ! end of loop over fld_types

      END DO ! end of loop over times

      write(6,*) '==================================='
      write(6,*) fieldn,' PP fields have been read in'
      write(6,*) '==================================='

!L 2.  If flddepc=t or compress = t:  Read in fields of constant and
!L     compressed field indices from levels dataset

!L 2.1 For ocean dumps, create compressed field indices and
!L     fields_const

      icode = 0

      if ((compress .or. flddepc) .and. .not. wave) then

! DEPENDS ON: calc_cfi_and_fld
        CALL calc_cfi_and_fld(ftin2,nlevels,len1_coldepc,               &
          cols_nowrap,len1_rowdepc,len1_flddepc,len2_flddepc,           &
          fields_const,fldsizelev,len_cfi,cfi1,cfi2,cfi3,compress,      &
          flddepc,ibm_to_cray,add_wrap_pts,imdi,l_bit_32,icode)

        if (icode  /=  0) then
         go to 9999
        endif

      end if

!L
!L 3. Over-ride elements in header arrays
!L
!      Arrays that can be over-ridden are fixed length, integer
!      constants, real constants and level dependent constants.

!     * set default lev_dep_consts for WAVE frequency using the
!       factor 1.1
!     * as set in real_const(15) by namelist
!     * need to set frmin as lev_dep_consts(1) in namelist input
!     * pick up the CO value from namelist input HEADER value

!L 3.1 Read in the header_data namelist which includes the
!L    lev_dep_consts,row_dep_consts,col_dep_consts,extra_consts

      rewind(5)
! DEPENDS ON: find_namelist
      I=FIND_NAMELIST(5,"HEADER_DATA")

      If(I == 0)then
        READ (UNIT=5, NML=HEADER_DATA, IOSTAT=ErrorStatus)
        CALL check_iostat(errorstatus, "namelist HEADER_DATA")
      Else
        write(6,*)'Cannot find namelist HEADER_DATA'
      End if
      write (6,*) ' '

!L For wave dumps calculate the lev_dep_consts again

      if (wave) then

        do m=2,len1_levdepc
!CC      FR(M) = CO*FR(M-1)
        lev_dep_consts(m) = real_const(15)*lev_dep_consts(m-1)
        enddo

      endif

!L 3.2 Amend the lookup tables

      do i = 1, len_look_user

        if ( ifld_int(i)  /=  imdi ) then

          if ( ifld_int(i)  ==  0 ) then
            do j = 1, len2_lookup
              lookup( item_int(i) , j ) = ival_int(i)
            end do
          else
            lookup( item_int(i) , ifld_int(i) ) = ival_int(i)
          end if

        end if    ! ifld_int(i)  /=  imdi

        if ( ifld_real(i)  /=  imdi ) then

          if ( ifld_real(i)  ==  0 ) then
            do j = 1, len2_lookup
              rlookup( item_real(i) , j ) = rval_real(i)
            end do
          else
            rlookup( item_real(i) , ifld_real(i) ) = rval_real(i)
          end if

        end if    ! ifld_real(i)  /=  imdi

      end do  ! i = 1, len_look_user

!L 3.3 Print out headers to screen

      if (pphead) then

      write(6,*) ' '
      write(6,*) 'fixhd'
      write(6,*) (fixhd(j),j=1,161)
      write(6,*) ' '
      write(6,*) 'int_const ; length ',len_intc
      write(6,*) (int_const(j),j=1,len_intc)
      write(6,*) ' '
      write(6,*) 'real_const ; length ',len_realc
      write(6,*) (real_const(j),j=1,len_realc)

        if (levdepc) then
         write(6,*) ' '
         write(6,*) 'level dependent constants '
          do j=1,len2_levdepc
           ipos=(j-1)*len1_levdepc
           write(6,*) 'variable ',j
           write(6,*) (lev_dep_consts(ipos+i),i=1,len1_levdepc)
          enddo
        endif
       write(6,*) ' '

        if (rowdepc) then
          write(6,*) ' '
          write(6,*) 'row dependent constants '
           do j=1,len2_rowdepc
            ipos=(j-1)*len1_rowdepc
            write(6,*) 'variable ',j
            write(6,*) (row_dep_consts(ipos+i),i=1,len1_rowdepc)
           enddo
        endif
       write(6,*) ' '

        if (coldepc) then
          write(6,*) ' '
          write(6,*) 'column dependent constants '
           do j=1,len2_coldepc
            ipos=(j-1)*len1_coldepc
            write(6,*) 'variable ',j
            write(6,*) (col_dep_consts(ipos+i),i=1,len1_coldepc)
           enddo
        endif
       write(6,*) ' '

        if (flddepc) then
          write(6,*) ' '
          write(6,*) 'fields constants'
          write(6,*)' len1_flddepc = ',len1_flddepc
          write(6,*)' len2_flddepc = ',len2_flddepc
        endif
       write(6,*) ' '

        if (extcnst) then
          write(6,*) ' '
          write(6,*) 'extra constants '
          write(6,*) (extra_const(i),i=1,len_extcnst)
        endif
       write(6,*) ' '

      endif
!L
!L 4. Write out header data to ancillary field file
!L

!L 4.1 Write out fixed, integer and real constants headers.
!L Set values for use in WRITHEAD and convert lookup and rlookup
!L into one array lookup_all using subroutine conv_real

      len_data=fixhd(161)

! DEPENDS ON: conv_real
      call conv_real(rlookup,lookup_all,len2_lookup)

      do i=1,len2_lookup
        lookup_all(1:45,i) = lookup(1:45,i)
      end do

! If logical lwfio (set in namelist LOGICALS) is true then set the
! LBEGIN and LBNREC fields in the LOOKUP Headers for VN 16 Type
! Dumpfiles.
      if (lwfio) then

! DEPENDS ON: set_dumpfile_address
         CALL set_dumpfile_address(fixhd,lfh,                           &
          lookup_all,len1_lookup_all,len2_lookup,                       &
          number_of_data_words_in_memory,                               &
          number_of_data_words_on_disk,                                 &
          disk_address)

      end if

!L 4.2 Use WRITHEAD to write the headers and constants

       CALL WRITHEAD(ftout,fixhd,lfh,int_const,len_intc,                &
       real_const,len_realc,lev_dep_consts,len1_levdepc,len2_levdepc,   &
       row_dep_consts,len1_rowdepc,len2_rowdepc,col_dep_consts,         &
       len1_coldepc,len2_coldepc,fields_const,len1_flddepc,             &
       len2_flddepc,extra_const,len_extcnst,dumphist,len_dumphist,      &
       cfi1,len_cfi(1),cfi2,len_cfi(2),cfi3,len_cfi(3),lookup_all,      &
       len1_lookup_all,len2_lookup,len_data,                            &
       umFortranIntegerSize()*8,                                        &
       start_block,icode,cmessage)

!L 5.0  Write out (rest of) data to ancillary file

!L 5.1  Return to start of pp input pp fields files

      write (6,*) ' '
      do n=1,n_pp_files
        rewind 29+n
        write (6,*) ' Rewinding PP file on Unit No ',29+n
      end do
      write (6,*) ' '

!L 5.2 Start loop over fields

      fieldn=0   ! field number counter

      do n=1,n_times

      do m=1,fld_types

!L 5.3 Do steps which are independant of level first

      if (n_pp_files == 1) then
        ftin1= unit_no(1)
      else
        if (field_order) then
          ftin1= unit_no(n)
        else
          ftin1= unit_no(m)
        endif
      endif

      fieldn = fieldn + 1

!L 5.4 Read pp header and determine length of field to
!L       be output to ancillary file

! DEPENDS ON: read_pp_header
      CALL read_pp_header (ftin1,pp_int,pp_real,ibm_to_cray,l_bit_32)

! (* extract the data dimensions and tracer/velocity grid type *)

      rows        = pp_int(lbrow)
      columns     = pp_int(lbnpt)
      len_extra   = MAX(pp_int(lbext), 0)  ! extra data in pp-field

      do i = 1, number_of_codes
        if ( pp_int(lbfc)  ==  field_code(i) ) then
          tracer_grid = grid_of_tracer(i)
          nlevs_this_code=nlevs_code(i)
          go to 30
        end if
      end do

30    continue

!L 5.5 Start loop over levels

      do levn = 1, nlevs_this_code

      if (levn  /=  1) then
        fieldn = fieldn + 1
! DEPENDS ON: read_pp_header
        CALL read_pp_header(ftin1,pp_int,pp_real,ibm_to_cray,l_bit_32)
      end if

!L 5.6 Set t_compress which depends on wave and nlevs_this_code

      if (compress .and. nlevs_this_code  /=  1) then

         t_compress = .true.

      elseif (wave .and. pp_int(lbfc) == 38) then

         print*,'re-setting compress false for lsmask'
         print*,'in ppheader section of anc_fld'
         t_compress = .false.

      else

         t_compress = .false.

      endif

!L 5.7 Calculate fieldsize

      if (add_wrap_pts) then
       if (t_compress) then
        if(.not.wave) then
          fieldsize= fldsizelev(levn)
        else
          fieldsize=n_sea_points
          if(pp_int(lbfc) == 38) then  ! LS mask data field
           print*,'fieldsize set uncomp for lsmask'
           fieldsize=rows*columns
          endif
        endif
       else
         fieldsize=rows*(columns+2)
       endif
      else
        if (t_compress) then
         if (.not. wave) then
           fieldsize= fldsizelev(levn)
         else
           fieldsize=n_sea_points
           if(pp_int(lbfc) == 38) then  ! LS mask data field
            print*,'fieldsize set uncomp for lsmask'
            fieldsize=rows*columns
           endif
         endif
        else
         fieldsize=rows*columns
        endif

      endif

!L 5.8 Call subroutine data to write the fields to the dump/ancillary
!L     file.

! DEPENDS ON: dataw
      CALL dataw(rows,columns,fieldsize,nlevels,levn,len_extra,         &
       fieldn,len1_lookup_all,lookup_all,fixhd,                         &
       len_cfi, cfi1, cfi2, cfi3,fldsizelev,ftin1,ftout,                &
       tracer_grid,add_wrap_pts,ibm_to_cray,t_compress,rmdi_input,      &
       wave,lsmask, l_bit_32,                                           &
       icode)

    END DO
    END DO
    END DO

      write (6,*) '========================================'
      write (6,*) fieldn,' fields written to Ancillary File'
      write (6,*) '========================================'

9999  continue
      return
      END SUBROUTINE anc_fld
