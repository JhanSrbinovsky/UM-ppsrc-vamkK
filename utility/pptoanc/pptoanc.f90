! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      PROGRAM PPTOANC
      USE IO
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE check_iostat_mod    
      USE ereport_mod, ONLY : ereport, ereport_finalise
      USE UM_Config, ONLY : &
          appInit, &
          exe_pptoanc
      USE lookup_addresses
      USE Submodel_Mod

      IMPLICIT NONE
!
! Routine: pptoanc -------------------------------------------------
!
! Description:
!  To create ancillary fields from pp fields.
!  pp fields are output in same order as they are input
!
! Method:
!
!     unit                    Description
!    ftin1=30 onwards  INPUT  pp files (unit #s provided by user)
!    ftin2=11          INPUT  levels dataset (only used if compress=t
!                             or flddepc=t)
!    ftout=10          OUTPUT ancillary file
!
! Use  namelists to set:
!
!
!      sizes          field_types,n_times,n_levels,n_pp_files,
!                     n_freq_waves,n_dir_waves,stash_code,field_code,
!                     nlevs_code,unit_no,len_intc,len_realc,
!                     len1_levdepc,len2_levdepc,len1_rowdepc,
!                     len2_rowdepc,len1_coldepc,len2_coldepc,
!                     len1_flddepc,len2_flddepc,len_extcnst,rmdi_input
!
!      logicals       add_wrap_pts,periodic,single_time,ibm_to_cray,
!                     compress,wave,levdepc,rowdepc,coldepc,
!                     flddepc,extcnst,pack32,pphead,grid_of_tracer,
!                     field_order
!
!      first_vt       fvhh,fvdd,fvmm,fvyy
!                     (first validity time)
!
!      last_vt        lvhh,lvdd,lvmm,lvyy
!                     (last validity time)
!
!      interval       year360,ivhh,ivdd,ivmm,ivyy
!                     (time interval between validity times)
!
!      header_data    fixhd, int_const,real_const,lev_dep_consts,
!                     row_dep_consts,col_dep_consts,extra_const,
!                     ifld_int, item_int, ival_int,
!                     ifld_real, item_real, rval_real
!
!  The last 6 variables are arrays allowing up to n= len_look_user
!  changes to the lookup tables. All arrays are initiated as missing
!  data. Changes are made in the order from n=1 to len_look_user. If
!  ifld_int(n) or ifld_real(n) is not a missing data value, changes
!  are made to the integer or real lookup tables.
!
!                for integer parts of lookup tables
!                     ifld_int(n)   field number to change
!                                   0 => all lookup tables
!                     item_int(n)   item number to change
!                     ival_int(n)   integer value to use
!
!                for real parts of lookup tables
!                     ifld_real(n)   field number to change
!                                   0 => all lookup tables
!                     item_real(n)   item number to change
!                     rval_real(n)   real value to use
!
!    The following elements are particularly worth checking
!                     fixhd(4)    grid type code (global is default)
!                     fixhd(8)  360 day calendar is default
!                     fixhd(12)   UM version number
!
!                     fixhd and lookup for dates of validity
!
!  The stash_code and field_code in the namelist sizes must be in the
!  same order as they are in the input pp fields
!
!  Do not put field types requiring different compression indices
!  in one ancillary file.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
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

! Routine arguments
!   Scalar arguments

      integer n_stash_codes    ,                                        &
                                 ! counter for number of stash codes
     &        n_unit_no        ,                                        &
                                 ! counter for number of unit numbers
     &        len2_lookup_max  ,                                        &
                                 ! 2nd dimension for lookup array
                                 ! in ancfld(maximum)
     &        cols_nowrap      ,                                        &
                                 ! no. of columns east-west without
                                 ! wrap_points
     &        n,i              ,                                        &
                                 ! loop counter
     &        icode
                                 ! error exit condition code

      !  Define variables from SIZES namelist

      integer field_types   ,                                           &
                               ! number of field types in I/O files
     &        n_times       ,                                           &
                               ! number of time periods in I/O files
     &        nlevels       ,                                           &
                               ! number of levels (default = 1)
     &        n_pp_files    ,                                           &
                               ! number of input pp files
     &        n_freq_waves  ,                                           &
                               ! number of wave frequencies
     &        n_dir_waves   ,                                           &
                               ! number of wave directions
     &        len_intc      ,                                           &
                               ! dimension for integer constants
     &        len_realc     ,                                           &
                               ! dimension for real constants
     &        len_extra     ,                                           &
                               ! dimension for extra data
     &        len1_levdepc  ,                                           &
                               ! dimension for levdepc array
     &        len2_levdepc  ,                                           &
                               ! 2nd dimension for levdepc array
     &        len1_rowdepc  ,                                           &
                               ! dimension for rowdepc array
     &        len2_rowdepc  ,                                           &
                               ! 2nd dimension for rowdepc array
     &        len1_coldepc  ,                                           &
                               ! dimension for coldepc array
     &        len2_coldepc  ,                                           &
                               ! 2nd dimension for coldepc array
     &        len1_flddepc  ,                                           &
                               ! dimension for flddepc array
     &        len2_flddepc  ,                                           &
                               ! 2nd dimension for flddepc array
     &        len_extcnst      ! dimension for extcnst array

      real rmdi_input          ! real missing data indicator
                               ! in input pp field

      !  Define variables from LOGICALS namelist

      logical add_wrap_pts ,                                            &
                             ! T => adds wrapping columns
                             !      e.g. for global grid
     &        periodic     ,                                            &
                             ! T => periodic in time
                             !      e.g. climate field
     &        single_time  ,                                            &
                             ! T => all fields input valid at one time
     &        ibm_to_cray  ,                                            &
                             ! T => input pp data is in IBM number
                             !      format and needs to be converted to
                             !      run on the Cray.
                             !      (Only use if running on Cray)
     &        compress     ,                                            &
                             ! T => fields are packed into ancillary
                             !      field compressed field indices are
                             !      calculated
     &        wave         ,                                            &
                             ! T => a wave dump is to be created
     &        levdepc      ,                                            &
                             ! T => if level dependent constants array
                             !      required
     &        rowdepc      ,                                            &
                             ! T => if row dependant constants are
                             !      required
     &        coldepc      ,                                            &
                             ! T => if column dependant constants are
                             !      required
     &        flddepc      ,                                            &
                             ! T => if fields of constants are
                             !      required
     &        extcnst      ,                                            &
                             ! T => if fields of constants are
                             !      required

     &        pack32       ,                                            &
                             ! T => use 32 bit Cray numbers
     &        pphead       ,                                            &
                             ! T => print out pp headers read in

     &        field_order ,                                             &
                             ! T => input pp fields ordered by time.
                             !      i.e. different months in input
                             !        files, same fields in all files
                             ! F => inout pp fields ordered by fields.
                             !      i.e. different fields in input
                             !        files, all months in all files

     &        lwfio          ! T => set the LBEGIN and LBNREC fields
                             !      in the LOOKUP Headers for VN 16
                             !      Type Dumpfiles.
                             ! F => Old dumpfiles

      CHARACTER(LEN=filenamelength) :: namelst
      CHARACTER(LEN=80) cmessage

! Parameters:
      integer ftin2                  ! input unit for mask file used
      parameter (ftin2=11)           ! for level dependent consts and
                                     ! compression indices.
                                     ! Only used when compress is T
      integer ftout
      parameter (ftout=10)           ! unit number for output ancillary
                                     ! file
      integer nolevsmax
      parameter (nolevsmax=200)      ! max number of levels; dimensions
                                     !  fldsizelev array
      integer number_of_codes
      parameter (number_of_codes=200)! max number of stash/field codes

      integer max_n_pp_files
      parameter (max_n_pp_files=50)  ! max number of input pp files

      INTEGER, PARAMETER :: Max_Filename_Len = 256

! Array arguments:

      integer len_cfi(3)  ,                                             &
                                 ! lengths of compressed field indices
     &        fldsizelev(nolevsmax)  ! size of packed field
                                     ! on each level

      logical grid_of_tracer(number_of_codes) ! T => fields are on a
                                              ! tracer grid

      ! Define variables from SIZES namelist

      integer stash_code(number_of_codes),                              &
                                          ! array of stash codes
     &        field_code(number_of_codes),                              &
                                          ! array of field codes
     &        nlevs_code(number_of_codes),                              &
                                          ! array of levels depending
                                          ! on field code
     &        unit_no(number_of_codes)    ! array of unit numbers for
                                          ! input

      CHARACTER(LEN=Max_Filename_Len) :: pp_file_name
      CHARACTER(LEN=1)                :: c_bit_32
      INTEGER                         :: err
      INTEGER                         :: errorstatus
      INTEGER                         :: bit_32
      LOGICAL                         :: L_bit_32   ! is this 32 bit?
      INTEGER ME_GC,NPROC_GC
! Function & Subroutine calls:
      integer FIND_NAMELIST

!- End of header

      namelist /sizes/ field_types,n_times,nlevels,n_pp_files,          &
     & n_freq_waves,n_dir_waves,stash_code,field_code,nlevs_code,       &
     & unit_no,len_intc,len_realc,len1_levdepc,len2_levdepc,            &
     & len1_rowdepc,len2_rowdepc,len1_coldepc,len2_coldepc,             &
     & len1_flddepc,len2_flddepc,len_extcnst,rmdi_input

      namelist /logicals/ add_wrap_pts,periodic,single_time,            &
     &  ibm_to_cray,compress,wave,levdepc,rowdepc,coldepc,flddepc,      &
     &  extcnst,pack32,pphead,grid_of_tracer,field_order,lwfio


!L 0 Initialise GC variables
      CALL GC_INIT(' ',ME_GC,NPROC_GC)
      CALL appInit(exe_pptoanc)
      CALL ioInit()

!L 1 Set values

!L 1.0 Set default values for SIZES NAMELIST


      field_types  = 2
      n_times      = 12
      nlevels      = 1
      n_pp_files   = 1
      n_freq_waves = 1
      n_dir_waves  = 1
      len_intc     = 40
      len_realc    = 40
      len1_levdepc = 1
      len2_levdepc = 1
      len1_rowdepc = 1
      len2_rowdepc = 1
      len1_coldepc = 1
      len2_coldepc = 1
      len1_flddepc = 1
      len2_flddepc = 1
      len_extcnst  = 1

      rmdi_input   = rmdi

!L 1.1 Initialise arrays in SIZES NAMELIST

      do n=1,number_of_codes
        field_code(n)=-99
        stash_code(n)=-99
        nlevs_code(n)=1
        unit_no(n)=-99
      enddo

!L 1.2 Open UNIT05 containing namelists and read in SIZES NAMELIST

      call get_file(5,namelst,filenamelength,icode)
      open (unit=5,file=namelst)

      rewind(5)
! DEPENDS ON: find_namelist
      I=FIND_NAMELIST(5,"SIZES")

      If(I == 0)then
        READ (UNIT=5, NML=SIZES, IOSTAT=ErrorStatus)
        CALL check_iostat(errorstatus, "namelist SIZES")
      Else
        write(6,*)'Cannot find namelist SIZES'
      End if

      write (6,*) ' '
      write (6,*) 'SIZES namelist is set up as follows:-'
      write (6,*) ' '
      write (6,sizes)

!L 1.3 Check that n_pp_files is not greater than max_n_pp_files

      if (n_pp_files <= 0 .or. n_pp_files >  max_n_pp_files) then
        write (6,*) ' '
        write (6,*) ' N_PP_FILES must in range 1-',MAX_N_PP_FILES
        write (6,*) ' N_PP_FILES must in range 1-',number_of_codes
        write (6,*) ' Resubmit job with new value for N_PP_FILES'
        go to 9999   !  Return
      endif

!L 1.4 Check that n_times and number of field_types is not greater
!L     than number_of_codes

      if (n_times  >   number_of_codes   .or.                           &
     &   field_types  >   number_of_codes ) then
        write (6,*) ' '
        write (6,*) ' ** WARNING ** WARNING ** '
        write (6,*) ' N_TIMES = ',n_times,' or FIELD_TYPES = ',         &
     &  field_types,' greater than NUMBER_OF_CODES = ',number_of_codes
        write (6,*) ' Dimension of UNIT_NO may be too small if used.'
      endif

!L 1.5 Count the number of stash codes and check they are not
!L     greater than number of field_types

      n_stash_codes = 0
      do n=1,number_of_codes
        if (stash_code(n) >= 0) then
          n_stash_codes = n_stash_codes + 1
        endif
      enddo

      if (n_stash_codes /= field_types) then
        write (6,*) ' '
        write (6,*) ' Wrong number of stash codes provided.'
        write (6,*) n_stash_codes,' stash codes in namelist.'
        write (6,*) field_types  ,' stash codes expected.'
        write (6,*) ' Rerun with correct no of stash codes'
        go to 9999   !  Return
      else
        write (6,*) ' '
        write (6,*) n_stash_codes,' stash codes in SIZES namelist.'
      endif

      if (nlevels  >   nolevsmax) then
         write(6,*) 'parameter nolevsmax is smaller than nlevels'
         write(6,*) 'increase nolevsmax in program create'
         go to 9999   !  Jump out
      end if

      if (rmdi_input  ==  rmdi) then
         write(6,*) 'rmdi_input should equal rmdi in input pp field'
         write(6,*) 'WARNING !!! '
         write(6,*) 'if not, RESUBMIT with the correct rmdi_input in    &
     & SIZES namelist.'
      end if

!L
!L 1.6 Set default values for LOGICALS NAMELIST
!L
      add_wrap_pts    = .false.
      periodic        = .false.
      single_time     = .false.
      ibm_to_cray     = .false.
      compress        = .false.
      wave            = .false.
      levdepc         = .false.
      rowdepc         = .false.
      coldepc         = .false.
      flddepc         = .false.
      extcnst         = .false.
      pack32          = .false.
      pphead          = .false.
      field_order     = .true.
      lwfio           = .true.

!L 1.7 Initialise array in LOGICAL NAMELIST

      do n = 1, number_of_codes
        grid_of_tracer(n)=.true.
      enddo

!L 1.8 Read in the LOGICALS NAMELIST

      rewind(5)
! DEPENDS ON: find_namelist
      I=FIND_NAMELIST(5,"LOGICALS")

      If(I == 0)then
         READ (UNIT=5, NML=LOGICALS, IOSTAT=ErrorStatus)
         CALL check_iostat(errorstatus, "namelist LOGICALS")
      Else
        write(6,*)'Cannot find namelist LOGICALS'
      End if

      write (6,*) ' '
      write (6,*) 'LOGICALS namelist is set up as follows:-'
      write (6,*) ' '
      write (6,logicals)

!L 1.9 Count number of unit numbers needed which depends on field_order,
!L     n_times and field_types

      n_unit_no = 0
      do n=1,number_of_codes
        if (unit_no(n) >  0) then
          n_unit_no = n_unit_no + 1
        endif
      enddo


      if (n_unit_no >  0) then
        do n=1,n_unit_no
          if (unit_no(n) <  30 .or. unit_no(n) >  29+n_pp_files) then
            write (6,*) ' '
            write (6,*) ' Unit no out of range in UNIT_NO :',unit_no(n)
            write (6,*) ' Range is 30-',29+n_pp_files
            write (6,*) ' Rerun with correct unit numbers'
            go to 9999   !  Return
          endif
        enddo
      else              ! n_unit_no
        do n=1,max_n_pp_files
          unit_no(n)=29+n
        enddo
      endif

!L 2 Set dimensions

!L 2.0 If data are to be compressed calculate the lengths of compression
!L    indices and number of points in field on each level using the
!L    levels dataset

      if (add_wrap_pts) then
        cols_nowrap = len1_coldepc-2
      else
        cols_nowrap = len1_coldepc
      endif

      if (compress .and. .not. wave)  then

! DEPENDS ON: calc_len_cfi
         call calc_len_cfi(ftin2,cols_nowrap,len1_rowdepc,              &
     &          nlevels,len_cfi,fldsizelev,ibm_to_cray,add_wrap_pts,    &
     &          l_bit_32,icode)

         if (icode  /=  0) then
         go to 9999        ! jump out
         end if

      else                 ! .not. compress

        len_cfi(1) = 1
        len_cfi(2) = 1
        len_cfi(3) = 1

      end if               ! compress

!L
!L 2.1 Calculate len2_lookup_max which depends on wave dimensions
!L
      icode = 0

      if (wave) then
        len2_lookup_max = field_types*n_times*nlevels                   &
     &   + (n_freq_waves*n_dir_waves -1)*n_times

      else
        len2_lookup_max = field_types*n_times*nlevels
      end if

      print*,'len2_lookup_max set to ',len2_lookup_max
      print*,' '
!L
!L 3 Read STASHmaster files

!L 3.1 Initialise N_INTERNAL_MODEL/INTERNAL_MODEL_INDEX

      N_INTERNAL_MODEL=2
      INTERNAL_MODEL_INDEX(1)=1    !  Atmos
      INTERNAL_MODEL_INDEX(2)=2    !  Ocean

!L 3.2 Determine ppxRecs from Stashmaster files

! DEPENDS ON: hdppxrf
      CALL HDPPXRF(22,'STASHmaster_A',ICODE,CMESSAGE)

!L 3.3 Read Ocean file and obtain number of records

! DEPENDS ON: hdppxrf
      CALL HDPPXRF(22,'STASHmaster_O',ICODE,CMESSAGE)

!  3.6 Find and open all pp files - check if 32 bit
!      data is required
      DO i = 1, n_pp_files
        CALL Get_File(29+i, pp_file_name, Max_Filename_Len, err)
        pp_file_name = TRIM(pp_file_name)
        OPEN(unit=29+i, file=pp_file_name, form="unformatted")
      END DO

      ! Get the "32 bit" flag from environment
      CALL Fort_Get_Env('BIT32', 5, c_bit_32, 1, err)
      c_bit_32 = TRIM(c_bit_32)
      READ(c_bit_32,'(i1)') bit_32
      IF (bit_32 == 1) THEN
        L_bit_32 = .TRUE.
      ELSE
        L_bit_32 = .FALSE.
      END IF


!L 4 Call main subroutine

! DEPENDS ON: anc_fld
      call anc_fld(ftin2,ftout,nolevsmax,number_of_codes,               &
     &  max_n_pp_files,len_cfi,fldsizelev,                              &
     &  field_types,n_times,nlevels,n_pp_files,stash_code,field_code,   &
     &  nlevs_code,unit_no,n_freq_waves,n_dir_waves,len_intc,len_realc, &
     &  len_extra,len1_levdepc,len2_levdepc,len1_rowdepc,len2_rowdepc,  &
     &  len1_coldepc,len2_coldepc,len1_flddepc,len2_flddepc,            &
     &  len_extcnst,rmdi_input,                                         &
     &  add_wrap_pts,periodic,single_time,ibm_to_cray,compress,wave,    &
     &  levdepc,rowdepc,coldepc,flddepc,extcnst,pack32,pphead,          &
     &  grid_of_tracer,field_order,lwfio,L_bit_32,                      &
!    #  len2_lookup_max,cols_nowrap,icode)
     &  len2_lookup_max,cols_nowrap,icode)

        if (icode  >   0) then
          go to 9999   !  Jump out
        end if

!L 5 Tidy up at end of program

!L 5.1 Close ancillary file and pp files

      call file_close (ftout,'ANCFILE',7,1,0,icode)
        if (icode  >   0) then
          write (6,*) ' Problem with FILE_CLOSE for unit no ',ftout
          go to 9999   !  Jump out
        end if

      DO i = i, n_pp_files
        CLOSE(i)
      END DO

      
!     ===========================================================
!L 5.2  NORMAL and ABNORMAL COMPLETION.
!     ===========================================================
9999  CONTINUE

      CALL ioShutdown()

      IF (icode == 0) THEN
        WRITE (6,'(A)') ' '
        WRITE (6,'(/,A)') 'Program completed normally'
        WRITE (6,'(A,I6,/)') 'Return code = ',icode
      ELSE
        WRITE (6,'(A)') 'PPTOANC Program'
        WRITE (6,'(A)') 'Return code has been set in program'
        WRITE (6,'(A,I6,/)') 'Return code = ',icode
        WRITE (6,'(A)') 'Program aborted'
        CALL EREPORT('PPTOANC', icode, 'Return code error')
      END IF


      CALL ereport_finalise( )
      END PROGRAM PPTOANC
