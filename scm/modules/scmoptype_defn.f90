! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Defines the derived-type, SCMop_type, necessary for declaring
!  SCMop - the structure that carries all the SCM diagnostic
!  information.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

MODULE scmoptype_defn

! Description:
! The NetCDF library routines were written to work in 32 bit and
! so we may have to call 32 bit routines from the 64 bit
! UM. Using the Intel compiler on the Linux platform this can
! easily be done by explicitly declaring the integers variables
! that will be inputs to or outputs from the NetCDF routines to
! be 32 bit, and the UM routines that call them
! (DUMP_STREAMS_INIT, DUMP_STREAMS and DUMP_STREAMS_END) can
! still be compiled in 64 bits.

! The NEC compiler however, will ignore any explicit requests for
! 32 bit integers in a routine that is being compiled at 64 bits
! - everything will be promoted to 64 bits. Thus on the NEC the UM
! routines that call the NetCDF library routines have to be compiled
! at 32 bits and have the inputs that they are recieving from the
! rest of the 64 bit UM explicitly declared as 64 bit. The upshot is
! that, in order to support both these methods, we will explictly
! declare the precision of every variable being passed to or
! from DUMP_STREAMS_INIT, DUMP_STREAMS, DUMP_STREAMS_END or any
! NetCDF library routine. This includes the contents of the
! SCMop structure, defined here.

USE um_types
IMPLICIT NONE

! Short-hands for integer64, real64 and logical64
  INTEGER, PARAMETER :: i64 = integer64
  INTEGER, PARAMETER :: r64 = real64
  INTEGER, PARAMETER :: l64 = logical64

! Precision of integers passed in to & out of NetCDF routines
  INTEGER, PARAMETER :: incdf = integer32

! Some limits
  INTEGER, PARAMETER :: maxndomprof = 99

! The lengths of certain character variables
  INTEGER, PARAMETER :: &
    lfilename = 100     &
  , lsname    = 30      &
  , llname    = 50      &
  , lunits    = 15

! If the filename of a stream is set to this the name will be
! ignored and the stream will be opened with the default name
! given the unit number.
  CHARACTER(len=lfilename), PARAMETER :: Default='<default>'

! The maximum no. of streams must be such that
! Stream(maxnstreams+1) gives an integer which is not
! too high for the machine precision
  INTEGER, PARAMETER :: maxnstreams  = 29
  INTEGER, PARAMETER :: inot_written = maxnstreams+1
  INTEGER, PARAMETER :: listlength   = 5000

!-----------------------------------------------------------------------
! A diagnostic's dump array has its own type so that we can
! declare an array of arrays and allocate memory to each
! separately
!-----------------------------------------------------------------------
  TYPE allocatable_array
    REAL(r64), POINTER :: dump(:)
  END TYPE allocatable_array

!-----------------------------------------------------------------------
! A stream has its own type to clearly separate its entries in
! SCMop from those corresponding to individual diagnostics
!-----------------------------------------------------------------------
  TYPE astream

   ! The unit the stream will write to
    INTEGER(i64) :: op_unit

   ! The name of the file it will create (if not set to default)
    CHARACTER (len=lfilename) :: filename
         ! The dumping period of the diagnostics sent to this stream
         ! (i.e. if a diagnostic is an average it is the number of
         ! timesteps it is averaged over, if it is an accumulation it
         ! is the number of timesteps it is accumulated over, if it is
         ! a maximum it is the number of timesteps it is the maximum
         ! over, etc.)

    INTEGER(i64) :: dump_step
         ! Flags whether diagnostics sent to this stream by the
         ! respective input to routine SCMoutput will actually be
         ! sent to this stream.

    INTEGER(i64) :: heed_hardwired
         ! Flags whether diagnostics sent to this stream by namelist
         ! request will actually be sent to this stream.

    INTEGER(i64) :: heed_acceptlist
         ! Flags whether diagnostics prevented from going to this
         ! stream by namelist request will actually be prevented.

    INTEGER(i64) :: heed_rejectlist

    CHARACTER (len=lsname), POINTER :: accept_list(:)
    CHARACTER (len=lsname), POINTER :: reject_list(:)

    INTEGER(i64) :: switch ! If zero, stream is not active.

    INTEGER(i64) :: FORMAT    ! Determines format of output file.
                              ! 0 = format intended for subsequent reading
                              !     by PV-wave routine scmread2.pro
                              !     (used by scmoutput.pro)
                              ! 1 = format geared to easy perusal by eye
                              ! 2 = format suitable for FSSI database
                              ! 3 = new PV-wave format designed to replace
                              !     format 0. Can be read by same routines.

    INTEGER(i64) :: n_output  ! The number of diagnostics that
                              ! will be output to this stream
  END TYPE astream


!-----------------------------------------------------------------------
! Define the derived type for SCMop. SCMop carries all
! necessary diagnostic information from the top(ish) level
! down to wherever any diagnostic is actually calculated,
! and then back up to the top for output.
!-----------------------------------------------------------------------
  TYPE SCMop_type

         ! Flags whether diagnostic system is "switched on"
    LOGICAL(l64) :: on

         ! first_pass will be true during all calls to SCMoutput in
         ! the first timestep (a formative stage for the list of
         ! diagnostics), and false thereafter (when the creation of
         ! new diagnostics will not be allowed).
    LOGICAL(l64) :: first_pass

         ! Pointers to daycount and stepcount, and knowledge of
         ! full_daysteps (all declared and defined in scm_main) are
         ! required so SCMop can tell the time.
    INTEGER(i64), POINTER :: daycount,stepcount

    INTEGER(i64) :: full_daysteps

         ! Knowledge of ntrad1 and ntrad (the first timestep on which
         ! radiation is called and the number of timesteps between calls
         ! thereafter) is required by SCMop to make sense of diagnostics
         ! only calculated on radiation timesteps.
    INTEGER(i64) :: ntrad1,ntrad

         ! An encoded integer representing which output streams are
         ! open (i.e. which streams will be output to file).
    INTEGER(i64) :: openstreams

         ! We want a certain number of streams, we will allocate
         ! exactly how many at runtime.
    TYPE(astream), POINTER :: strm(:)

         ! maxnentries will be the size of all the arrays associated
         ! with the diagnostic entries (sname, etc.) once allocated.
         ! It can be increased at runtime in routine expand_SCMop.
    INTEGER(i64) :: maxnentries

         ! nentries, n_output and nSCMoutput will be the total number
         ! of diagnostic entries in SCMop, the number of those entries
         ! being output to any stream (with no multiple counting of
         ! entries resulting from the same call to SCMoutput), and
         ! the number of calls made to SCMoutput so far this timestep
         ! respectively.
    INTEGER(i64) :: nentries,n_output,nSCMoutput

         ! Will hold the total number of expected/observed sub-steps
         ! and the current sub-step number
    INTEGER(i64) :: num_substeps,substep_number

         ! The diagnostic entries, 1:nentries have been set by newdiag
         ! via a call to SCMoutput. IF A NEW ARRAY IS ADDED HERE, THERE
         ! MUST BE AN ASSOCIATED SECTION OF CODE IN ROUTINE EXPAND_SCMOP.
    CHARACTER(len=lsname), POINTER :: sname(:) ! Short name
    CHARACTER(len=llname), POINTER :: lname(:) ! Long name
    CHARACTER(len=lunits), POINTER :: units(:) ! Units

    INTEGER(i64), POINTER :: &
      domprof(:)       &! Domain profile
    , timprof(:)       &! Time profile
    , streams(:)       &! List of streams to write to
    , dump_step(:)     &! Dumping period
    , nadd2dump(:)     &! Number of calls to add2dump this period
    , ncols(:)         &! No. of columns
    , nrows(:)         &! No. of rows
    , nlevs(:)         &! No. of levels
    , nelements(:)     &! elements in total (given the domain profile)
    , sname_id(:)      &! An integer unique to each sname
    , wd(:)            &! Index of another entry upon which this entry depends
    , lastencounter(:) &! Last timestep this entry was seen by SCMoutput
    , substep(:)        ! Substep this entry pertains to

    LOGICAL(l64), POINTER :: &
      only_radsteps(:)  ! Only defined on radiation timesteps?

    INTEGER(incdf), POINTER :: &
      netcdf_id(:)      ! A NetCDF integer id

    TYPE(allocatable_array), POINTER :: &
      diag(:)           ! Dump array

    ! The domain profiles, set by define_domprof.
    CHARACTER (len=15) ::   &
      d_name(maxndomprof)

    INTEGER(i64) ::         &
      d_rowa1 (maxndomprof) &
    , d_rowa2 (maxndomprof) &
    , d_rowb1 (maxndomprof) &
    , d_rowb2 (maxndomprof) &
    , d_lev1  (maxndomprof) &
    , d_lev2  (maxndomprof)

         ! This array will be used to memorise the order in which the
         ! calls to SCMoutput take place, which saves time in SCMoutput
         ! by avoiding the translation between the inputs and the
         ! corresponding diagnostic entries in this structure
    INTEGER(i64), POINTER :: diag_mem(:,:)

  END TYPE SCMop_type

END MODULE scmoptype_defn

