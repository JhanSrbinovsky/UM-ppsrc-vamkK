! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Initialise the SCM output diagnostic system

SUBROUTINE setup_diags                                                        &
  ( row_length, rows, model_levels, wet_levels, bl_levels, sm_levels          &
  , st_levels, land_points, ntiles, n_vis_thresh, cloud_levels, total_nsteps  &
  , timestep, full_daysteps, a_sw_radstep_prog, a_sw_radstep_diag             &
  , ntrad1, daycount, stepcount, nscmdpkgs, l_scmdiags, scmop )


  USE scm_cntl_mod, ONLY:                                                     &
      scm_nml, scmop_type, i64, maxndomprof, maxnstreams, main_diag_switch    &
    , netcdf_chunksize, diags, listlength, inot_written, lsname, default      &
    , strm_switch, strm_unit, strm_format, strm_filename, strm_dumpstep       &
    , strm_rejectlist, strm_acceptlist                                        &
    , strm_heed_hardwired, strm_heed_acceptlist, strm_heed_rejectlist         &
    , l_scmdiag_gen,   l_scmdiag_rad, l_scmdiag_bl,  l_scmdiag_surf           &
    , l_scmdiag_land,  l_scmdiag_sea, l_scmdiag_lsp, l_scmdiag_conv           &
    , l_scmdiag_lscld, l_scmdiag_pc2, l_scmdiag_gwd, l_scmdiag_forc           &
    , l_scmdiag_incs

  USE scm_utils, ONLY:                                                        &
      zhook_in, zhook_out, jprb, lhook, dr_hook

  IMPLICIT NONE

! Description:
!   Perform all pre-timestepping initialisation for the output
!   diagnostic system, i.e. initialise SCMop. In particular, the
!   DIAGS namelist is read in, and the output streams and domain
!   profiles are set up. The output files are -not- initialised
!   here since this requires information which is not available
!   until the end of at least one timestep (see dump_streams_init).
!   All arguments are INTENT In, apart from main_diag_switch,
!   the l_SCMdiags logicals and SCMop which are INTENT Out.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran90

  TYPE(SCMop_type) :: SCMop ! Out The derived-type structure
                            ! containing all the diagnostic
                            ! information

  INTEGER ::          &
    row_length        &
  , rows              &
  , model_levels      &
  , wet_levels        &
  , bl_levels         &
  , sm_levels         &
  , st_levels         &
  , land_points       &
  , ntiles            &
  , n_vis_thresh      &
  , cloud_levels      &
  , total_nsteps      &
  , full_daysteps     &
  , ntrad1            &
  , a_sw_radstep_prog &
  , a_sw_radstep_diag


  REAL :: timestep

  ! SCMop has pointers which will point to stepcount and daycount
  ! so that it can know which step it's on
  INTEGER, TARGET ::  &
    daycount          &
  , stepcount

! Parameters and stuff used around and about the internal
! workings of the diagnostic system
! Start of include file: s_scmop_internal.h
! Description:
!  Defines/declares some parameters and statement functions used
!  around and about the internal workings of the SCM diagnostic system.


! Include the parameters that are used in the calls
! to SCMoutput...
! Start of include file: s_scmop.h
! Description:
!  Declares and defines some parameters necessary for calling SCMoutput
!
!

! Integers to represent the different time profiles. All must
! be non-negative and less than "only_radsteps".

  INTEGER, PARAMETER :: &
    t_inst        = 1   &! Give the instantaneous value
  , t_avg         = 2   &! Construct the average value
  , t_max         = 3   &! " maximum value
  , t_min         = 4   &! " minimum value
  , t_acc         = 5   &! " accumulated value
  , t_div         = 7   &! " average value divided
                         !   by another diagnostic
  , t_mult        = 8   &! " average value multiplied
                         !   by another diagnostic
  , t_acc_div     = 9   &! " accumulated value divided
                         !   by another diagnostic
  , t_acc_mult    = 10  &! " accumulated value multiplied
                         !   by another diagnostic
  , t_const       = 11  &! The value is constant.
  , only_radsteps = 100  ! When added to one of the above parameters,
                         ! flags that the diagnostic is only available
                         ! on radiation timesteps

! Integers to represent the different domain profiles
  INTEGER, PARAMETER :: &
    d_sl      = 1       &
  , d_soilt   = 2       &
  , d_bl      = 3       &
  , d_wet     = 4       &
  , d_all     = 5       &
  , d_soilm   = 6       &
  , d_tile    = 7       &
  , d_vis     = 9       &
  , d_point   = 13      &
  , d_allxtra = 14      &
  , d_land    = 15      &
  , d_cloud   = 16

! Statement function to encode a stream number into an integer
  INTEGER :: &
    Stream   &
  , strm

  Stream(strm) = 2**(strm-1)

! The default streams for diagnostics to go to will be 1,2,3,4,5 and 6.
! The following should thus be equal to:
!
! Stream(1) [2^0=1] + Stream(2) [2^1=2]  + Stream(3) [2^2=4]
! Stream(4) [2^3=8] + Stream(5) [2^4=16] + Stream(6) [2^5=32]
! Total = 63
!
! where Stream() is the statement function defined above.
! default is 63 (all)

  INTEGER, PARAMETER :: &
    default_streams = 63

! Integers to represent the different diagnostics packages
  INTEGER, PARAMETER :: &
    SCMDiag_gen   = 1   & ! General diagnostics
  , SCMDiag_rad   = 2   & ! Radiation
  , SCMDiag_bl    = 3   & ! Boundary layer
  , SCMDiag_surf  = 4   & ! Surface
  , SCMDiag_land  = 5   & ! Land points only
  , SCMDiag_sea   = 6   & ! Sea points only
  , SCMDiag_lsp   = 7   & ! Large scale precip
  , SCMDiag_conv  = 8   & ! Convection
  , SCMDiag_lscld = 9   & ! Large scale cloud
  , SCMDiag_pc2   = 10  & ! PC2
  , SCMDiag_forc  = 11  & ! Forcing
  , SCMDiag_incs  = 12  & ! Increments
  , SCMDiag_gwd   = 13    ! Gravity Wave Drag

! End of include file: s_scmop.h

! Statement function to translate daycount and stepcount into
! an integer representing the number of timesteps since the
! start of the run

  INTEGER :: stepnmbr
  stepnmbr(SCMop) = (SCMop%daycount-1)*SCMop%full_daysteps &
                  +  SCMop%stepcount

!-------- Stream inquiry/manipulation statement functions -----------

! These need to be declared for the statement functions below
  integer(i64) :: streamlist

! Statement function to modify an encoded stream list
! to flag that the streams should not be written to file
  INTEGER :: DoNotWrite
  DoNotWrite(streamlist) = streamlist+2**(inot_written-1)

! Statement function to test whether a particular stream is
! switched on in an encoded stream list
  LOGICAL :: StreamIsOn
  StreamIsOn(streamlist,strm)=                                      &
       (streamlist-INT(streamlist/Stream(strm+1))*Stream(strm+1))/  &
       Stream(strm).GE.1

! Statement function to test whether the "do not write" flag
! has been switched on in an encoded stream list
  LOGICAL :: NotWritten
  NotWritten(streamlist)=StreamIsOn(streamlist,inot_written)

! Statement function to test whether any stream is switched
! on in an encoded stream list
  LOGICAL :: AnyStreamOn
  AnyStreamOn(streamlist)=                                          &
       (.NOT.NotWritten(streamlist).AND.streamlist.NE.0)            &
       .OR.                                                         &
       (NotWritten(streamlist).AND.streamlist.NE.DoNotWrite(0))

! End of include file: s_scmop_internal.h

  INTEGER ::          &
    i,j               &! Counters.
  , pos               &! Current position in string.
  , last_comma        &! Position of last comma.
  , n_words            ! No. of comma separated words found.

  CHARACTER (len=5) :: suffix        ! Holds frmt-dependent filename suffix
  CHARACTER (len=listlength) :: list ! Temporarily holds either an accept or
                                     ! reject list

  CHARACTER (len=lsname) :: name     ! Temporarily holds a (short) name of a
                                     ! diagnostic

  CHARACTER(LEN=*), PARAMETER ::  routinename = 'setup_diags'

  INTEGER :: istatus                 ! Error code

  ! Some constants which will be used to scale diagnostics
  REAL, PARAMETER :: rsec_day = 86400.0
  REAL :: ntimestepsperday,oneKsecday

  INTEGER :: nSCMdpkgs             ! No of SCM diagnostics packages
                                   ! (set in Scm_Main)
  LOGICAL :: l_SCMdiags(nSCMdpkgs) ! Logicals array for SCM
                                   ! diagnostics packages

  ! Dr Hook
  !=============================================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('SETUP_DIAGS',zhook_in,zhook_handle)

!
!     Note that s_scmop.h is included in s_scmop_internal.h
!

  SCMop%first_pass     =  .TRUE.       ! Diagnostic list not yet finalised
  SCMop%on             =  .FALSE.      ! Timestepping not yet started
  SCMop%daycount       => daycount
  SCMop%stepcount      => stepcount
  SCMop%full_daysteps  =  full_daysteps
  SCMop%ntrad          =  a_sw_radstep_diag
  SCMop%ntrad1         =  ntrad1
  SCMop%maxnentries    =  0            ! No memory allocated for diagnostics yet
  SCMop%nentries       =  0            ! No diagnostics requested
                                       ! to be stored yet
  SCMop%n_output       =  0            ! No diagnostics requested
                                       ! to be output yet
  SCMop%nSCMoutput     =  0            ! No calls to SCMoutput yet
  SCMop%num_substeps   =  1            ! Expected number of sub-steps
  SCMop%substep_number =  0            ! No substepping started yet

  ! Initialise one element of the domain profiles to flag
  ! undefined profiles
  DO i=1, maxndomprof
    SCMop%d_lev1(i) = -1
  END DO

  ! Allocate memory for streams
  ALLOCATE(SCMop%strm(maxnstreams))

  ! Allocate some token space to the diag_mem array (ensures the
  ! SIZE(diag_mem) operation in SCMoutput doesn't fall over the
  ! first time)
  ALLOCATE(SCMop%diag_mem(1,2))


  ! Modify default values of some namelist variables...

  ! Diagnostic system on by default (cannot be set by data
  ! statement because it is intent out).
  main_diag_switch = 1

  ! Default NetCDF file chunk size
  netcdf_chunksize = 8192

  ! l_SCMdiags logicals (intent Out) are also set to default values
  l_SCMdiag_gen   = .TRUE.
  l_SCMdiag_rad   = .FALSE.
  l_SCMdiag_bl    = .FALSE.
  l_SCMdiag_surf  = .FALSE.
  l_SCMdiag_land  = .FALSE.
  l_SCMdiag_sea   = .FALSE.
  l_SCMdiag_lsp   = .FALSE.
  l_SCMdiag_conv  = .FALSE.
  l_SCMdiag_lscld = .FALSE.
  l_SCMdiag_PC2   = .FALSE.
  l_SCMdiag_gwd   = .FALSE.
  l_SCMdiag_forc  = .FALSE.
  l_SCMdiag_incs  = .FALSE.

  ! Stream 1 will be on by default
  strm_switch(1) = 1

  ! The default output units for each stream (37 upwards)...
  DO i=1, maxnstreams
    strm_unit(i) = 36+i
  END DO

  ! Read the DIAGS namelist scm_nml
  OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus, STATUS='old')
  IF (istatus == 0) THEN
    READ(10, DIAGS, iostat=istatus)
    CLOSE(10)
    IF (istatus /= 0) THEN
      WRITE(6,*)'------------------------------------------------'
      WRITE(6,*)'SETUP_DIAGS WARNING: DIAGS namelist not found in'
      WRITE(6,*)'namelist file. Diagnostic system will operate on'
      WRITE(6,*)'defaults.'
      WRITE(6,*)'------------------------------------------------'
    END IF
  ELSE
    WRITE(6,*) "-------------------------------------------------"
    WRITE(6,*) " Error opening " // TRIM(ADJUSTL(scm_nml))
    WRITE(6,*) " IO status = ", istatus
    WRITE(6,*) " Using default diagnostic settings"
    WRITE(6,*) "-------------------------------------------------"
  END IF

  ! No point doing any more if the diagnostic system has been
  ! switched off
  IF (main_diag_switch == 0) RETURN

  ! NetCDF I/O is done using C routines which don't utilise Fortran
  ! unit numbers. Set to -1 so that this is reflected in the stream
  ! summary print-out.
  DO i=1, maxnstreams
    IF (strm_format(i) == 4) THEN
      strm_unit(i)=-1
    END IF
  END DO

  ! Set up the output streams...
  DO i=1, maxnstreams
    SCMop%strm(i)%switch  = strm_switch(i)     ! Stream on/off
    SCMop%strm(i)%op_unit = strm_unit(i)       ! Output unit
    SCMop%strm(i)%format  = strm_format(i)     ! Output format

    ! Output filename...
    IF (strm_filename(i) /= Default) THEN
      SCMop%strm(i)%filename = strm_filename(i)
    ELSE

      ! The filename suffix will be format dependent
      suffix='.out'     ! Unrecognised format
      IF (SCMop%strm(i)%format == 0) THEN
        suffix = '.dat'   ! Old Wave format
      ELSE IF (SCMop%strm(i)%format == 1) THEN
        suffix = '.txt'   ! Easy to read format
      ELSE IF (SCMop%strm(i)%format == 2) THEN
        suffix = '.fssi'  ! Format for FSSI
      ELSE IF (SCMop%strm(i)%format == 3) THEN
        suffix = '.dat'   ! New Wave format
      ELSE IF (SCMop%strm(i)%format == 4) THEN
        suffix = '.nc'    ! NetCDF format
      END IF
      WRITE(SCMop%strm(i)%filename,'(A,I2.2,A)')'stream',i,suffix
    END IF

    ! If the dumping period is positive, take this as a number
    ! of timesteps. If it is negative take it as a number of
    ! seconds. If it is zero, replace it with the number of
    ! timesteps in the whole run (i.e. one dump on the
    ! last step).
    IF (strm_dumpstep(i) == 0) strm_dumpstep(i)=total_nsteps

    IF (strm_dumpstep(i) >  0) THEN

      ! The dumping period has been specified in timesteps
      SCMop%strm(i)%dump_step = strm_dumpstep(i)
    ELSE IF (strm_dumpstep(i) <  0) THEN

      ! The dumping period has been specified in seconds, convert
      ! to nearest non-zero number of timesteps
      SCMop%strm(i)%dump_step = MAX(NINT(-strm_dumpstep(i)/timestep),1)
    END IF

    ! Accept diagnostics sent to stream by respective SCMoutput
    ! parameter? (0=no)
    SCMop%strm(i)%heed_hardwired = strm_heed_hardwired(i)

    ! Accept diagnostics sent to stream by namelist? (0=no)
    SCMop%strm(i)%heed_acceptlist = strm_heed_acceptlist(i)

    ! Reject diagnostics sent to stream by any method? (0=no)
    SCMop%strm(i)%heed_rejectlist = strm_heed_rejectlist(i)

    ! No. of diagnostics sent to this stream (none until at least
    ! first call to SCMoutput)
    SCMop%strm(i)%n_output = 0

    ! We now want to extract the individual diagnostic names from
    ! the comma-separated lists in strm_acceptlist(i) and
    ! strm_rejectlist(i). We will scan through both strings
    ! twice, first to count the number of names so the correct
    ! space can be allocated (j=1 and j=3) and then again to put
    ! the names into the allocated arrays (j=2 and j=4).
    DO j=1, 4

      IF (j <= 2) THEN
        ! First two loops work on the acceptlist
        list = strm_acceptlist(i)
      ELSE
        ! Next two are indepedent of the first and work
        ! on the reject list
        list = strm_rejectlist(i)
      END IF

      ! We're going to scan through the list of
      ! comma-separated words one character at a time
      ! starting from the beginning.
      ! Make sure there's a trailing comma on the list
      list = TRIM(list)//','
      pos = 1               ! Current position
      last_comma = 0        ! Pos'n of last comma found
      n_words = 0           ! No. of words found so far

      ! Loop through the characters of the string
      DO WHILE (pos <= LEN_TRIM(list))

        IF &
            (list(pos:pos) == '\\') &
            THEN
          ! This character is a backslash - it escapes the
          ! next character causing us to jump forward one.
          pos = pos+1

        ELSE IF (list(pos:pos) == ' '.AND.                        &
                 last_comma == pos-1) THEN
          ! This is a space trailing a comma - ignore by
          ! "moving" the last comma forward one
          last_comma = pos
        ELSE IF (list(pos:pos) == ',') THEN

          ! We've found a (non-escaped) comma. Extract the name
          ! enclosed by this comma and the last.
          name = list(last_comma+1:pos-1)

          ! Update the position of the last comma found
          last_comma = pos

          ! Count the number of names we've found
          n_words = n_words+1

          ! On the second pass for each list, fill the allocated
          ! arrays with the names.
          IF (j == 2) SCMop%strm(i)%accept_list(n_words) = name
          IF (j == 4) SCMop%strm(i)%reject_list(n_words) = name
        END IF

          ! Update the position
          pos = pos+1

      END DO ! over characters in a list

      ! After the first pass for each list, allocate the arrays
      ! into which the names will be put on the second pass.
      IF (j == 1) ALLOCATE(SCMop%strm(i)%accept_list(n_words))
      IF (j == 3) ALLOCATE(SCMop%strm(i)%reject_list(n_words))

    END DO ! j=1,4

     ! FSSI-format output files require a dump for the first
     ! timestep as well as at the end of every dumping period. The
     ! way this has been coded in routine dump_streams means that
     ! if the dumping period is equal to one then this additional
     ! dump will not be produced and the file will be one line
     ! shorter than expected. I am assured that the SSFM will never
     ! be run with a dumping period of one, but it's worth printing
     ! a warning in case it ever is.
     IF (SCMop%strm(i)%format == 2.AND.                             &
          SCMop%strm(i)%dump_step == 1) THEN
        WRITE(6,*)'---'
        WRITE(6,*)'setup_diags WARNING: a FSSI-format output file'
        WRITE(6,*)'has been requested with a dumping period of'
        WRITE(6,*)'one. This file will not have an additional dump'
        WRITE(6,*)'for the first timestep and so may be one line'
        WRITE(6,*)'shorter than expected.'
        WRITE(6,*)'---'
     END IF

  END DO                     ! i=1,maxnstreams


! Set up the domain profiles
! DEPENDS ON: define_domprof
  CALL define_domprof(d_all,'all_levels', &
       1,row_length,                      &! Every column
       1,rows,                            &! Every row
       1,model_levels,SCMop)               ! Every level

! DEPENDS ON: define_domprof
  CALL define_domprof(d_sl,'single_level',&
       1,row_length,                      &! Every column
       1,rows,                            &! Every row
       1,1,SCMop)                          ! Only the lowest level

! DEPENDS ON: define_domprof
  CALL define_domprof(d_wet,'wet_levels',            &
       1,row_length,                                 &
       1,rows,                                       &
       1,wet_levels,SCMop)

! DEPENDS ON: define_domprof
  CALL define_domprof(d_bl,'bl_levels',              &
       1,row_length,                                 &
       1,rows,                                       &
       1,bl_levels,SCMop)

! DEPENDS ON: define_domprof
  CALL define_domprof(d_soilm,'soil_moist_levels',   &
       1,row_length,                                 &
       1,rows,                                       &
       1,sm_levels,SCMop)

! DEPENDS ON: define_domprof
  CALL define_domprof(d_soilt,'soil_temp_levels',    &
       1,row_length,                                 &
       1,rows,                                       &
       1,st_levels,SCMop)

! DEPENDS ON: define_domprof
  CALL define_domprof(d_land,'land_points',          &
       1,land_points,                                &! This is a fudge
       1,1,                                          &
       1,1,SCMop)

! DEPENDS ON: define_domprof
  CALL define_domprof(d_allxtra,'all_levs_plus1',    &
       1,row_length,                                 &
       1,rows,                                       &
       1,model_levels+1,SCMop)

! DEPENDS ON: define_domprof
  CALL define_domprof(d_tile,'tile_types',           &
       1,row_length,                                 &
       1,rows,                                       &
       1,ntiles,SCMop)

! DEPENDS ON: define_domprof
  CALL define_domprof(d_point,'single_point',1,1,1,1,1,1,SCMop)

! DEPENDS ON: define_domprof
  CALL define_domprof(d_vis,'vis_thresholds',        &
       1,row_length,                                 &
       1,rows,                                       &
       1,n_vis_thresh,SCMop)

! DEPENDS ON: define_domprof
  CALL define_domprof(d_cloud,'cloud_levels',        &
       1,row_length,                                 &
       1,rows,                                       &
       1,cloud_levels,SCMop)

  ! Some diagnostics are just constants which will be combined
  ! with other, variable diagnostics. May as well just call
  ! SCMoutput for them once then, rather than every timestep.

! DEPENDS ON: scmoutput
  CALL scmoutput(rsec_day, 'sec_day'                                 &
    , 'No. of seconds per day','-'                                   &
    , t_const, d_point, DoNotWrite(default_streams), '', RoutineName)

  oneKsecday=1000.0*rsec_day

! DEPENDS ON: scmoutput
  CALL scmoutput(oneKsecday, 'oneKsecday'                            &
    , '1000.*rsec_day','-'                                           &
    , t_const, d_point, DoNotWrite(default_streams), '', RoutineName)

  ntimestepsperday = rsec_day/timestep

! DEPENDS ON: scmoutput
  CALL scmoutput(ntimestepsperday, 'ntspday'                         &
    , 'No. of timesteps per day', '-'                                &
    , t_const, d_point, DoNotWrite(default_streams), '', RoutineName)

  ! Place details from read in diagnostics packages logicals
  ! into logicals array
  ! Note: general is default true, the rest default false
  l_SCMdiags(SCMdiag_gen)   = l_SCMdiag_gen
  l_SCMdiags(SCMdiag_rad)   = l_SCMdiag_rad
  l_SCMdiags(SCMdiag_bl)    = l_SCMdiag_bl
  l_SCMdiags(SCMdiag_surf)  = l_SCMdiag_surf
  l_SCMdiags(SCMdiag_land)  = l_SCMdiag_land
  l_SCMdiags(SCMdiag_sea)   = l_SCMdiag_sea
  l_SCMdiags(SCMdiag_lsp)   = l_SCMdiag_lsp
  l_SCMdiags(SCMdiag_conv)  = l_SCMdiag_conv
  l_SCMdiags(SCMdiag_lscld) = l_SCMdiag_lscld
  l_SCMdiags(SCMdiag_PC2)   = l_SCMdiag_PC2
  l_SCMdiags(SCMdiag_gwd)   = l_SCMdiag_gwd
  l_SCMdiags(SCMdiag_forc)  = l_SCMdiag_forc
  l_SCMdiags(SCMdiag_incs)  = l_SCMdiag_incs

  IF (lhook) CALL dr_hook('SETUP_DIAGS',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE setup_diags
