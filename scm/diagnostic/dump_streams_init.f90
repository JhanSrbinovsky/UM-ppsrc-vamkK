! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Initialise SCM diagnostic output files

SUBROUTINE dump_streams_init                                                  &
  ( SCMop, row_length, rows, model_levels, wet_levels, bl_levels              &
  , cloud_levels, ozone_levels, st_levels, sm_levels, ntiles, year_init       &
  , month_init, day_init, hour_init, min_init, sec_init, timestep, ndayin     &
  , nminin, nsecin, sec_day, tot_nsteps, a_sw_radstep_prog                    &
  , a_sw_radstep_diag, z_top_of_model, first_constant_r_rho_level             &
  , eta_theta, eta_rho, orog, r_theta_levels, r_rho_levels, netcdf_chunksize )

  USE netcdf
  USE UM_types

  USE scm_utils, ONLY:                                                        &
      scm_timestep, scm_timestep_count, time_info, scm_nc_time_id             &
    , scm_nc_date_id, scm_nc_hour_id, scm_nc_local_id, zhook_in, zhook_out    &
    , jprb, lhook, dr_hook

  ! SCMop_type is defined in here...
  USE scmoptype_defn

  IMPLICIT NONE

! Description:
!   Perform all steps necessary for subsequent writing of output stream data :
!   write the headers of the data files and any auxilliary files that may
!   accompany them

! Method:

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran77 with some bits of Fortran90

! NOTE: on the NEC this routine may be compiled in 32 bits because requests
!       to create 32 bit integers (necessary to communicate with the NetCDF
!       library routines - see incdf) are ignored by the compiler when
!       compiling at 64 bits. Thus all input and output variables are
!       explicitly declared as 64 bit using the KIND types i64, r64 and l64
!       delcared in scmoptype_defn module

  TYPE(SCMop_type) :: SCMop ! InOut The derived-type structure containing
                            !       all the diagnostic information

! In Horizontal and vertical model domain...
  INTEGER(i64) :: &
    row_length    &
  , rows          &
  , model_levels  &
  , wet_levels    &
  , bl_levels     &
  , cloud_levels  &
  , ozone_levels  &
  , st_levels     &
  , sm_levels     &
  , ntiles

  INTEGER(i64) :: &
    year_init     &! In Initial year
  , month_init    &! In Initial month
  , day_init      &! In Initial day
  , hour_init     &! In Initial hour
  , min_init      &! In Initial minute
  , sec_init       ! In Initial second

  REAL(r64) ::    &
    timestep       ! In Timestep

  INTEGER(i64) ::     &
    ndayin            &! In No. of days requested in run
  , nminin            &! In No. of minutes requested in run
  , nsecin            &! In No. of seconds requested in run
  , sec_day           &! In No. of seconds in a day (Why this
                       !    isn't in an include file I don't know)
  , tot_nsteps        &! In Total no. of steps to be done
  , a_sw_radstep_prog &! In No. of timesteps between prognostic
                       !    calls to radiation (3C/3Z)
  , a_sw_radstep_diag  ! In No. of timesteps between diagnostic
                       !    calls to radiation (3C/3Z)

  REAL(r64) ::        &
     z_top_of_model    ! In The height of the "top of the
                       !    atmosphere"

  INTEGER(i64) ::     &! In Lowest rho level that has constant height
     first_constant_r_rho_level

  REAL(r64) ::                  &
     eta_theta(model_levels+1)  &! In The etas of the theta and rho levels
   , eta_rho(model_levels)      &!
   , orog(row_length,rows)       ! In Orography height

  ! IN The physical heights of the theta and rho levels...
  REAL(r64) ::                                      &
     r_theta_levels(row_length,rows,0:model_levels) &
   , r_rho_levels(row_length,rows,model_levels)

  ! IN The ChunkSize input to NF90_Create(), the routine that
  ! creates a NetCDF file. Controls a space versus time trade-off:
  ! memory allocated in the netcdf library versus number of system
  ! calls.
  INTEGER(i64) :: netcdf_chunksize

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

  INTEGER :: get_um_version ! Function subroutine which returns
                            ! the UM version as defined by
                            ! environment variable $VN.
  INTEGER :: um_version_env ! Stores UM version from $VN.

  CHARACTER (len=3) :: um_version_c ! The UM version as a string.
  CHARACTER (len=100) :: fmt        ! Holds format specifiers

  INTEGER ::               &
    i, n, m                &! Counters
  , nhourin                 ! No. of hours requested in run

  INTEGER(i64) ::          &! 64-bit counter for passing to other
    n64                     ! routines compiled in 64 bit

  INTEGER :: unit

  INTEGER(i64) ::          &
    diags(SCMop%nentries)  &
  , ndiags

  REAL :: ndump

  ! Copy of netcdf_chunksize at the correct precision for passing
  ! into NetCCDF routines. May be altered by NF90_Create.
  INTEGER(incdf) :: netcdf_chunksize_io

  ! NetCDF identifiers for each dimension of each domain
  ! profile. Integers that are passed into NetCDF library routines
  ! should have the same precision as those used internally within the
  ! library (which may have been compiled at 32 bit), hence the use of
  ! incdf.
  INTEGER (incdf) ::                 &
    netcdf_dimension_time            &
  , netcdf_dimension(3,maxndomprof)

  ! Other NetCDF-format related variables
  INTEGER (incdf) :: Status,NcID
  CHARACTER (len=5) :: c_unit

  ! Character function that incorporates the substep number into
  ! the short name of a diagnostic. Somewhat longer than a normal
  ! short name.
  CHARACTER (len=lsname+10) :: add_substep_to_sname

  ! Temporary variables to hold output from add_substep_to_sname
  CHARACTER (len=lsname+10) :: short_name
  CHARACTER (len=lsname+10), ALLOCATABLE :: short_names(:)
  CHARACTER (len=100) :: start_date_str


  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('DUMP_STREAMS_INIT',zhook_in,zhook_handle)

  ! Get the UM version from the environment variable $VN
! DEPENDS ON: get_um_version
  um_version_env=get_um_version()

  ! Check for a sensible version number
  IF (um_version_env < 0) THEN
    WRITE(*,'(A)') '---'
    WRITE(6,'(A)') 'Dump_Streams_Init: WARNING, the '           //  &
                   'environment variable $VN does not seem'
    WRITE(6,'(A)') 'to be set to a valid UM version number. I ' //  &
                   'will assume this is UM VN7.4,'
    WRITE(6,'(A)') 'but if incorrect it could lead to '         //  &
                   'unreadable output files.'
    WRITE(*,'(A)') '---'
    um_version_env = 704
  END IF

  ! Convert the UM version number to character form
  WRITE(um_version_c,'(I1,A1,I1)')                                  &     
    INT(um_version_env/100.), '.', MOD(INT(um_version_env),10)

  ! Write out the UM version number
  WRITE(*,'(A)') '==================='
  WRITE(*,'(A)') 'UM Version='//TRIM(um_version_c)
  WRITE(*,'(A)') '==================='
  WRITE(*,'(A)') ' '

  ! Give a summary of the output streams
  IF (SCMop%n_output /= 0) THEN
    WRITE(*,'(65("-"))')
    WRITE(*,'(A)')'Summary of open streams:'
    WRITE(*,'(65("-"))')
    WRITE(*,'(A)')'No.| Unit | Time between dumps | No. dgnstcs '//&
            '| Format | Filename '
    WRITE(*,'(A)')'   |      | (steps) |(minutes) |             '//&
            '|        |          '
    WRITE(*,'(A)')'---+------+---------+----------+-------------'//&
            '+--------+----------'
    DO i=1, maxnstreams
      IF (SCMop%strm(i)%switch /= 0) THEN

        ! Write the stream output unit into a character variable,
        ! unless this is a NetCDF stream, in which case the unit
        ! is not used.
        IF (SCMop%strm(i)%format /= 4) THEN
          WRITE(c_unit,'(I5)') SCMop%strm(i)%op_unit
        ELSE
          c_unit='n/a'
          c_unit=ADJUSTR(c_unit)
        END IF

        WRITE(*,'(I2," |",A5," |",I8," | ",G8.3," |",I12," |",I7," | ",A)') &
          i, c_unit,                                                        &
          SCMop%strm(i)%dump_step,                                          & 
          SCMop%strm(i)%dump_step*timestep/60.0,                            &
          SCMop%strm(i)%n_output,SCMop%strm(i)%format,                      &
          TRIM(SCMop%strm(i)%filename)
      END IF
    END DO

    WRITE(*,'(65("-"))')
    WRITE(*,'(A)')' '

  ELSE
    PRINT*,'Dump_Streams_Init: No diagnostics to be written.'
    ! Nothing to do if there's no diagnostics to be written
    GOTO 9999
  END IF

  ! Write the scumlist file
  n64 = 0
! DEPENDS ON: write_scumlist
  CALL write_scumlist (SCMop, n64)

  ! write the no. of levels, etc into a file
! DEPENDS ON: write_domain_data
  CALL write_domain_data                                           &
    ( row_length, rows, model_levels, z_top_of_model               &
    , first_constant_r_rho_level, orog, eta_theta,eta_rho )

  ! Calculate nhourin from nminin and nsecin
  nhourin = INT(nminin/60.0+nsecin/3600.0)

  ! Write the headers of the output data files
  DO n=1, maxnstreams
    IF (SCMop%strm(n)%switch /= 0.AND.                             &
        SCMop%strm(n)%n_output >  0) THEN
      ! This stream is open and there are diagnostics
      ! to be written out from it
      unit = SCMop%strm(n)%op_unit

      ! Check this unit is free (not applicable for NetCDF files)
      DO m=1, n-1
        IF (SCMop%strm(m)%switch /= 0.AND.                         &
            SCMop%strm(m)%format /= 4.AND.                         &
            SCMop%strm(m)%n_output > 0.AND.                        &
            SCMop%strm(m)%op_unit == unit) THEN
          PRINT*,'Dump_Streams_Init ERROR: two output files '//    &
                 'with the same unit numbers: ',n,m,unit
          PRINT*,'Output redirected to unit 103!'

          unit = 103
          SCMop%strm(n)%op_unit = unit
        END IF
      END DO

! The header depends on the format

!-------------------------------------------------------------------
! Format=0: designed for reading by PV-wave routine
!           scmread2.pro
!-------------------------------------------------------------------
      IF (SCMop%strm(n)%format == 0) THEN

        ! Open the stream's output file
        OPEN(unit=unit,file=SCMop%strm(n)%filename)

        ! Write the UM version and the stream format
        WRITE(unit,'(A,I1)')'vn'//um_version_c//',old_format'

        ! Write the total no. of days, hours, dumps per day and
        ! steps, plus the timestep, the dumping period of this
        ! stream and the no. of diagnostics per dump
        ndump = sec_day/(SCMop%strm(n)%dump_step*timestep)
        WRITE(unit,                                              &
             '(I4,1X,I4,1X,F7.3,1X,I7,1X,F9.2,1X,I9,1X,I4)')     &
              ndayin, nhourin, ndump, tot_nsteps, timestep,      &
              SCMop%strm(n)%dump_step,                           &
              SCMop%strm(n)%n_output

        ! Make a list the diagnostics we're going to write
        ndiags = 0
        DO i=1, SCMop%nentries
          IF (StreamIsOn(SCMop%streams(i),n).AND..NOT.           &
              NotWritten(SCMop%streams(i))) THEN

            ndiags = ndiags+1
            diags(ndiags) = i
          END IF
        END DO

        IF (ndiags /= SCMop%strm(n)%n_output) THEN
          PRINT*,'Dump_Streams_Init ERROR: an inconsistency '//  &
                 'has ocurred !',ndiags, SCMop%strm(n)%n_output, &
                  n,SCMop%nentries

          ! Switch the diagnostic system off so dodgy data is
          ! not mistakenly used in good faith.
          PRINT*,'Switching diagnostic system off!'
          SCMop%on = .FALSE.
        END IF

        ! Write that list to the file
        WRITE(fmt,'(A,I5,A)') '(',ndiags,'(I5,1X))'
        WRITE(unit,fmt)(SCMop%sname_id(diags(i)),i=1,ndiags)

!-------------------------------------------------------------------
! Format=1 : designed for easy visual perusal.
!-------------------------------------------------------------------
      ELSE IF (SCMop%strm(n)%format == 1) THEN

        ! Open the stream's output file
        OPEN(unit=unit,file=SCMop%strm(n)%filename)

        WRITE(unit,'(A,A)')'UM version=',um_version_c

        ! Write the total no. of days, hours, dumps per day and
        ! steps, plus the timestep, the dumping period of this
        ! stream and the no. of diagnostics per dump
        ndump = sec_day/(SCMop%strm(n)%dump_step*timestep)
        WRITE(unit,'(A)')'General run information...'
        WRITE(unit,'(A,I4,1X,A,I4,1X,A,F7.3,1X,A,I7)')                  &
             'No. of days=',ndayin,'No. of hrs=',nhourin,               &
             'No. of dumps=',ndump,'Total number of steps=',            &
              tot_nsteps
        WRITE(unit,'(A,F9.2,1X,A,I9,1X,A,I4)')                          &
             'Timestep=',timestep,                                      &
             'No. of steps between dumps=',                             &
             SCMop%strm(n)%dump_step,                                   &
             'No. of diags outputting to this stream=',                 &
             SCMop%strm(n)%n_output

        ! list the diagnostics we're going to write
        ndiags = 0
        DO i=1, SCMop%nentries
          IF (StreamIsOn(SCMop%streams(i),n).AND..NOT.                  &
            NotWritten(SCMop%streams(i))) THEN

            ndiags = ndiags+1
            diags(ndiags) = i
          END IF
        END DO

        IF (ndiags /= SCMop%strm(n)%n_output) THEN
          PRINT*,'Dump_Streams_Init ERROR: an inconsistency has ocurred !' &
                , ndiags, SCMop%strm(n)%n_output, n,SCMop%nentries
 
          ! Switch the diagnostic system off so dodgy data is
          ! not mistakenly used in good faith.
          PRINT*,'Switching diagnostic system off!'
          SCMop%on = .FALSE.
        END IF

        WRITE(unit,'(A)')'Diagnostics in this file... '//               &
             '(subset of entries in scumlist file)'

        ! Ensure that the sub-step is appended to the short
        ! name of each diagnostic entry, if required.

        ALLOCATE(short_names(ndiags))

        DO i=1, ndiags
! DEPENDS ON: add_substep_to_sname
          short_names(i) = add_substep_to_sname(SCMop,diags(i))
        END DO ! i

        ! Write the information
        WRITE(unit,'(I3,1X,A,1X,A)')                                    &
             (SCMop%sname_id(diags(i)),short_names(i),                  &
              SCMop%lname(diags(i)),i=1,ndiags)

        DEALLOCATE(short_names)

!-------------------------------------------------------------------
! Format=2 : format required for FSSI database
!-------------------------------------------------------------------
      ELSE IF (SCMop%strm(n)%format == 2) THEN

        ! Initialisation done in DUMP_STREAMS for the time being.

!-------------------------------------------------------------------
! Format=3 : new format designed for reading by PV-wave
!            routine scmread2.pro
!-------------------------------------------------------------------
      ELSE IF (SCMop%strm(n)%format == 3) THEN

        ! Open the stream's output file
        OPEN(unit=unit,file=SCMop%strm(n)%filename)

        ! Write the UM version, plus any extra strings necessary
        ! for PV-Wave to know what to expect in the file.
        WRITE(unit,'(A,I1)')'vn'//um_version_c//',start_time_present'

        ! Write information about the size and shape of the domain
        WRITE(unit,'(I4,I4,8I4)') row_length,rows,                      &
             model_levels,wet_levels,bl_levels,cloud_levels,            &
             ozone_levels,st_levels,sm_levels,ntiles

        ! Write the starting time of the run as specified in
        ! the INDATA namelist
        WRITE(unit,'(I6,1x,I4,1x,I5,1x,I4,1x,I4,1x,I4,1x)')             &
             year_init,month_init,day_init,hour_init,min_init,          &
             sec_init

        ! Write information pertaining to the length of
        ! the model run and time in general
        Write(unit,'(I9,1x,F9.2,1x,I9,1x,I4,1x,I4,1x,I5,1x,I5,1x,I5)')  &
             tot_nsteps,timestep,SCMop%strm(n)%dump_step,               &
             a_sw_radstep_prog,a_sw_radstep_diag,ndayin,nminin,nsecin

        ! Write information about the levels and orography...
        WRITE(fmt,'(A,I3,A)') '(',model_levels+1 ,'(1PE18.10E3," "))'
        WRITE(unit,fmt) eta_theta

        WRITE(fmt,'(A,I3,A)') '(',model_levels   ,'(1PE18.10E3," "))'
        WRITE(unit,fmt) eta_rho

        WRITE(fmt,'(A,I6,A)') '(',row_length*rows,'(1PE18.10E3," "))'
        WRITE(unit,fmt) orog

        WRITE(fmt,'(A,I6,A)') '(',row_length*rows*(model_levels+1)      &
                             ,'(1PE18.10E3," "))'
        WRITE(unit,fmt) r_theta_levels

        WRITE(fmt,'(A,I6,A)') '(',row_length*rows*(model_levels  )      &
                             ,'(1PE18.10E3," "))'

        WRITE(unit,fmt) r_rho_levels

        ! Write data into the file pertaining to the
        ! individual diagnostics that it will contain.

        n64 = n
! DEPENDS ON: write_scumlist
        CALL write_scumlist(SCMop,n64)

!-------------------------------------------------------------------
! Format=4 : NetCDF
!-------------------------------------------------------------------
      ELSE IF (SCMop%strm(n)%format == 4) THEN

        ! Create the NetCDF file. Chunksize is a parameter,
        ! the size of which may greatly affect the I/O speed
        ! of the NetCDF file writing.
        netcdf_chunksize_io = netcdf_chunksize

        Status = Nf90_Create(SCMop%strm(n)%filename,                    &
        Nf90_Clobber,NcID,Chunksize = netcdf_chunksize_io)

        ! Set global attributes
        status = nf90_put_att                  &
               ( ncid   = ncid                 &
               , varid  = nf90_global          &
               , name   = "Start date"         &
               , values = TRIM(start_date_str) )
        
        ! Declare the three spatial dimensions of each domain
        ! profile
        DO i=1, maxndomprof

          ! Check this domain has been defined.
          IF (SCMop%d_lev1(i) /= -1) THEN
            Status = Nf90_Def_Dim(                                      &
              NcID  = NcID,                                             &
              Name  = TRIM(ADJUSTL(SCMop%d_name(i)))//'_i',             &
              Len   = INT(SCMop%d_rowa2(i)-SCMop%d_rowa1(i)+1,incdf),   &
              DimID = netcdf_dimension(1,i))

            Status = Nf90_Def_Dim(                                      &
              NcID  = NcID,                                             &
              Name  = TRIM(ADJUSTL(SCMop%d_name(i)))//'_j',             &
              Len   = INT(SCMop%d_rowb2(i)-SCMop%d_rowb1(i)+1,incdf),   &
              DimID = netcdf_dimension(2,i))

            Status = Nf90_Def_Dim(                                      &
              NcID  = NcID,                                             &
              Name  = TRIM(ADJUSTL(SCMop%d_name(i)))//'_k',             &
              Len   = INT(SCMop%d_lev2(i)-SCMop%d_lev1(i)+1,incdf),     &
              DimID = netcdf_dimension(3,i))
          END IF
        END DO

        ! Declare the time dimension
        Status = Nf90_Def_Dim(                                          &
          NcID  = NcID,                                                 &
          Name  = "time",                                               &
          Len   = INT(tot_nsteps/SCMop%strm(n)%dump_step,incdf),        &
          DimID = netcdf_dimension_time)

        ! Define a seconds variable
        status = nf90_def_var               &
          ( ncid   = ncid                   &
          , name   = "seconds"              &
          , xtype  = nf90_float             &
          , dimids = netcdf_dimension_time  &
          , varid  = scm_nc_time_id         )

        ! Define seconds variable attributes
        status = nf90_put_att               &
          ( ncid   = ncid                   &
          , name   = "units"                &
          , varid  = scm_nc_time_id         &
          , values = "s"                    )

        status = nf90_put_att               &
          ( ncid   = ncid                   &
          , name   = "long_name"            &
          , varid  = scm_nc_time_id         &
          , values = "Number of seconds since start of simulations" )


        ! Define a date variable
        status = nf90_def_var               &
          ( ncid   = ncid                   &
          , name   = "date"                 &
          , xtype  = nf90_int               &
          , dimids = netcdf_dimension_time  &
          , varid  = scm_nc_date_id         )

        ! Define date variable attributes
        status = nf90_put_att               &
          ( ncid   = ncid                   &
          , name   = "units"                &
          , varid  = scm_nc_date_id         &
          , values = "yyyymmdd"             )

        status = nf90_put_att               &
          ( ncid   = ncid                   &
          , name   = "long_name"            &
          , varid  = scm_nc_date_id         &
          , values = "Date string"          )


        ! Define each diagnostic that will go to this stream
        DO i=1, SCMop%nentries

          ! Is this entry to be sent to this stream, and
          ! not been flagged as one not to output?
          IF (StreamIsOn(SCMop%streams(i),n).AND..NOT.                  &
              NotWritten(SCMop%streams(i))) THEN

            ! Get the short name with the substep number
            ! appended if necessary.
            n64 = i
            short_name = add_substep_to_sname(SCMop,n64)

            ! Define this variable
            Status = Nf90_Def_Var(                                      &
              NcID   = NcID,                                            &
              Name   = short_name,                                      &
              XType  = Nf90_Float,                                      &
              DimIDs = (/                                               &
                netcdf_dimension(1,SCMop%domprof(i)),                   &
                netcdf_dimension(2,SCMop%domprof(i)),                   &
                netcdf_dimension(3,SCMop%domprof(i)),                   &
                netcdf_dimension_time/),                                &
              VarID  = SCMop%netcdf_id(i))

            IF (Status /= nf90_NoErr) THEN
              PRINT*,'Error defining NetCDF variable: ', short_name
            END IF

            ! Define its attributes
            Status = Nf90_Put_Att(                                      &
              NcID   = NcID,                                            &
              VarID  = SCMop%netcdf_id(i),                              &
              Name   = "units",                                         &
              Values = SCMop%units(i))

            Status = Nf90_Put_Att(                                      &
              NcID   = NcID,                                            &
              VarID  = SCMop%netcdf_id(i),                              &
              Name   = "long_name",                                     &
              Values = SCMop%lname(i))
          END IF
        END DO

        ! End of file definition
        Status = Nf90_EndDef(NcID = NcID)

        ! Record the file ID so it can be used later for output.
        SCMop%strm(n)%op_unit = NcID

      ELSE
        PRINT*,'Dump_Streams_Init ERROR: unknown format '//             &
               'for stream ',n
      END IF

    END IF                  ! SCMop%strm(n)%switch /= 0
  END DO                     ! do n=1,maxnstreams

9999 CONTINUE

  IF (lhook) CALL dr_hook('DUMP_STREAMS_INIT',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE dump_streams_init

