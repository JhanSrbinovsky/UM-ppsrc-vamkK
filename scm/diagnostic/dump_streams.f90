! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Write SCM diagnostic data to output files

SUBROUTINE dump_streams                                                 &
  ( SCMop, day, time, row_length, rows, model_levels, dayno_init        &
  , steptime, site, lat, lon, timeinit, year, lcal360 )

  USE netcdf
  USE UM_types
  USE scm_utils, ONLY:                                                  &
      scm_timestep, scm_timestep_count, time_info, scm_nc_time_id       &
    , scm_nc_date_id, scm_nc_hour_id, scm_nc_local_id, zhook_in         &
    , zhook_out, jprb, lhook, dr_hook

! SCMop_type is defined in here...
  USE scmoptype_defn

  IMPLICIT NONE

! Description:
!   Write the diagnostic output data to file for all streams
!   with a dumping period ending now.

! Method:

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran77 with some bits of Fortran90

! NOTE: on the NEC this routine may be compiled in 32 bits
!       because requests to create 32 bit integers (necessary to
!       communicate with the NetCDF library routines - see incdf) are
!       ignored by the compiler when compiling at 64 bits. Thus all
!       input and output variables are explicitly declared as 64 bit
!       using the KIND types i64, r64 and l64 delcared in scmoptype_defn module

  TYPE(SCMop_type) :: SCMop ! In The derived-type structure
                            ! containing all the diagnostic
                            ! information

  INTEGER(i64) :: model_levels ! In The number of model levels
  INTEGER(i64) :: day          ! In Day
  INTEGER(i64) :: time         ! In Time in seconds


  INTEGER :: unit

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

  INTEGER :: n,m,i,j,k,l       ! Counters
  INTEGER(i64) :: n64          ! 64-bit copy of N
  INTEGER :: dump              ! The dump number (1:ndumps)
  INTEGER :: ncols,nrows,nlevs ! No. of columns,rows,levels of a
                               ! diagnostic
  INTEGER :: nlines            ! Number of lines written to file
  INTEGER :: dump_step         ! A dumping period

  ! Local variables for NetCDF-format output. Integers that are
  ! passed into NetCDF library routines should have the same
  ! precision as those used internally within the library (which
  ! may have been compiled at 32 bit), hence the use of incdf.
  INTEGER(incdf) :: Status,shape(4)

  ! Local variables for FSSI-format output
  TYPE fssi_dump_type
     INTEGER :: ndiags
     INTEGER, POINTER :: diags(:)
     INTEGER :: dumpsize
  END TYPE fssi_dump_type

  TYPE(fssi_dump_type) :: fssi_dump(maxnstreams)
  SAVE fssi_dump

  INTEGER, ALLOCATABLE :: diags(:)
  INTEGER :: ndiags

  INTEGER(i64) :: row_length,rows
  INTEGER(i64) :: dayno_init            ! Initial value of DAYNUMBER
  INTEGER(i64) :: steptime              ! Elapsed time (=stepcount*timestep)
  INTEGER(i64) :: site(row_length,rows) ! SSFM site WMO number
  INTEGER(i64) :: timeinit              ! Start time of integration
  INTEGER(i64) :: year                  ! Year (4 digits!!!)

  REAL(r64) :: lat(row_length,rows)     ! SSFM site latitude
  REAL(r64) :: lon(row_length,rows)     ! SSFM site longitude

  LOGICAL(l64) :: lcal360               ! Flag for 360 day calendar

  INTEGER(i64) ::  &
    ElapsedDays    &! Whole days elapsed from start of run
  , ElapsedSecs    &! Whole seconds elapsed from start of run
  , BasisTimeDays  &! Whole days to basis time from start of calendar
  , BasisTimeSecs  &! Seconds in day at basis time
  , CurrentYear    &! Year at current time in model run
  , Month          &! Month at current time in model run
  , Day2           &! Day at current time in model run
  , DayNumber      &! Day number at current time in model run
  , Hour           &! Hour at current time in model run
  , Minute         &! Minute at current time in model run
  , Second          ! Second at current time in model run

  INTEGER ::       &
    RecLen         &! File record length
  , isec           &! Second at current time in model run
  , Forecastime    &! Forecast time in `hundred' hours
  , DumpSize        ! Size of dump containing diagnostics required for output
                    ! to file

  INTEGER, PARAMETER ::   &
    HdrOffset = 11*12      ! Offset for heading at start of each output line,
                           ! containing site, date, time, etc.

  INTEGER ::              &
    datestamp             &! Date stamp for netcdf output file
  , localstamp             ! Local time stampfor netcdf output file

  CHARACTER(LEN=200) ::       &
    fmt605                &
  , fmt620                &
  , fmt650

  CHARACTER (len=100) ::  &
    fmt                   &! Hold format specifiers
  , fmt2                   !

  CHARACTER (len=3) ::    &
    clsname                ! Holds the short-name length as string

  ! Character function that incorporates the substep number into
  ! the short name of a diagnostic. Somewhat longer than a normal
  ! short name.
  CHARACTER (len=lsname+10) :: add_substep_to_sname,short_name

  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('DUMP_STREAMS',zhook_in,zhook_handle)

  ! Reset no. of calls to SCMoutput (this should be done at the
  ! end of every timestep since it acts as the number of calls
  ! to SCMoutput within the current timestep)
  SCMop%nSCMoutput = 0

  ! Loop over the streams
  DO m=1, maxnstreams

    ! If this stream is not switched on then skip to the next
    IF (SCMop%strm(m)%switch == 0) THEN
      CYCLE
    END IF

    ! If there are no diagnostics in this stream then skip to
    ! the next
    IF (SCMop%strm(m)%n_output <= 0) THEN
      CYCLE
    END IF

    ! If we are not at the end of the dumping period of this stream
    ! (or, in the case of FSSI-format streams, it's not the first
    ! timestep) then skip to the next.
    dump_step = SCMop%strm(m)%dump_step

    IF (MOD(stepnmbr(SCMop)-1,dump_step)+1 /= dump_step                 &
      .AND..NOT.                                                        &
      (SCMop%strm(m)%format == 2.AND.stepnmbr(SCMop) == 1)) THEN
      CYCLE
    END IF

    ! The output unit for this stream
    unit = SCMop%strm(m)%op_unit

    ! What we write to the unit depends on the stream format

!-----------------------------------------------------------------------
! Formats 0 and 3: designed for reading by PV-wave routine scmread2.pro
!-----------------------------------------------------------------------
     IF (SCMop%strm(m)%format == 0.OR.SCMop%strm(m)%format == 3) THEN

       ! Write the day and time
       WRITE(unit,'(I9,1PE18.10E3)') day,time

       ! Write the diagnostic data
       DO N=1, SCMop%nentries

         ! Is this entry to be sent to this stream, and
         ! not been flagged as one not to output?
         IF (StreamIsOn(SCMop%streams(n),m).AND..NOT.                   &
             NotWritten(SCMop%streams(n))) THEN
           ! Yes, write it.
           WRITE(unit,'(5(1PE18.10E3,","))')                            &
                (SCMop%diag(N)%dump(i),                                 &
                i=1,SCMop%nelements(n))
         END IF
       END DO

!-------------------------------------------------------------------
! Format=1 : designed for easy visual perusal.
!-------------------------------------------------------------------
     ELSE IF (SCMop%strm(m)%format == 1) THEN

       ! Which dump is this?
       dump = stepnmbr(SCMop)/SCMop%strm(m)%dump_step

       ! Write the dump and step numbers
       WRITE(unit,*) 'Dump= ', dump, ' ; Timestep=', stepnmbr(SCMop)

       ! Write the day and time
       WRITE(unit,*) 'Day= ', day, '; Time= ', time

       ! We will keep tabs on the number of lines we've
       ! written in this dump
       nlines = 0

       ! Write the diagnostic data
       DO N=1, SCMop%nentries

         ! Is this entry to be written to this stream's output file?
         IF (StreamIsOn(SCMop%streams(n),m).AND..NOT.                   &
             NotWritten(SCMop%streams(n))) THEN

           ! Yes.

           ! Get the number of columns, rows and levels
           ! this diagnostic is defined on
           ncols = SCMop%ncols(n)
           nrows = SCMop%nrows(n)
           nlevs = SCMop%nlevs(n)

           ! Ensure that the sub-step is appended to the short
           ! name of each diagnostic entry, if required.
           n64 = n

! DEPENDS ON: add_substep_to_sname
           short_name = add_substep_to_sname(SCMop,n64)

           ! Create a format to get all the numbers
           ! corresponding to one point on one line
           ! (i.e. all levels on one line)
           WRITE(fmt ,'(A,I3,A)') '(I3,1X,A,',                          &
                 nlevs, '(1PE18.10E3,","))'
           WRITE(clsname,'(I3)')lsname-6
           WRITE(fmt2,'(A,A,A,I3,A)')                                   &
                '(',clsname,'X,A,I4,A,',model_levels,                   &
                '(" Level=",I3,8X," "))'

           DO i=1, ncols
             DO j=1, nrows
               IF (MOD(nlines,15) == 0) THEN
                 WRITE(unit,fmt2)'Dump=',dump,':',                      &
                      (k,k=1,model_levels)
               END IF

               nlines = nlines+1
               WRITE(unit,fmt)SCmop%sname_id(n),                        &
                     short_name,(SCMop%diag(n)%dump                     &
                     ((k-1)*nrows*ncols+(j-1)*ncols+i),                 &
                     k=1,nlevs)
             END DO
           END DO
         END IF
       END DO

!-------------------------------------------------------------------
! Format=2 : format required for FSSI database
!-------------------------------------------------------------------
     ELSE IF (SCMop%strm(m)%format == 2) THEN

       ! Calculate times and dates required for output
       ElapsedDays = dayno_init + SCMop%daycount-1                      &
             + INT((TimeInit + StepTime) / 86400 )
       ElapsedSecs = MOD(TimeInit + StepTime, INT(86400,i64))

       ! Calculate Model Basis Time (i.e. begining of current year)
       ! relative to calender zero
! DEPENDS ON: time2sec
       CALL time2sec                                                    &
         ( Year, 1, 0, 0, 0, 0, 0, 0, BasisTimeDays, BasisTimeSecs      &
         , LCal360 )

        ! Calculate full date at current time in model integration
! DEPENDS ON: sec2time
       CALL sec2time                                                    &
         ( ElapsedDays, ElapsedSecs, BasisTimeDays, BasisTimeSecs       &
         , CurrentYear, Month, Day2, Hour, Minute, Second               &
         , DayNumber, LCal360 )

       IF (SCMop%stepcount == 1 .AND. SCMop%daycount == 1) THEN
         ISec = 0
       END IF

       ! Calculate forecast time in `hundred' hours
       Forecastime = ((SCMop%daycount-1)*24+INT(StepTime/3600))*100     &
             +MOD(INT(StepTime/60.0),60)

       ! If this is the first timestep then make a list of the
       ! diagnostics for output to this stream, diags(1:ndiags),
       ! together with the total size of one dump for one site,
       ! DumpSize. Save this info in the structure fssi_dump(m)
       ! so we don't have to re-calculate it next time.
       IF (stepnmbr(SCMop) == 1) THEN
         ndiags   = 0
         DumpSize = 0
         ALLOCATE(diags(SCMop%nentries))

         ! Unless heed_acceptlist=2 and the acceptlist has
         ! non-zero size then output all diagnostics in the same
         ! order as the calls to SCMoutput.
         IF (SCMop%strm(m)%heed_hardwired /= 0.OR.                      &
             SCMop%strm(m)%heed_acceptlist == 0) THEN

           DO n=1, SCMop%nentries

             ! Is this entry to be written to this stream's output file?
             IF (StreamIsOn(SCMop%streams(n),m).AND..NOT.               &
                 NotWritten(SCMop%streams(n))) THEN

               ! Yes.
               ndiags = ndiags+1
               diags(ndiags) = N

               ! This diagnostic will contribute a number of
               ! elements to the dump for one site equal to
               ! the number of levels it's defined on.
               DumpSize = DumpSize+SCMop%nlevs(n)

             END IF
           END DO
         ELSE

           ! Output only those diagnostics specified by the
           ! accept_list and do it in the order given in the
           ! accept_list.
           DO i=1, SIZE(SCMop%strm(m)%accept_list)
             DO n=1, SCMop%nentries

               ! Does the entry have the right short name?
               IF (TRIM(SCMop%sname(n)) ==                              &
                   TRIM(SCMop%strm(m)%accept_list(i))) THEN

                 ! Is this entry to be written to this
                 ! stream's output file?
                 IF (StreamIsOn(SCMop%streams(n),m).AND..NOT. &
                     NotWritten(SCMop%streams(n))) THEN

                   ! Yes.
                   ndiags = ndiags+1
                   diags(ndiags) = N

                   ! This diagnostic will contribute a
                   ! number of elements to the dump for
                   ! one site equal to the number of levels
                   ! it's defined on.
                   DumpSize = DumpSize+SCMop%nlevs(n)

                 END IF
               END IF
             END DO
           END DO
         END IF

        ! Store ndiags, diags(1:ndiags) and DumpSize in
        ! fssi_dump(m) so we don't have to execute this section
        ! of code again.
        fssi_dump(m)%ndiags = ndiags
        ALLOCATE(fssi_dump(m)%diags(ndiags))
        fssi_dump(m)%diags = diags(1:ndiags)
        fssi_dump(m)%dumpsize = dumpsize
        DEALLOCATE(diags)
      END IF               ! if (stepnmbr(SCMop) == 1)

      ! Format statements for output
      fmt605=                                                           &
         '(2a12,2(i2.2,''/''),i4.4,1x,'//                               &
         'i2.2,'':'',i2.2,''.'',i2.2,'' Diags'',i4,'' RL'',i8)'
      WRITE(fmt620,                                                     &
         fmt='("(2x,''Namelist'',7x,''Day'',8x,''Month'',"'//           &
             '"7x,''Year'',8x,''Hour'',"'                  //           &
             '"7x,''Minute'',6x,''Second'',9x,''FT'',8x,'  //           &
             '''Site'',7x,''Lat'',8x,''Long'',4x,",'       //           &
             'I4,"(5x,i3,'':'',i3))")') fssi_dump(m)%dumpsize
      WRITE(fmt650,                                                     &
         fmt='("(a12,8(i12),2F12.5,",I4.4,"(g12.5))")')                 &
         fssi_dump(m)%dumpsize

      ! Write output to file, one line for each site at each
      ! time, all diagnostics (single and multiple level) on
      ! the same line
      DO i=1, row_length
        DO j=1, rows

          ! Perform initialisation of the file if this is the
          ! first (i,j) cycle on the first timestep
          IF (SCMop%stepcount == 1.AND.SCMop%daycount == 1              &
             .AND.i == 1.AND.j == 1) THEN

            RecLen = fssi_dump(m)%dumpsize*12 + HdrOffset
            OPEN(Unit=unit,Recl=RecLen,Access='Sequential',             &
                 File=SCMop%strm(M)%filename)

            ! Write out initial headings for outputs to file
            ! (a dummy string is substituted for DiagRefFile at
            ! this stage)
            WRITE(unit,fmt=fmt605)'SSM         ',                       &
                 'XXXXXXXXXXXX',day2,month,currentyear,                 &
                 Hour,Minute,ISec,SCMop%strm(m)%n_output,reclen

            WRITE(unit,fmt=fmt620)                                      &
                 ((SCMop%sname_id(fssi_dump(m)%diags(l)),k,             &
                 k=1,SCMop%nlevs(fssi_dump(m)%diags(l))),               &
                 l=1,fssi_dump(m)%ndiags)
          END IF

          ! Write out diagnostics data to file
          ! (a dummy string is substituted for DiagRefFile at
          ! this stage).
          WRITE(unit,fmt=fmt650)'XXXXXXXXXXXX',Day2,Month,              &
               CurrentYear,Hour,Minute,Second,Forecastime,              &
               Site(i,j),Lat(i,j),Lon(i,j),                             &

          ! What follows is a little bit complicated: two
          ! nested implied do-loops. The inner index is k,
          ! representing vertical level, and the outer is l,
          ! representing a given output diagnostic. For a given
          ! k and l, what is written out is the element of the
          ! dump array of lth output diagnostic corresponding
          ! to point(i,j) and level k. Get it?
               (                                                        &
                (                                                       &
                 SCMop%diag(fssi_dump(m)%diags(l))%dump                 &
                  (                                                     &
                   (k-1)*SCMop%nrows(fssi_dump(m)%diags(l))*            &
                    SCMop%ncols(fssi_dump(m)%diags(l)) +                &
                   (j-1)*SCMop%ncols(fssi_dump(m)%diags(l))+i           &
                  ),                                                    &
                 k=1,SCMop%nlevs(fssi_dump(m)%diags(l))                 &
                ),                                                      &
                l=1,fssi_dump(m)%ndiags                                 &
               )
        END DO            ! j=1,rows
      END DO               ! i=1,row_length

!-------------------------------------------------------------------
! Format=4 : NetCDF
!-------------------------------------------------------------------
    ELSE IF (SCMop%strm(m)%format == 4) THEN

      ! Which dump is this?
      dump = stepnmbr(SCMop)/SCMop%strm(m)%dump_step

      datestamp = time_info(1)*10000  &
        + time_info(2)*100            &
        + time_info(3)

      localstamp = time_info(4)*10000 &
        + time_info(5)*100            &
        + time_info(6)

      ! Write out time variable
      status = nf90_put_var                           &
        ( ncid   = INT(scmop%strm(m)%op_unit, incdf)  &
        , varid  = scm_nc_time_id                     &
        , start  = (/INT(dump, incdf)/)               &
        , count  = (/INT(1,incdf)/)                   & 
        , values = (/scm_timestep_count*scm_timestep/))

      ! Write out date variable
      status = nf90_put_var                           &
        ( ncid   = INT(scmop%strm(m)%op_unit, incdf)  &
        , varid  = scm_nc_date_id                     &
        , start  = (/INT(dump, incdf)/)               &
        , count  = (/INT(1,incdf)/)                   & 
        , values = (/INT(datestamp, incdf)/)          )

      ! Write out hour variable
      status = nf90_put_var                           &
        ( ncid   = INT(scmop%strm(m)%op_unit, incdf)  &
        , varid  = scm_nc_hour_id                     &
        , start  = (/INT(dump, incdf)/)               &
        , count  = (/INT(1,incdf)/)                   & 
        , values = (/FLOAT(time_info(4))/)            )

      ! Write out local time variable
      status = nf90_put_var                           &
        ( ncid   = INT(scmop%strm(m)%op_unit, incdf)  &
        , varid  = scm_nc_local_id                    &
        , start  = (/INT(dump, incdf)/)               &
        , count  = (/INT(1,incdf)/)                   & 
        , values = (/INT(localstamp, incdf)/)         )

      ! Write out each diagnostic
      DO N=1,SCMop%nentries

        ! Is this entry to be sent to this stream, and
        ! not been flagged as one not to output?
        IF (StreamIsOn(SCMop%streams(n),m).AND..NOT.                    &
            NotWritten(SCMop%streams(n))) THEN

          ! Yes. Get the proper shape of the variable (at the
          ! moment it's being held as a 1D array).
          shape = INT((/SCMop%ncols(n),                                 &
                        SCMop%nrows(n),                                 &
                        SCMop%nlevs(n),                                 &
                        INT(1,i64)/),incdf)

          Status = Nf90_Put_Var(                                        &
            NcID   = INT(SCMop%strm(m)%op_unit,incdf),                  &
            VarID  = SCMop%netcdf_id(n),                                &
            start  = INT((/1,1,1,dump/),incdf),                         &
            count  = shape,                                             &
            Values = RESHAPE(SCMop%diag(n)%dump,shape))
        END IF
      END DO

    ELSE
      PRINT*,'Dump_Streams ERROR: unknown format for stream ', m

    END IF                  ! if (SCMop%strm(m) == ...)

  END DO ! M=1,maxnstreams

  IF (lhook) CALL dr_hook('DUMP_STREAMS',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE dump_streams

