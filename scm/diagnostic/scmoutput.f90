! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Create a SCM output diagnostic

SUBROUTINE scmoutput                                                          &
  ( x, sname, lname, units, timprof, domprof, streams, sname2                 &
  , calling_routine )

  USE global_scmop
  USE scm_utils, ONLY: scm_trap_nan

  IMPLICIT NONE

! Description:
!   Create output diagnostic based on given inputs. This routine
!   cannot be called more than once per timestep with the same
!   "sname" unless the call is inside a sub-stepped part of the
!   model and the sub-stepping has been delimited by calls to
!   scm_substep_start and scm_substepping_end. The order of calls
!   to SCMoutput should not change between timesteps.

! Method:

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran77 with some bits of Fortran90

  INTEGER :: domprof  ! In Domain profile for the diagnostic
  REAL :: x  &        ! In Variable from which diagnostic will be constructed
      ((SCMop%d_rowa2(domprof)-SCMop%d_rowa1(domprof)+1)*           &
       (SCMop%d_rowb2(domprof)-SCMop%d_rowb1(domprof)+1)*           &
       (SCMop%d_lev2 (domprof)-SCMop%d_lev1 (domprof)+1))

  CHARACTER(len=*) :: &
    sname             &! In Short name for the diagnostic,
                       !    this should be unique
  , lname             &! In Long name for the diagnostic
  , units             &! In Units of the diagnostic
  , sname2             ! In Short name of another, previously
                       !    defined diagnostic which will be used
                       !    in the construction of this one
                       !    according to the time profile

  CHARACTER(len=*) :: &
    calling_routine    ! In Routine that has called scmoutput

  INTEGER ::          &
    timprof           &! In The time profile for the diagnostic
  , streams            ! In An encoded integer specifying
                       !    which output streams the diagnostic
                       !    is to go to

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


  INTEGER :: d,i,j     ! General use

  CHARACTER(len=30) :: sdum0

  LOGICAL ::          &
    startperiod       &! Will be used to flag if we are at
  , endperiod          ! the start or end of a dumping period

  INTEGER ::                     &
    ndistinct_diags              &
  , distinct_diags(maxnstreams)  &
  , nThroughPeriod               &
  , ntrad                        &
  , ntrad1

  LOGICAL ::       &
    call_add2dump  &
  , order_changing

  ! Will hold the contents of SCMop%diag_mem while being
  ! re-allocated
  INTEGER, ALLOCATABLE :: itemp(:,:)

  ! Perform no action if SCMop is not turned on, except in the
  ! case of a constant diagnostic
  IF (.NOT. SCMop%on .AND. timprof /= t_const) GOTO 9999

  ! Enforce a rule that sname must be at least one character long
  IF (LEN(sname) == 0) THEN
    IF (SCMop%first_pass) THEN
      PRINT*,'SCMoutput ERROR: sname is a null string, this '//   &
           'is not allowed, diagnostic ignored: ',                &
           lname(1:LEN(lname))
    END IF
    GOTO 9999
  END IF

  ! Increment recorded no. of calls to this routine this timestep
  IF (SCMop%on) SCMop%nSCMoutput = SCMop%nSCMoutput+1

  ! If the requested list of streams span a range of dumping
  ! periods then the given inputs will correspond to more than one
  ! diagnostic entry in SCMop. Thus we need to get
  ! distinct_diags(1:ndistinct_diags) - the indices of those
  ! entries. The routine getdistinctdiags can do this but, if this
  ! is not the first timestep, and assuming the order of the calls
  ! to this routine doesn't change between timesteps, we can just
  ! use our memory of a previous timestep instead.

  order_changing = .FALSE.
  IF (.NOT.SCMop%first_pass) THEN
    ! Use memory of previous timesteps to know which diagnostic
    ! entries this call pertains to.
    ndistinct_diags = SCMop%diag_mem(SCMop%nSCMoutput,1)
    DO i=1, ndistinct_diags
      distinct_diags(i) = SCMop%diag_mem(SCMop%nSCMoutput,2)+i-1

      ! Check we're right...
      IF (TRIM(sname(1:MIN(lsname,LEN(sname)))) /=                  &
          TRIM(SCMop%sname(distinct_diags(i)))) THEN
        WRITE(6,*)'*******************************************'
        WRITE(*,*) 'SCMoutput warning: the order of the calls'//    &
                   ' to SCMoutput seems to be changing:'
        WRITE(*,'(A,I4,A)') '  On step ',stepnmbr(SCMop),           &
                  ' expected '//                                    &
        TRIM(SCMop%sname(distinct_diags(i)))//                      &
                  ', but encountered '//sname(1:MIN(lsname,LEN(sname)))
        WRITE(6,*)'*******************************************'
        order_changing = .TRUE.

      ELSE IF (SCMop%substep_number /=                              &
               SCMop%substep(distinct_diags(i))) THEN
        WRITE(6,*)'*******************************************'
        WRITE(6,'(A)') ' SCMoutput error: the order of the '  //    &
                       'sub-steps seems to be changing.'
        WRITE(6,'(A,I5,A,I3,A,I3)')                                 &
                ' On step ',stepnmbr(SCMop),' diagnostic '    //    &
                  sname(1:MIN(lsname,LEN(sname)))             //    &
                ' expected sub-step ',                              &
                  SCMop%substep(distinct_diags(i)),                 &
                ' but got sub-step ',SCMop%substep_number
        WRITE(6,'(A)') ' This does not make sense and is a '  //    &
                          'sign of a potentially serious problem.'
        WRITE(6,*)'*******************************************'
        order_changing = .TRUE.
      END IF
    END DO
  END IF

  IF (SCMop%first_pass.OR.order_changing) THEN

    ! Either this is the first timestep or the order of the calls
    ! to SCMoutput is changing

! DEPENDS ON: getdistinctdiags
    CALL getdistinctdiags                                       &
      ( sname, lname, units, timprof, domprof, streams, sname2  &! In
      , ndistinct_diags, distinct_diags, SCMop )      ! Out,Out,InOut

    ! Don't want to do this next bit for constant diagnostics
    ! (which should be declared when the system is off)
    IF (SCMop%on) THEN

      ! Store the values of ndistinct_diags and distinct_diags(1)
      ! for future reference (to avoid unnecessary calls to
      ! getdistinctdiags). But is there enough space in the
      ! SCMop%diag_mem array?
      IF (SIZE(SCMop%diag_mem,1) == SCMop%nSCMoutput-1) THEN
        ! No. Make the array bigger...
        ! Allocate a temporary array with the same size and
        ! shape as diag_mem
        ALLOCATE(itemp(SCMop%nSCMoutput-1,2))

        ! Copy the contents of diag_mem into it
        itemp = SCMop%diag_mem

        ! Re-allocate diag_mem with a larger size
        DEALLOCATE(SCMop%diag_mem)
        ALLOCATE(SCMop%diag_mem(SCMop%nSCMoutput+49,2))

        ! Copy the original contents back in
        SCMop%diag_mem(1:SCMop%nSCMoutput-1,1:2) = itemp
      END IF

      SCMop%diag_mem(SCMop%nSCMoutput,1) = ndistinct_diags
      SCMop%diag_mem(SCMop%nSCMoutput,2) = distinct_diags(1)
    END IF
  END IF

  ! From here on in, none of the input parameters are
  ! referred to at all, they have been distilled to
  ! distinct_diags(1:ndistinct_diags) and information
  ! in SCMop

  ntrad  = SCMop%ntrad       ! No. of timesteps between calls to rad'n
  ntrad1 = SCMop%ntrad1      ! Timestep containing 1st call to rad'n

  DO i=1, ndistinct_diags
    d = distinct_diags(i)

    ! If this is not the first time we've seen this diagnostic,
    ! check the last time was the previous timestep
    IF (SCMop%lastencounter(d) >= 0.AND.                            & 
        SCMop%lastencounter(d) /= stepnmbr(SCMop)-1) THEN
      PRINT*,'SCMoutput ERROR: last encounter with this '//         &
             'diagnostic was not last timestep: ',sname,            &
             SCMop%lastencounter(d),stepnmbr(SCMop),                &
             (SCMop%daycount-1),SCMop%full_daysteps,                &
             SCMop%stepcount
    END IF

    ! Record the fact that this diagnostic was seen by
    ! this routine on this timestep.
    SCMop%lastencounter(d) = stepnmbr(SCMop)

    ! Calculate how many timesteps we are through the current
    ! dumping period. 1=first time step, dump_step=last timestep
    ! (and so a dump will occur this timestep).
    nThroughPeriod = MOD(stepnmbr(SCMop)-1,SCMop%dump_step(d))+1

    ! Decide whether we are at the start of a dumping period, at
    ! the end of of a dumping period, and whether we need to call
    ! add2dump (using nThroughPeriod this is trivial for most
    ! diagnostics, but has added complications in the case of
    ! diagnostics only calculated on radiation timesteps)
    startperiod   = .FALSE.
    endperiod     = .FALSE.
    call_add2dump = .TRUE.

    IF (.NOT.SCMop%only_radsteps(d)) THEN
      ! Diagnostic d is a normal diagnostic, valid at every timestep
      IF (nThroughPeriod == 1) startperiod = .TRUE.
      IF (nThroughPeriod == SCMop%dump_step(d)) endperiod = .TRUE.

    ELSE
      ! Diagnostic d is based on a variable which only has
      ! valid values on radiation timesteps.

      IF (MOD(SCMop%stepcount-ntrad1,ntrad) /= 0) THEN

        ! This is not a radiation timestep, assume input
        ! array x contains nonsense information - do not
        ! call add2dump.
        call_add2dump = .FALSE.
      ELSE

        ! The criteria for startperiod and endperiod are now
        ! altered slightly, since startperiod must be true if
        ! this is the first radiation time step during this
        ! dumping period, and endperiod must be true if this is
        ! the last radiation time step during this dumping
        ! period.
        IF (nThroughPeriod-1 <  ntrad) startperiod = .TRUE.
        IF (SCMop%dump_step(d)-nThroughPeriod <  ntrad) endperiod = .TRUE.

        ! This is a radiation timestep and so we can call
        ! add2dump, but the no. of timesteps by which to divide
        ! in order to calculate the average (or whatever) is
        ! not dump_step, but the no. of times this part of the
        ! code has been reached during this dumping period,
        ! stored in SCMop%nadd2dump(d) (which for normal
        ! diagnostics is set to dump_step in newdiag)
   
        IF (startperiod) THEN
          SCMop%nadd2dump(d) = 1
        ELSE
          SCMop%nadd2dump(d) = SCMop%nadd2dump(d)+1
        END IF

      END IF            ! (mod(SCMop%stepcount-ntrad1,ntrad) /= 0)
    END IF             ! (.not.SCMop%only_radsteps(d))

    IF (call_add2dump) THEN
! DEPENDS ON: add2dump
      CALL add2dump (x, SCMop%nelements(d), d, SCMop, startperiod, endperiod)

      DO j=1, SCMop%nelements(d)
        WRITE(sdum0,*) x(j)

        IF (INDEX(sdum0, 'NaN') /= 0) THEN
          CALL scm_trap_nan (sname, calling_routine)
          EXIT
        END IF

        IF (INDEX(sdum0, 'nan') /= 0) THEN
          CALL scm_trap_nan (sname, calling_routine)
          EXIT
        END IF

        IF (INDEX(sdum0, 'NAN') /= 0) THEN
          CALL scm_trap_nan (sname, calling_routine)
          EXIT
        END IF

      END DO

    END IF

  END DO

9999 CONTINUE

  RETURN

END SUBROUTINE scmoutput
