! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Create a new SCM diagnostic entry in SCMop.

SUBROUTINE newdiag                                                            &
  ( sname, lname, units, timprof, domprof, istrm, lnot_written, sname2, d     &
  , SCMop )

! SCMop_type is defined in here...
  USE scmoptype_defn

  IMPLICIT NONE

! Description:
!   Create a new diagnostic entry in SCMop. Returns index of
!   newly created diagnostic, which is unchanged from its input
!   value if an error ocurrs and the entry could not be created.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran 90

  TYPE(SCMop_type) :: SCMop       ! InOut The derived-type structure
                                  !       containing all the diagnostic
                                  !       information
  CHARACTER (len=lsname) :: sname ! In Short name for the diagnostic
  CHARACTER (len=*) :: lname      ! In Long name for the diagnostic
  CHARACTER (len=*) :: units      ! In Units of the diagnostic
  CHARACTER (len=*) :: sname2     ! In Short name of a previously
                                  !    defined diagnostic used in the
                                  !    construction of this one
  INTEGER ::  &
    timprof   &! In Time profile for the diagnostic
  , domprof   &! In Domain profile for the diagnostic
  , istrm      ! In Stream to which diagnostic is to be sent

  LOGICAL ::  &
    lnot_written ! In If true, diagnostic will not be written out

  INTEGER ::  &
    d            ! InOut Index of the newly created
                 !       diagnostic entry in SCMop. Unchanged from
                 !       input value if entry could not be created.

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

  INTEGER :: d_input,i,d1,d2,d3

  LOGICAL ::                     &
    sname2_found                 &
  , right_dumping_period

  d_input = d ! record input value of d

  IF (SCMop%nentries == SCMop%maxnentries) THEN

    ! No room at the inn, make the inn bigger.
! DEPENDS ON: expand_scmop
    CALL expand_scmop(SCMop)
  END IF

  ! Increment the recorded number of entries
  SCMop%nentries = SCMop%nentries+1

  ! d is the index of this new entry
  d = SCMop%nentries

  ! Record the details of this entry in the respective arrays...
  SCMop%sname(d)   = sname
  SCMop%lname(d)   = lname
  SCMop%units(d)   = units
  SCMop%domprof(d) = domprof

  IF (timprof <  only_radsteps) THEN

    ! This is a normal diagnostic - calculated on every timestep
    SCMop%timprof(d) = timprof
    SCMop%only_radsteps(d) = .FALSE.
  ELSE

    ! This diagnostic is based on an array only valid on
    ! radiation timesteps
    SCMop%timprof(d) = timprof-only_radsteps
    SCMop%only_radsteps(d) = .TRUE.
  END IF

  SCMop%streams(d) = Stream(istrm)
  IF (lnot_written) THEN
    SCMop%streams(d) = DoNotWrite(SCMop%streams(d))
  END IF

  ! The dumping period of the diagnostic is that of the stream
  ! it's being sent to
  SCMop%dump_step(d) = SCMop%strm(istrm)%dump_step
  SCMop%nadd2dump(d) = SCMop%strm(istrm)%dump_step

  ! The dimensions of the diagnostic array can be obtained from
  ! the domain profile
  d1 = SCMop%d_rowa2(domprof)-SCMop%d_rowa1(domprof)+1
  d2 = SCMop%d_rowb2(domprof)-SCMop%d_rowb1(domprof)+1
  d3 = SCMop%d_lev2 (domprof)-SCMop%d_lev1 (domprof)+1
  SCMop%ncols(d) = d1
  SCMop%nrows(d) = d2
  SCMop%nlevs(d) = d3
  SCMop%nelements(d) = d1*d2*d3

  ! Allocate the space for the dump array
  ALLOCATE(SCMop%diag(d)%dump(d1*d2*d3))
  ! Initialise it for initialisation's sake
  SCMop%diag(d)%dump = -999.0

  ! lastencounter will be set in SCMoutput
  SCMop%lastencounter(d) = -1

  ! Set the substep number that this entry is being created for
  SCMop%substep(d) = SCMop%substep_number

  ! wd will be the index of a diagnostic entry upon which this
  ! diagnostic depends (set to zero if sname2 is a null string)
  SCMop%wd(d) = 0

  IF (LEN(sname2) >  0) THEN
    ! Find the index of the weighting diagnostic
    sname2_found = .FALSE. ! (flags if at least found the right name)
    right_dumping_period = .FALSE. ! (" " " " right dumping period)

    DO i=1, d-1
      IF (SCMop%sname(i) == sname2) THEN

        ! We have found a diagnostic of the correct name
        sname2_found = .TRUE.

        ! But it must have the same dumping period or be constant
        IF (SCMop%dump_step(i) == SCMop%strm(istrm)%dump_step    &
          .OR.                                                   &
            SCMop%timprof(i) == t_const) THEN

          right_dumping_period = .TRUE.

          ! But it must also be defined on the same substep as
          ! the current entry, or have been defined outside of
          ! a sub-stepped part of the model.

          IF (SCMop%substep(i) == SCMop%substep_number .OR.     &
              SCMop%substep(i) == 0) THEN

            ! We have found an entry with the correct name & substep
            SCMop%wd(d) = i
            EXIT

          END IF
        END IF
      END IF
    END DO

    IF (SCMop%wd(d) == 0) THEN
      WRITE(*,*)' '
      WRITE(*,*)'************ERROR IN ROUTINE NEWDIAG************'
      WRITE(*,*)'* The following error stems from a call to '//   &
                'SCMoutput...'
      WRITE(*,*)'* You have requested that diagnostic "',         &
                TRIM(sname),'" be dependent '
      WRITE(*,*)'* on diagnostic "',TRIM(sname2),                 &
                '" (which is non-constant),'
      IF (.NOT.sname2_found) THEN
        WRITE(*,*)'* but the latter diagnostic has not yet '//    &
                  'been defined.'
        WRITE(*,*)'* Diagnostic "',TRIM(sname),'" will '//        &
                  'therefore not be calculated or output.'
        WRITE(*,*)'* Note: this message will be '//               &
                  'repeated for every stream you have '
        WRITE(*,*)'* requested for this diagnostic'

      ELSE IF (.NOT. right_dumping_period) THEN
        WRITE(*,*)'* but you have requested that the former '//   &
                  'be sent to a stream with a'
        WRITE(*,*)'* dumping period of ',                         &
                  SCMop%strm(istrm)%dump_step,                    &
                  ', while the latter is not. This is not'
        WRITE(*,*)'* permitted: if diagnostic "A" is to '//       &
                  'depend on diagnostic "B", and'
        WRITE(*,*)'* diagnostic "B" is not constant, then '//     &
                  'diagnostic "B" must be calculated'
        WRITE(*,*)'* for every dumping period that is '//         &
                  'requested for diagnostic "A". This'
        WRITE(*,*)'* can be guaranteed by sending "B" to all '//  &
                  'the streams which you have'
        WRITE(*,*)'* requested for "A".'
        WRITE(*,*)'* Diagnostic "',TRIM(sname),'" will '//        &
                  'therefore not be calculated or'
        WRITE(*,*)'* output with a dumping period of ',           &
                  SCMop%strm(istrm)%dump_step
      ELSE IF (SCMop%substep_number /= 0) THEN
        WRITE(6,*)'* but while that diagnostic appears to be ' // &
                  'defined within a sub-stepped'
        WRITE(6,*)'* part of the code, no entry could be '     // &
                  'found for it on the current substep '
        WRITE(6,*)'* of ',SCMop%substep_number
        WRITE(*,*)'**************************************'//      &
                  '**********'
        WRITE(*,*)' '
      ELSE
        WRITE(6,*)'* but the latter is defined within a '      // &
                  'sub-stepped part of the '
        WRITE(6,*)'* code while the former is not. This does ' // &
                  'not make sense.'
      END IF

      WRITE(6,*)'**************************************'       // &
                '**********'
      WRITE(6,*)' '
      d = d_input
      SCMop%nentries = SCMop%nentries-1
      GOTO 9999
    END IF
  ELSE
    IF (timprof == t_div.OR.timprof == t_mult.OR.                 &
        timprof == t_acc_div.OR.                                  &
        timprof == t_acc_mult) THEN
      WRITE(*,'(A,1X,A,1X,A,1X,I3)')                              &
           'newdiag ERROR: you have requested this '//            &
           'diag to be dependent on another but have not '//      &
           'specified which.',TRIM(sname),timprof
      d = d_input
      SCMop%nentries = SCMop%nentries-1
      GOTO 9999
    END IF
  END IF

  ! Assign this diagnostic an integer unique to its sname
  ! (do this by searching for a previously defined diagnostic
  ! with the same sname, if you find it give it that number,
  ! if you don't give it the highest number you came across
  ! plus one)
  SCMop%sname_id(d) = 0

  DO i=1, SCMop%nentries-1
    IF (SCMop%sname(i) == sname) THEN

      ! i has the same sname as d, give it the same sname_id
      SCMop%sname_id(d) = SCMop%sname_id(i)
      EXIT
    ELSE

      ! record the largest sname_id, but store as -ve to
      ! indicate that we have not found a match for sname
      SCMop%sname_id(d) = MIN(SCMop%sname_id(d),-SCMop%sname_id(i))
    END IF
  END DO

  IF (SCMop%sname_id(d) <= 0) THEN

    ! A diagnostic with the same sname was not found: assign a
    ! value one larger than the current largest value of sname_id
    SCMop%sname_id(d) = -SCMop%sname_id(d)+1

    ! Record the number of diagnostics that will be output (with
    ! no double counting from the same diagnostic being
    ! calculated with different dumping periods)
    IF (.NOT.lnot_written) SCMop%n_output = SCMop%n_output+1
  END IF

  ! Record the number of diagnostics that will be sent to
  ! this stream
  IF (.NOT.lnot_written) THEN 
    SCMop%strm(istrm)%n_output = SCMop%strm(istrm)%n_output+1
  END IF

9999 CONTINUE

  RETURN

END SUBROUTINE newdiag

