! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Obtain indices of diagnostic entries in SCMop given SCMoutput inputs

SUBROUTINE getdistinctdiags                                                   &
  ( sname_i, lname, units, timprof, domprof, streams, sname2_i                &
  , ndistinct_diags, distinct_diags, SCMop )

! SCMop_type is defined in here...
  USE scmoptype_defn

  IMPLICIT NONE

! Description:
!   On the first call to this routine with a given sname_i, create
!   the corresponding diagnostic entries in SCMop and return their
!   indices in distinct_diags(1:ndistinct_diags). On subsequent
!   calls, simply look up the previously created entries and return
!   in distinct_diags(1:ndistinct_diags).

! Method:
!   For all streams to which the diagnostic is to be sent, look
!   for a corresponding entry in SCMop with the correct dump_step.
!   If one does not exist create it by calling NEWDIAG. Record its
!   index if it's different to the index found in any previous
!   cycle of the loop over streams. Thus build up a list of the
!   entries in SCMop to which this diagnostic is associated, and
!   return.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran90

  ! See SCMoutput for a description of the input variables
  CHARACTER(len=*) :: &! In
    sname_i           &! In
  , lname             &! In
  , units             &! In
  , sname2_i           ! In

  INTEGER ::  &
    timprof   &! In
  , domprof   &! In
  , streams    ! In

  INTEGER ::                    &
    ndistinct_diags             &! Out
  , distinct_diags(maxnstreams)  ! Out

  TYPE(SCMop_type) :: SCMop ! InOut The derived-type structure
                            !       containing all the diagnostic
                            !       information

! Parameters and stuff used around and about the internal
! workings of the diagnostic system...
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

  ! d will equal this if the entry could not be created
  INTEGER, PARAMETER :: unset = 0

  INTEGER :: i,j,d,initial_nentries

  CHARACTER(len=lsname) :: sname,sname2

  LOGICAL :: distinct

  ndistinct_diags = 0

  ! Elsewhere we are going to rely on the "sname" for any
  ! diagnostic being exactly lsname characters long, so here we
  ! make them the right length.
  sname  = sname_i
  sname2 = sname2_i

  IF (LEN(sname_i) >  lsname.AND.SCMop%first_pass) THEN

    ! Warn the user that we have done this
    PRINT*,'GetDistinctDiags WARNING: diagnostic name truncated: '            &
          ,sname_i(1:LEN(sname_i)),' -> ', sname
  END IF

  ! If this is the first timestep, check this "sname" has not been
  ! used in a previous call to SCMoutput
  IF (SCMop%first_pass) THEN
    DO i=1, SCMop%nentries
      IF (SCMop%sname(i) == sname .AND.                                       &
          SCMop%substep(i) == SCMop%substep_number) THEN
        PRINT*,'GetDistinctDiags ERROR: same sname used in >1 '//             &
               'call to SCMoutput on same sub-step:',sname
        GOTO 9999
      END IF
    END DO
  END IF

  ! Get list of streams to which the diagnostic is to be
  ! sent. Nominally this is determined by the variable, streams
  ! (the so-called "hard-wired" list of streams). But namelist
  ! information can alter this by choosing to ignore the
  ! hard-wired list and/or to add/remove diagnostics to/from
  ! the list. The two lists of extra diagnostics to be added to,
  ! and rejected from, stream X are in SCMop%strm(X)%accept_list and
  ! SCMop%strm(X)%reject_list respectively.

  streamlist = 0 ! This will be the final list of streams encoded into
                 ! one integer. For now it is set to "no streams".

  ! Loop over all streams.
  DO j=1, maxnstreams

    ! If this stream is closed do nothing.
    IF (SCMop%strm(j)%switch == 0) THEN
      CYCLE ! Go to next value of j.
    END IF

    ! Shall we pay attention to the list of streams as
    ! specified in the call to SCMoutput?
    IF (SCMop%strm(j)%heed_hardwired /= 0) THEN

      ! Yes. Add stream j if requested in the call.
      IF (StreamIsOn(streams,j).AND..NOT.StreamIsOn(streamlist,j)) THEN
        streamlist = streamlist+Stream(j)
      END IF

      ! If it has been requested in the call that this
      ! diagnostic should not be written out, make sure this
      ! is so.
      IF (NotWritten(streams).AND..NOT.NotWritten(streamlist)) THEN
        streamlist = DoNotWrite(streamlist)
      END IF
    END IF

    ! Shall we pay attention to the list of diagnostics
    ! requested by namelist to be sent to this stream?
    IF (SCMop%strm(j)%heed_acceptlist /= 0) THEN
      ! Yes.

      ! If the diagnostic is not already being sent to this stream
      IF (.NOT.StreamIsOn(streamlist,j)) THEN

        ! then check if it's in the list of extra diagnostics
        ! to suck into this stream.

        ! Loop over the the names of diagnostics to be accepted
        DO i=1, SIZE(SCMop%strm(j)%accept_list)
          ! Do we have a name match?
          IF (TRIM(SCMop%strm(j)%accept_list(i)) == TRIM(sname)) THEN

            ! The name was in the list for stream j. Send
            ! this diagnostic to stream j.
            streamlist = streamlist+Stream(j)
            EXIT   ! this inner loop
          END IF
        END DO

      END IF

    END IF ! SCMop%strm(j)%heed_acceptlist /= 0

    ! Shall we pay attention to list of diagnostics requested
    ! by namelist to be prevented from going to this stream?
    IF (SCMop%strm(j)%heed_rejectlist /= 0) THEN
      ! Yes.

      ! If the diagnostic is being sent to this stream
      IF (StreamIsOn(streamlist,j)) THEN

        ! then check if it's in the list of diagnostics
        ! which are not to be sent to this stream
        DO i=1, SIZE(SCMop%strm(j)%reject_list)
          IF (TRIM(SCMop%strm(j)%reject_list(i)) == TRIM(sname)) THEN

            ! This diagnostic is not to be sent to this stream.
            streamlist = streamlist-Stream(j)
            EXIT
          END IF
        END DO

      END IF ! StreamIsOn(streamlist,j)

    END IF    ! SCMop%strm(j)%heed_rejectlist /= 0

  END DO       ! j=1,maxnstreams

  ! From now on we will use streamlist instead of streams as the
  ! integer representing the list of streams to which the
  ! diagnostic is to be sent.

  initial_nentries = SCMop%nentries

  ! Loop over all output streams
  DO j=1, maxnstreams
    ! j represents a stream, is the diagnostic to go
    ! to this stream?
    IF (StreamIsOn(streamlist,j).AND.SCMop%strm(j)%switch /= 0) THEN

      ! Yes, look to see if the entry already exists in SCMop
      d = unset ! (label d as unset for now)

      DO i=1, SCMop%nentries
        IF (SCMop%sname(i) == sname.AND.                                      &
            SCMop%substep(i) == SCMop%substep_number .AND.                    &
            SCMop%dump_step(i) == SCMop%strm(j)%dump_step) THEN
          ! This diagnostic exists.
          d = i
          EXIT
        END IF
      END DO

      IF (d == unset) THEN
        IF (SCMop%first_pass) THEN

          ! We didn't find an existing diagnostic fitting the
          ! inputs, create a new diagnostic entry in SCMop
! DEPENDS ON: newdiag
          CALL newdiag                                                        &
            ( sname, lname, units, timprof, domprof, j                        &
            , NotWritten(streamlist), sname2(1:lsname*MIN(LEN(sname2_i),1))   &
            , d, SCMop )
        ELSE

          ! Should not be having new diagnostics beyond
          ! the first timestep
          PRINT*,'GetDistinctDiags ERROR: new diag after first timestep:',    &
                 sname(1:LEN(sname)),                                         &
                 SCMop%stepcount,SCMop%daycount,d
          GOTO 9999
        END IF
      ELSE
        IF (SCMop%first_pass) THEN

          ! d is the index of a diagnostic with the same sname
          ! as given in the input parameter list and the same
          ! dump_step as stream j. Send diagnostic d to stream
          ! j as well as it's existing streams then.
          IF (.NOT.StreamIsOn(SCMop%streams(d),j)) THEN

            SCMop%streams(d) = SCMop%streams(d) + Stream(j)

            ! Increment n_output to count the number we're
            ! going to output to this stream
            IF (.NOT.NotWritten(streamlist)) THEN
              SCMop%strm(j)%n_output = SCMop%strm(j)%n_output+1
            END IF

          ELSE

            ! Diagnostic d is already going to this stream,
            ! this should not ocurr.
            PRINT*,'GetDistinctDiags ERROR: Same diag. sent '//               &
                   'to same stream twice',d,sname(1:LEN(sname)),              &
                   j,SCMop%streams(d)
          END IF
        ELSE
          ! Check that diagnostic d is set up to go to
          ! this stream
          IF (.NOT.StreamIsOn(SCMop%streams(d),j)) THEN
            PRINT*,'GetDistinctDiags ERROR: the requested '//                 &
                   'streams for this diagnostic have changed',                &
                   d,sname(1:LEN(sname)),j,SCMop%streams(d)
          END IF
        END IF
      END IF

      ! We should now have a value for d, but it may still be
      ! unset if an error ocurred in newdiag
      IF (d /= unset) THEN

        ! Is this diagnostic distinct? i.e. is the value of d at
        ! this point different than for any previous value of j
        ! in this loop?
        distinct = .TRUE.

        DO i=1, ndistinct_diags
          IF (distinct_diags(i) == d) THEN
            distinct = .FALSE.
            EXIT
          END IF
        END DO

        IF (distinct) THEN
          ndistinct_diags = ndistinct_diags+1
          distinct_diags(ndistinct_diags) = d

          ! Make a requirement that diagnostic entries must be
          ! created in order (see use of diag_mem in SCMoutput)
          IF (ndistinct_diags >  1) THEN
            IF (d /= distinct_diags(ndistinct_diags-1)+1) THEN
              PRINT*,'GetDistinctDiags ERROR: non-consecutive diags',         &
                     distinct_diags(1:ndistinct_diags)
            END IF
          END IF
        END IF
      END IF

    END IF                  ! (StreamIsOn(streams,j))
  END DO                     ! j=1,maxnstreams

9999 CONTINUE

  RETURN

END SUBROUTINE getdistinctdiags
