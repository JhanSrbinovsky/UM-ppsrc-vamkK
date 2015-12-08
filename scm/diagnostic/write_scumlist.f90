! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Write out info about the output diagnostics for one or all streams

SUBROUTINE write_scumlist (SCMop, istrm)

! SCMop_type is defined in here...
  USE scmoptype_defn

  IMPLICIT NONE

! Description:
!   For a given stream, writes information about each domain profile and
!   diagnostic that is active in that stream to the stream's output file
!   as part of its header. PV-wave can then read this info to make sense
!   of the rest of the data file. If the requested stream number is zero,
!   a new file will be opened and all the same information will be written
!   to it, but now pertaining to all streams. Such a file will list all
!   domain profiles and diagnostics, and will have additional comments.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran 90

  TYPE(SCMop_type) :: &
    SCMop              ! In The derived-type structure containing
                       !    all the diagnostic information

  INTEGER ::                    &
    istrm                       &! In The stream in question. Zero=all streams
  , ndiags                      &! The number of diagnostic entries found in
                                 ! SCMop which we will write about.
  , diags(SCMop%nentries)       &! Contains in 1:niags the indices of the
                                 ! entries in SCMop
  , domprof_sparse(maxndomprof)  ! A sparse array indicating the domain
                                 ! profiles possessed by the ndiags entries.

  INTEGER :: i,j    ! Counters

  ! Character function that incorporates the substep number into
  ! the short name of a diagnostic. Somewhat longer than a normal
  ! short name.
  CHARACTER (len=lsname+10) :: add_substep_to_sname,short_name

  ! The longest possible sname after the substep has been appended
  INTEGER :: sname_length

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

  INTEGER :: unit ! The unit to which we'll write
  CHARACTER (len=100) :: fmt
  CHARACTER (len=lsname) :: sname
  INTEGER :: substep

  ! Make a list the diagnostics we're going to write and
  ! which domain profiles they use. ndiags will hold the total
  ! number of diagnostic entries in SCMop fitting the bill, and
  ! domprof_sparse will be a sparse array indicating which domain
  ! profiles are used by these diagnostics.
  ndiags = 0
  domprof_sparse = 0

  ! If a specific stream has been specified...
  IF (istrm /= 0) THEN

    ! ... select only those diagnostics going to that stream
    DO i=1, SCMop%nentries
      IF (StreamIsOn(SCMop%streams(i),istrm).AND..NOT.              &
            NotWritten(SCMop%streams(i))) THEN
        ndiags = ndiags+1
        diags(ndiags) = i
        domprof_sparse(SCMop%domprof(i)) = 1
      END IF
    END DO

    ! Make a check just for the sake of it
    IF (ndiags /= SCMop%strm(istrm)%n_output) THEN
      PRINT*,'dump_streams_init ERROR: an inconsistency '//         &
             'has ocurred !',ndiags,SCMop%strm(istrm)%n_output,     &
             istrm,SCMop%nentries

      ! Switch the diagnostic system off so dodgy data is
      ! not mistakenly used in good faith.
      PRINT*,' -> switching diagnostic system off!'
      SCMop%on = .FALSE.
    END IF
  ELSE
    ! If strm=0 then use all diagnostics which are going to
    ! at least one stream. But since the same diagnostic can
    ! go to several streams, and thus have several entries in
    ! SCMop, then make a condition that if the i'th entry has
    ! the same short name as the (i-1)'th, then skip the i'th
    ! to avoid multiple identical lines. We -do- want multiple
    ! lines for different sub-steps though, so allow identical
    ! consecutive snames if their sub-step numbers are different.

    sname = ''
    substep = -999

    DO i=1, SCMop%nentries
      IF (AnyStreamOn(SCMop%streams(i)).AND.                        &
           (SCMop%sname(i) /= sname .OR.                            &
            SCMop%substep(i) /= substep)) THEN

        ndiags = ndiags+1
        diags(ndiags) = i
        domprof_sparse(SCMop%domprof(i)) = 1
        sname = SCMop%sname(i)
        substep = SCMop%substep(i)

      END IF
    END DO
  END IF

  ! Determine the unit number to write to
  IF (istrm /= 0) THEN

    ! Write to the unit of the specified stream. The file
    ! should already be open.
    unit = SCMop%strm(istrm)%op_unit
  ELSE

    ! Write to a new file we'll attach to unit 10
    unit = 10
    OPEN (unit=unit,file='scumlist.dat')
  END IF

!-------------------------------------------------------------
!     Write about the domain profiles
!-------------------------------------------------------------

  ! Write the number of domain profiles used by
  ! diagnostics in this stream
  IF (istrm == 0) WRITE(unit,*)'No. of domain profiles:'
                  WRITE(unit,'(I3)') SUM(domprof_sparse)

  ! Write the format with which we will write some of
  ! the upcoming lines
  fmt = '(I3,1X,A15,I3,1X,I3,1X,I3,1X,I3,1X,I3,1X,I3)'
  IF (istrm == 0) WRITE(unit,'(A)')'Line format:'
                  WRITE(unit,'(A)') TRIM(fmt)

  ! Write info about each domain profile
  IF (istrm == 0) WRITE(unit,*)'List of domain profiles:'

  DO i=1, maxndomprof
    IF (domprof_sparse(i) == 1) THEN
      WRITE(unit,fmt)                                          &
        i,SCMop%d_name(i),                                     &
        SCMop%d_rowa1(i),SCMop%d_rowa2(i),                     &
        SCMop%d_rowb1(i),SCMop%d_rowb2(i),                     &
        SCMop%d_lev1(i),SCMop%d_lev2(i)
    END IF
  END DO

!-------------------------------------------------------------
!     Write about the diagnostics themselves
!-------------------------------------------------------------

  ! Write the number
  IF (istrm == 0) WRITE(unit,*)'No. of diagnostics:'
                  WRITE(unit,'(I3)') ndiags

  ! If there are multiple sub-steps then the short-name will
  ! have "_#N" tagged onto it, where N is the sub-step number.
  ! In this case we need to increase the amount of space made
  ! available for the short name.

  IF (SCMop%num_substeps <= 1) THEN
    sname_length = lsname
  ELSE
    sname_length = lsname                                      &
                 + 2+(INT(LOG10(float(SCMop%num_substeps)))+1)
  END IF

  ! Compose the format of the line that will describe
  ! each diagnostic
  WRITE(fmt,'(A,I2, A,I2, A,I2, A)')'(I3' //                   &
       ',1X,A',sname_length,                                   &
       ',1X,A',llname,                                         &
       ',1X,A',lunits,                                         &
       ',I2)'
  ! Write the format to the file
  IF (istrm == 0) WRITE(unit,'(A)')'Line format:'
                  WRITE(unit,'(A)') TRIM(fmt)

  ! Write a 1-line description of each diagnostic
  ! consisting of a unique integer ID, its short name,
  ! its long name, its units and its domain profile.
  IF (istrm == 0) WRITE(unit,*)'List of diagnostics '//        &
                  '(i,sname,lname,units,domprof):'
  DO i=1, ndiags
    j = diags(i)
! DEPENDS ON: add_substep_to_sname
    short_name = add_substep_to_sname(SCMop,j)
    WRITE(unit,fmt)                                            &
         SCMop%sname_id(j),short_name,                         &
         SCMop%lname(j),SCMop%units(j),SCMop%domprof(j)
  END DO

  IF (istrm == 0) THEN
    ! We need to close the new file we opened above
    CLOSE(unit)
  END IF

  RETURN
END SUBROUTINE write_scumlist

