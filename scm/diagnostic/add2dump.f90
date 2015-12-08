! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Update the dump array of a SCM diagnostic entry in SCMop

SUBROUTINE add2dump (x, nelements, d, SCMop, startperiod, endperiod)

! SCMop_type is defined in here...
  USE scmoptype_defn

  IMPLICIT NONE

! Description:
!   Take input array x and add to the dump array in the manner
!   specified by the time profile. This performs this timestep's
!   role in constructing the diagnostic.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran90

  TYPE(SCMop_type) :: SCMop ! InOut The derived-type structure
                            !       containing all the diagnostic
                            !       information

  INTEGER :: nelements      ! The no. of elements in array x

  REAL :: x(nelements)      ! In The array from which the diagnostic
                            !    is constructed

  INTEGER :: d              ! In The diagnostic index in SCMop

  LOGICAL ::               &
    startperiod            &! In These flag the start and end
  , endperiod               !    of the dumping period, when special
                            !    actions may need to be taken                 


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

  REAL, POINTER :: dump(:)  ! The dump array in which the
                            ! diagnostic is constructed

  INTEGER :: timprof        ! The time profile of the diagnostic

  INTEGER :: nadd2dump      ! The number by which to divide
                            ! accumulated diag's in order to get the mean
  INTEGER :: wnelements     ! The no. of elements in array wdump

  REAL, POINTER :: wdump(:) ! A dump array of another
                            ! diagnostic which may be used in the
                            ! construction of this diagnostic

  INTEGER :: wtimprof       ! The timeprofile of the other diag.
  INTEGER :: wnadd2dump     ! As nadd2dump, but for the other diagnostic

  INTEGER :: i              ! Counter

  ! We don't want lots of '%'s around, make some abbreviations...
  dump      => SCMop%diag(d)%dump
  nelements =  SCMop%nelements(d)
  timprof   =  SCMop%timprof(d)
  nadd2dump =  SCMop%nadd2dump(d)

  IF (SCMop%wd(d) >  0) THEN
    wdump      => SCMop%diag(SCMop%wd(d))%dump
    wnelements =  SCMop%nelements(SCMop%wd(d))
    wtimprof   =  SCMop%timprof(SCMop%wd(d))
    wnadd2dump =  SCMop%nadd2dump(SCMop%wd(d))
  ELSE
    NULLIFY(wdump)
    wnelements = 0
    wtimprof   = 0
    wnadd2dump = 0
  END IF

!-T_CONST--A constant value---------------------------------------------
  IF (timprof == t_const) THEN
    DO i=1, nelements
      dump(i) = x(i)
    END DO

!-T_INST--Instaneous value----------------------------------------------
  ELSE IF (timprof == t_inst) THEN
    IF (startperiod) THEN

      ! Set the dump to a default value. Shouldn't really be
      ! necessary, but handy for checking everything's working OK.
      DO i=1, nelements
        dump(i) = -999.0
      END DO
    END IF

    IF (endperiod) THEN

      ! The dump will ocurr at this timestep, so we want the value now.
      DO i=1, nelements
        dump(i) = x(i)
      END DO
    END IF

!-T_ACC--Accumulated value----------------------------------------------
  ELSE IF (timprof == t_acc) THEN
    IF (startperiod) THEN

      ! Take this timestep's values
      DO i=1, nelements
        dump(i) = x(i)
      END DO
    ELSE

      ! Add to the current values
      DO i=1, nelements
        dump(i) = dump(i) + x(i)
      END DO
    END IF

!-T_AVG--Value averaged over a period----------------------------------
  ELSE IF (timprof == t_avg) THEN
    IF (startperiod) THEN

      ! Take this timestep's values
      DO i=1, nelements
        dump(i) = x(i)
      END DO
    ELSE

      ! Add to the current values
      DO i=1, nelements
        dump(i) = dump(i) + x(i)
      END DO
    END IF

    ! And form average if it's the end of the period
    IF (endperiod) THEN
      DO i=1, nelements
        dump(i) = dump(i)/nadd2dump
      END DO
    END IF

!-T_MIN--The minimum value over a time period---------------------------
  ELSE IF (timprof == t_min) THEN
    IF (startperiod) THEN

      ! Take this timestep's values
      DO i=1, nelements
        dump(i) = x(i)
      END DO
    ELSE

      ! Overwrite the current values if smaller
      DO i=1, nelements
        IF (x(i) <  dump(i)) THEN
          dump(i) = x(i)
        END IF
      END DO
    END IF

!-T_MAX--The maximum value over a time period---------------------------
  ELSE IF (timprof == t_max) THEN
    IF (startperiod) THEN

      ! Take this timestep's values
      DO i=1, nelements
        dump(i) = x(i)
      END DO
    ELSE

      ! Overwrite the current value if bigger
      DO i=1, nelements
        IF (x(i) >  dump(i)) THEN
          dump(i) = x(i)
        END IF
      END DO
    END IF

!-T_DIV-T_MULT--The accumulated value of one diagnostic multiplied or---
!---------------divided by the value of another-------------------------
  ELSE IF (timprof == t_div    .OR.timprof == t_mult .OR.               &
           timprof == t_acc_div.OR.timprof == t_acc_mult) THEN

    IF (startperiod) THEN
      DO i=1, nelements
        dump(i) = x(i)
      END DO
    ELSE

      ! Add to the current value
      DO i=1, nelements
        dump(i) = dump(i) + x(i)
      END DO
    END IF

    IF (endperiod) THEN
      IF (timprof == t_mult.OR.timprof == t_acc_mult) THEN

        ! Multiply the two diagnostics
        DO i=1, nelements
          dump(i) = dump(i)*wdump(MOD(i-1,wnelements)+1)
        END DO
      ELSE IF (timprof == t_div.OR.timprof == t_acc_div) THEN

        ! Divide the two diagnostics
        DO i=1, nelements
          IF (wdump(MOD(i-1,wnelements)+1) /= 0) THEN
            dump(i) = dump(i)/                                          &
                      wdump(MOD(i-1,wnelements)+1)
          ELSE IF (dump(i) /= 0) THEN
             WRITE(*,'(A,1X,A,1X,I3,1X,A,1X,I6)')                       &
                  'Add2Dump WARNING:divide by zero avoided:'
          END IF
        END DO
      END IF
 
      ! Divide by dumping period?
      IF (timprof == t_div.OR.timprof == t_mult) THEN
        DO i=1, nelements
          dump(i) = dump(i)/nadd2dump
        END DO
      END IF
    END IF

!-----------------------------------------------------------------------
  ELSE
    PRINT*,'Add2Dump ERROR: I do not know about this ',                 &
           'temporal type:',timprof
  END IF
!-----------------------------------------------------------------------

  RETURN

END SUBROUTINE add2dump

