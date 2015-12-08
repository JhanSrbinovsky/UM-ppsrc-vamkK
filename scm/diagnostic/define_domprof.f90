! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Set up a domain profile for SCM diagnostics

SUBROUTINE define_domprof                                                     &
  ( tindex, tname, rowa1, rowa2, rowb1, rowb2, lev1, lev2, SCMop )
  USE UM_types
! SCMop_type is defined in here...
  USE scmoptype_defn

  IMPLICIT NONE

! Description:
! Set up a domain profile, i.e. a spatial region over which
! diagnostics may be defined. This will be stored in SCMop
! and in future can be referred to using tindex.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
! Language: Fortran 90

  INTEGER ::           &
    tindex              ! In An index for the profile

  CHARACTER (len=*) :: &
    tname               ! In A name for the profile

  INTEGER ::           &
    rowa1              &! In The horizontal area over which
  , rowa2              &! diagnostics of this type are defined
  , rowb1              &
  , rowb2

  INTEGER ::           &
    lev1               &! In The vertical range over which " " "
  , lev2

  TYPE(SCMop_type) :: SCMop ! InOut The derived-type structure
                            !       containing all the diagnostic
                            !       information

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

  ! Check that tindex is not greater than maxndomprof (i.e.
  ! avoid array out-of-bounds errors)
  IF (tindex >  maxndomprof) THEN
    PRINT*,'Define_Domprof ERROR: trying to define a domain'
    PRINT*,'profile with an index above the maximum. Increase'
    PRINT*,'maxndomprof or use a lower index. Domain profile'
    PRINT*,'not defined:',TRIM(tname)
    GOTO 9999
  END IF

  ! For safety, enforce a rule that no diagnostic may be defined
  ! unless *every* domain has been defined. This is to ensure
  ! SCMop%nelements(x) can never become incorrect by a domain being
  ! re-defined.
  IF (SCMop%nentries >  0) THEN
    PRINT*,'Define_Domprof ERROR: nentries is non-zero, setting ',   &
           'to zero now. Some diagnostics may be discarded',         &
    SCMop%nentries
    SCMop%nentries = 0
  END IF

  ! Fill the relevant section of the SCMop structure
  SCMop%d_name(tindex)  = tname
  SCMop%d_rowa1(tindex) = rowa1
  SCMop%d_rowa2(tindex) = rowa2
  SCMop%d_rowb1(tindex) = rowb1
  SCMop%d_rowb2(tindex) = rowb2
  SCMop%d_lev1(tindex)  = lev1
  SCMop%d_lev2(tindex)  = lev2

9999 CONTINUE

  RETURN

END SUBROUTINE define_domprof

