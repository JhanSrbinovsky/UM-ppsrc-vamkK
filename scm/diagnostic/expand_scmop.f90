! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Increase the size of the SCM diagnostic arrays in SCMop

SUBROUTINE expand_scmop (SCMop)
  USE UM_types
! SCMop_type is defined in here...
  USE scmoptype_defn

  IMPLICIT NONE

! Description:
!   Increases the maximum number of diagnostics entries
!   allowed in SCMop by reallocating the relevant arrays.

! Method:
!   On the first timestep, as calls to SCMoutput are being made,
!   memory has to be allocated to arrays to hold the resulting
!   information. Since it is not known at the start of the run how
!   many calls to SCMoutput there will be and what their input
!   parameters are, no memory is allocated at the outset and the
!   variable SCMop%maxnentries (the maximum no. of diagnostic
!   "entries" that the arrays in SCMop can handle before they run
!   out of space) is zero. In this case most of the
!   statements in this routine are ignored (since they start with
!   "if (maxnentries >  0)"), and the arrays are simply allocated
!   with a size equal to "chunk". On subsequent calls the contents
!   of the allocatable arrays are copied into temporary arrays,
!   de-allocated, re-allocated with their original size plus
!   "chunk", and then the data is copied back in. i.e. the arrays
!   are re-allocated with large sizes without losing their
!   information.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran 90

  TYPE(SCMop_type) :: SCMop ! INOUT The derived-type structure
                            ! containing all the diagnostic
                            ! information

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

  ! Temporary arrays to hold data while SCMop arrays are being reallocated
  INTEGER, ALLOCATABLE :: Ixxx(:)
  LOGICAL, ALLOCATABLE :: Lxxx(:)
  CHARACTER (len=llname), ALLOCATABLE :: Cxxx(:)
  TYPE(allocatable_array), ALLOCATABLE :: Dxxx(:)

  ! The number by which to increment SCMop%maxnentries
  INTEGER, PARAMETER :: chunk=100

  ! Holds SCMop%maxnentries
  INTEGER :: maxnentries

  IF (SCMop%nentries /= SCMop%maxnentries) THEN
    PRINT*,'Expand_SCMop WARNING: SCMop is being expanded '//      &
           'before it is full, this could be dangerous'
  END IF

  maxnentries = SCMop%maxnentries

  ! Allocate the space required for the temporary arrays
  IF (maxnentries >  0) ALLOCATE(Cxxx(maxnentries))
  IF (maxnentries >  0) ALLOCATE(Ixxx(maxnentries))
  IF (maxnentries >  0) ALLOCATE(Lxxx(maxnentries))
  IF (maxnentries >  0) ALLOCATE(Dxxx(maxnentries))

  ! Increase the size of all the arrays in SCMop associated
  ! to specific diagnostics...

  IF (maxnentries >  0) Cxxx = SCMop%sname
  IF (maxnentries >  0) DEALLOCATE(SCMop%sname)
  ALLOCATE(SCMop%sname(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%sname(1:maxnentries) = Cxxx

  IF (maxnentries >  0) Cxxx = SCMop%lname
  IF (maxnentries >  0) DEALLOCATE(SCMop%lname)
  ALLOCATE(SCMop%lname(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%lname(1:maxnentries) = Cxxx

  IF (maxnentries >  0) Cxxx = SCMop%units
  IF (maxnentries >  0) DEALLOCATE(SCMop%units)
  ALLOCATE(SCMop%units(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%units(1:maxnentries) = Cxxx

  IF (maxnentries >  0) Ixxx = SCMop%domprof
  IF (maxnentries >  0) DEALLOCATE(SCMop%domprof)
  ALLOCATE(SCMop%domprof(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%domprof(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Ixxx = SCMop%timprof
  IF (maxnentries >  0) DEALLOCATE(SCMop%timprof)
  ALLOCATE(SCMop%timprof(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%timprof(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Ixxx = SCMop%streams
  IF (maxnentries >  0) DEALLOCATE(SCMop%streams)
  ALLOCATE(SCMop%streams(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%streams(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Ixxx = SCMop%dump_step
  IF (maxnentries >  0) DEALLOCATE(SCMop%dump_step)
  ALLOCATE(SCMop%dump_step(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%dump_step(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Ixxx = SCMop%nadd2dump
  IF (maxnentries >  0) DEALLOCATE(SCMop%nadd2dump)
  ALLOCATE(SCMop%nadd2dump(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%nadd2dump(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Lxxx = SCMop%only_radsteps
  IF (maxnentries >  0) DEALLOCATE(SCMop%only_radsteps)
  ALLOCATE(SCMop%only_radsteps(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%only_radsteps(1:maxnentries) = Lxxx

  IF (maxnentries >  0) Ixxx = SCMop%ncols
  IF (maxnentries >  0) DEALLOCATE(SCMop%ncols)
  ALLOCATE(SCMop%ncols(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%ncols(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Ixxx = SCMop%nrows
  IF (maxnentries >  0) DEALLOCATE(SCMop%nrows)
  ALLOCATE(SCMop%nrows(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%nrows(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Ixxx = SCMop%nlevs
  IF (maxnentries >  0) DEALLOCATE(SCMop%nlevs)
  ALLOCATE(SCMop%nlevs(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%nlevs(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Ixxx = SCMop%nelements
  IF (maxnentries >  0) DEALLOCATE(SCMop%nelements)
  ALLOCATE(SCMop%nelements(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%nelements(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Ixxx = SCMop%sname_id
  IF (maxnentries >  0) DEALLOCATE(SCMop%sname_id)
  ALLOCATE(SCMop%sname_id(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%sname_id(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Ixxx = SCMop%wd
  IF (maxnentries >  0) DEALLOCATE(SCMop%wd)
  ALLOCATE(SCMop%wd(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%wd(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Ixxx = SCMop%lastencounter
  IF (maxnentries >  0) DEALLOCATE(SCMop%lastencounter)
  ALLOCATE(SCMop%lastencounter(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%lastencounter(1:maxnentries) = Ixxx

  IF (maxnentries > 0) Ixxx = SCMop%substep
  IF (maxnentries > 0) DEALLOCATE(SCMop%substep)
  ALLOCATE(SCMop%substep(maxnentries+chunk))
  IF (maxnentries > 0) SCMop%substep(1:maxnentries) = Ixxx

  IF (maxnentries >  0) Dxxx = SCMop%diag
  IF (maxnentries >  0) DEALLOCATE(SCMop%diag)
  ALLOCATE(SCMop%diag(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%diag(1:maxnentries) = Dxxx

  IF (maxnentries >  0) Ixxx = SCMop%netcdf_id
  IF (maxnentries >  0) DEALLOCATE(SCMop%netcdf_id)
  ALLOCATE(SCMop%netcdf_id(maxnentries+chunk))
  IF (maxnentries >  0) SCMop%netcdf_id(1:maxnentries) = Ixxx

  IF (maxnentries >  0) DEALLOCATE(Cxxx)
  IF (maxnentries >  0) DEALLOCATE(Ixxx)
  IF (maxnentries >  0) DEALLOCATE(Lxxx)
  IF (maxnentries >  0) DEALLOCATE(Dxxx)

  SCMop%maxnentries = SCMop%maxnentries+chunk

  RETURN

END SUBROUTINE expand_scmop

