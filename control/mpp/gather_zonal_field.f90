! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Gathers a zonal field from many processors onto one

! Subroutine interface:

SUBROUTINE gather_zonal_field ( local_field, global_field,             &
                                local_size,  global_size,              &
                                levels,      grid_code,                &
                                grid_type,   halo_type,     gather_pe)
USE mpl, ONLY : &
    mpl_real

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE UM_ParVars
IMPLICIT NONE

! Description:
! Takes a decomposed zonal field on many processors, and gathers
! it to just one processor.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine arguments:
INTEGER, INTENT(IN) :: local_size    ! size of level of LOCAL_FIELD
INTEGER, INTENT(IN) :: global_size   ! size of level of GLOBAL FIELD
INTEGER, INTENT(IN) :: levels        ! number of levels in field
INTEGER, INTENT(IN) :: grid_code     ! ppx grid code of field
INTEGER, INTENT(IN) :: grid_type     ! P,U or V
INTEGER, INTENT(IN) :: halo_type     ! halo type of field
INTEGER, INTENT(IN) :: gather_pe     ! PE to gather GLOBAL_FIELD

REAL, INTENT(IN) ::   local_field(local_size,levels)   ! field to gather
REAL, INTENT(OUT) ::  global_field(global_size,levels) ! gathered field


! Local variables
INTEGER  :: info      ! MPL return code
INTEGER  :: idx       ! index
INTEGER  :: iproc     ! loop counter

INTEGER  :: counts(0:nproc-1)           ! Counts for gather
INTEGER  :: displacements(0:nproc-1)    ! Counts for gather

! Parameters
INTEGER, PARAMETER :: unset = -1234     ! for initialising saved vars

! Saved variables
INTEGER, SAVE :: old_local_size  = unset
INTEGER, SAVE :: old_global_size = unset
INTEGER, SAVE :: old_levels      = unset
INTEGER, SAVE :: old_grid_type   = unset
INTEGER, SAVE :: old_halo_type   = unset
INTEGER, SAVE :: old_gather_pe   = unset
INTEGER, SAVE :: send_type       = unset

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------

IF (lhook) CALL dr_hook('GATHER_ZONAL_FIELD',zhook_in,zhook_handle)

! Only processors at the West need to take part in this communication
! i.e. a column of processors at the west using a column communicator
IF (at_extremity(pwest)) THEN
  ! 1.0 Setup MPI datatypes - can we reuse?
  IF ((local_size   /= old_local_size  )  .OR.                  &
      (global_size  /= old_global_size )  .OR.                  &
      (levels       /= old_levels      )  .OR.                  &
      (grid_type    /= old_grid_type   )  .OR.                  &
      (halo_type    /= old_halo_type   )  .OR.                  &
      (gather_pe    /= old_gather_pe   )) THEN

    ! Need to create new send type
    ! Must get rid of old one if this isn't the first call
    IF (send_type /= unset) THEN
      CALL mpl_type_free(send_type, info)
    END IF

    ! New type is effectively a 1D vector with a halo. This halo is
    ! actually unlikely to be there, but we need to account for it
    ! just in case.
    CALL mpl_type_vector(local_size, 1, 1 + 2*halosize(2,halo_type), &
                         mpl_real, send_type, info)
    CALL mpl_type_commit(send_type, info)
  END IF ! New MPI type

  ! 2.0 Setup counts and displacements for recv end
  idx = 0
  DO iproc = 0, nproc-1
    IF (g_gridpos(1,iproc) == 0 ) THEN  ! at west of LPG
      counts(idx) = g_blsize(2, grid_type, iproc)

      IF (idx == 0) THEN
        displacements(idx) = 0
      ELSE
        displacements(idx) = SUM(counts(0:idx-1))
      END IF

      idx = idx + 1
    END IF
  END DO

  ! 3.0 Do the exhange of data
  CALL mpl_gatherv(local_field(1+halosize(2,halo_type),1), levels,     &
                   send_type,                                          &
                   global_field, counts, displacements, mpl_real,      &
                   gather_pe, gc_proc_col_group, info)
END IF

IF (lhook) CALL dr_hook('GATHER_ZONAL_FIELD',zhook_out,zhook_handle)

RETURN
END SUBROUTINE gather_zonal_field
