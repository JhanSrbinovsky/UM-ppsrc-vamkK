! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Scatters zonal field from one processor to many processors

! Subroutine interface:

SUBROUTINE scatter_zonal_field (                                  &
  local_field , global_field ,                                    &
  local_size  , global_size  ,                                    &
  levels, grid_code, grid_type ,halo_type,scatter_pe)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE mpl, ONLY : mpl_real
USE ereport_mod, ONLY: ereport
USE UM_ParVars, ONLY: mype, nproc, g_blsize, blsize,              &
                      halo_type_no_halo, nproc_x, nproc_y

IMPLICIT NONE

! Description:
! Takes a zonal field on a single processor, and decomposes it over
! many processors.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine arguments:

INTEGER, INTENT(IN) :: local_size     ! IN: size of level of LOCAL_FIELD
INTEGER, INTENT(IN) :: global_size    ! IN: size of level of GLOBAL FIELD
INTEGER, INTENT(IN) :: levels         ! IN: number of levels in field
INTEGER, INTENT(IN) :: grid_code      ! IN: ppx grid code of field
INTEGER, INTENT(IN) :: grid_type      ! IN: P,U or V
INTEGER, INTENT(IN) :: halo_type      ! IN: halo type of field
INTEGER, INTENT(IN) :: scatter_pe     ! IN:  PE on which GLOBAL_FIELD resides

REAL, INTENT(IN) :: global_field(global_size,levels) ! IN : field to scatter
REAL, INTENT(OUT) :: local_field(local_size,levels)  ! OUT : local part of field

! Local variables

INTEGER :: i, j, k                    ! loopers
INTEGER :: icode                      ! error code
INTEGER :: ierr                       ! mpl return code
INTEGER :: my_comm                    ! communicator
INTEGER :: disp                       ! cumulative displacement

INTEGER :: counts(0:nproc-1)          ! numbers of elements to send
INTEGER :: displacements(0:nproc-1)   ! offset from start of array

CHARACTER (LEN=80) :: cmessage        ! error message
CHARACTER (LEN=*), PARAMETER :: routinename='scatter_zonal_field'


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!====================================================================

IF (lhook) CALL dr_hook('SCATTER_ZONAL_FIELD',zhook_in,zhook_handle)

! All zonal fields to scatter should have no halos. Throw an error
! if this isn't the case
IF (halo_type /= halo_type_no_halo) THEN
  icode = 1
  WRITE(cmessage, '(a)') 'Cannot scatter zonal fields with a halo'
  CALL ereport(routinename, icode, cmessage)
END IF

! Set up counts and displacements for data being scattered out
IF (mype == scatter_pe) THEN
  DO i = 0, nproc-1
    counts(i) = g_blsize(2, grid_type, i)
  END DO

  disp = 0
  DO j = 0, nproc_y-1
    IF (j > 0) THEN
      disp = disp + counts((j-1)*nproc_x)
    END IF

    DO i = 0, nproc_x-1
      displacements(j*nproc_x + i ) = disp
    END DO
  END DO
END IF

CALL gc_get_communicator(my_comm, ierr)

! Do scatters, 1 per level. Not the most efficient, but in practise
! only a single level gets to this routine at a time anyway and doing
! it this way simplifies things a lot.
DO k = 1, levels
  CALL mpl_scatterv(global_field(1,k), counts, displacements,     &
                    mpl_real, local_field(1,k),                   &
                    blsize(2, grid_type), mpl_real, scatter_pe,   &
                    my_comm, ierr)
END DO

IF (lhook) CALL dr_hook('SCATTER_ZONAL_FIELD',zhook_out,zhook_handle)
RETURN

END SUBROUTINE scatter_zonal_field
