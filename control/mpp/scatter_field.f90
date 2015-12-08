! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Scatters a field from one processor to many processors

! Subroutine Interface:
SUBROUTINE scatter_field(local_field,    global_field,     &
                         local_row_len,  local_rows,       &
                         global_row_len, global_rows,      &
                         grid_type,      halo_type,        &
                         scatter_pe,     proc_group,       &
                         icode,          cmessage)

USE mpp_conf_mod, ONLY : gcom_coll_limit

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE UM_ParVars
  IMPLICIT NONE

!
! Description:
! Interface to potentially 2 methods of taking a field that is on a single
! processor and distributing it over a group of processors using the
! standard UM decomposition.

! Method:
! For c96_1c GCOM and MPL versions are available. The choice is made
! on a user definable threshold that is given in the UMUI.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine Arguments:
INTEGER :: local_row_len   ! length of rows in local part of field
INTEGER :: local_rows      ! number of rows in local part of field
INTEGER :: global_row_len  ! length of rows in global field
INTEGER :: global_rows     ! number of rows in global field
INTEGER :: grid_type       ! type (p,u or v) of grid
INTEGER :: halo_type       ! halo type (hence width) of grid
INTEGER :: scatter_pe      ! processor to scatter global field from
INTEGER :: proc_group      ! group id of processors involved here
INTEGER :: icode           ! out return code

REAL    :: local_field(local_row_len*local_rows)
                                         ! local part of field
REAL    :: global_field(global_row_len*global_rows)
                                         ! (on pe scatter_pe) global field

CHARACTER (LEN=80)  :: cmessage  ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Parameters and Common Blocks

IF (lhook) CALL dr_hook('SCATTER_FIELD',zhook_in,zhook_handle)


! Use GCOM version if number of processors less than or equal to
! threshold
IF (nproc <= gcom_coll_limit) THEN

! DEPENDS ON: scatter_field_gcom
  CALL scatter_field_gcom(local_field,    global_field,     &
                          local_row_len,  local_rows,       &
                          global_row_len, global_rows,      &
                          grid_type,      halo_type,        &
                          scatter_pe,     proc_group,       &
                          icode,          cmessage)


ELSE

! DEPENDS ON: scatter_field_mpl
  CALL scatter_field_mpl(local_field,    global_field,      &
                         local_row_len,  local_rows,        &
                         global_row_len, global_rows,       &
                         grid_type,      halo_type,         &
                         scatter_pe,     proc_group,        &
                         icode,          cmessage)


END IF

IF (lhook) CALL dr_hook('SCATTER_FIELD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE scatter_field
