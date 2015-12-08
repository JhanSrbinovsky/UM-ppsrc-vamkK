! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Gathers a field from many processors to one processor

! Subroutine Interface:
SUBROUTINE gather_field(local_field,    global_field,     &
                        local_row_len,  local_rows,       &
                        global_row_len, global_rows,      &
                        grid_type,      halo_type,        &
                        gather_pe,      proc_group,       &
                        icode,          cmessage)


USE mpp_conf_mod, ONLY: gcom_coll_limit
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
IMPLICIT NONE

!
! Description:
! Interface to potentially 2 methods of taking a field that is decomposed
! over a group of processors and gathering the data so a single processor
! contains the entire global field.

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
INTEGER :: gather_pe       ! processor to gather global field to
INTEGER :: proc_group      ! group id of processors involved here
INTEGER :: icode           ! out return code

REAL    :: local_field(local_row_len*local_rows)
                                         ! local part of field
REAL    :: global_field(global_row_len*global_rows)
                                         ! (on pe gather_pe) global field

CHARACTER (LEN=80)  :: cmessage  ! error message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



! Parameters and Common Blocks

IF (lhook) CALL dr_hook('GATHER_FIELD',zhook_in,zhook_handle)


! Use GCOM version if number of processors less than or equal to
! threshold
IF (nproc <= gcom_coll_limit) THEN

! DEPENDS ON: gather_field_gcom
  CALL gather_field_gcom(local_field,    global_field,     &
                         local_row_len,  local_rows,       &
                         global_row_len, global_rows,      &
                         grid_type,      halo_type,        &
                         gather_pe,      proc_group,       &
                         icode,          cmessage)


ELSE

! DEPENDS ON: gather_field_mpl
  CALL gather_field_mpl(local_field,    global_field,      &
                        local_row_len,  local_rows,        &
                        global_row_len, global_rows,       &
                        grid_type,      halo_type,         &
                        gather_pe,      proc_group,        &
                        icode,          cmessage)


END IF

IF (lhook) CALL dr_hook('GATHER_FIELD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE gather_field
