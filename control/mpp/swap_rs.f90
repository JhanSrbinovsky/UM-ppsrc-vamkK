! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine to call swap_bounds_mv for R_u and R_v to avoid a supected
! compiler bug

SUBROUTINE swap_rs(r_u, r_v, row_length, u_rows, v_rows, off_x, off_y,  &
                   model_levels, fld_type_u, fld_type_v)

USE swapable_field_mod, ONLY:                                           &
    swapable_field_pointer_type

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
! Purpose:
!   This subroutine is to call swap_bounds_mv for R_u and R_v from
!   ni_conv_ctl. The sole reason for this to be in a subroutine is
!   to avoid a suspected compiler bug (Intel).
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Arguments
INTEGER, INTENT(IN) :: fld_type_u
INTEGER, INTENT(IN) :: fld_type_v
INTEGER, INTENT(IN) :: off_x
INTEGER, INTENT(IN) :: off_y
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: u_rows
INTEGER, INTENT(IN) :: v_rows

REAL, TARGET, INTENT(INOUT) :: r_u(1-off_x:row_length+off_x,            &
                                   1-off_y:u_rows+off_y,    model_levels)
REAL, TARGET, INTENT(INOUT) :: r_v(1-off_x:row_length+off_x,            &
                                   1-off_y:v_rows+off_y,    model_levels)

! Local variables
INTEGER :: i_field
TYPE (swapable_field_pointer_type) :: fields_to_swap(2) ! for swapbounds

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SWAP_RS',zhook_in,zhook_handle)

i_field = 0
i_field = i_field + 1
fields_to_swap(i_field) % field       => r_u(:,:,:)
fields_to_swap(i_field) % field_type  =  fld_type_u
fields_to_swap(i_field) % levels      =  model_levels
fields_to_swap(i_field) % rows        =  u_rows
fields_to_swap(i_field) % vector      =  .TRUE.

i_field = i_field + 1
fields_to_swap(i_field) % field       => r_v(:,:,:)
fields_to_swap(i_field) % field_type  =  fld_type_v
fields_to_swap(i_field) % levels      =  model_levels
fields_to_swap(i_field) % rows        =  v_rows
fields_to_swap(i_field) % vector      =  .TRUE.

! DEPENDS ON: swap_bounds_mv
CALL swap_bounds_mv(fields_to_swap, i_field, row_length,        &
                    off_x,          off_y)

IF (lhook) CALL dr_hook('SWAP_RS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE swap_rs
