! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Calculate diagnostic quantities from the initial atmosphere dump
MODULE uc_to_ub_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE uc_to_ub                                             &
 (uc, row_length, rows, n_rows, levels, halo_x,halo_y, ub)

USE dynamics_grid_mod, ONLY: l_vatpoles

USE atm_fields_bounds_mod, ONLY:   udims, udims_s, vdims, vdims_s

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE Field_Types
USE UM_ParVars
USE domain_params

USE eg_v_at_poles_mod
USE eg_parameters_mod, ONLY: pole_consts

IMPLICIT NONE

! Description:
!   uC_to_uB performs a simple horizontal interpolation of a field uC
!   (without halos) at the u position for Arakawa 'C' grid staggering
!   onto the uv position of the 'B' grid.
! Method:
!   1. Copy field (without halos) to an array with halos.
!   2. Populate halos using swapbounds.
!   3. Make a simple average of adjacent rows and copy to output field
!     (without halos).

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

! Declarations

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows            ! horizontal dimensions
INTEGER, INTENT(IN) :: n_rows          ! rows for last (N) row of pes
INTEGER, INTENT(IN) :: levels          ! vertical levels
INTEGER, INTENT(IN) :: halo_x
INTEGER, INTENT(IN) :: halo_y          ! fixed halo sizes

!   Array  arguments with intent(in):
REAL, INTENT(IN) :: uc(udims%i_start:udims%i_end,           & ! 0:row_len -1
                       udims%j_start:udims%j_end,           & ! 1:rows
                       levels) ! Field on u points on C grid

! On B grid, u & v occupy same point
! In the i-direction they match the C grid u-points
! In the j-direction they match the C grid v-points

!   Array  arguments with intent(inout):
! This is needed since call to eg_uv_at_poles_Bgrid uses ub
REAL, INTENT(INOUT) :: ub(udims%i_start:udims%i_end,        &  ! 0:row_len -1
                        vdims%j_start:vdims%j_end,          &  ! 0:n_rows -1
                        levels) ! Field on uv points on B grid


! Local parameters:
CHARACTER(LEN=*), PARAMETER ::  routinename='uC_to_uB'


! Local scalars:
INTEGER  ::  i,j,k     !  loop indices

! Local dynamic arrays:
REAL  :: uc_halo(udims_s%i_start:udims_s%i_end,               &
                 udims_s%j_start:udims_s%j_end,levels)
                 
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! End of header

!   1. Copy field (without halos) to an array with halos

IF (lhook) CALL dr_hook('UC_TO_UB',zhook_in,zhook_handle)
DO k=1,levels
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end
      uc_halo(i,j,k)=uc(i,j,k)
    END DO ! i
  END DO ! j
END DO ! k 


! The following lines are required for LAM where the top and bottom halo
! at the extremity of the model is not filled in by SWAP_BOUNDS
IF (sb_model_domain /= mt_global) THEN
  DO k=1,levels
    DO i=udims%i_start,udims%i_end
      uc_halo(i,udims%j_end+1,k)=uc(i,udims%j_end,k)
      uc_halo(i,udims%j_start-1,k)=uc(i,udims%j_start,k)
    END DO
  END DO ! k
END IF

!   2. Populate halos using swapbounds.

! DEPENDS ON: swap_bounds
CALL swap_bounds(uc_halo,row_length,rows,levels,                  &
                 halo_x,halo_y,fld_type_u,.TRUE.)

!   3. Make a simple average of adjacent rows and copy to output field
!     (without halos). Set to zero on polar rows.

! Here we use udims to define row_length of B grid
!  and vdims to define number of rows

DO k=1,levels
  DO j=vdims%j_start,vdims%j_end
    DO i=udims%i_start,udims%i_end
      ub(i,j,k)=(uc_halo(i,j,k)+uc_halo(i,j+1,k)) * 0.5
    END DO ! i
  END DO ! j
END DO ! k

IF (l_vatpoles) THEN
IF (sb_model_domain == mt_global) THEN

  DO k=1,levels

    ! S polar row
    IF (at_extremity(psouth)) THEN 
      CALL eg_uv_at_poles_Bgrid &
                  (ub(:,:,k), row_length, rows, n_rows,                 &
                   halo_x, halo_y, pole_consts, -1.0,                   &
                   vdims%j_start+1, vdims%j_start)
    END IF

    ! N polar row
    IF (at_extremity(pnorth)) THEN
      CALL eg_uv_at_poles_Bgrid &
                  (ub(:,:,k), row_length, rows, n_rows,                 &
                   halo_x, halo_y, pole_consts, 1.0,                    &
                   vdims%j_end-1, vdims%j_end)
    END IF

  END DO

END IF
END IF ! vatpoles

IF (lhook) CALL dr_hook('UC_TO_UB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE uc_to_ub
END MODULE uc_to_ub_mod
