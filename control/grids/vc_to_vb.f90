! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Calculate diagnostic quantities from the initial atmosphere dump
MODULE vc_to_vb_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE vc_to_vb &
(vc, rows, row_length, n_rows, levels, halo_x, halo_y,               &
                               global_row_length, vb                 &
, ub_arg                                                             &
)

USE dynamics_grid_mod, ONLY: l_vatpoles

USE atm_fields_bounds_mod, ONLY:    udims, udims_s, vdims, vdims_s

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE Field_Types
USE UM_ParVars
USE domain_params

USE eg_v_at_poles_mod
USE eg_parameters_mod, ONLY: pole_consts
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! Description:
!   vC_to_vB performs a simple horizontal interpolation of a field vC
!   (without halos) at the v position for Arakawa 'C' grid staggering
!   onto the uv position of the 'B' grid.
! Method:
!   1. Copy field (without halos) to an array with halos.
!   2. Populate halos using swapbounds.
!   3. Make a simple average of adjacent columns and copy to output
!      field (without halos).
!   4. In the vatpoles case, compute v at the poles using routine
!                                              eg_v_at_poles_Bgrid

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


! Subroutine arguments
!   Scalar arguments with intent(in):

INTEGER, INTENT(IN) :: rows            ! rows
INTEGER, INTENT(IN) :: row_length      ! horizontal dimensions
INTEGER, INTENT(IN) :: n_rows          ! rows for last (N) row of pes
INTEGER, INTENT(IN) :: levels          ! vertical levels (=1 in dyn_diag)
INTEGER, INTENT(IN) :: halo_x
INTEGER, INTENT(IN) :: halo_y            ! fixed halo sizes
INTEGER, INTENT(IN) :: global_row_length

!   Array  arguments with intent(in):
REAL, INTENT(IN) ::  vc(vdims%i_start:vdims%i_end,            &
                        vdims%j_start:vdims%j_end,          &
                        levels) ! v-Field on v points on C grid

REAL, INTENT(INOUT),OPTIONAL :: ub_arg(udims%i_start:udims%i_end, &  
                          vdims%j_start:vdims%j_end,          &  
                          levels) ! u-Field on uv points on B grid

! On B grid, u & v occupy same point
! In the i-direction they match the C grid u-points
! In the j-direction they match the C grid v-points

!   Array  arguments with intent(out):
REAL, INTENT(OUT) :: vb(udims%i_start:udims%i_end,            &
                        vdims%j_start:vdims%j_end,          &
                        levels) ! Field on uv points on B grid

! Local parameters:
CHARACTER(LEN=*), PARAMETER :: routinename='vC_to_vB'
CHARACTER(LEN=80)           :: cmessage


! Local scalars:
INTEGER :: i,j,k     !  loop indices
INTEGER :: istat     !  return code

! Local dynamic arrays:
REAL  ::  vc_halo(vdims_s%i_start:vdims_s%i_end,               &
                  vdims_s%j_start:vdims_s%j_end,levels)
REAL, ALLOCATABLE  :: ub(:,:) ! u-Field on uv points on B grid

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header


!   1. Copy field (without halos) to an array with halos

IF (lhook) CALL dr_hook('VC_TO_VB',zhook_in,zhook_handle)

IF (l_vatpoles .NEQV. PRESENT(ub_arg)) THEN
  istat = 1
  cmessage = "Incompatible options for V at poles."
  CALL ereport(RoutineName,istat,cmessage)
END IF


DO k=1,levels
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end
      vc_halo(i,j,k)=vc(i,j,k)
    END DO ! i
  END DO ! j
END DO ! k     
        
! The following line is required for LAM where the RHS halo
! at the extremity of the model is not filled in by SWAP_BOUNDS
IF (sb_model_domain /= mt_global) THEN
  DO k=1,levels
    DO j=vdims%j_start,vdims%j_end
      vc_halo(vdims%i_end+1,j,k)=vc(vdims%i_end,j,k)
      vc_halo(vdims%i_start-1,j,k)=vc(vdims%i_start,j,k)
    END DO ! j
  END DO ! k
END IF

!   2. Populate halos using swapbounds.

! DEPENDS ON: swap_bounds
CALL swap_bounds(vc_halo,row_length,n_rows,levels,                &
                 halo_x,halo_y,fld_type_v,.TRUE.)

!   3. Make a simple average of adjacent columns and copy to output
!      field (without halos). 
!      On polar rows set every point ot the absolute average for the row.

! here we use udims to define row_length of B grid
! and vdims to define number of rows.

DO k=1,levels
  DO j=vdims%j_start,vdims%j_end
    DO i=udims%i_start,udims%i_end
      vb(i,j,k)=(vc_halo(i,j,k)+vc_halo(i+1,j,k)) * 0.5
    END DO ! i
  END DO ! j
END DO ! k

! The presence of the optional argument ub signifies l_vatpoles
IF (PRESENT(ub_arg)) THEN
IF (sb_model_domain == mt_global) THEN
  IF (at_extremity(psouth) .OR. at_extremity(pnorth)) THEN
    ALLOCATE(ub(udims%i_start:udims%i_end, &  
                vdims%j_start:vdims%j_end))
  ELSE
    ALLOCATE(ub(1,1))
  END IF

  ! S polar row
  IF (at_extremity(psouth)) THEN
    DO k=1,levels
      ub = ub_arg(:,:,k)
      CALL eg_uv_at_poles_Bgrid &
                  (ub, row_length, rows, n_rows,                        &
                   halo_x, halo_y, pole_consts, -1.0,                   &
                   vdims%j_start+1, vdims%j_start, vb(:,:,k) )
    END DO
  END IF

  ! N polar row
  IF (at_extremity(pnorth)) THEN
    DO k=1,levels
      ub = ub_arg(:,:,k)
      CALL eg_uv_at_poles_Bgrid &
                  (ub, row_length, rows, n_rows,                        &
                   halo_x, halo_y, pole_consts, 1.0,                    &
                   vdims%j_end-1, vdims%j_end, vb(:,:,k) )
    END DO
  END IF
  IF (ALLOCATED(ub)) THEN
    DEALLOCATE(ub)
  END IF
END IF
END IF


IF (lhook) CALL dr_hook('VC_TO_VB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE vc_to_vb
END MODULE vc_to_vb_mod
