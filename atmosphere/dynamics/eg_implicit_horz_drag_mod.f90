! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_implicit_horz_drag_mod

! Subroutine: eg_init_horz_drag
!
! Description: An "implicit" implementations of a simple boundary layer
!              drag scheme. It is implicit in a sense as fast physics 
!              is. It computes the same source term as the explicit
!              version but uses the latest available estimate of the
!              velocity field.
!
! Method: (see eg_explicit_horz_drag_mod)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: DYNAMICS ADVECTION
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

CONTAINS

SUBROUTINE eg_implicit_horz_drag(r_u,r_v,u,v)


USE atm_fields_bounds_mod, ONLY : udims,vdims,udims_s,vdims_s
USE timestep_mod,          ONLY : timestep
! use the f1sp predictors (coming in through argument list now) rather than the latest np1 value!
!USE atm_step_local,        ONLY : u_np1,v_np1
USE eg_horz_drag_mod,      ONLY : cd_u,cd_v,r_u_store,r_v_store,  &
                                  l_expl_horz_drag

USE proc_info_mod,         ONLY : me,n_proc

USE parkind1,              ONLY : jprb, jpim
USE yomhook,               ONLY : lhook, dr_hook

USE PrintStatus_mod

IMPLICIT NONE

REAL, INTENT (IN) ::                                                  &
    u(udims_s%i_start:udims_s%i_end,                                  &
      udims_s%j_start:udims_s%j_end,                                  &
      udims_s%k_start:udims_s%k_end)                                  &
  , v(vdims_s%i_start:vdims_s%i_end,                                  &
      vdims_s%j_start:vdims_s%j_end,                                  &
      vdims_s%k_start:vdims_s%k_end) 

REAL, INTENT (INOUT) ::                                               &
  r_u(udims_s%i_start:udims_s%i_end,                                  &
      udims_s%j_start:udims_s%j_end,                                  &
      udims_s%k_start:udims_s%k_end)                                  &
, r_v(vdims_s%i_start:vdims_s%i_end,                                  &
      vdims_s%j_start:vdims_s%j_end,                                  &
      vdims_s%k_start:vdims_s%k_end)    


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER i
INTEGER j
INTEGER k
INTEGER ierr

REAL, ALLOCATABLE :: diagu(:,:,:)   ! only used for soure term diagnostic print
REAL, ALLOCATABLE :: diagv(:,:,:)   ! at the end of this routine if 
                                    ! PrintStatus == PrStatus_Diag (stores source
                                    ! terms at entry in roder to subtract at the end)

REAL max_r_u   ! diagnostic, min/max of r_u/r_v
REAL max_r_v
REAL min_r_u
REAL min_r_v


IF (lhook) CALL dr_hook('EG_IMPLICIT_HORZ_DRAG',zhook_in,zhook_handle)

IF ( PrintStatus == PrStatus_Diag) THEN

  ALLOCATE(diagu(udims_s%i_start:udims_s%i_end,                       &
                 udims_s%j_start:udims_s%j_end,                       &
                 udims_s%k_start:udims_s%k_end))
  ALLOCATE(diagv(vdims_s%i_start:vdims_s%i_end,                       &
                 vdims_s%j_start:vdims_s%j_end,                       &
                 vdims_s%k_start:vdims_s%k_end))

  diagu(:,:,:) = r_u(:,:,:)
  diagv(:,:,:) = r_v(:,:,:)

END IF

DO k=udims%k_start,udims%k_end
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end

      r_u(i,j,k) =  r_u(i,j,k) -                                      &
                    cd_u(i,j,k)* SQRT(                                &
                                 (.25*(v(i,j  ,k)+v(i+1,j  ,k)        &
                                      +v(i,j-1,k)+v(i+1,j-1,k)))**2   &
                                 +u(i,j,k)**2                         &
                                 ) * u(i,j,k)

    END DO
  END DO
END DO

DO k=vdims%k_start,vdims%k_end
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end

      r_v(i,j,k) =  r_v(i,j,k) -                                      &
                    cd_v(i,j,k)* SQRT(                                &
                                 (.25*(u(i,j  ,k)+u(i-1,j  ,k)        &
                                      +u(i,j+1,k)+u(i-1,j+1,k)))**2   &
                                 +v(i,j,k)**2                         &
                                 ) * v(i,j,k)

    END DO
  END DO
END DO

IF(l_expl_horz_drag) THEN
! the explict source contribution has to be subtracted, to avoid double counting
  r_u(:,:,:) = r_u(:,:,:) - r_u_store(:,:,:)
  r_v(:,:,:) = r_v(:,:,:) - r_v_store(:,:,:)
END IF


IF ( PrintStatus == PrStatus_Diag) THEN

     max_r_u = MAXVAL(R_u-diagu)
     min_r_u = MINVAL(R_u-diagu)
     max_r_v = MAXVAL(R_v-diagv)
     min_r_v = MINVAL(R_v-diagv)

     DEALLOCATE (diagv)
     DEALLOCATE (diagu)

     CALL gc_rmax(1,n_proc,ierr,max_r_u)
     CALL gc_rmax(1,n_proc,ierr,max_r_v)
     CALL gc_rmax(1,n_proc,ierr,min_r_u)
     CALL gc_rmax(1,n_proc,ierr,min_r_v)

  IF ( me == 0 ) THEN

    write(6,fmt='(A)') '=============================================='
    WRITE(6,fmt='(A)') 'EG_IMPLICIT_HORZ_DRAG, min / max of S_u/S_v'
    WRITE(6,fmt='(A,2E15.5)') 'S_U:',max_r_u,min_r_u
    WRITE(6,fmt='(A,2E15.5)') 'S_V:',max_r_v,min_r_v
    write(6,fmt='(A)') '=============================================='

  END IF

END IF


IF (lhook) CALL dr_hook('EG_IMPLICIT_HORZ_DRAG',zhook_out,zhook_handle)
RETURN

END SUBROUTINE eg_implicit_horz_drag

END MODULE eg_implicit_horz_drag_mod
