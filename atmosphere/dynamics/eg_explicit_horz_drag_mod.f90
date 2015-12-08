! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_explicit_horz_drag_mod

! Subroutine: eg_explicit_horz_drag
!
! Description: Computes the source terms arising from a simple
!              implementation of boundary layer drag.
!
! Method: 
!              Du/dt= -Cd |u|* u /h
!
!              where
!                      |u| = sqrt(u^2 + v^2),
!                       h  = BL depth (set to 1 km)
!              and the drag coefficient Cd is a function of height
!              above the ground, z. According to Stull's boundary layer
!              book (p266) for a neutrally stratified atmosphere:
!
!              Cd= K^2 [ ln (z/z0) ] ^(-2)
!
!              where
!                     z0 is the roughness length (0.1 m) and
!                     K  is the von Karman constant.
!
!              NOTE: in the implementation Cd already contains
!                    the factor dt/h to save on the multiplication.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: DYNAMICS ADVECTION
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

CONTAINS

SUBROUTINE eg_explicit_horz_drag(r_u,r_v,u,v)


USE atm_fields_bounds_mod, ONLY : udims,vdims,udims_s,vdims_s
USE timestep_mod,          ONLY : timestep

USE eg_horz_drag_mod,      ONLY : cd_u,cd_v,l_impl_horz_drag,         &
                                  r_u_store,r_v_store

USE proc_info_mod,         ONLY : me,n_proc
USE parkind1,              ONLY : jprb, jpim
USE yomhook,               ONLY : lhook, dr_hook

USE PrintStatus_mod

IMPLICIT NONE

REAL, INTENT (INOUT) ::                                               &
  r_u(udims_s%i_start:udims_s%i_end,                                  &
      udims_s%j_start:udims_s%j_end,                                  &
      udims_s%k_start:udims_s%k_end)                                  &
, r_v(vdims_s%i_start:vdims_s%i_end,                                  &
      vdims_s%j_start:vdims_s%j_end,                                  &
      vdims_s%k_start:vdims_s%k_end)    

REAL, INTENT (IN) ::                                                  &
    u(udims_s%i_start:udims_s%i_end,                                  &
      udims_s%j_start:udims_s%j_end,                                  &
      udims_s%k_start:udims_s%k_end)                                  &
,   v(vdims_s%i_start:vdims_s%i_end,                                  &
      vdims_s%j_start:vdims_s%j_end,                                  &
      vdims_s%k_start:vdims_s%k_end)


!     local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*),PARAMETER :: RoutineName='EG_EXPLICIT_HORZ_DRAG' 
INTEGER                    ErrorStatus                        

INTEGER i
INTEGER j
INTEGER k
INTEGER ierr

REAL usource(udims_s%i_start:udims_s%i_end,                           &
             udims_s%j_start:udims_s%j_end) ! temporary storage for u
                                            ! increment (source term) to
                                            ! avoid if test in inner most
                                            ! do loop
REAL vsource(vdims_s%i_start:vdims_s%i_end,                           &
             vdims_s%j_start:vdims_s%j_end) ! temporary storage for v
                                            ! increment (source term)
  
REAL max_r_u   ! diagnostic, min/max of r_u/r_v
REAL max_r_v
REAL min_r_u
REAL min_r_v


ErrorStatus=0 

IF (lhook) CALL dr_hook('EG_EXPLICIT_HORZ_DRAG',zhook_in,zhook_handle)

usource = 0.
vsource = 0.

DO k=udims%k_start,udims%k_end
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end

      usource(i,j) =                                                  &
                    cd_u(i,j,k)* SQRT(                                &
                                 (.25*(v(i,j  ,k)+v(i+1,j  ,k)        &
                                      +v(i,j-1,k)+v(i+1,j-1,k)))**2   &
                                 +u(i,j,k)**2                         &
                                 ) * u(i,j,k)

    END DO
  END DO

  r_u(:,:,k) =  r_u(:,:,k) - usource(:,:)

  IF (l_impl_horz_drag) r_u_store(:,:,k) = - usource(:,:)
 
END DO

DO k=vdims%k_start,vdims%k_end
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end

      vsource(i,j) =                                                  &
                    cd_v(i,j,k)* SQRT(                                &
                                 (.25*(u(i,j  ,k)+u(i-1,j  ,k)        &
                                      +u(i,j+1,k)+u(i-1,j+1,k)))**2   &
                                 +v(i,j,k)**2                         &
                                 ) * v(i,j,k)

    END DO
  END DO

  r_v(:,:,k) =  r_v(:,:,k) - vsource(:,:)

  IF (l_impl_horz_drag) r_v_store(:,:,k) = - vsource(:,:)

END DO


IF ( PrintStatus == PrStatus_Diag) THEN

     max_r_u = MAXVAL(R_u)
     min_r_u = MINVAL(R_u)
     max_r_v = MAXVAL(R_v)
     min_r_v = MINVAL(R_v)

     CALL gc_rmax(1,n_proc,ierr,max_r_u)
     CALL gc_rmax(1,n_proc,ierr,max_r_v)
     CALL gc_rmax(1,n_proc,ierr,min_r_u)
     CALL gc_rmax(1,n_proc,ierr,min_r_v)

  IF ( me == 0 ) THEN

    write(6,fmt='(A)') '=============================================='
    WRITE(6,fmt='(A)') 'EG_EXPLICIT_HORZ_DRAG, min / max of R_u/R_v'
    WRITE(6,fmt='(A,2E15.5)') 'R_U:',max_r_u,min_r_u
    WRITE(6,fmt='(A,2E15.5)') 'R_V:',max_r_v,min_r_v
    write(6,fmt='(A)') '=============================================='

  END IF

END IF

IF (lhook) CALL dr_hook('EG_EXPLICIT_HORZ_DRAG',zhook_out,zhook_handle)
RETURN

END SUBROUTINE eg_explicit_horz_drag

END MODULE eg_explicit_horz_drag_mod
