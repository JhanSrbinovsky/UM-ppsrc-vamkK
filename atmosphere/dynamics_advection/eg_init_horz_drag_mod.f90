! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_init_horz_drag_mod
! Subroutine: eg_init_horz_drag
!
! Description: Initialises the height dependent horizontal drag
!              coefficient
!
! Method: 
!              Cd = K^2 [ ln (z/z0) ] ^(-2)
!
!              where
!                     z0 is the roughness length (0.1 m) and
!                     K  is the von Karman constant.
!
!              NOTE: in the implementation Cd contains
!                    the factor dt/h to save on the multiplication.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: DYNAMICS ADVECTION
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE
CONTAINS

SUBROUTINE eg_init_horz_drag(l_impl_horz_drag_in,l_expl_horz_drag_in)

USE atmos_constants_mod,   ONLY : vkman
USE atm_fields_bounds_mod, ONLY : udims_s,vdims_s,udims,vdims
USE earth_constants_mod,   ONLY : earth_radius
USE level_heights_mod,     ONLY : xi3_at_u => r_at_u,                 &
                                  xi3_at_v => r_at_v,                 &
                                  xi3_at_theta => r_theta_levels
USE timestep_mod,          ONLY : timestep
USE conversions_mod,       ONLY : pi

USE eg_horz_drag_mod,      ONLY : cd_u,cd_v,r_u_store,r_v_store,      &
                                  l_impl_horz_drag,l_expl_horz_drag

USE proc_info_mod,         ONLY : me
USE parkind1,              ONLY : jprb, jpim
USE yomhook,               ONLY : lhook, dr_hook
USE Ereport_mod
USE PrintStatus_mod
IMPLICIT NONE

LOGICAL, INTENT (IN) :: l_impl_horz_drag_in
LOGICAL, INTENT (IN) :: l_expl_horz_drag_in

!     local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*),PARAMETER :: RoutineName='EG_INIT_HORZ_DRAG' ! Name of 
                                                           ! routine

INTEGER                    ErrorStatus                     ! Error code

CHARACTER(LEN=80)             Message                         ! Text for 
                                                           ! output.

INTEGER i
INTEGER j
INTEGER k
INTEGER ierr

REAL eta

REAL heightaboveground
REAL, PARAMETER :: boundarylayerheight = 2000.
REAL, PARAMETER :: transitionheight    = 5000.


REAL, PARAMETER :: z0 = 0.1    ! roughness length
REAL, PARAMETER :: h  = 1000.  ! boundary layer depth


ErrorStatus=0 

l_impl_horz_drag=l_impl_horz_drag_in
l_expl_horz_drag=l_expl_horz_drag_in


IF (lhook) CALL dr_hook('EG_INIT_HORZ_DRAG',zhook_in,zhook_handle)

ALLOCATE (cd_u(  udims%i_start:udims%i_end,                           & 
                 udims%j_start:udims%j_end,                           &
                 udims%k_start:udims%k_end),STAT=ierr)

IF(ierr.ne.0) THEN
  WRITE(message,*) 'allocation of cd_u failed'
  CALL Ereport(RoutineName,ierr,message)
END IF

ALLOCATE (cd_v(  vdims%i_start:vdims%i_end,                           & 
                 vdims%j_start:vdims%j_end,                           &
                 vdims%k_start:vdims%k_end),STAT=ierr)

IF(ierr.ne.0) THEN
  WRITE(message,*) 'allocation of cd_v failed'
  CALL Ereport(RoutineName,ierr,message)
END IF


IF(l_impl_horz_drag .AND. l_expl_horz_drag) THEN

  ALLOCATE (r_u_store(udims_s%i_start:udims_s%i_end,                  & 
                      udims_s%j_start:udims_s%j_end,                  &
                      udims_s%k_start:udims_s%k_end),STAT=ierr)

  IF(ierr.ne.0) THEN
    WRITE(message,*) 'allocation of r_u_store failed'
    CALL Ereport(RoutineName,ierr,message)
  END IF

  r_u_store = 0.

  ALLOCATE (r_v_store(vdims_s%i_start:vdims_s%i_end,                  & 
                      vdims_s%j_start:vdims_s%j_end,                  &
                      vdims_s%k_start:vdims_s%k_end),STAT=ierr)

  IF(ierr.ne.0) THEN
    WRITE(message,*) 'allocation of r_v_store failed'
    CALL Ereport(RoutineName,ierr,message)
  END IF

  r_v_store = 0.

END IF



DO k=udims%k_start,udims%k_end
  DO j=udims%j_start,udims%j_end
    DO i=udims%i_start,udims%i_end


      heightaboveground = (xi3_at_u(i,j,k)-0.5*(xi3_at_theta(i,j,0)+  &
                                                xi3_at_theta(i+1,j,0)))

      cd_u(i,j,k) =  vkman**2  * 1./( LOG (heightaboveground/z0))**2  &
                                   / h * timestep

      IF(heightaboveground>transitionheight) THEN

        cd_u(i,j,k) = 0.

      ELSE

        IF(heightaboveground>boundarylayerheight) THEN

           eta = 1.-(heightaboveground-boundarylayerheight)           &
                    /transitionheight

           cd_u(i,j,k) = cd_u(i,j,k) * SIN( 0.5*pi*(eta))**2

        END IF

      END IF

    END DO
  END DO
END DO



DO k=vdims%k_start,vdims%k_end
  DO j=vdims%j_start,vdims%j_end
    DO i=vdims%i_start,vdims%i_end

      heightaboveground = (xi3_at_v(i,j,k)-0.5*(xi3_at_theta(i,j,0)+  &
                                                xi3_at_theta(i,j+1,0)))


      cd_v(i,j,k) =  vkman**2  * 1./( LOG (heightaboveground/z0))**2  &
                                   / h * timestep

      IF(heightaboveground>transitionheight) THEN

        cd_v(i,j,k) = 0.

      ELSE

        IF(heightaboveground>boundarylayerheight) THEN

           eta = 1.-(heightaboveground-boundarylayerheight)           &
                    /transitionheight

           cd_v(i,j,k) = cd_v(i,j,k) * SIN( 0.5*pi*(eta))**2

        END IF

      END IF

    END DO
  END DO
END DO


IF ( PrintStatus == PrStatus_Diag .AND. me == 0 ) THEN

  write(6,fmt='(A)') '================================================'
  write(6,fmt='(A)') 'level  | horizontal drag coefficient at i=1,j=1'
  write(6,fmt='(A)') '------------------------------------------------'

  DO k=udims%k_start,udims%k_end

    write(6,fmt='(I4,A,E15.5)') k,'   |', cd_u(1,1,k)*h/ timestep

  END DO

  write(6,fmt='(A)') '================================================'

END IF

IF (lhook) CALL dr_hook('EG_INIT_HORZ_DRAG',zhook_out,zhook_handle)
RETURN

END SUBROUTINE eg_init_horz_drag

END MODULE eg_init_horz_drag_mod
