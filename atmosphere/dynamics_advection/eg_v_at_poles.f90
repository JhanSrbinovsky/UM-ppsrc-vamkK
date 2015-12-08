! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_v_at_poles_mod

USE global_2d_sums_mod, ONLY: &
    global_2d_sums

IMPLICIT NONE
CONTAINS
SUBROUTINE eg_v_at_poles(u,v, pole_sgn_flip,                                 &
                   jj_u_pole, jj_v_pole, udim_in,vdim_in,u_pole_v,u_pole_u,  &
                   v_pole_out,xi1_pole_out)


USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE conversions_mod
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE proc_info_mod, ONLY : global_row_length, gc_proc_row_group
USE eg_parameters_mod, ONLY : pole_consts
USE ereport_mod, ONLY : ereport

IMPLICIT NONE
!
! Description: Calcualte v at the north and south poles
!  
!
! Method: Chapter 8, ENDGame formulation version 1.01
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


REAL,    INTENT(IN) :: pole_sgn_flip
INTEGER, INTENT(IN) :: jj_v_pole, jj_u_pole

TYPE (array_dims)   :: udim_in,vdim_in

REAL, INTENT(INOUT) :: v(vdim_in%i_start:vdim_in%i_end,                        &
                         vdim_in%j_start:vdim_in%j_end,                        &
                         vdim_in%k_start:vdim_in%k_end)

REAL, INTENT(IN)    :: u(udim_in%i_start:udim_in%i_end,                        &
                         udim_in%j_start:udim_in%j_end,                        &
                         udim_in%k_start:udim_in%k_end)

! u-component on pole on v-staggered (or p-like) point
REAL, OPTIONAL, INTENT(INOUT)  :: u_pole_v(vdims%i_start:vdims%i_end,1,        &
                                           vdims%k_start:vdims%k_end)

! u-component at pole on u-staggered point
REAL, OPTIONAL, INTENT(INOUT)  :: u_pole_u(udims%i_start:udims%i_end,1,        &
                                           udims%k_start:udims%k_end)

REAL, OPTIONAL, INTENT(OUT) ::   v_pole_out(udim_in%k_start:udim_in%k_end)
REAL, OPTIONAL, INTENT(OUT) :: xi1_pole_out(udim_in%k_start:udim_in%k_end)

! Local variables

INTEGER :: i, k, info
REAL    :: c, d, e, f, dx, pi_fac
REAL    :: a, b, g, xi1_pole, v_pole
REAL    :: glob_sum(udims%k_start:udims%k_end,3)
REAL    :: temp1, temp2, cs, sn



REAL    :: tmp_loc(udims%i_start:udims%i_end,udims%k_start:udims%k_end)


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_V_AT_POLES',zhook_in,zhook_handle)

pi_fac = 1.0/pi

! First do south pole

c = pole_consts(1)
d = pole_consts(2)
e = pole_consts(3)
f = pole_consts(4)
DO k = udims%k_start, udims%k_end
   DO i = udims%i_start, udims%i_end
      dx = pi_fac*( xi1_p(i+1) - xi1_p(i) )





      tmp_loc(i,k) = dx*u(i,jj_u_pole,k)*SIN(xi1_u(i))

   END DO
END DO






CALL global_2d_sums(tmp_loc(:,:),udims%i_end-udims%i_start+1,        & 
                    1, 0, 0, (udims%k_end-udims%k_start+1),          &
                    glob_sum(:,1),  gc_proc_row_group) 


DO k = udims%k_start, udims%k_end
   DO i = udims%i_start, udims%i_end
      dx = pi_fac*( xi1_p(i+1) - xi1_p(i) )
      tmp_loc(i,k) = dx*u(i,jj_u_pole,k)*COS(xi1_u(i))
   END DO
END DO
 


CALL global_2d_sums(tmp_loc(:,:),udims%i_end-udims%i_start+1,        & 
                    1, 0, 0, (udims%k_end-udims%k_start+1),          &
                    glob_sum(:,2),  gc_proc_row_group) 

DO k = udims%k_start, udims%k_end
   DO i = udims%i_start, udims%i_end
      dx = pi_fac*( xi1_p(i+1) - xi1_p(i) )
     tmp_loc(i,k) = dx*u(i,jj_u_pole,k)
   END DO
END DO

CALL global_2d_sums(tmp_loc(:,:),udims%i_end-udims%i_start+1,        & 
                    1, 0, 0, (udims%k_end-udims%k_start+1),          &
                    glob_sum(:,3),  gc_proc_row_group) 




DO k = udims%k_start, udims%k_end
   a     = glob_sum(k,1)
   b     = glob_sum(k,2)
   g     = glob_sum(k,3)

   temp1 = (1.0-c-2.0*e**2)*(b-f*g) - (a-e*g)*(d-2.0*e*f)
   temp2 =-(1.0+c-2.0*f**2)*(a-e*g) + (b-f*g)*(d-2.0*e*f)

   temp1 = pole_sgn_flip*temp1
   temp2 = pole_sgn_flip*temp2

   IF( temp1 == 0.0 .AND. temp2 == 0.0 ) THEN
         xi1_pole = 0.0
   ELSE
         xi1_pole = ATAN2(temp1,temp2)
   END IF

   cs     = COS(xi1_pole)
   sn     = SIN(xi1_pole)
   dx     = SIN(xi2_p(jj_u_pole))
   temp1  = dx*(1.0 - c*(1.0-2.0*sn**2) - 2.0*d*sn*cs                 &
                    - 2.0*(e*cs - f*sn)**2 )
   v_pole = ( (a-e*g)*cs - (b-f*g)*sn )/temp1

   DO i = vdims%i_start, vdims%i_end
      v(i,jj_v_pole,k) = v_pole*COS(xi1_p(i)-xi1_pole)
   END DO

   IF (PRESENT(u_pole_v)) THEN
      DO i = vdims%i_start, vdims%i_end
         u_pole_v(i,1,k) =  -1.0*pole_sgn_flip*v_pole*SIN(xi1_p(i)-xi1_pole)
      END DO

   END IF

   IF (PRESENT(u_pole_u)) THEN

      DO i = udims%i_start, udims%i_end
         u_pole_u(i,1,k) = -1.0*pole_sgn_flip*v_pole*SIN(xi1_u(i)-xi1_pole)
      END DO

   END IF

   IF (PRESENT(  v_pole_out))   v_pole_out(k)=   v_pole
   IF (PRESENT(xi1_pole_out)) xi1_pole_out(k)= xi1_pole


END DO


 


IF (lhook) CALL dr_hook('EG_V_AT_POLES',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_v_at_poles


! =====================================================================

SUBROUTINE eg_uv_at_poles_Bgrid &
                  (u, row_length, rows, n_rows,                       &
                   offx, offy, pole_consts, pole_sgn_flip,            &
                   j_pole_src, j_pole_dest, v)


USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE conversions_mod
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE proc_info_mod, ONLY :  global_row_length, gc_proc_row_group

IMPLICIT NONE
!
! Description: Calculate u or v at the north and south poles on 
!              B grid pressure levels, one level at a time
!  
!
! Method: Chapter 8, ENDGame formulation version 1.01
!  
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


REAL,    INTENT(IN) :: pole_sgn_flip
INTEGER, INTENT(IN) :: j_pole_src, j_pole_dest

! Array dimensions
INTEGER :: offx, offy
INTEGER :: row_length, rows, n_rows

REAL  pole_consts(4)

! u and v on B grid
! Note: if v present, v is computed; otherwise u is computed
REAL, INTENT(INOUT) :: u(udims%i_start:udims%i_end, vdims%j_start:vdims%j_end)
REAL, INTENT(OUT), OPTIONAL :: &
                       v(udims%i_start:udims%i_end, vdims%j_start:vdims%j_end)

! Local variables

INTEGER :: i, info
REAL    :: c, d, e, f, dx, pi_fac
REAL    :: a, b, g, xi1_pole
REAL    :: mag_pole_wind  ! Magnitude of polar wind
REAL    :: glob_sum(3)
REAL    :: temp1, temp2, cs, sn
REAL, DIMENSION(udims%i_start:udims%i_end,3) :: tmp_loc

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_V_AT_POLES_Bgrid',zhook_in,zhook_handle)

pi_fac = 1.0/pi

c = pole_consts(1)
d = pole_consts(2)
e = pole_consts(3)
f = pole_consts(4)

tmp_loc(:,:) = 0.0

DO i = udims%i_start, udims%i_end
  dx = pi_fac*( xi1_p(i+1) - xi1_p(i) )
  tmp_loc(i,1) = dx *u(i,j_pole_src) *SIN(xi1_u(i))
  tmp_loc(i,2) = dx *u(i,j_pole_src) *COS(xi1_u(i))
  tmp_loc(i,3) = dx *u(i,j_pole_src)
END DO


CALL global_2d_sums(tmp_loc,row_length,        & 
                    1, 0, 0, 3,                &
                    glob_sum,  gc_proc_row_group) 

a     = glob_sum(1)
b     = glob_sum(2)
g     = glob_sum(3)

temp1 = (1.0-c-2.0*e**2)*(b-f*g) - (a-e*g)*(d-2.0*e*f)
temp2 =-(1.0+c-2.0*f**2)*(a-e*g) + (b-f*g)*(d-2.0*e*f)

temp1 = pole_sgn_flip*temp1
temp2 = pole_sgn_flip*temp2

IF( temp1 == 0.0 .AND. temp2 == 0.0 ) THEN
      xi1_pole = 0.0
ELSE
      xi1_pole = ATAN2(temp1,temp2)
END IF

cs     = COS(xi1_pole)
sn     = SIN(xi1_pole)
dx     = SIN(xi2_p(j_pole_src))
temp1  = dx*(1.0 - c*(1.0-2.0*sn**2) - 2.0*d*sn*cs                 &
                 - 2.0*(e*cs - f*sn)**2 )
mag_pole_wind = ( (a-e*g)*cs - (b-f*g)*sn )/temp1

DO i = udims%i_start, udims%i_end
  IF (PRESENT(v)) THEN
    v(i,j_pole_dest) = mag_pole_wind*COS(xi1_u(i)-xi1_pole)
  ELSE 
    u(i,j_pole_dest) = pole_sgn_flip*mag_pole_wind*SIN(xi1_u(i)-xi1_pole)
  END IF
END DO


IF (lhook) CALL dr_hook('EG_V_AT_POLES_Bgrid',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_uv_at_poles_Bgrid

END MODULE eg_v_at_poles_mod
