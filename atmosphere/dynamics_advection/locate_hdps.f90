! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE locate_hdps_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE locate_hdps(i_out, j_out, wght_xi1, wght_xi2,                    &
                       xi1_out, xi2_out, dsm1,                              &
                       pnt_type, dim_i, dim_j, dim_k) 

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE um_types, ONLY: integer32

USE horiz_grid_mod
USE lookup_table_mod
USE UM_ParParams
USE proc_info_mod,         ONLY : mype=>me
USE dynamics_input_mod,    ONLY : l_regular

IMPLICIT NONE

!
! Description:
!
!          Locate the position of the departure points
!          and return the integer location of the pivot point
!          and offset/weight for use in interpolation.
!
! Method:
!
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

INTEGER, INTENT(IN)  :: dim_i, dim_j, dim_k, pnt_type, dsm1

INTEGER(KIND=integer32), INTENT(OUT) ::                                     &
                        i_out(dim_i, dim_j, dim_k),                         &
                        j_out(dim_i, dim_j, dim_k)

REAL,    INTENT(IN)  :: xi1_out(dim_i, dim_j, dim_k),                       &
                        xi2_out(dim_i, dim_j, dim_k)

REAL,    INTENT(OUT) :: wght_xi1(dim_i, dim_j, dim_k),                      &
                        wght_xi2(dim_i, dim_j, dim_k)

INTEGER              :: i, j, k, id
REAL                 :: x_off, y_off, temp
REAL                 :: rdxi1, rdxi2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

REAL :: real_temp_var
INTEGER :: my_floor
my_floor(real_temp_var) = INT(real_temp_var - 0.5 + SIGN(0.5,real_temp_var))


IF (lhook) CALL dr_hook('LOCATE_HDPS',zhook_in,zhook_handle)


IF ( l_regular ) THEN

  rdxi1 = 1.0/delta_xi1
  rdxi2 = 1.0/delta_xi2

  SELECT CASE(pnt_type)
     CASE(fld_type_u)
        x_off = 1.0
        y_off = 0.5
     CASE(fld_type_v)
        x_off = 0.5
        y_off = 1.0
     CASE(fld_type_w, fld_type_p)
        x_off = 0.5
        y_off = 0.5
  END SELECT

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(static) PRIVATE(k,j,i, temp)    &
!$OMP&             SHARED(dim_k, dim_j, dim_i, dsm1,                      &
!$OMP&                    xi1_out,base_xi1,rdxi1,i_out,wght_xi1,x_off,    &
!$OMP&                    xi2_out,base_xi2,rdxi2,j_out,wght_xi2,y_off)
!cdir collapse
  DO k = 1, dim_k
     DO j = 1, dim_j
        DO i = 1, dim_i
           temp            = (xi1_out(i,j,k)-base_xi1)*rdxi1
           i_out(i,j,k)    = my_floor(x_off + temp)
           wght_xi1(i,j,k) = temp - ( i_out(i,j,k) - x_off )
        END DO
     END DO

     DO j = 1, dim_j
        DO i = 1, dim_i
           temp            = (xi2_out(i,j,k) - base_xi2)*rdxi2
           j_out(i,j,k)    = my_floor(y_off + temp)
           wght_xi2(i,j,k) = temp - ( j_out(i,j,k) - y_off)
           j_out(i,j,k)    = j_out(i,j,k) - dsm1
        END DO
     END DO
  END DO
!$OMP END PARALLEL DO

ELSE  ! variable resolution

   IF( pnt_type == fld_type_u ) THEN
      rdxi1 = 1.0/delta_xi1_u
   ELSE
      rdxi1 = 1.0/delta_xi1_p
   END IF

   IF( pnt_type == fld_type_v ) THEN
      rdxi2 = 1.0/delta_xi2_v
   ELSE
      rdxi2 = 1.0/delta_xi2_p
   END IF

! DO xi1 direction first

!$OMP  PARALLEL DO  DEFAULT(NONE) SCHEDULE(static) PRIVATE(k,j,i,temp,id)   &
!$OMP&              SHARED(dim_k, dim_j, dim_i, dsm1,                       &
!$OMP&                     xi1_out,base_xi1_u,base_xi1_p, rdxi1,            &
!$OMP&                     xi2_out,base_xi2_v,base_xi2_p, rdxi2,            &
!$OMP&                     i_out,wght_xi1, j_out,wght_xi2,                  &
!$OMP&                     glob_xi1_u, glob_xi1_p, glob_xi2_v, glob_xi2_p,  &
!$OMP&                     glob_rdxi1_u, glob_rdxi1_p, glob_rdxi2_v,        &
!$OMP&                     glob_rdxi2_p, i_lkup_u, i_lkup_p,                &
!$OMP&                     j_lkup_v, j_lkup_p, pnt_type)
   DO k = 1, dim_k
     IF(pnt_type == fld_type_u) THEN
       DO j = 1, dim_j
         DO i = 1, dim_i
           temp            = (xi1_out(i,j,k)-base_xi1_u)*rdxi1
           id              = i_lkup_u( INT(temp) )
           wght_xi1(i,j,k) = (xi1_out(i,j,k) - glob_xi1_u(id))               &
                             *glob_rdxi1_u(id)
           IF( wght_xi1(i,j,k) > 1.0 ) THEN
              id              = id + 1
              wght_xi1(i,j,k) = (xi1_out(i,j,k) - glob_xi1_u(id))            &
                                  *glob_rdxi1_u(id)
!          ELSE IF( wght_xi1(i,j,k) < 0.0 ) THEN
!             id              = id - 1
!             wght_xi1(i,j,k) = (xi1_out(i,j,k) - glob_xi1_u(id))            &
!                               *glob_rdxi1_u(id)
             END IF
             i_out(i,j,k) = id + 1 ! offset for u -> p
           END DO
         END DO
     ELSE
       DO j = 1, dim_j
         DO i = 1, dim_i
           temp            = (xi1_out(i,j,k)-base_xi1_p)*rdxi1
           id              = i_lkup_p( INT(temp) )
           wght_xi1(i,j,k) = (xi1_out(i,j,k) - glob_xi1_p(id))               &
                             *glob_rdxi1_p(id)
           IF( wght_xi1(i,j,k) > 1.0 ) THEN
              id              = id + 1
              wght_xi1(i,j,k) = (xi1_out(i,j,k) - glob_xi1_p(id))          &
                                  *glob_rdxi1_p(id)
!          ELSE IF( wght_xi1(i,j,k) < 0.0 ) THEN
!             id              = id - 1
!             wght_xi1(i,j,k) = (xi1_out(i,j,k) - glob_xi1_p(id))          &
!                                 *glob_rdxi1_p(id)
           END IF
           i_out(i,j,k) = id
         END DO
       END DO
     END IF

! Now xi2

     IF(pnt_type == fld_type_v) THEN
       DO j = 1, dim_j
         DO i = 1, dim_i
           temp            = (xi2_out(i,j,k)-base_xi2_v)*rdxi2
           id              = j_lkup_v( INT(temp) )
           wght_xi2(i,j,k) = (xi2_out(i,j,k) - glob_xi2_v(id))               &
                             *glob_rdxi2_v(id)
           IF( wght_xi2(i,j,k) > 1.0 ) THEN
              id              = id + 1
              wght_xi2(i,j,k) = (xi2_out(i,j,k) - glob_xi2_v(id))            &
                                *glob_rdxi2_v(id)
!          ELSE IF( wght_xi2(i,j,k) < 0.0 ) THEN
!             id              = id - 1
!             wght_xi2(i,j,k) = (xi2_out(i,j,k) - glob_xi2_v(id))            &
!                               *glob_rdxi2_v(id)
           END IF
           j_out(i,j,k) = id + 1 - dsm1 ! offset for v -> p
         END DO
       END DO
     ELSE
       DO j = 1, dim_j
         DO i = 1, dim_i
           temp            = (xi2_out(i,j,k)-base_xi2_p)*rdxi2
           id              = j_lkup_p( INT(temp) )
           wght_xi2(i,j,k) = (xi2_out(i,j,k) - glob_xi2_p(id))               &
                             *glob_rdxi2_p(id)
           IF( wght_xi2(i,j,k) > 1.0 ) THEN
              id              = id + 1
              wght_xi2(i,j,k) = (xi2_out(i,j,k) - glob_xi2_p(id))            &
                                *glob_rdxi2_p(id)
!          ELSE IF( wght_xi2(i,j,k) < 0.0 ) THEN
!             id              = id - 1
!             wght_xi2(i,j,k) = (xi2_out(i,j,k) - glob_xi2_p(id))            &
!                               *glob_rdxi2_p(id)
           END IF
           j_out(i,j,k) = id - dsm1
         END DO
       END DO
     END IF
   END DO
!$OMP END PARALLEL DO

END IF ! L_regular

IF (lhook) CALL dr_hook('LOCATE_HDPS',zhook_out,zhook_handle)

  END SUBROUTINE locate_hdps
END MODULE locate_hdps_mod
