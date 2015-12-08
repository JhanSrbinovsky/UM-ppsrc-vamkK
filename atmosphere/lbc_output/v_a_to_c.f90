! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolate wind components from P to U/V grids.
!
! Description:
!   Interpolates v-component from the p-grid to
!   v-grid.
!
! Method:
!   u-comp and v-comp are both at all p-grid points after being rotated.
!   The u_comp is interpolated in the x_direction to the u-grid.
!   The v_comp is interpolated in the y_direction to the v-grid.
!   Linear interpolation is used. The u-grid and v-grid are within
!   the p-grid so there is no approximation round the edges.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

      SUBROUTINE v_a_to_c (                                             &
                 v_p                                                    &
      ,          v                                                      &
      ,          lbc_row_len_p                                          &
      ,          lbc_rows_p                                             &
      ,          lbc_row_len_v                                          &
      ,          lbc_rows_v                                             &
      ,          offset                                                 &
      ,          l_var_lbc                                              &
      ,          l_same_rot                                             &
      ,          offset_vr                                              &
      ,          lbc_rows                                               &
      ,          halo_y                                                 &
      ,          phi_p_in                                               &
      ,          phi_v_in                                               &
       )

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      INTEGER ::  lbc_row_len_p
      INTEGER ::  lbc_row_len_v
      INTEGER ::  lbc_rows_p
      INTEGER ::  lbc_rows_v
      INTEGER ::  offset(2)
      INTEGER ::  lbc_rows  
      INTEGER ::  halo_y      
      INTEGER ::  offset_vr 
      
      REAL    ::  v_p (lbc_row_len_p, lbc_rows_p)
      REAL    ::  v   (lbc_row_len_v, lbc_rows_v)

      LOGICAL l_var_lbc
      LOGICAL l_same_rot

! Input VarRes grid info in degrees           
      REAL    ::  phi_p_in ( 1-halo_y: lbc_rows + halo_y )   
      REAL    ::  phi_v_in ( 1-halo_y: lbc_rows + halo_y ) 
     
      INTEGER :: i,j
      REAL  weight1, weight2 

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('V_A_TO_C',zhook_in,zhook_handle)

      IF (L_same_rot) THEN  ! v's do not need averaging

        DO j = 1, lbc_rows_v
          DO i = 1, lbc_row_len_v
            v(i,j) = v_p(offset(1)+i,offset(2)+j)
          END DO
        END DO

      ELSE  ! grids have different rotations    

       IF (l_var_lbc) THEN
      
        DO j = 1, lbc_rows_v 
          weight1=(phi_v_in(offset_vr+j) - phi_p_in(offset_vr+j))/      &
                 (phi_p_in(offset_vr+j+1) - phi_p_in(offset_vr+j))
          weight2= 1.0 - weight1
          DO i = 1, lbc_row_len_v
            v(i,j) = weight1 * v_p(offset(1)+i,j)                          &
                    + weight2 * v_p(offset(1)+i,j+1) 
          END DO
        END DO

       ELSE  ! regular grid

        DO j = 1, lbc_rows_v
          DO i = 1, lbc_row_len_v
            v(i,j) = 0.5 * ( v_p(offset(1)+i,j) + v_p(offset(1)+i,j+1) )
          END DO
        END DO

       END IF ! l_var_lbc

      END IF ! L_same_rot

      IF (lhook) CALL dr_hook('V_A_TO_C',zhook_out,zhook_handle)

      RETURN
      END SUBROUTINE v_a_to_c
