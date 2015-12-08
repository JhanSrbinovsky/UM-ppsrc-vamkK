! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolate wind components from P to U/V grids.
!
! Description:
!   Interpolates u-component from the p-grid to
!   u-grid respectively.
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

      SUBROUTINE u_a_to_c (                                             &
                 u_p                                                    &
      ,          u                                                      &
      ,          lbc_row_len_p                                          &
      ,          lbc_rows_p                                             &
      ,          lbc_row_len_u                                          &
      ,          lbc_rows_u                                             &
      ,          offset                                                 &
      ,          l_var_lbc                                              &
      ,          l_same_rot                                             &
      ,          offset_vr                                              &
      ,          lbc_row_len                                            &
      ,          halo_x                                                 &
      ,          lambda_p_in                                            &
      ,          lambda_u_in                                            &
       )    

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      INTEGER ::  lbc_row_len_p
      INTEGER ::  lbc_row_len_u
      INTEGER ::  lbc_rows_p
      INTEGER ::  lbc_rows_u
      INTEGER ::  offset(2)
      INTEGER ::  lbc_row_len
      INTEGER ::  halo_x      
      INTEGER ::  offset_vr  
      
      REAL    ::  u_p (lbc_row_len_p, lbc_rows_p)
      REAL    ::  u   (lbc_row_len_u, lbc_rows_u)
      
      LOGICAL  l_var_lbc
      LOGICAL  l_same_rot

! Input VarRes grid info in degrees      
      REAL  :: Lambda_p_in ( 1-halo_x: lbc_row_len + halo_x )  
      REAL  :: Lambda_u_in ( 1-halo_x: lbc_row_len + halo_x )
      
      INTEGER :: i,j
      REAL  weight1, weight2 
     
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('U_A_TO_C',zhook_in,zhook_handle)

      IF (L_same_rot) THEN  ! u's do not need averaging

        DO j = 1, lbc_rows_u
          DO i = 1, lbc_row_len_u
            u(i,j) = u_p(offset(1)+i,offset(2)+j)
          END DO
        END DO

      ELSE ! for grid with different rotation     

       IF (l_var_lbc) THEN
                     
        DO i = 1, lbc_row_len_u
          weight1=(lambda_u_in(offset_vr+i) - lambda_p_in(offset_vr+i)) &
               /(lambda_p_in(offset_vr+i+1) - lambda_p_in(offset_vr+i))  
          weight2= 1.0 - weight1
         
          DO j = 1, lbc_rows_u        
            u(i,j) = weight1 *  u_p(i,offset(2)+j)                         &
                    + weight2 * u_p(i+1,offset(2)+j)
          END DO
        END DO

      ELSE ! standard regular grid     

        DO j = 1, lbc_rows_u
          DO i = 1, lbc_row_len_u
            u(i,j) = 0.5 * ( u_p(i,offset(2)+j) + u_p(i+1,offset(2)+j) )
          END DO
        END DO

       END IF !  l_var_lbc

      END IF !  L_same_rot

      IF (lhook) CALL dr_hook('U_A_TO_C',zhook_out,zhook_handle)

      RETURN
      END SUBROUTINE u_a_to_c
