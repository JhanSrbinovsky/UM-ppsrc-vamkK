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
!   u-grid.
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

      Subroutine lbc_u_a_to_c (                                         &
                 u_p                                                    &
      ,          u                                                      &
      ,          lbc_size_p                                             &
      ,          lbc_size_u                                             &
      ,          lbc_levels                                             &
      ,          lbc_row_len                                            &
      ,          lbc_rows                                               &
      ,          rimwidth                                               &
      ,          halo_x                                                 &
      ,          halo_y                                                 &
      ,          l_var_lbc                                              &
      ,          l_same_rot                                             &
      ,          l_eg_grid                                              &
      ,          lambda_p_in                                            &
      ,          lambda_u_in                                            &
       )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

      Integer ::  lbc_size_p
      Integer ::  lbc_size_u
      Integer ::  lbc_levels
      Integer ::  lbc_row_len
      Integer ::  lbc_rows
      Integer ::  rimwidth
      Integer ::  halo_x
      Integer ::  halo_y

      Real    ::  u_p (lbc_size_p, lbc_levels)
      Real    ::  u   (lbc_size_u, lbc_levels)

      LOGICAL ::  l_var_lbc
      LOGICAL ::  l_same_rot
      LOGICAL ::  l_eg_grid
      
! Input VarRes grid info in degrees      
      REAL :: lambda_p_in ( 1-halo_x: lbc_row_len + halo_x ) 
      REAL :: lambda_u_in ( 1-halo_x: lbc_row_len + halo_x )
      
      INTEGER ::  ipt_p, ipt_u
      INTEGER ::  lbc_row_len_p
      INTEGER ::  lbc_row_len_u
      INTEGER ::  lbc_rows_p
      INTEGER ::  lbc_rows_u
      INTEGER ::  level
      INTEGER ::  offset(2) ! Offset in x, offset in y
      INTEGER ::  offset_vr 

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
          
! North

      IF (lhook) CALL dr_hook('LBC_U_A_TO_C',zhook_in,zhook_handle)
      ipt_p = 1
      ipt_u = 1
      
      IF (l_eg_grid) THEN
        lbc_row_len_p = lbc_row_len + 2*halo_x + 1
        lbc_rows_p    = rimwidth + halo_y + 1
        lbc_row_len_u = lbc_row_len_p - 1
        lbc_rows_u    = lbc_rows_p - 1
        offset(1)     = 1
        offset(2)     = 0
        offset_vr  = - halo_x
      ELSE
        lbc_row_len_p = lbc_row_len + 2*halo_x
        lbc_rows_p    = rimwidth + halo_y + 1
        lbc_row_len_u = lbc_row_len_p - 1
        lbc_rows_u    = lbc_rows_p - 1
        offset(1)     = 0
        offset(2)     = 1
        offset_vr  = - halo_x
      END IF
      
      do level =1, lbc_levels

! DEPENDS ON: u_a_to_c
       call u_a_to_c (                                                  &
            u_p(ipt_p,level),                                           &
            u(ipt_u,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_u,                                              &
            lbc_rows_u,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            l_same_rot,                                                 &
            offset_vr,                                                  &
            lbc_row_len,                                                &
            halo_x,                                                     &
            lambda_p_in,                                                &
            lambda_u_in )

      enddo

! East

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_u = ipt_u + lbc_row_len_u * lbc_rows_u
      
      IF (l_eg_grid) THEN
        lbc_row_len_p = rimwidth + halo_x + 1
        lbc_rows_p    = lbc_rows - 2 * rimwidth + 2
        lbc_row_len_u = lbc_row_len_p - 1
        lbc_rows_u    = lbc_rows_p - 2
        offset(1)     = 1
        offset(2)     = 1
        offset_vr  = lbc_row_len - rimwidth - 1
      ELSE
        lbc_row_len_p = rimwidth + halo_x + 1
        lbc_rows_p    = lbc_rows - 2 * rimwidth
        lbc_row_len_u = lbc_row_len_p - 1
        lbc_rows_u    = lbc_rows_p
        offset(1)     = 0
        offset(2)     = 0
        offset_vr  = lbc_row_len - rimwidth - 1
      END IF
      
      do level =1, lbc_levels

! DEPENDS ON: u_a_to_c
       call u_a_to_c (                                                  &
            u_p(ipt_p,level),                                           &
            u(ipt_u,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_u,                                              &
            lbc_rows_u,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            l_same_rot,                                                 &
            offset_vr,                                                  &
            lbc_row_len,                                                &
            halo_x,                                                     &
            lambda_p_in,                                                &
            lambda_u_in )

      enddo

! South

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_u = ipt_u + lbc_row_len_u * lbc_rows_u
      
      IF (l_eg_grid) THEN
        lbc_row_len_p = lbc_row_len + 2*halo_x + 1
        lbc_rows_p    = rimwidth + halo_y + 1
        lbc_row_len_u = lbc_row_len_p - 1
        lbc_rows_u    = lbc_rows_p - 1
        offset(1)     = 1
        offset(2)     = 1
        offset_vr = - halo_x
      ELSE
        lbc_row_len_p = lbc_row_len + 2*halo_x
        lbc_rows_p    = rimwidth + halo_y + 1
        lbc_row_len_u = lbc_row_len_p - 1
        lbc_rows_u    = lbc_rows_p - 1
        offset(1)     = 0
        offset(2)     = 0
        offset_vr = - halo_x
      END IF
      
       
      do level =1, lbc_levels

! DEPENDS ON: u_a_to_c
       call u_a_to_c (                                                  &
            u_p(ipt_p,level),                                           &
            u(ipt_u,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_u,                                              &
            lbc_rows_u,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            l_same_rot,                                                 &
            offset_vr,                                                  &
            lbc_row_len,                                                &
            halo_x,                                                     &
            lambda_p_in,                                                &
            lambda_u_in )

      enddo

! West

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_u = ipt_u + lbc_row_len_u * lbc_rows_u
      
      IF (l_eg_grid) THEN
        lbc_row_len_p = rimwidth + halo_x + 1
        lbc_rows_p    = lbc_rows - 2 * rimwidth + 2
        lbc_row_len_u = lbc_row_len_p - 1
        lbc_rows_u    = lbc_rows_p - 2
        offset(1)     = 1
        offset(2)     = 1
        offset_vr = - halo_x 
      ELSE
        lbc_row_len_p = rimwidth + halo_x + 1
        lbc_rows_p    = lbc_rows - 2 * rimwidth
        lbc_row_len_u = lbc_row_len_p - 1
        lbc_rows_u    = lbc_rows_p
        offset(1)     = 0
        offset(2)     = 0
        offset_vr = - halo_x 
      END IF
      
      do level =1, lbc_levels

! DEPENDS ON: u_a_to_c
       call u_a_to_c (                                                  &
            u_p(ipt_p,level),                                           &
            u(ipt_u,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_u,                                              &
            lbc_rows_u,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            l_same_rot,                                                 &
            offset_vr,                                                  &
            lbc_row_len,                                                &
            halo_x,                                                     &
            lambda_p_in,                                                &
            lambda_u_in )

      enddo

      IF (lhook) CALL dr_hook('LBC_U_A_TO_C',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE lbc_u_a_to_c
