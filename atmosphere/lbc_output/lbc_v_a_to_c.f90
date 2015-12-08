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

! *****************************************************************
! Routine to control v-comp interpolation from p-grid to v-grid
! *****************************************************************

      Subroutine lbc_v_a_to_c (                                         &
                 v_p                                                    &
      ,          v                                                      &
      ,          lbc_size_p                                             &
      ,          lbc_size_v                                             &
      ,          lbc_levels                                             &
      ,          lbc_row_len                                            &
      ,          lbc_rows                                               &
      ,          rimwidth                                               &
      ,          halo_x                                                 &
      ,          halo_y                                                 &
      ,          l_var_lbc                                              &
      ,          l_same_rot                                             &
      ,          l_eg_grid                                              &
      ,          phi_p_in                                               &
      ,          phi_v_in                                               &
       )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

      Integer ::  lbc_size_p    ! grid size on p grid
      Integer ::  lbc_size_v    ! grid size on v grid
      Integer ::  lbc_levels    ! no of lbc levels
      Integer ::  lbc_row_len   ! row length of lbc grid
      Integer ::  lbc_rows      ! rows on lbc grid
      Integer ::  rimwidth      !
      Integer ::  halo_x
      Integer ::  halo_y
      
      Real    ::  v_p (lbc_size_p, lbc_levels)  ! lbc v on p grid
      Real    ::  v   (lbc_size_v, lbc_levels)  ! lbc v on v grid

      Logical l_var_lbc
      LOGICAL l_same_rot
      LOGICAL ::  l_eg_grid
      
! Input VarRes grid info in degrees
      Real    ::  phi_p_in ( 1-halo_y: lbc_rows + halo_y )
      Real    ::  phi_v_in ( 1-halo_y: lbc_rows + halo_y )

! Local variables
      Integer ::  ipt_p, ipt_v
      Integer ::  lbc_row_len_p
      Integer ::  lbc_row_len_v
      Integer ::  lbc_rows_p
      Integer ::  lbc_rows_v
      Integer ::  level
      Integer ::  offset(2)
      Integer ::  offset_vr 

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      
! To be replaced with loop over iside (5.3)
! Set up row length, rows and start addresses in parlbcs

! North

      IF (lhook) CALL dr_hook('LBC_V_A_TO_C',zhook_in,zhook_handle)

      ! This is related to grid in lbc_interp_coeffs which sets up the grid
      ! which the p field is on.
      ipt_p = 1
      ipt_v = 1
      
      IF (l_eg_grid) THEN
        lbc_row_len_p = lbc_row_len + 2*halo_x + 1
        lbc_rows_p    = rimwidth + halo_y + 1
        lbc_row_len_v = lbc_row_len_p - 1
        lbc_rows_v    = lbc_rows_p - 1
        offset(1)     = 1
        offset(2)     = 1
        offset_vr = lbc_rows - rimwidth - 1      
      ELSE
        lbc_row_len_p = lbc_row_len + 2*halo_x
        lbc_rows_p    = rimwidth + halo_y + 1
        lbc_row_len_v = lbc_row_len_p
        lbc_rows_v    = lbc_rows_p - 1
        offset(1)     = 0
        offset(2)     = 0
        offset_vr = lbc_rows - rimwidth - 1
      END IF
       
      do level =1, lbc_levels

! DEPENDS ON: v_a_to_c
       call v_a_to_c (                                                  &
            v_p(ipt_p,level),                                           &
            v(ipt_v,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_v,                                              &
            lbc_rows_v,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            l_same_rot,                                                 &
            offset_vr,                                                  &
            lbc_rows,                                                   &
            halo_y,                                                     &
            phi_p_in,                                                   &
            phi_v_in )

      enddo

! East

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_v = ipt_v + lbc_row_len_v * lbc_rows_v
      
      IF (l_eg_grid) THEN
        lbc_row_len_p = rimwidth + halo_x + 1
        lbc_rows_p    = lbc_rows - 2 * rimwidth + 2
        lbc_row_len_v = lbc_row_len_p - 1
        lbc_rows_v    = lbc_rows_p - 1
        offset(1)     = 1
        offset(2)     = 1
        offset_vr  = rimwidth
      ELSE
        lbc_row_len_p = rimwidth + halo_x + 1
        lbc_rows_p    = lbc_rows - 2 * rimwidth
        lbc_row_len_v = lbc_row_len_p - 1
        lbc_rows_v    = lbc_rows_p - 1
        offset(1)     = 1
        offset(2)     = 0
        offset_vr = rimwidth
      END IF
      
      do level =1, lbc_levels

! DEPENDS ON: v_a_to_c
       call v_a_to_c (                                                  &
            v_p(ipt_p,level),                                           &
            v(ipt_v,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_v,                                              &
            lbc_rows_v,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            l_same_rot,                                                 &
            offset_vr,                                                  &
            lbc_rows,                                                   &
            halo_y,                                                     &
            phi_p_in,                                                   &
            phi_v_in )

      enddo

! South

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_v = ipt_v + lbc_row_len_v * lbc_rows_v
      
      IF (l_eg_grid) THEN
        lbc_row_len_p = lbc_row_len + 2*halo_x + 1
        lbc_rows_p    = rimwidth + halo_y + 1
        lbc_row_len_v = lbc_row_len_p - 1
        lbc_rows_v    = lbc_rows_p - 1
        offset(1)     = 1
        offset(2)     = 1
        offset_vr =  - halo_y          
      ELSE
        lbc_row_len_p = lbc_row_len + 2*halo_x
        lbc_rows_p    = rimwidth + halo_y + 1
        lbc_row_len_v = lbc_row_len_p
        lbc_rows_v    = lbc_rows_p - 1
        offset(1)     = 0
        offset(2)     = 0
        offset_vr =  - halo_y
      END IF 
       
      do level =1, lbc_levels

! DEPENDS ON: v_a_to_c
       call v_a_to_c (                                                  &
            v_p(ipt_p,level),                                           &
            v(ipt_v,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_v,                                              &
            lbc_rows_v,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            l_same_rot,                                                 &
            offset_vr,                                                  &
            lbc_rows,                                                   &
            halo_y,                                                     &
            phi_p_in,                                                   &
            phi_v_in )

      enddo

! West

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_v = ipt_v + lbc_row_len_v * lbc_rows_v
      
      IF (l_eg_grid) THEN
        lbc_row_len_p = rimwidth + halo_x + 1
        lbc_rows_p    = lbc_rows - 2 * rimwidth + 2
        lbc_row_len_v = lbc_row_len_p - 1
        lbc_rows_v    = lbc_rows_p - 1
        offset(1)     = 1
        offset(2)     = 1
        offset_vr = rimwidth
      ELSE
        lbc_row_len_p = rimwidth + halo_x + 1
        lbc_rows_p    = lbc_rows - 2 * rimwidth
        lbc_row_len_v = lbc_row_len_p - 1
        lbc_rows_v    = lbc_rows_p - 1
        offset(1)     = 0
        offset(2)     = 0
        offset_vr = rimwidth
      END IF
       
      do level =1, lbc_levels

! DEPENDS ON: v_a_to_c
       call v_a_to_c (                                                  &
            v_p(ipt_p,level),                                           &
            v(ipt_v,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_v,                                              &
            lbc_rows_v,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            l_same_rot,                                                 &
            offset_vr,                                                  &
            lbc_rows,                                                   &
            halo_y,                                                     &
            phi_p_in,                                                   &
            phi_v_in )

      enddo

      IF (lhook) CALL dr_hook('LBC_V_A_TO_C',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE lbc_v_a_to_c
