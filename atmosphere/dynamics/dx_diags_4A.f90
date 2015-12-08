! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_dx_diags_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_dx_diags       (dx_frame_cnt, mype,                       &
                              row_length, rows, n_rows, model_levels,   &
                              offx, offy, halo_i, halo_j,               &
                              u, v, w, rho, thetav, exner, exner_star,  &
                              l_dry, m_v, m_cl, m_cf, m_r, m_gr, m_cf2)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE earth_constants_mod
USE level_heights_mod
USE eg_dxout_mod
USE horiz_grid_mod

USE eg_parameters_mod, ONLY : l_slice

IMPLICIT NONE
!
! Description: Writes a set of dynamics diagnostics in dx format 
!  
!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

       LOGICAL l_dry

       INTEGER dx_frame_cnt, mype
       INTEGER row_length, rows, n_rows, model_levels
       INTEGER offx, offy, halo_i, halo_j

       REAL                                                             &
        u(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),       &
        v(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),     &
        w(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
        thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels), &
        exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1),  &
        exner_star(1-offx:row_length+offx, 1-offy:rows+offy),           &
        rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels),      &
        m_v(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
        m_cl(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
        m_cf(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
        m_r(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),    &
        m_gr(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
        m_cf2(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

       REAL temp(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

       REAL tlid(1-offx:row_length+offx,1-offy:rows+offy)

       INTEGER i, j, k

       INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
       INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
       REAL(KIND=jprb)               :: zhook_handle



! 1.0 Start of subroutine code: perform the calculation.

       IF (lhook) CALL dr_hook('EG_DX_DIAGS',zhook_in,zhook_handle)




       IF( dx_frame_cnt == 0 ) THEN
! DEPENDS ON: eg_dx_out_grid
         CALL eg_dx_out_grid(1,xi1_u,xi2_p,r_at_u,                      &
                              row_length, rows, model_levels,           &
                              halo_i, halo_j, offx, offy, mype,         &
                              earth_radius)
! DEPENDS ON: eg_dx_out_grid
         CALL eg_dx_out_grid(2,xi1_p,xi2_v,r_at_v,                      &
                              row_length, n_rows, model_levels,         &
                              halo_i, halo_j, offx, offy, mype,         &
                              earth_radius)
! DEPENDS ON: eg_dx_out_grid
         CALL eg_dx_out_grid(3,xi1_p,xi2_p,r_theta_levels,              &
                              row_length, rows, model_levels+1,         &
                              halo_i, halo_j, offx, offy, mype,         &
                              earth_radius)
! DEPENDS ON: eg_dx_out_grid
         CALL eg_dx_out_grid(4,xi1_p,xi2_p, r_rho_levels,               &
                              row_length, rows, model_levels,           &
                              halo_i, halo_j, offx, offy, mype,         &
                              earth_radius)

       END IF

       CALL eg_dxout('Uvel',dx_frame_cnt,u, row_length, rows,           &
                      model_levels,offx, offy,1,mype)

       CALL eg_dxout('Vvel',dx_frame_cnt, v, row_length, n_rows,        &
                      model_levels,offx, offy, 2,mype)

       CALL eg_dxout('Wvel',dx_frame_cnt, w, row_length, rows,          &
                      model_levels+1,offx, offy,3,mype)

       CALL eg_dxout('theta',dx_frame_cnt, thetav, row_length, rows,    &
                      model_levels+1, offx, offy,3,mype)

       CALL eg_dxout('Exner',dx_frame_cnt, exner,                       &
                      row_length, rows, model_levels, offx, offy,4,mype)

       CALL eg_dxout('Exner_star',dx_frame_cnt, exner_star,             &
                      row_length, rows, 1, offx, offy,3,mype)

       CALL eg_dxout('rho',dx_frame_cnt, rho, row_length, rows,         &
                      model_levels,offx, offy,4,mype)

       If (.Not. l_dry) Then
         CALL EG_DXout('m_v',dx_frame_cnt, m_v, row_length, rows,         &
                        model_levels+1, offx, offy,3,mype)

         CALL EG_DXout('m_cl',dx_frame_cnt, m_cl, row_length, rows,       &
                        model_levels+1, offx, offy,3,mype)

         CALL EG_DXout('m_cf',dx_frame_cnt, m_cf, row_length, rows,       &
                        model_levels+1, offx, offy,3,mype)

         CALL EG_DXout('m_r',dx_frame_cnt, m_r, row_length, rows,        &
                        model_levels+1, offx, offy,3,mype)

         CALL EG_DXout('m_gr',dx_frame_cnt, m_gr, row_length, rows,       &
                        model_levels+1, offx, offy,3,mype)

         CALL EG_DXout('m_cf2',dx_frame_cnt, m_cf2, row_length, rows,       &
                        model_levels+1, offx, offy,3,mype)

       Endif

! Calculate temperature
       DO k = 1, model_levels
         DO j = 1, rows
           DO i = 1, row_length
             temp(i,j,k) = 0.5*exner(i,j,k)                             &
                                *(thetav(i,j,k)+thetav(i,j,k-1))
           END DO
         END DO
       END DO

!      Calculate lid temperature
       k = model_levels
       DO j = 1, rows
         DO i = 1, row_length
           tlid(i,j) = 0.5*(exner(i,j,k) + exner(i,j,k+1))              &
                        *thetav(i,j,k)
         END DO
       END DO


       CALL eg_dxout('Temp',dx_frame_cnt,  temp,                        &
                     row_length, rows, model_levels, offx, offy,4,mype)

       CALL eg_dxout('Ttop',dx_frame_cnt,  tlid,                        &
                     row_length, rows, 1, offx, offy,3,mype)


       dx_frame_cnt = dx_frame_cnt + 1


       IF (lhook) CALL dr_hook('EG_DX_DIAGS',zhook_out,zhook_handle)
       RETURN
       END SUBROUTINE eg_dx_diags
       END MODULE eg_dx_diags_mod
