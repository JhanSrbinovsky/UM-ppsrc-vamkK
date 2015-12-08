! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE stph_diag_mod

IMPLICIT NONE

!  ---------------------------------------------------------------------
!  Sructure containing stochastic physics diagnostics.
!  This permits easier addition diagnostics without additional
!  passing of arguments through code tree.
!  It also does not require the addition
!  of extra subroutine arguments when adding a new diagnostic.
!  Module is used atm_step and stph_skeb2
!  ---------------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Stochastic Physics

 TYPE strstphdiag

! Need to create a flag and a pointer

 LOGICAL ::  l_skeb2_u
 LOGICAL ::  l_skeb2_v
 LOGICAL ::  l_skeb2_u_incr
 LOGICAL ::  l_skeb2_v_incr
 LOGICAL ::  l_skeb2_u_rot
 LOGICAL ::  l_skeb2_v_rot
 LOGICAL ::  l_skeb2_u_div
 LOGICAL ::  l_skeb2_v_div
 LOGICAL ::  l_skeb2_disp_smag
 LOGICAL ::  l_skeb2_disp_conv
 LOGICAL ::  l_skeb2_disp_skeb1
 LOGICAL ::  l_skeb2_smodfield
 LOGICAL ::  l_skeb2_streamfunction
 LOGICAL ::  l_skeb2_random_pattern
 LOGICAL ::  l_skeb2_ke_psif
 LOGICAL ::  l_skeb2_ke_sdisp
 LOGICAL ::  l_skeb2_ke_cdisp
 LOGICAL ::  l_skeb2_ke_kdisp
 LOGICAL ::  l_skeb2_ke_m_psif
 LOGICAL ::  l_skeb2_ke_prewindincr
 LOGICAL ::  l_skeb2_ke_windincr
 LOGICAL ::  l_skeb2_ke_postwindincr

 REAL, POINTER :: skeb2_u(:, :, :)
!         1     U Wind after SKEB2
 REAL, POINTER :: skeb2_v(:, :, :)
!         2     V Wind after SKEB2
 REAL, POINTER :: skeb2_u_incr(:, :, :)
!         3     SKEB2 Full u increment
 REAL, POINTER :: skeb2_v_incr(:, :, :)
!         4     SKEB2 Full v increment
 REAL, POINTER :: skeb2_u_rot(:, :, :)
!         5     SKEB2: rotational u incr
 REAL, POINTER :: skeb2_v_rot(:, :, :)
!         6     SKEB2: rotational v incr
 REAL, POINTER :: skeb2_u_div(:, :, :)
!         7     SKEB2: divergent u incr
 REAL, POINTER :: skeb2_v_div(:, :, :)
!         8     SKEB2: divergent v incr
 REAL, POINTER :: skeb2_disp_smag(:, :, :)
!         9     SKEB2: dissipation field from smagorinsky code
 REAL, POINTER :: skeb2_disp_conv(:, :, :)
!         10     SKEB2: dissipation field from convection
 REAL, POINTER :: skeb2_disp_skeb1(:, :, :)
!         11     SKEB2: dissipation field from SKEB1 code
 REAL, POINTER :: skeb2_smodfield(:, :, :)
!         12     SKEB2: Smoothed Modulating field
 REAL, POINTER :: skeb2_streamfunction(:, :, :)
!         13     SKEB2: final stream function field
 REAL, POINTER :: skeb2_random_pattern(:, :, :)
!         14     SKEB2: intial random pattern
 REAL, POINTER :: skeb2_ke_psif(:, :)
!         15     SKEB2: Vert Integ. KE of initial SF forcing
 REAL, POINTER :: skeb2_ke_sdisp(:, :)
!         16     SKEB2: Vert Integ. KE of numerical diss
 REAL, POINTER :: skeb2_ke_cdisp(:, :)
!         17     SKEB2: Vert Integ. KE of convection diss
 REAL, POINTER :: skeb2_ke_kdisp(:, :)
!         18     SKEB2: Vert Integ. KE of SKEB1 dissipation
 REAL, POINTER :: skeb2_ke_m_psif(:, :)
!         19     SKEB2: Vert Integ. KE of smoothed dissipation
 REAL, POINTER :: skeb2_ke_prewindincr(:, :)
!         20     SKEB2: Vert Integ. KE of total wind incr before SKEB2
 REAL, POINTER :: skeb2_ke_windincr(:, :)
!         21     SKEB2: Vert Integ. KE of wind incr from SKEB2
 REAL, POINTER :: skeb2_ke_postwindincr(:, :)
!         22     SKEB2: Vert Integ. KE of total wind incr after SKEB2
 END TYPE strstphdiag

END MODULE stph_diag_mod
