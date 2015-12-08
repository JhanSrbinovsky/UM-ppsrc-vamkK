! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_calc_p_star_mod
IMPLICIT NONE
CONTAINS

SUBROUTINE eg_calc_p_star(                                            &
                 model_levels, row_length, rows,   p,                 &
                 thetav,  m_v, m_cl, m_cf, m_r, m_gr,                 &
                 g,  p_star, w_surf, w_lid )

USE um_parvars,          ONLY : offx, offy, halo_i, halo_j
USE level_heights_mod,   ONLY : xi3_at_theta=>r_theta_levels,            &
                                xi3_at_rho=>r_rho_levels
USE metric_terms_mod,    ONLY :  h3_p_eta
USE parkind1,            ONLY : jpim, jprb       !DrHook
USE yomhook,             ONLY : lhook, dr_hook   !DrHook
USE atmos_constants_mod
USE ref_pro_mod,         ONLY : p_ref_pro => exner_ref_pro
USE atm_fields_bounds_mod
USE Field_Types
USE proc_info_mod,       ONLY : model_domain
USE domain_params

IMPLICIT NONE
!
! Description:
!  
!         Calculates surface pressure.
!
! Method: Adaptation of New Dynamics code.
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


INTEGER, INTENT(IN) ::                                                &
  row_length,                                                         &
                   ! number of point on a row.
  rows,                                                               &
                   ! number of rows.
  model_levels

REAL, INTENT(IN) ::                                                   &
  thetav(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
   g(1-offx:row_length+offx, 1-offy:rows+offy,0:model_levels),         &
  w_surf(1-offx:row_length+offx,1-offy:rows+offy),                    &
  w_lid(1-offx:row_length+offx,1-offy:rows+offy)

REAL p(1-offx:row_length+offx, 1-offy:rows+offy, model_levels+1)

REAL,                                                                 &
 DIMENSION(1-offx:row_length+offx,1-offy:rows+offy,                   &
           0:model_levels),                                           &
 INTENT(IN) :: m_v, m_cl, m_cf, m_r, m_gr

! Arguments with Intent OUT. ie: Output variables.

REAL, INTENT(INOUT) ::                                                &
  p_star (1-offx:row_length+offx, 1-offy:rows+offy)


! Local Variables.

INTEGER  :: i, j,k       ! Loop indices
REAL     :: a, b

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_CALC_P_STAR',zhook_in,zhook_handle)


k = model_levels
DO j = pdims%j_start, pdims%j_end
   DO i = pdims%i_start, pdims%i_end

      b = 1.0 + m_v(i,j,0) + m_cl(i,j,0) + m_cf(i,j,0)                &
              + m_r(i,j,0) + m_gr(i,j,0)
      a = h3_p_eta(i,j,0)                                             &
               *(xi3_at_rho(i,j,1)-xi3_at_theta(i,j,0))               &
               *b/( cp*thetav(i,j,0) )

      p_star(i,j) = p(i,j,1) + a*( g(i,j,0) - w_surf(i,j) )


      b = 1.0 + m_v(i,j,k) + m_cl(i,j,k) + m_cf(i,j,k)                &
              + m_r(i,j,k) + m_gr(i,j,k)
      a = -2.0*h3_p_eta(i,j,k)                                        &
               *(xi3_at_rho(i,j,k)-xi3_at_theta(i,j,k))               &
               *b/( cp*thetav(i,j,k) )

      p(i,j,k+1) = p(i,j,k) - a*( g(i,j,k) - w_lid(i,j) )

   END DO
END DO

! p_star is interpolated, so we do need the halo fill (apart from 
! the fact that all fields with halos should have them filled or have no
! halos to start with!)
 
! DEPENDS ON: swap_bounds
    CALL swap_bounds(p_star,row_length,rows,1,offx,offy,                &
                       fld_type_p,.false.)
                       
IF (lhook) CALL dr_hook('EG_CALC_P_STAR',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_calc_p_star
END MODULE eg_calc_p_star_mod
