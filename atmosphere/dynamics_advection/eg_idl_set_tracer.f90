! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
SUBROUTINE eg_idl_set_tracer(                                         &
             row_length, rows, n_rows, model_levels,                  &
             halo_i, halo_j, offx, offy, earth_radius,                &
             delta_xi1, delta_xi2, base_xi1, base_xi2,                &
             xi1_p, xi2_p,                                            &
             dt_bubble, idl_bubble_width, idl_bubble_height,          &
             idl_bubble_depth, idl_bubble_xoffset,                    &
             idl_bubble_yoffset, xi3_at_theta, xi3_at_rho, tracer )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod
IMPLICIT NONE
!
! Description: create a tracer field
!  
!
! Method:
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
                     ! number of points on a processor row
  rows,                                                               &
                     ! number of rows in a processor theta field
  n_rows,                                                             &
                     ! number of rows in a processor v field
  model_levels,                                                       &
                     ! number of model levels
  halo_i,                                                             &
                     ! Size of halo in i direction.
  halo_j,                                                             &
                     ! Size of halo in j direction.
  offx,                                                               &
                     ! small i direction halo
  offy
                     ! small j direction halo

! Bubble idealised options

! Idealised options
REAL, INTENT(IN)  :: idl_bubble_width
                    ! Essentially a scaling factor (m)
REAL, INTENT(IN)  :: idl_bubble_height
                    ! Height of bubble centre (m)
REAL, INTENT(IN)  :: idl_bubble_depth
                    ! Bubble's flat radius (m)
REAL, INTENT(IN)  :: idl_bubble_xoffset
         ! Bubble x-offset (normalised units: 0.5 = domain centre)
REAL, INTENT(IN)  :: idl_bubble_yoffset
         ! Bubble y-offset (normalised units: 0.5 = domain centre)
REAL, INTENT(IN) :: earth_radius,                                     &
                       delta_xi1,                                     &
                       delta_xi2,                                     &
                        base_xi1,                                     &
                        base_xi2,                                     &
                       dt_bubble

REAL, INTENT(IN) ::                                                   &
  xi1_p(1-halo_i:row_length+halo_i),                                  &
  xi2_p(1-halo_j:rows+halo_j)

REAL, INTENT(IN) ::                                                   &
  xi3_at_theta(1-halo_i:row_length+halo_i,                            &
                 1-halo_j:rows+halo_j,0:model_levels),                &
  xi3_at_rho(1-halo_i:row_length+halo_i,                              &
               1-halo_j:rows+halo_j, model_levels)

! Output Arrays from this routine
REAL, INTENT(OUT) ::                                                  &
  tracer(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Local Variables

REAL ::                x0,                                            &
                       y0,                                            &
                       z0,                                            &
                       dist1,                                         &
                       dist2,                                         &
                       a1=3.0,                                        &
                       b1=1.0,                                        &
                       a2=1.0,                                        &
                       b2=4.0
INTEGER ::              i,                                            &
                        j,                                            &
                        k

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_IDL_SET_TRACER',zhook_in,zhook_handle)




x0 = base_xi1 + pdims%i_end*idl_bubble_xoffset*delta_xi1
y0 = base_xi2 + pdims%j_end*idl_bubble_yoffset*delta_xi2
z0 = earth_radius + idl_bubble_height


tracer(:,:,:) = 1.0

DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      dist1 = SQRT( ((xi1_p(i)-x0)/a1)**2 +                           &
                    ((xi3_at_theta(i,1,k)-z0)/b1)**2 )
      dist2 = SQRT( ((xi1_p(i)-x0)/a2)**2 +                           &
                    ((xi3_at_theta(i,1,k)-z0)/b2)**2 )

      IF ( dist1 <= idl_bubble_depth .OR. dist2 <= idl_bubble_depth)  &
             tracer(i,j,k) = tracer(i,j,k) + dt_bubble

! Add bubble temperature perturbation to the initial state

    END DO
  END DO
END DO


! End EG_IDL_Set_tracer

IF (lhook) CALL dr_hook('EG_IDL_SET_TRACER',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_set_tracer
