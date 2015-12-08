! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************
!
! Subroutine Bottom_w_Calc.

      Subroutine Bottom_w_Calc(                                         &
     &                         r_theta_levels, u, v,                    &
     &                         w, w_adv, sec_theta_latitude,            &
     &                         delta_lambda, delta_phi,                 &
     &                         lambda_p, phi_p, lambda_u, phi_v,        &
     &                         rows, n_rows, row_length, model_levels,  &
     &                         model_domain, L_regular,                 &
     &                         off_x, off_y, halo_i, halo_j,            &
     &                         proc_row_group, at_extremity,            &
     &                         g_row_length)

! Purpose:
!          Calculates etadot from given u, v and w fields.
!
! Method:
!          Is described in ;
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE global_2d_sums_mod, ONLY : global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Integer                                                           &
              ! model dimensions
     &  row_length                                                      &
                         ! number of points on a row
     &, g_row_length                                                    &
                         ! global number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels
     &, off_x                                                           &
                   ! Size of small halo in i direction
     &, off_y                                                           &
                   ! Size of small halo in j direction
     &, halo_i                                                          &
                ! Size of halo in i direction
     &, halo_j                                                          &
                ! Size of halo in j direction
     &, proc_row_group ! Group id for processors on the same row

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south,east or west of the processor grid
     &, L_regular

      Real                                                              &
           ! horizontal co-ordinate information
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_p    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : rows + halo_j )                           &
     &, phi_v    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : n_rows+halo_j )

      Real                                                              &
           ! trigonometric functions
     &  sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)

! Primary Arrays

      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      0:model_levels)                                             &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          0:model_levels)

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k                                                         &
                          ! Loop indices
     &, j0                                                              &
     &, j1

      Real                                                              &
     &  Recip_delta_lambda                                              &
     &, Recip_delta_phi                                                 &
     &, weight1                                                         &
     &, weight2                                                         &
     &, sum_row

      Real :: sum_tmp(1)

      Integer info

! arrays

      Real                                                              &
     &  u_on_theta(1-off_x:row_length+off_x, rows)                                  &
     &, v_on_theta(row_length, 1-off_y:n_rows+off_y)                                &
     &, term(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! No External Routines:

! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Interpolate u and v to w points.
!              Not done at top level.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('BOTTOM_W_CALC',zhook_in,zhook_handle)
      j0 = 1
      j1 = rows
      IF (.NOT. l_vatpoles) THEN
      If (model_domain  /=  mt_bi_cyclic_LAM) then
        If (at_extremity(PSouth)) j0 = 2
        If (at_extremity(PNorth)) j1 = rows-1
      Endif
      END IF ! vatpoles
      recip_delta_lambda = 1./ delta_lambda
      recip_delta_phi = 1./ delta_phi

!  Level 1 u assumed to be same as surface

      Do j = 1, rows
        Do i = 0, row_length
          u_on_theta (i,j) = u(i,j,1)
        End Do
      End Do

!  Level 1 v assumed to be same as surface
      Do j = 0, n_rows
        Do i = 1, row_length
          v_on_theta (i,j) =  v(i,j,1)
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 2.   Calculate u.grad(r)
! ----------------------------------------------------------------------
      k = 0
      
      IF (L_regular) then

! Calculate u dr/d lambda
      IF (l_vatpoles) THEN
      DO j = j0, j1
        DO i = 1, row_length+1
          term(i,j) = u_on_theta(i,j) *                                 &
                     (r_theta_levels(i,j,k) -                           &
                      r_theta_levels(i-1,j,k) ) * recip_delta_lambda    &
                      * 2.0 / (r_theta_levels(i,j,k) +                  &
                               r_theta_levels(i-1,j,k) )
        END DO
      END DO
      ELSE
      DO j = j0, j1
        DO i = 0, row_length
          term(i,j) = u_on_theta(i,j) *                                 &
                     (r_theta_levels(i+1,j,k) -                         &
                      r_theta_levels(i,j,k) ) * recip_delta_lambda      &
                      * 2.0 / (r_theta_levels(i+1,j,k) +                &
                               r_theta_levels(i,j,k) )
        END DO
      END DO
      END IF ! vatpoles

! average quantity and store in w(i,j,0)
      IF (l_vatpoles) THEN
      DO j = j0, j1
        DO i = 1, row_length
          w(i,j,k) = (term(i+1,j) + term(i,j) ) * 0.5 *                 &
                              sec_theta_latitude(i,j)
        END DO
      END DO
      ELSE
      DO j = j0, j1
        DO i = 1, row_length
          w(i,j,k) = (term(i,j) + term(i-1,j) ) * 0.5 *                 &
                              sec_theta_latitude(i,j)
        END DO
      END DO
      END IF ! vatpoles

! Calculate v dr/d phi
      IF (l_vatpoles) THEN
      Do j = j0, j1+1
        Do i = 1, row_length
          term(i,j) = v_on_theta(i,j) *                                 &
                     (r_theta_levels(i,j,k) -                           &
                      r_theta_levels(i,j-1,k) ) * recip_delta_phi       &
                        * 2.0 / (r_theta_levels(i,j,k) +                &
                                 r_theta_levels(i,j-1,k) )
        End Do
      End Do
      ELSE
      Do j = j0-1, n_rows
        Do i = 1, row_length
          term(i,j) = v_on_theta(i,j) *                                 &
                     (r_theta_levels(i,j+1,k) -                         &
                      r_theta_levels(i,j,k) ) * recip_delta_phi         &
                        * 2.0 / (r_theta_levels(i,j+1,k) +              &
                                 r_theta_levels(i,j,k) )
        End Do
      End Do
      END IF ! vatpoles

! average quantity and add on to u term
      IF (l_vatpoles) THEN
      DO j = j0, j1
        DO i = 1, row_length
          w(i,j,k) = (term(i,j+1) + term(i,j)) * 0.5 + w(i,j,k)
          w_adv(i,j,k) = w(i,j,k)
        END DO
      END DO
      ELSE
      DO j = j0, j1
        DO i = 1, row_length
          w_adv(i,j,k) = w(i,j,k)
        END DO
      END DO
      END IF ! vatpoles

      ELSE  ! variable resolution

! Calculate u dr/d lambda
      Do j = j0, j1
        Do i = 0, row_length
        weight1 =  lambda_p(i+1)-lambda_u(i)
        weight2 =  lambda_u(i) - lambda_p(i)
        term(i,j) = u_on_theta(i,j)                                     &
     &            * (r_theta_levels(i+1,j,k) - r_theta_levels(i,j,k))   &
     &            / ( weight1 * r_theta_levels(i,j,k) +                 &
     &                weight2 * r_theta_levels(i+1,j,k) )
        End Do
      End Do

! average quantity and store in w(i,j,0)
      Do j = j0, j1
        Do i = 1, row_length
        recip_delta_lambda = 1.0/(lambda_u(i)-lambda_u(i-1))
        weight1 = (lambda_u(i)-lambda_p(i)) * recip_delta_lambda
        weight2 = 1.0-weight1
        w(i,j,k) = ( weight1 * term(i-1,j) + weight2 * term(i,j) )      &
     &                   *  sec_theta_latitude(i,j)
        End Do
      End Do

! Calculate v dr/d phi
      Do j = j0-1, n_rows
        Do i = 1, row_length
          weight1 =  phi_p(i,j+1) - phi_v(i,j)
          weight2 =  phi_v(i,j) - phi_p(i,j)
          term(i,j) = v_on_theta(i,j) *                                 &
     &               (r_theta_levels(i,j+1,k) - r_theta_levels(i,j,k))/ &
     &                              ( weight1 * r_theta_levels(i,j,k) + &
     &                                weight2 * r_theta_levels(i,j+1,k))
        End Do
      End Do

! average quantity and add on to u term
      Do j = j0, j1
        Do i = 1, row_length
          recip_delta_phi  =  1.0 / (phi_v(i,j) - phi_v(i,j-1))
          weight1 =  (phi_v(i,j) - phi_p(i,j)) * recip_delta_phi
          w(i,j,k) = ( weight1 * term(i,j-1) +                          &
     &                (1.0 - weight1) * term(i,j) ) +  w(i,j,k)
          w_adv(i,j,k) = w(i,j,k)
        End Do
      End Do

      End If   ! L_regular

IF (.NOT. l_vatpoles) THEN
! northern and southern boundary updates.
! save values for summing.

      If (model_domain  ==  mt_global) Then    ! global model domain

        If(at_extremity(PSouth))then

          CALL global_2d_sums(term(1:row_length,1:1), row_length,       &
                              1, 0, 0, 1, sum_tmp, proc_row_group)

          sum_row = sum_tmp(1) / g_row_length
          Do i = 1, row_length
            w(i,1,k) = sum_row
          End Do
        End If

        If(at_extremity(PNorth))then
          CALL global_2d_sums(term(1:row_length,n_rows:n_rows),         &
                              row_length, 1, 0, 0, 1,  sum_tmp,         &
                              proc_row_group)

          sum_row = sum_tmp(1) / g_row_length
          Do i = 1, row_length
            w(i,rows,k) = sum_row
          End Do
        End If

      Else     ! limited area model.
!  Assume values outside boundary are zero.

       If(at_extremity(PSouth) .and.                                    &
     &                 model_domain  /=  mt_bi_cyclic_LAM) then
          Do i= 1,row_length
            w(i,1,k) = 0.
          End Do
        End If

       If(at_extremity(PNorth) .and.                                    &
     &                 model_domain  /=  mt_bi_cyclic_LAM) then
          j = rows
          Do i=1,row_length
            w(i,j,k) = 0.
          End Do
        End If

      End If   !( model_domain  ==  mt_global )
END IF ! vatpoles

! End of routine  Bottom_w_Calc.
      IF (lhook) CALL dr_hook('BOTTOM_W_CALC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Bottom_w_Calc
      
