! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Quintic_Lagrange_niv

      Subroutine Quintic_Lagrange_niv(                                  &
     &                          Ext_Data,                               &
     &                          dim_i_in, dim_j_in, dim_k_in,           &
     &                          dim_i_out, dim_j_out, dim_k_out,        &
     &                          halo_i, halo_j,                         &
     &                          weight_lambda, weight_phi,              &
     &                          i_out, j_out, depart_level,             &
     &                          halo_i_out, halo_j_out,                 &
     &                          Data_out)

! Purpose:
!          Performs quintic Lagrange interpolation of the input field to
!          a set of points defined by i_out, j_out, depart_level, and
!          weight_lambda,weight_phi,r_out, but using no interpolation
!          in the vertical.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.


      Integer                                                           &
     &  dim_i_in                                                        &
                    ! Dimension of Data_in in i direction.
     &, dim_j_in                                                        &
                    ! Dimension of Data_in in j direction.
     &, dim_k_in                                                        &
                    ! Dimension of Data_in in k direction.
     &, dim_i_out                                                       &
                    ! Dimension of Data_out in i direction.
     &, dim_j_out                                                       &
                    ! Dimension of Data_out in j direction.
     &, dim_k_out                                                       &
                    ! Dimension of Data_out in k direction.
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, halo_i_out                                                      &
                    ! Size of data out halo in i direction.
     &, halo_j_out  ! Size of data out halo in j direction.

      Real                                                              &
     &  Ext_Data (1-halo_i:dim_i_in+halo_i,                             &
     &            1-halo_j:dim_j_in+halo_j,-1:dim_k_in+2)               &
                                                          !data to be
                                                          ! interpolated
     &, weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
                                                        ! a number
                                                      ! between 0 & 1
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)  ! a number between
                                                      ! 0 & 1

      Integer                                                           &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! point such that
     &, j_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! the desired
                                                      ! output point
             ! lies between it and it+1
     &, depart_level (dim_i_out, dim_j_out, dim_k_out)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Data_out (1-halo_i_out:dim_i_out+halo_i_out,                    &
     &            1-halo_j_out:dim_j_out+halo_j_out, dim_k_out)
                  ! data interpolated to desired locations.

! Local Variables.

      Integer                                                           &
     &  i, j, k  ! Loop indices

      Real                                                              &
     &  one_sixth                                                       &
     &, one_twelve                                                      &
     &, one_twentyfour                                                  &
     &, one_hundredtwenty                                               &
     &, val_i_minus                                                     &
     &, val_i_minus2                                                    &
     &, val_i                                                           &
     &, val_i_plus                                                      &
     &, val_i_plus2                                                     &
     &, val_i_plus3                                                     &
     &, phi                                                             &
     &, phi_cubed                                                       &
     &, phi_sq                                                          &
     &, phi_fourth                                                      &
     &, phi_fifth                                                       &
     &, lambda_cubed                                                    &
     &, lambda_sq                                                       &
     &, lambda                                                          &
     &, lambda_fourth                                                   &
     &, lambda_fifth                                                    &
     &, coeff_minus2                                                    &
     &, coeff_minus                                                     &
     &, coeff_zero                                                      &
     &, coeff_plus                                                      &
     &, coeff_plus2                                                     &
     &, coeff_plus3                                                     &
     &, coeffl_minus2                                                   &
     &, coeffl_minus                                                    &
     &, coeffl_zero                                                     &
     &, coeffl_plus                                                     &
     &, coeffl_plus2                                                    &
     &, coeffl_plus3

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines: None.

! ----------------------------------------------------------------------
! Section 1.   Perform quintic Lagrange Interpolation in j direction and
!              then in i direction.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('QUINTIC_LAGRANGE_NIV',zhook_in,zhook_handle)
      one_sixth = 1./6.
      one_twelve = 1./12.
      one_twentyfour = 1./24.
      one_hundredtwenty = 1./120.

! begin loop over levels which ends at end of subroutine

      Do k = 1, dim_k_out

        Do j = 1, dim_j_out
          Do i = 1, dim_i_out
! interpolate in j at each index level for the six points needed for
! the i interpolation.
              phi = weight_phi (i,j,k)
              phi_sq = weight_phi (i,j,k) * weight_phi (i,j,k)
              phi_cubed = weight_phi (i,j,k) * phi_sq
              phi_fourth = phi_sq * phi_sq
              phi_fifth = phi_sq * phi_cubed

              coeff_plus3 = one_hundredtwenty * (phi_fifth -            &
     &                      5. * phi_cubed + 4. * phi)
              coeff_plus2 = - one_twentyfour * (phi_fifth -             &
     &                      phi_fourth - 7. * phi_cubed + phi_sq        &
     &                       + 6. * phi)
              coeff_plus  = one_twelve * (phi_fifth - 2. * phi_fourth   &
     &                       - 7. * phi_cubed + 8. * phi_sq +           &
     &                      12. * phi)
              coeff_zero  = - one_twelve * (phi_fifth -                 &
     &                      3. * phi_fourth - 5. * phi_cubed +          &
     &                      15. * phi_sq + 4. * phi -12. )
              coeff_minus =  one_twentyfour * (phi_fifth -              &
     &                       4. * phi_fourth - phi_cubed +              &
     &                       16. * phi_sq - 12. * phi )
              coeff_minus2 = - one_hundredtwenty * (phi_fifth -         &
     &                        5. * phi_fourth + 5. * phi_cubed +        &
     &                        5. * phi_sq - 6. * phi )

              val_i_minus2 = coeff_plus3 *                              &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)+3,depart_level(i,j,k)) &
     &                         + coeff_plus2 *                          &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         + coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k),depart_level(i,j,k))   &
     &                         + coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)-1,depart_level(i,j,k)) &
     &                         + coeff_minus2 *                         &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)-2,depart_level(i,j,k))

              val_i_minus = coeff_plus3 *                               &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+3,depart_level(i,j,k)) &
     &                         + coeff_plus2 *                          &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         + coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k),depart_level(i,j,k))   &
     &                         + coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-1,depart_level(i,j,k)) &
     &                         + coeff_minus2 *                         &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-2,depart_level(i,j,k))

              val_i       = coeff_plus3 *                               &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+3,depart_level(i,j,k))   &
     &                         + coeff_plus2 *                          &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+2,depart_level(i,j,k))   &
     &                         + coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k))   &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k))     &
     &                         + coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)-1,depart_level(i,j,k))   &
     &                         + coeff_minus2 *                         &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)-2,depart_level(i,j,k))

              val_i_plus  = coeff_plus3 *                               &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+3,depart_level(i,j,k)) &
     &                         + coeff_plus2 *                          &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         + coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k))   &
     &                         + coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-1,depart_level(i,j,k)) &
     &                         + coeff_minus2 *                         &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-2,depart_level(i,j,k))

              val_i_plus2 = coeff_plus3 *                               &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+3,depart_level(i,j,k)) &
     &                         + coeff_plus2 *                          &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         + coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k),depart_level(i,j,k))   &
     &                         + coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-1,depart_level(i,j,k)) &
     &                         + coeff_minus2 *                         &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-2,depart_level(i,j,k))


              val_i_plus3 = coeff_plus3 *                               &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)+3,depart_level(i,j,k)) &
     &                         + coeff_plus2 *                          &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         + coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k),depart_level(i,j,k))   &
     &                         + coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)-1,depart_level(i,j,k)) &
     &                         + coeff_minus2 *                         &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)-2,depart_level(i,j,k))


! interpolate in i and store.

              lambda = weight_lambda (i,j,k)
              lambda_sq = weight_lambda (i,j,k) * weight_lambda (i,j,k)
              lambda_cubed = weight_lambda (i,j,k) * lambda_sq
              lambda_fourth = lambda_sq * lambda_sq
              lambda_fifth = lambda_sq * lambda_cubed

              coeffl_plus3 = one_hundredtwenty * (lambda_fifth -        &
     &                       5. * lambda_cubed + 4. * lambda)
              coeffl_plus2 = - one_twentyfour * (lambda_fifth -         &
     &                       lambda_fourth - 7. * lambda_cubed +        &
     &                       lambda_sq + 6. * lambda )
              coeffl_plus  = one_twelve * (lambda_fifth -               &
     &                       2. * lambda_fourth - 7. * lambda_cubed     &
     &                       +  8. * lambda_sq + 12. * lambda )
              coeffl_zero  = - one_twelve * (lambda_fifth -             &
     &                       3. * lambda_fourth - 5. * lambda_cubed     &
     &                       + 15. * lambda_sq + 4. * lambda - 12. )
              coeffl_minus =  one_twentyfour * (lambda_fifth -          &
     &                        4. * lambda_fourth - lambda_cubed         &
     &                        + 16. * lambda_sq  - 12. * lambda )
              coeffl_minus2 = -1. * one_hundredtwenty * (lambda_fifth - &
     &                        5. * lambda_fourth + 5. * lambda_cubed    &
     &                        + 5. *lambda_sq - 6. * lambda )

              data_out (i,j,k) = coeffl_plus3 * val_i_plus3             &
     &                             + coeffl_plus2 * val_i_plus2         &
     &                             + coeffl_plus * val_i_plus           &
     &                             + coeffl_zero * val_i                &
     &                             + coeffl_minus * val_i_minus         &
     &                             + coeffl_minus2 * val_i_minus2

          End Do
        End Do

! End loop over levels.
      End Do

! End of routine.
      IF (lhook) CALL dr_hook('QUINTIC_LAGRANGE_NIV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Quintic_Lagrange_niv

