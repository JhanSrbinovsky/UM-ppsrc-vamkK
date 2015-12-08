! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Quintic_Lagrange

      Subroutine eg_Quintic_Lagrange(                                   &
     &                          Ext_Data,                               &
     &                          dim_i_in, dim_j_in, dim_k_in,           &
     &                          dim_i_out, dim_j_out, dim_k_out,        &
     &                          halo_i, halo_j, number_of_inputs,       &
     &                          weight_lambda, weight_phi,              &
     &                          i_out, j_out, k_out,                    &
     &                          coeff_z,                                &
     &                          Data_out)

! Purpose:
!          Performs quintic Lagrange interpolation of the input field to
!          a set of points defined by i_out, j_out, k_out, and
!          weight_lambda,weight_phi,r_out.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Dynamics Advection
!
      USE um_types, ONLY: integer32
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
     &, number_of_inputs ! number of fields to interpolate

      Real                                                              &
     &  Ext_Data (1-halo_i:dim_i_in+halo_i+1,                           &
     &            1-halo_j:dim_j_in+halo_j,                             &
                                                   !data to be
     &           -1:dim_k_in+2, number_of_inputs)                       &
                                                   ! interpolated
     &, weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
                                                        ! a number
                                                      ! between 0 & 1
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)  ! a number between
                                                      ! 0 & 1

      INTEGER(KIND=integer32) ::                                        &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! point such that
     &, j_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! the desired
                                                      ! output point
     &, k_out (dim_i_out, dim_j_out, dim_k_out)       ! lies between it
                                                      ! and it+1

! Arguments with Intent IN
      Real                                                              &
     &  coeff_z(-2:3, dim_i_out, dim_j_out, dim_k_out)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Data_out (dim_i_out, dim_j_out, dim_k_out, number_of_inputs)
                      ! data interpolated to desired locations.

! Local Variables.

      Integer                                                           &
     &  i, j, k, n, index ! Loop indices

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
     &, coeff_minus2(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_minus(dim_i_out, dim_j_out, dim_k_out)                    &
     &, coeff_zero(dim_i_out, dim_j_out, dim_k_out)                     &
     &, coeff_plus(dim_i_out, dim_j_out, dim_k_out)                     &
     &, coeff_plus2(dim_i_out, dim_j_out, dim_k_out)                    &
     &, coeff_plus3(dim_i_out, dim_j_out, dim_k_out)                    &
     &, coeffl_minus2(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeffl_minus(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeffl_zero(dim_i_out, dim_j_out, dim_k_out)                    &
     &, coeffl_plus(dim_i_out, dim_j_out, dim_k_out)                    &
     &, coeffl_plus2(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeffl_plus3(dim_i_out, dim_j_out, dim_k_out)

      Real                                                              &
     &  col_data (dim_i_out, dim_j_out, -2:3) ! Horizontally interpolatd
                                              ! data ready for vertical
                                              ! interpolation.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! External Routines: None.

! NB: This routine has not been optimised unlike the other SL
!     interpolation routines.

! ----------------------------------------------------------------------
! Section 1.   Perform quintic Lagrange Interpolation in j direction and
!              then in i direction.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('EG_QUINTIC_LAGRANGE',zhook_in,zhook_handle)
      one_sixth = 1./6.
      one_twelve = 1./12.
      one_twentyfour = 1./24.
      one_hundredtwenty = 1./120.

      Do k = 1, dim_k_out

        Do j = 1, dim_j_out
        Do i = 1, dim_i_out




          phi = weight_phi (i,j,k)
          phi_sq = weight_phi (i,j,k) * weight_phi (i,j,k)
          phi_cubed = weight_phi (i,j,k) * phi_sq
          phi_fourth = phi_sq * phi_sq
          phi_fifth = phi_sq * phi_cubed

          coeff_plus3(i,j,k) = one_hundredtwenty * (phi_fifth -         &
     &                      5. * phi_cubed + 4. * phi)
          coeff_plus2(i,j,k) = - one_twentyfour * (phi_fifth -          &
     &                      phi_fourth - 7. * phi_cubed + phi_sq        &
     &                       + 6. * phi)
          coeff_plus(i,j,k)  = one_twelve * (phi_fifth                  &
     &                       - 2. * phi_fourth                          &
     &                       - 7. * phi_cubed + 8. * phi_sq +           &
     &                      12. * phi)
          coeff_zero(i,j,k)  = - one_twelve * (phi_fifth -              &
     &                      3. * phi_fourth - 5. * phi_cubed +          &
     &                      15. * phi_sq + 4. * phi -12. )
          coeff_minus(i,j,k) =  one_twentyfour * (phi_fifth -           &
     &                       4. * phi_fourth - phi_cubed +              &
     &                       16. * phi_sq - 12. * phi )
          coeff_minus2(i,j,k) = - one_hundredtwenty * (phi_fifth -      &
     &                        5. * phi_fourth + 5. * phi_cubed +        &
     &                        5. * phi_sq - 6. * phi )

          lambda = weight_lambda (i,j,k)
          lambda_sq = weight_lambda (i,j,k) * weight_lambda (i,j,k)
          lambda_cubed = weight_lambda (i,j,k) * lambda_sq
          lambda_fourth = lambda_sq * lambda_sq
          lambda_fifth = lambda_sq * lambda_cubed

          coeffl_plus3(i,j,k) = one_hundredtwenty * (lambda_fifth -     &
     &                       5. * lambda_cubed + 4. * lambda)
          coeffl_plus2(i,j,k) = - one_twentyfour * (lambda_fifth -      &
     &                       lambda_fourth - 7. * lambda_cubed +        &
     &                       lambda_sq + 6. * lambda )
          coeffl_plus(i,j,k)  = one_twelve * (lambda_fifth -            &
     &                       2. * lambda_fourth - 7. * lambda_cubed     &
     &                       +  8. * lambda_sq + 12. * lambda )
          coeffl_zero(i,j,k)  = - one_twelve * (lambda_fifth -          &
     &                       3. * lambda_fourth - 5. * lambda_cubed     &
     &                       + 15. * lambda_sq + 4. * lambda - 12. )
          coeffl_minus(i,j,k) =  one_twentyfour * (lambda_fifth -       &
     &                        4. * lambda_fourth - lambda_cubed         &
     &                        + 16. * lambda_sq  - 12. * lambda )
          coeffl_minus2(i,j,k) = -1. * one_hundredtwenty *              &
     &                        (lambda_fifth -                           &
     &                        5. * lambda_fourth + 5. * lambda_cubed    &
     &                        + 5. *lambda_sq - 6. * lambda )

        End Do ! i loop

        End Do ! j loop

      End Do

! begin loop over number of inputs which ends at end of subroutine

      Do n = 1, number_of_inputs

! begin loop over levels which ends at end of subroutine

      Do k = 1, dim_k_out

        Do index = -2, 3

            Do j = 1, dim_j_out
            Do i = 1, dim_i_out




!CDIR UNROLL=6
! interpolate in j at each index level for the six points needed for
! the i interpolation.


              val_i_minus2 = coeff_plus3(i,j,k) *                       &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)+3,k_out(i,j,k)+index,n)&
     &                         + coeff_plus2(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)+2,k_out(i,j,k)+index,n)&
     &                         + coeff_plus(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)+1,k_out(i,j,k)+index,n)&
     &                         + coeff_zero(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k),k_out(i,j,k)+index,n)  &
     &                         + coeff_minus(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)-1,k_out(i,j,k)+index,n)&
     &                         + coeff_minus2(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)-2,j_out(i,j,k)-2,k_out(i,j,k)+index,n)

              val_i_minus = coeff_plus3(i,j,k) *                        &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+3,k_out(i,j,k)+index,n)&
     &                         + coeff_plus2(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+2,k_out(i,j,k)+index,n)&
     &                         + coeff_plus(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,k_out(i,j,k)+index,n)&
     &                         + coeff_zero(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k),k_out(i,j,k)+index,n)  &
     &                         + coeff_minus(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-1,k_out(i,j,k)+index,n)&
     &                         + coeff_minus2(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-2,k_out(i,j,k)+index,n)

              val_i       = coeff_plus3(i,j,k) *                        &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+3,k_out(i,j,k)+index,n)  &
     &                         + coeff_plus2(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+2,k_out(i,j,k)+index,n)  &
     &                         + coeff_plus(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k)+index,n)  &
     &                         + coeff_zero(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+index,n)    &
     &                         + coeff_minus(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)-1,k_out(i,j,k)+index,n)  &
     &                         + coeff_minus2(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)-2,k_out(i,j,k)+index,n)

              val_i_plus  = coeff_plus3(i,j,k) *                        &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+3,k_out(i,j,k)+index,n)&
     &                         + coeff_plus2(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+2,k_out(i,j,k)+index,n)&
     &                         + coeff_plus(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k)+index,n)&
     &                         + coeff_zero(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+index,n)  &
     &                         + coeff_minus(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-1,k_out(i,j,k)+index,n)&
     &                         + coeff_minus2(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-2,k_out(i,j,k)+index,n)

              val_i_plus2 = coeff_plus3(i,j,k) *                        &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+3,k_out(i,j,k)+index,n)&
     &                         + coeff_plus2(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+2,k_out(i,j,k)+index,n)&
     &                         + coeff_plus(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,k_out(i,j,k)+index,n)&
     &                         + coeff_zero(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k),k_out(i,j,k)+index,n)  &
     &                         + coeff_minus(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-1,k_out(i,j,k)+index,n)&
     &                         + coeff_minus2(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-2,k_out(i,j,k)+index,n)


              val_i_plus3 = coeff_plus3(i,j,k) *                        &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)+3,k_out(i,j,k)+index,n)&
     &                         + coeff_plus2(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)+2,k_out(i,j,k)+index,n)&
     &                         + coeff_plus(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)+1,k_out(i,j,k)+index,n)&
     &                         + coeff_zero(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k),k_out(i,j,k)+index,n)  &
     &                         + coeff_minus(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)-1,k_out(i,j,k)+index,n)&
     &                         + coeff_minus2(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+3,j_out(i,j,k)-2,k_out(i,j,k)+index,n)


! interpolate in i and store.

              col_data (i,j,index) = coeffl_plus3(i,j,k) * val_i_plus3  &
     &                             + coeffl_plus2(i,j,k) * val_i_plus2  &
     &                             + coeffl_plus(i,j,k) * val_i_plus    &
     &                             + coeffl_zero(i,j,k) * val_i         &
     &                             + coeffl_minus(i,j,k) * val_i_minus  &
     &                             + coeffl_minus2(i,j,k) * val_i_minus2

            End Do     ! i loop

            End Do     ! j loop

        End Do    ! index loop

! ----------------------------------------------------------------------
! Section 2.   Perform quintic Lagrange Interpolation in k direction.
! ----------------------------------------------------------------------

! Interpolate data
        Do j = 1, dim_j_out
          Do i = 1, dim_i_out
            Data_out (i,j,k,n) = coeff_z(-2,i,j,k) *                    &
     &                           col_data (i,j,-2) +                    &
     &                           coeff_z(-1,i,j,k) *                    &
     &                           col_data (i,j,-1) +                    &
     &                           coeff_z(0,i,j,k) *                     &
     &                           col_data (i,j,0) +                     &
     &                           coeff_z(1,i,j,k) *                     &
     &                           col_data (i,j,1) +                     &
     &                           coeff_z(2,i,j,k) *                     &
     &                           col_data (i,j,2) +                     &
     &                           coeff_z(3,i,j,k) *                     &
     &                           col_data (i,j,3)

          End Do
        End Do

! End loop over levels.
      End Do

! End loop over number of inputs.
      End Do

! End of routine.
      IF (lhook) CALL dr_hook('EG_QUINTIC_LAGRANGE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE eg_Quintic_Lagrange
