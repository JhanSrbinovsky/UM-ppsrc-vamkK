! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine ECMWF_mono_quasi_cubic_niv

      Subroutine ECMWF_mono_quasi_cubic_niv(                            &
     &                          Ext_Data,                               &
     &                          dim_i_in, dim_j_in, dim_k_in,           &
     &                          dim_i_out, dim_j_out, dim_k_out,        &
     &                          halo_i, halo_j,                         &
     &                          weight_lambda, weight_phi,              &
     &                          i_out, j_out, depart_level,             &
     &                          halo_i_out, halo_j_out,                 &
     &                          Data_out)

! Purpose:
!          Performs monotone quasi-cubic interpolation, as developed
!          by Mariano Hortal at ECMWF, of the input field to a
!          set of points defined by i_out, j_out, k_out, and
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
! History:
! Date     Version     Comment
! ----     -------     -------
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
     &, depart_level (dim_i_out, dim_j_out, dim_k_out)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Data_out (1-halo_i_out:dim_i_out+halo_i_out,                    &
     &            1-halo_j_out:dim_j_out+halo_j_out, dim_k_out)
                  ! data interpolated to desired locations.

! Local Variables.

      Integer                                                           &
     &  i, j, k ! Loop indices

      Real                                                              &
     &  one_sixth                                                       &
     &, val_i_minus                                                     &
     &, val_i                                                           &
     &, val_i_plus                                                      &
     &, val_i_plus2                                                     &
     &, phi_cubed                                                       &
     &, phi_sq                                                          &
     &, lambda_cubed                                                    &
     &, lambda_sq                                                       &
     &, coeff_minus                                                     &
     &, coeff_zero                                                      &
     &, coeff_plus                                                      &
     &, coeff_plus2                                                     &
     &, max_mono                                                        &
     &, min_mono

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines: None.

! ----------------------------------------------------------------------
! Section 1.  Perform cubic Lagrange Interpolation in j direction and
!             then in i direction, for central lines, linear for outside
!             lines, and ensure all pieces are monotone.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('ECMWF_MONO_QUASI_CUBIC_NIV',zhook_in,zhook_handle)
      one_sixth = 1./6.
! begin loop over levels which ends at end of subroutine

      Do k = 1, dim_k_out

          Do j = 1, dim_j_out
!dir$ suppress
            Do i = 1, dim_i_out
! interpolate in j at each index level for the four points needed for
! the i interpolation.
              phi_sq = weight_phi (i,j,k) * weight_phi (i,j,k)
              phi_cubed = weight_phi (i,j,k) * phi_sq

              coeff_plus2 = one_sixth * (phi_cubed - weight_phi(i,j,k))
              coeff_plus  = 0.5 * (phi_cubed - phi_sq                   &
     &                             - 2.*weight_phi(i,j,k))
              coeff_zero  = 0.5 * (phi_cubed - 2. * phi_sq -            &
     &                             weight_phi(i,j,k) + 2.)
              coeff_minus = one_sixth * (phi_cubed - 3. * phi_sq +      &
     &                                      2. * weight_phi(i,j,k))

              val_i_minus = weight_phi(i,j,k) *                         &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + (1.- weight_phi(i,j,k) ) *             &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k),depart_level(i,j,k))

              val_i       = coeff_plus2 *                               &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+2,depart_level(i,j,k))   &
     &                         - coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k))   &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k))     &
     &                         - coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)-1,depart_level(i,j,k))

              max_mono = max (                                          &
     &     Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k)),   &
     &     Ext_Data(i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k)) )
              min_mono = min (                                          &
     &     Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k)),   &
     &     Ext_Data(i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k)) )
              If (val_i  >   max_mono ) val_i = max_mono
              If (val_i  <   min_mono ) val_i = min_mono

              val_i_plus  = coeff_plus2 *                               &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+2,depart_level(i,j,k)) &
     &                         - coeff_plus *                           &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + coeff_zero *                           &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k))   &
     &                         - coeff_minus *                          &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-1,depart_level(i,j,k))

              max_mono = max (                                          &
     &     Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k)), &
     &     Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k)) )
              min_mono = min (                                          &
     &     Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k)), &
     &     Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k)) )
              If (val_i_plus  >   max_mono ) val_i_plus = max_mono
              If (val_i_plus  <   min_mono ) val_i_plus = min_mono

              val_i_plus2 = weight_phi(i,j,k) *                         &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                         + (1. - weight_phi(i,j,k) )*             &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k),depart_level(i,j,k))

! interpolate in i and store.
              lambda_sq = weight_lambda (i,j,k) * weight_lambda (i,j,k)
              lambda_cubed = weight_lambda (i,j,k) * lambda_sq

              data_out (i,j,k) = one_sixth * (lambda_cubed -            &
     &                                            weight_lambda(i,j,k)) &
     &                                         * val_i_plus2            &
     &                               - 0.5     * (lambda_cubed -        &
     &                                            lambda_sq -           &
     &                                         2.*weight_lambda(i,j,k)) &
     &                                         * val_i_plus             &
     &                               + 0.5     * (lambda_cubed -        &
     &                                            2. * lambda_sq -      &
     &                                            weight_lambda(i,j,k)  &
     &                                            + 2.)                 &
     &                                         * val_i                  &
     &                             - one_sixth * (lambda_cubed - 3. *   &
     &                                            lambda_sq + 2. *      &
     &                                            weight_lambda(i,j,k)) &
     &                                         * val_i_minus

              max_mono = max (val_i, val_i_plus)
              min_mono = min (val_i, val_i_plus)
              If (data_out  (i,j,k)  >   max_mono )                     &
     &            data_out  (i,j,k) = max_mono
              If (data_out  (i,j,k)  <   min_mono )                     &
     &            data_out  (i,j,k) = min_mono
          End Do
        End Do

! End loop over levels.
      End Do

! End of routine.
      IF (lhook) CALL dr_hook('ECMWF_MONO_QUASI_CUBIC_NIV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ECMWF_mono_quasi_cubic_niv
