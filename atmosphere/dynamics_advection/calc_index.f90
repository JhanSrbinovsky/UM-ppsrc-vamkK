! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Calc_Index
!

      Subroutine Calc_Index(                                            &
     &                      dim_i_out, dim_j_out, dim_k_out,            &
     &                      delta_lambda_in, delta_phi_in,              &
     &                      base_lambda, base_phi,                      &
     &                      glambda_p, phi_p, grecip_dlamp, recip_dphip,&
     &                      recip_dlam, recip_dphi, max_look,           &
     &                      look_lam, look_phi, halo_lam, halo_phi,     &
     &                      L_regular, lambda_out, phi_out,             &
     &                      halo_i, halo_j,                             &
     &                      g_row_length, row_length, rows, datastart,  &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi)

! Purpose:
!          Calculates indices of lambda_out and phi_out.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 90
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer, Intent(IN) ::                                            &
     &  dim_i_out                                                       &
                    ! Dimension of Data_out in i direction.
     &, dim_j_out                                                       &
                    ! Dimension of Data_out in j direction.
     &, dim_k_out                                                       &
                    ! Dimension of Data_out in k direction.
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, max_look                                                        &
                        ! max size of look-up arrays for searches
     &, row_length                                                      &
                        ! row_length for dynamic arrays
     &, rows                                                            &
                        ! rows for p-field dynamic arrays
     &, g_row_length                                                    &
                        ! global number of points on a row
     &, datastart(3)
                     ! First gridpoints held by this processor.

      Logical, Intent(IN) ::                                            &
     &  L_regular      ! False if variable resolution

      Real, Intent(IN) ::                                               &
     &  delta_lambda_in                                                 &
                         ! holds spacing between points in the i
                         ! direction for the input data field.
     &, delta_phi_in                                                    &
                         ! holds spacing between points in the j
                         ! direction for the input data field.
     &, recip_dlam                                                      &
     &, recip_dphi                                                      &
     &, base_lambda                                                     &
     &, base_phi

! look-up table halos
       Integer, Intent(IN) ::                                           &
     &  halo_lam                                                        &
     &, halo_phi

!VarRes horizontal co-ordinate look-up table
       Integer, Intent(IN) ::                                           &
     &  look_lam(1-halo_lam:max_look-halo_lam)                          &
     &, look_phi(1-halo_phi:max_look-halo_phi)

!  VarRes horizontal co-ordinate spacing etc
       Real, Intent(IN) ::                                              &
     &  glambda_p(1-halo_i : g_row_length + halo_i)                     &
     &, phi_p( 1-halo_i : row_length + halo_i                           &
     &,        1-halo_j : rows + halo_j )                               &
     &, grecip_dlamp(1-halo_i : g_row_length + halo_i)                  &
     &, recip_dphip( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows + halo_j )

      Real, Intent(IN) ::                                               &
     &  lambda_out (dim_i_out, dim_j_out, dim_k_out)                    &
                                                      ! Lambda
                                                      ! co-ordinate of
                                                      ! output data on
                                                      ! input.
     &, phi_out (dim_i_out, dim_j_out, dim_k_out)
                                                      ! Phi Co-ordinate
                                                      ! of output data
                                                      ! on input.

! Arguments with Intent OUT. ie: Output variables.

      Integer, Intent(Out) ::                                           &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
     &, j_out (dim_i_out, dim_j_out, dim_k_out)

      Real, Intent(Out) ::                                              &
     &  weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
     &, weight_phi    (dim_i_out, dim_j_out, dim_k_out)

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k                                                         &
                   ! Loop indices
     &, dsm1y                                                           &
     &, index

      Real                                                              &
     &  Recip_delta_lambda_in                                           &
     &, Recip_delta_phi_in                                              &
     &, temp

      Real :: out_p
      Real :: out_l

! Statement function for optimising floor
! Can only be used when the result is cast into an integer
      REAL :: real_temp_var
      REAL :: my_floor
      my_floor(real_temp_var) = real_temp_var - 0.5 + sign(0.5,real_temp_var)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
!  Section 1.   Find i and j point indices
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CALC_INDEX',zhook_in,zhook_handle)
      If ( L_regular ) Then

        Recip_delta_lambda_in = 1./delta_lambda_in
        Recip_delta_phi_in = 1./delta_phi_in
        dsm1y = datastart(2)-1

!cdir collapse
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i,         &
!$OMP& out_l,out_p) SHARED(dim_k_out,dim_j_out,dim_i_out,               &
!$OMP& lambda_out,recip_delta_lambda_in,i_out,weight_lambda,            &
!$OMP& phi_out, base_phi ,recip_delta_phi_in, dsm1y,j_out,weight_phi) 
        Do k = 1, dim_k_out
          Do j = 1, dim_j_out
            Do i = 1, dim_i_out
              out_l = lambda_out(i,j,k) * recip_delta_lambda_in
              i_out(i,j,k) = my_floor(1.0 + out_l)
              weight_lambda(i,j,k) = out_l - ( i_out(i,j,k) - 1.0 )
              out_p = (phi_out(i,j,k) - base_phi) * recip_delta_phi_in
              j_out(i,j,k) = my_floor(1.0 + out_p)
              weight_phi(i,j,k) = out_p - ( j_out(i,j,k) - 1.0)
              j_out(i,j,k) = j_out(i,j,k) - dsm1y
            End Do
          End Do
        End Do
!$OMP END PARALLEL DO

      Else  ! variable resolution
        dsm1y = datastart(2)-1

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(k,j,i,index)   &
!$OMP& SHARED(dim_k_out,dim_j_out,dim_i_out,                            &
!$OMP& lambda_out,recip_dlam,i_out,look_lam,weight_lambda,base_lambda,  &
!$OMP& glambda_p,phi_p,grecip_dlamp,weight_phi,recip_dphip,             &
!$OMP& phi_out,base_phi,recip_dphi,dsm1y,j_out,look_phi) 
        Do k = 1, dim_k_out
          Do j = 1, dim_j_out
            Do i = 1, dim_i_out
            index = 1 + my_floor(lambda_out(i,j,k) * recip_dlam)
            i_out(i,j,k) = look_lam(index)
            If ( lambda_out(i,j,k) + base_lambda                        &
     &                >= glambda_p(i_out(i,j,k)+1) )                    &
     &                i_out(i,j,k) = i_out(i,j,k) + 1
            weight_lambda(i,j,k) = ( lambda_out(i,j,k) + base_lambda -  &
     &                                     glambda_p(i_out(i,j,k)) ) *  &
     &                                  grecip_dlamp(i_out(i,j,k) )

            index = 1 + my_floor((phi_out(i,j,k)-base_phi) * recip_dphi)
            j_out(i,j,k) = look_phi(index) - dsm1y
            If ( phi_out(i,j,k) >= phi_p(i, j_out(i,j,k)+1) )           &
     &                     j_out(i,j,k) = j_out(i,j,k) + 1
            weight_phi(i,j,k) = ( phi_out(i,j,k) -                      &
     &                             phi_p(i, j_out(i,j,k)) ) *           &
     &                       recip_dphip(i, j_out(i,j,k))
            End Do
          End Do
        End Do  !  k = 1, dim_k_out
!$OMP END PARALLEL DO

      End If ! L_regular

      IF (lhook) CALL dr_hook('CALC_INDEX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Calc_Index
