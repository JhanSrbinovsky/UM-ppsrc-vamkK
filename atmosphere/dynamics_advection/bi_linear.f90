! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Bi_Linear.
      subroutine Bi_Linear (dim_i_out, dim_j_out, dim_k_out,            &
     &                      dim_i_in, dim_j_in, dim_k_in,               &
     &                      halo_i, halo_j, Data_in,                    &
     &                      i_out, j_out, weight_lambda, weight_phi,    &
     &                      Data_out)

! Purpose:
!          Performs bi-linear horizontal interpolation of a field
!          defined on one grid to another grid.
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
                    ! Dimension of Data_out in j direction.
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j      ! Size of halo in j direction.

      INTEGER (KIND=integer32) ::                                       &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
     &, j_out (dim_i_out, dim_j_out, dim_k_out)

      Real                                                              &
     &  weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)

      Real                                                              &
     &  Data_in (1-halo_i:dim_i_in+halo_i, 1-halo_j:dim_j_in+halo_j,    &
     &           dim_k_in)                 ! data to be interpolated

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Data_out(dim_i_out, dim_j_out, dim_k_out)

! Local variables

      Integer i, j, k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines: None

! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Perform bi-linear interpolation.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('BI_LINEAR',zhook_in,zhook_handle)

!$OMP  PARALLEL DO DEFAULT(NONE) SHARED(data_out, data_in, weight_phi, &
!$OMP& j_out, i_out, weight_lambda, dim_k_out, dim_j_out, dim_i_out)   &
!$OMP& PRIVATE(i, j, k) SCHEDULE(STATIC)
      Do k = 1, dim_k_out
         Do j = 1, dim_j_out
            Do i = 1, dim_i_out

               Data_out (i,j,k) = (1.-weight_lambda(i,j,k)) *           &
     &              (1.-weight_phi(i,j,k)) *                            &
     &              Data_in (i_out(i,j,k),j_out(i,j,k),1)               &
     &              + weight_lambda(i,j,k) * (1.-weight_phi(i,j,k))     &
     &              * Data_in (i_out(i,j,k)+1,j_out(i,j,k),1)           &
     &              + (1.-weight_lambda(i,j,k)) * weight_phi(i,j,k)     &
     &              * Data_in (i_out(i,j,k),j_out(i,j,k)+1,1)           &
     &              + weight_lambda(i,j,k) * weight_phi(i,j,k) *        &
     &              Data_in (i_out(i,j,k)+1,j_out(i,j,k)+1,1)
            End Do
         End Do
      End Do
!$OMP END PARALLEL DO 

! end of routine

      IF (lhook) CALL dr_hook('BI_LINEAR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Bi_Linear
