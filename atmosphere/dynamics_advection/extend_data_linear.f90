! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Extend_data_linear.

      Subroutine Extend_data_linear(                                    &
     &                             Data_in1, Data_in2, Data_in3,        &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             number_of_inputs, halo_i, halo_j,    &
     &                             Ext_data)

! Purpose:
!          Extends data arrays into bigger arrays to
!          allow use of more efficient linear interpolation.
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
!  6.0  18/08/03  NEC SX-6 optimisation - R Barnes & J-C Rioual.
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
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, number_of_inputs !number of fields to interpolate.

      Real                                                              &
     &  Data_in1 (1-halo_i:dim_i_in+halo_i,                             &
     &            1-halo_j:dim_j_in+halo_j, dim_k_in)                   &
                                                      ! data to be
                                                      ! interpolated
     &, Data_in2 (1-halo_i:dim_i_in+halo_i,                             &
     &            1-halo_j:dim_j_in+halo_j, dim_k_in)                   &
                                                      ! optional second
                                                      ! field of data to
                                                      ! be interpolated
     &, Data_in3 (1-halo_i:dim_i_in+halo_i,                             &
     &            1-halo_j:dim_j_in+halo_j, dim_k_in) ! optional third
                                                      ! field of data to
                                                      ! be interpolated

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Ext_Data (1-halo_i:dim_i_in+halo_i+1,1-halo_j:dim_j_in+halo_j,  &
     &            -1:dim_k_in+2,number_of_inputs)

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k ,n     ! Loop indices

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines: None

! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Extend input data to bigger area to
!              allow interpolation to be done without having to re-do
!              any end points.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('EXTEND_DATA_LINEAR',zhook_in,zhook_handle)
      Do k = 1, dim_k_in

        Do j = 1-halo_j, dim_j_in+halo_j
          Do i = 1-halo_i, dim_i_in+halo_i
            Ext_Data (i,j,k,1) = Data_in1 (i,j,k)
          End Do
        End Do
      End Do

      If (number_of_inputs  >=  2) Then
        Do k = 1, dim_k_in
          Do j = 1-halo_j, dim_j_in+halo_j
             Do i = 1-halo_i, dim_i_in+halo_i
              Ext_Data (i,j,k,2) = Data_in2 (i,j,k)
            End Do
          End Do
        End Do
      End If

      If (number_of_inputs  ==  3) Then
        Do k = 1, dim_k_in
          Do j = 1-halo_j, dim_j_in+halo_j
             Do i = 1-halo_i, dim_i_in+halo_i
              Ext_Data (i,j,k,3) = Data_in3 (i,j,k)
            End Do
          End Do
        End Do
      End If

! Extend in k direction not required.

! End of routine.
      IF (lhook) CALL dr_hook('EXTEND_DATA_LINEAR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Extend_data_linear
