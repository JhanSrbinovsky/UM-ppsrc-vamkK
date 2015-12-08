! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Mono_Enforce.

      Subroutine Mono_Enforce(                                          &
     &                        Ext_Data, number_of_inputs,               &
     &                        dim_i_in, dim_j_in, dim_k_in,             &
     &                        dim_i_out, dim_j_out, dim_k_out,          &
     &                        halo_i, halo_j,                           &
     &                        i_out, j_out, k_out,                      &
     &                        Data_out_high)

! Purpose:
!          Ensures scheme is monotone.
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
                    ! Dimension of input Data arrays in i direction.
     &, dim_j_in                                                        &
                    ! Dimension of input Data arrays in j direction.
     &, dim_k_in                                                        &
                    ! Dimension of input Data arrays in k direction.
     &, dim_i_out                                                       &
                    ! Dimension of output Data arrays in i direction.
     &, dim_j_out                                                       &
                    ! Dimension of output Data arrays in j direction.
     &, dim_k_out                                                       &
                    ! Dimension of output Data arrays in k direction.
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, number_of_inputs

      INTEGER (KIND=integer32) ::                                       &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
     &, j_out (dim_i_out, dim_j_out, dim_k_out)                         &
     &, k_out (dim_i_out, dim_j_out, dim_k_out)

      Real                                                              &
     &  Ext_Data (1-halo_i:dim_i_in+halo_i+1,                           &
     &            1-halo_j:dim_j_in+halo_j,-1:dim_k_in+2                &
     &            ,number_of_inputs) !data to be interpolated


! Arguments with Intent IN/OUT.

      Real                                                              &
     &  Data_out_high (dim_i_out, dim_j_out, dim_k_out,                 &
     &                 number_of_inputs)

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k, n        ! Loop indices

      Real                                                              &
     &  max_mono                                                        &
     &, min_mono

      Integer :: ii      ! Indexing
      Integer :: jj      ! Indexing
      Integer :: kk      ! Indexing
      Real :: t1         ! Temporary
      Real :: t2         ! Temporary
      Real :: t3         ! Temporary
      Real :: t4         ! Temporary
      Real :: t5         ! Temporary
      Real :: t6         ! Temporary
      Real :: t7         ! Temporary
      Real :: t8         ! Temporary
      Real :: tmax1      ! Temporary max
      Real :: tmax2      ! Temporary max
      Real :: tmax3      ! Temporary max
      Real :: tmax4      ! Temporary max
      Real :: tmin1      ! Temporary min
      Real :: tmin2      ! Temporary min
      Real :: tmin3      ! Temporary min
      Real :: tmin4      ! Temporary min

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines: None
! subroutines: None
! Functions: None

! ----------------------------------------------------------------------
!  Section 1.  Set alpha_max to 1 as both input fields are monotone.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('MONO_ENFORCE',zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,k,j,i,max_mono,min_mono)      &
!$OMP& PRIVATE(ii,jj,kk, t1,t2,t3,t4,t5,t6,t7,t8)                       &
!$OMP& PRIVATE(tmax1,tmax2,tmax3,tmax4, tmin1,tmin2,tmin3,tmin4)        &
!$OMP& SHARED(number_of_inputs,dim_k_out,dim_j_out,dim_i_out,Ext_Data,  &
!$OMP&  data_out_high,i_out,j_out,k_out) SCHEDULE(STATIC)
      Do n = 1, number_of_inputs
        Do k = 1, dim_k_out
          Do j = 1, dim_j_out
            Do i = 1, dim_i_out

! Find max and min monotone values for the point concerned.

              max_mono = max (                                          &
     &       Ext_Data(i_out(i,j,k),j_out(i,j,k),k_out(i,j,k),n),        &
     &       Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k),n),      &
     &       Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k),n),      &
     &       Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k),n),    &
     &       Ext_Data(i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+1,n),      &
     &       Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+1,n),    &
     &       Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k)+1,n),    &
     &       Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k)+1,n) )

              min_mono = min (                                          &
     &       Ext_Data(i_out(i,j,k),j_out(i,j,k),k_out(i,j,k),n),        &
     &       Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k),n),      &
     &       Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k),n),      &
     &       Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k),n),    &
     &       Ext_Data(i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+1,n),      &
     &       Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+1,n),    &
     &       Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k)+1,n),    &
     &       Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k)+1,n) )

              If ( data_out_high(i,j,k,n)  >   max_mono )               &
     &             data_out_high(i,j,k,n) = max_mono
              If ( data_out_high(i,j,k,n)  <   min_mono )               &
     &             data_out_high(i,j,k,n) = min_mono

            End Do
          End Do
        End Do
      End Do
!$OMP END PARALLEL DO

! End of routine.
      IF (lhook) CALL dr_hook('MONO_ENFORCE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Mono_Enforce
