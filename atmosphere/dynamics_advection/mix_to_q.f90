! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine mix_to_q

      Subroutine mix_to_q                                               &
     &                  (row_length, rows, model_levels,                  &
     &                   halo_i, halo_j,                                &
     &                   mix_v, mix_cl, mix_cf,                         &
     &                   mix_cf2, mix_rain, mix_graup,                  &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   q, qcl, qcf                                    &
     &                  ,qcf2,qrain,qgraup                              &
     &                   )

! Purpose:
!          Convert from mixing ratios to specific humidities
!
! Method:
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

      Integer                                                           &
     &  row_length                                                      &
                        ! number of points on a row
     &, rows                                                            &
                        ! number of rows of data
     &, model_levels                                                      &
     &, halo_i                                                          &
     &, halo_j

      Logical                                                           &
     &  L_mcr_qcf2                                                      &
                          ! true if second cloud ice active
     &, L_mcr_qrain                                                     &
                          ! true if rain active
     &, L_mcr_qgraup      ! true if graupel active

! Arguments with Intent IN. ie: Input

      Real                                                              &
     &  mix_v  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                             &
     &, mix_cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                             &
     &, mix_cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                             &
     &, mix_cf2   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &             model_levels)                                          &
     &, mix_rain  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &             model_levels)                                          &
     &, mix_graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &             model_levels)

! Arguments with Intent OUT. ie: Output

      Real                                                              &
     &  q    (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        model_levels)                                               &
     &, qcl  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        model_levels)                                               &
     &, qcf  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        model_levels)                                               &
     &, qcf2    (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           model_levels)                                            &
     &, qrain   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           model_levels)                                            &
     &, qgraup  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           model_levels)

! local variables
      REAL :: moist   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)

      Integer                                                           &
     & i, j, k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 1. convert mix_v, mix_cl,mix_cf to q, qcl,qcf
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('MIX_TO_Q',zhook_in,zhook_handle)
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k,moist)
      Do k = 1, model_levels

        Do j = 1, rows
          Do i = 1, row_length
            moist(i,j) =1.0+ mix_v (i,j,k)+mix_cl(i,j,k)+mix_cf(i,j,k)
          End Do
        End Do
        
        If (L_mcr_qcf2) Then
          
          Do j = 1, rows
            Do i = 1, row_length
              moist(i,j) = moist(i,j)+mix_cf2(i,j,k)
            End Do
          End Do
          
        End If
        If (L_mcr_qrain) Then
          
          Do j = 1, rows
            Do i = 1, row_length
              moist(i,j)= moist(i,j)+mix_rain(i,j,k)
            End Do
          End Do
          
        End If
        If (L_mcr_qgraup) Then
          
          Do j = 1, rows
            Do i = 1, row_length
              moist(i,j)= moist(i,j)+mix_graup(i,j,k)
            End Do
          End Do
          
        End If

        Do j = 1, rows
          Do i = 1, row_length
            moist(i,j)= 1./moist(i,j)
          End Do
        End Do

        Do j = 1, rows
          Do i = 1, row_length
            q  (i,j,k) = mix_v (i,j,k) * moist(i,j)
            qcl(i,j,k) = mix_cl(i,j,k) * moist(i,j)
            qcf(i,j,k) = mix_cf(i,j,k) * moist(i,j)
          End Do
        End Do
        
        If (L_mcr_qcf2) Then
          
          Do j = 1, rows
            Do i = 1, row_length
              qcf2(i,j,k)  = mix_cf2(i,j,k) * moist(i,j)
            End Do
          End Do
          
        End If
        If (L_mcr_qrain) Then
          
          Do j = 1, rows
            Do i = 1, row_length
              qrain(i,j,k)  = mix_rain(i,j,k) * moist(i,j)
            End Do
          End Do
          
        End If
        If (L_mcr_qgraup) Then
          
          Do j = 1, rows
            Do i = 1, row_length
              qgraup(i,j,k)  = mix_graup(i,j,k) * moist(i,j)
            End Do
          End Do
          
        End If
        
      End Do
!$OMP END PARALLEL DO

! end of routine

      IF (lhook) CALL dr_hook('MIX_TO_Q',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE mix_to_q
