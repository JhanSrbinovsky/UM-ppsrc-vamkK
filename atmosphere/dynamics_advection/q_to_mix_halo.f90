! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine q_to_mix


! Subroutine q_to_mix_halo

      SUBROUTINE q_to_mix_halo                                          &
     &                  (row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j,                                &
     &                   q, qcl, qcf,                                   &
     &                   qcf2, qrain, qgraup,                           &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   mix_v, mix_cl, mix_cf                          &
     &                  ,mix_cf2, mix_rain, mix_graup                   &
     &                   )

! Purpose:
!          Convert from specific humidities to mixing ratios
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER                                                           &
     &  row_length                                                      &
                        ! number of points on a row
     &, rows                                                            &
                        ! number of rows of data
     &, wet_levels                                                      &
     &, halo_i                                                          &
     &, halo_j

      LOGICAL                                                           &
     &  L_mcr_qcf2                                                      &
                          ! true if second cloud ice active
     &, L_mcr_qrain                                                     &
                          ! true if rain active
     &, L_mcr_qgraup      ! true if graupel active

! Arguments with Intent IN. ie: Input

      REAL                                                              &
     &  q    (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        wet_levels)                                               &
     &, qcl  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        wet_levels)                                               &
     &, qcf  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        wet_levels)                                               &
     &, qcf2    (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           wet_levels)                                            &
     &, qrain   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           wet_levels)                                            &
     &, qgraup  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           wet_levels)

! Arguments with Intent OUT. ie: Output

      REAL                                                              &
     &  mix_v  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_levels)                                             &
     &, mix_cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_levels)                                             &
     &, mix_cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_levels)                                             &
     &, mix_cf2   (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &             wet_levels)                                          &
     &, mix_rain  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &             wet_levels)                                          &
     &, mix_graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &             wet_levels)

! local variables
      REAL :: moist(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)

      INTEGER                                                           &
     & i, j, k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 1. convert q, qcl,qcf to mix_v, mix_cl,mix_cf
! ----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('Q_TO_MIX_HALO',zhook_in,zhook_handle)

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED)                     &
!$OMP& PRIVATE(i,j,k,moist) 
      DO k=1,wet_levels

        DO j=1-halo_j,rows+halo_j
          DO i=1-halo_i,row_length+halo_i
            moist(i,j) = 1 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k)
          END DO
        END DO

        IF (L_mcr_qcf2) THEN
          DO j=1-halo_j,rows+halo_j
            DO i=1-halo_i,row_length+halo_i
              moist(i,j)      = moist(i,j) - qcf2(i,j,k)
            END DO
          END DO
        END IF

        IF (L_mcr_qrain) THEN
          DO j=1-halo_j,rows+halo_j
            DO i=1-halo_i,row_length+halo_i
              moist(i,j)      = moist(i,j) - qrain(i,j,k)
            END DO
          END DO
        END IF
        
        IF (L_mcr_qgraup) THEN
          DO j=1-halo_j,rows+halo_j
            DO i=1-halo_i,row_length+halo_i
              moist(i,j)      = moist(i,j) - qgraup(i,j,k)
            END DO
          END DO
        END IF
        
        DO j=1-halo_j,rows+halo_j
          DO i=1-halo_i,row_length+halo_i       
            moist(i,j) = 1./moist(i,j)
          END DO
        END DO

        DO j=1-halo_j,rows+halo_j
          DO i=1-halo_i,row_length+halo_i
            mix_v (i,j,k) = q  (i,j,k) * moist(i,j)
            mix_cl(i,j,k) = qcl(i,j,k) * moist(i,j)
            mix_cf(i,j,k) = qcf(i,j,k) * moist(i,j)
          END DO
        END DO

        IF (L_mcr_qcf2) THEN
          DO j=1-halo_j,rows+halo_j
            DO i=1-halo_i,row_length+halo_i
              mix_cf2(i,j,k) = qcf2(i,j,k) * moist(i,j)
            END DO
          END DO
        END IF

        IF (L_mcr_qrain) THEN
          DO j=1-halo_j,rows+halo_j
            DO i=1-halo_i,row_length+halo_i
              mix_rain(i,j,k) = qrain(i,j,k) * moist(i,j)
            END DO
          END DO
        END IF
        
        IF (L_mcr_qgraup) THEN
          DO j=1-halo_j,rows+halo_j
            DO i=1-halo_i,row_length+halo_i
              mix_graup(i,j,k) = qgraup(i,j,k) * moist(i,j)
            END DO
          END DO
        END IF
               
      END DO
!$OMP END PARALLEL DO 

! end of routine

      IF (lhook) CALL dr_hook('Q_TO_MIX_HALO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE q_to_mix_halo
