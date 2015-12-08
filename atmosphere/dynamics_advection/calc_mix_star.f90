! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine calc_mix_star

      SUBROUTINE calc_mix_star                                          &
     &                  (row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, offx, offy,                    &
     &                   q, qcl, qcf,                                   &
     &                   qcf2, qrain, qgraup,                           &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   qcf2_star, qrain_star, qgraup_star,            &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   mix_v_star, mix_cl_star, mix_cf_star           &
     &                  ,mix_cf2_star, mix_rain_star, mix_graup_star    &
     &                   )

! Purpose:
!          calculate mixing ratio increments
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: Fortran 90 + extensions
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
     &, halo_j                                                          &
     &, offx                                                            &
     &, offy

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

      REAL                                                              &
     &  q_star    (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &        wet_levels)                                               &
     &, qcl_star  (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &        wet_levels)                                               &
     &, qcf_star  (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &        wet_levels)                                               &
     &, qcf2_star    (1-offx:row_length+offx, 1-offy:rows+offy,         &
     &                wet_levels)                                       &
     &, qrain_star   (1-offx:row_length+offx, 1-offy:rows+offy,         &
     &                wet_levels)                                       &
     &, qgraup_star  (1-offx:row_length+offx, 1-offy:rows+offy,         &
     &                wet_levels)

! Arguments with Intent OUT. ie: Output

      REAL                                                              &
     &  mix_v_star  (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &          wet_levels)                                             &
     &, mix_cl_star (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &          wet_levels)                                             &
     &, mix_cf_star (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &          wet_levels)                                             &
     &, mix_cf2_star   (1-offx:row_length+offx, 1-offy:rows+offy,       &
     &                  wet_levels)                                     &
     &, mix_rain_star  (1-offx:row_length+offx, 1-offy:rows+offy,       &
     &                  wet_levels)                                     &
     &, mix_graup_star (1-offx:row_length+offx, 1-offy:rows+offy,       &
     &                  wet_levels)

! local variables
      REAL                                                              &
     & sum_q,sum_q_star                                                 &
     &, conv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)          &
     &, moist(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, moist_star(1-offx:row_length+offx, 1-offy:rows+offy)

      INTEGER                                                           &
     & i, j, k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 1. calculate mix_v_star, mix_cl_star, mix_cf_star
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CALC_MIX_STAR',zhook_in,zhook_handle)

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, j, k,    & 
!$OMP& conv, moist, moist_star)
      DO k=1,wet_levels
        DO j=1,rows
          DO i=1,row_length
            moist(i,j) = 1 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k)
            moist_star(i,j) = q_star(i,j,k) + qcl_star(i,j,k) +         &
                 &                          qcf_star(i,j,k)
          END DO
        END DO

        IF (L_mcr_qcf2) THEN
          DO j=1,rows
            DO i=1,row_length
              moist(i,j)      = moist(i,j) - qcf2(i,j,k)
              moist_star(i,j) = moist_star(i,j) + qcf2_star(i,j,k)
            END DO
          END DO
        END IF

        IF (L_mcr_qrain) THEN
          DO j=1,rows
            DO i=1,row_length
              moist(i,j)      = moist(i,j) - qrain(i,j,k)
              moist_star(i,j) = moist_star(i,j) + qrain_star(i,j,k)
            END DO
          END DO
        END IF
        IF (L_mcr_qgraup) THEN
          DO j=1,rows
            DO i=1,row_length
              moist(i,j)      = moist(i,j) - qgraup(i,j,k)
              moist_star(i,j) = moist_star(i,j) + qgraup_star(i,j,k)
            END DO
          END DO
        END IF

        DO j = 1, rows
          DO i = 1, row_length
            conv(i,j)= 1./(  moist(i,j) *                               &
                 (moist(i,j) - moist_star(i,j)) )
          END DO
        END DO

        DO j = 1, rows
          DO i = 1, row_length
            mix_v_star(i,j,k)  =                                        &
                 ( q_star(i,j,k)*moist(i,j) +                           &
                 q(i,j,k)*moist_star(i,j) )  *conv(i,j)
            mix_cl_star(i,j,k) =                                        &
                 ( qcl_star(i,j,k)*moist(i,j) +                         &
                 qcl(i,j,k)*moist_star(i,j) )  *conv(i,j)
            mix_cf_star(i,j,k) =                                        &
                 ( qcf_star(i,j,k)*moist(i,j) +                         &    
                 qcf(i,j,k)*moist_star(i,j) )  *conv(i,j)
          END DO
        END DO
        
        IF (L_mcr_qcf2) THEN
          DO j = 1, rows
            DO i = 1, row_length
              mix_cf2_star(i,j,k)  =                                    &
     &               ( qcf2_star(i,j,k)*moist(i,j) +                    &
     &                 qcf2(i,j,k)*moist_star(i,j) )  *conv(i,j)
            END DO
          END DO
        END IF

        IF (L_mcr_qrain) THEN
          DO j = 1, rows
            DO i = 1, row_length              
              mix_rain_star(i,j,k)  =                                   &
     &               ( qrain_star(i,j,k)*moist(i,j) +                   &
     &                 qrain(i,j,k)*moist_star(i,j) )  *conv(i,j)
            END DO
          END DO
        END IF

      IF (L_mcr_qgraup) THEN
        DO j = 1, rows
          DO i = 1, row_length
            mix_graup_star(i,j,k)  =                                    &
     &               ( qgraup_star(i,j,k)*moist(i,j) +                  &
     &                 qgraup(i,j,k)*moist_star(i,j) )  *conv(i,j)
            END DO
          END DO
        END IF

      END DO
!$OMP END PARALLEL DO
! end of routine

      IF (lhook) CALL dr_hook('CALC_MIX_STAR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE calc_mix_star
