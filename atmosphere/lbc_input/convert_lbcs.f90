! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE convert_lbcs(theta, q, qcl, qcf, qcf2, qr, qgr,               &
                              l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,           &
                              rimwidth)

      USE atm_fields_bounds_mod, ONLY : tdims_s, qdims_l
      USE atmos_constants_mod,   ONLY : recip_epsilon

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE
!
! Description: Convert ND prgnotic lbc data for ENDGame
!              Temporary fix until lbc output changed to
!              ENDGame prognostics.
!
! Code Owner: See Unified Model Code Owners HTML page
!
! This file belongs in section: LBC Input

      INTEGER, INTENT(IN)    :: rimwidth 
!                                 ( = LENRIM(fld_type_p,halo_type_extended) )
      REAL,    INTENT(INOUT) :: theta(rimwidth,tdims_s%k_start:tdims_s%k_end)
      REAL,    INTENT(INOUT) :: q(rimwidth,qdims_l%k_start:qdims_l%k_end)
      REAL,    INTENT(INOUT) :: qcl(rimwidth,qdims_l%k_start:qdims_l%k_end)
      REAL,    INTENT(INOUT) :: qcf(rimwidth,qdims_l%k_start:qdims_l%k_end)
      REAL,    INTENT(INOUT) :: qcf2(rimwidth,qdims_l%k_start:qdims_l%k_end)
      REAL,    INTENT(INOUT) :: qr(rimwidth,qdims_l%k_start:qdims_l%k_end)
      REAL,    INTENT(INOUT) :: qgr(rimwidth,qdims_l%k_start:qdims_l%k_end)

      LOGICAL, INTENT(IN)    :: l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup

      INTEGER                :: i, k
      REAL                   :: q2m_conv(rimwidth,tdims_s%k_start:tdims_s%k_end)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('convert_lbcs',zhook_in,zhook_handle)

      DO k = qdims_l%k_start, qdims_l%k_end
         DO i = 1, rimwidth
            q2m_conv(i,k) = 1.0 - q(i,k) - qcf(i,k) - qcl(i,k)
         END DO
      END DO

      IF( l_mcr_qcf2 ) THEN
         DO k = qdims_l%k_start, qdims_l%k_end
            DO i = 1, rimwidth
               q2m_conv(i,k) = q2m_conv(i,k) - qcf2(i,k)
            END DO
         END DO
      END IF

      IF( l_mcr_qrain ) THEN
         DO k = qdims_l%k_start, qdims_l%k_end
            DO i = 1, rimwidth
               q2m_conv(i,k) = q2m_conv(i,k) - qr(i,k)
            END DO
         END DO
      END IF

      IF( l_mcr_qgraup ) THEN
         DO k = qdims_l%k_start, qdims_l%k_end
            DO i = 1, rimwidth
               q2m_conv(i,k) = q2m_conv(i,k) - qgr(i,k)
            END DO
         END DO
      END IF

      ! Assumed here that theta levels and moisture levels match up.
      DO k = qdims_l%k_start, qdims_l%k_end
         DO i = 1, rimwidth
            q2m_conv(i,k) = 1.0/q2m_conv(i,k)

            q(i,k)     = q(i,k)   * q2m_conv(i,k)
            qcl(i,k)   = qcl(i,k) * q2m_conv(i,k)
            qcf(i,k)   = qcf(i,k) * q2m_conv(i,k)
            theta(i,k) = ( 1.0 + recip_epsilon*q(i,k) ) *theta(i,k)
         END DO
      END DO

      IF( l_mcr_qcf2 ) THEN
         DO k = qdims_l%k_start, qdims_l%k_end
            DO i = 1, rimwidth
               qcf2(i,k) = qcf2(i,k)   * q2m_conv(i,k)
            END DO
         END DO
      END IF

      IF( l_mcr_qrain ) THEN
         DO k = qdims_l%k_start, qdims_l%k_end
            DO i = 1, rimwidth
               qr(i,k) = qr(i,k)   * q2m_conv(i,k)
            END DO
         END DO
      END IF

      IF( l_mcr_qgraup ) THEN
         DO k = qdims_l%k_start, qdims_l%k_end
            DO i = 1, rimwidth
               qgr(i,k) = qgr(i,k)   * q2m_conv(i,k)
            END DO
         END DO
      END IF

      IF (lhook) CALL dr_hook('convert_lbcs',zhook_out,zhook_handle)

      END SUBROUTINE convert_lbcs

