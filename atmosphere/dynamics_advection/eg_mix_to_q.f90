! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_mix_to_q_mod
CONTAINS
SUBROUTINE eg_mix_to_q                                                &
                  (qdims,mixdims,                                     &
                   mix_v, mix_cl, mix_cf,                             &
                   mix_cf2, mix_rain, mix_graup,                      &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,             &
                   q, qcl, qcf,                                       &
                   qcf2,qrain,qgraup                                  &
                   )


USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod, ONLY : array_dims

IMPLICIT NONE
!
! Description:
!  
!          Convert from mixing ratios to specific humidities
!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

TYPE (array_dims) , INTENT(IN) :: qdims,mixdims


LOGICAL                                                               &
  l_mcr_qcf2                                                          &
                    ! true if second cloud ice active
, l_mcr_qrain                                                         &
                    ! true if rain active
, l_mcr_qgraup      ! true if graupel active

! Arguments with Intent IN. ie: Input

REAL :: q (qdims%i_start:qdims%i_end,                                 &
           qdims%j_start:qdims%j_end,                                 &
           qdims%k_start:qdims%k_end)                                 

REAL ::  qcl  (qdims%i_start:qdims%i_end,                             &
           qdims%j_start:qdims%j_end,                                 &
           qdims%k_start:qdims%k_end) 

REAL ::  qcf (qdims%i_start:qdims%i_end,                              &
           qdims%j_start:qdims%j_end,                                 &
           qdims%k_start:qdims%k_end) 

REAL ::  qcf2  (qdims%i_start:qdims%i_end,                            &
           qdims%j_start:qdims%j_end,                                 &
           qdims%k_start:qdims%k_end) 

REAL ::  qrain  (qdims%i_start:qdims%i_end,                           &
           qdims%j_start:qdims%j_end,                                 &
           qdims%k_start:qdims%k_end) 

REAL ::  qgraup  (qdims%i_start:qdims%i_end,                          &
           qdims%j_start:qdims%j_end,                                 &
           qdims%k_start:qdims%k_end) 


! Arguments with Intent OUT. ie: Output

REAL ::  mix_v  (mixdims%i_start:mixdims%i_end,                       &
                 mixdims%j_start:mixdims%j_end,                       &
                 mixdims%k_start:mixdims%k_end)

REAL ::  mix_cl (mixdims%i_start:mixdims%i_end,                       &
                 mixdims%j_start:mixdims%j_end,                       &
                 mixdims%k_start:mixdims%k_end)

REAL ::  mix_cf  (mixdims%i_start:mixdims%i_end,                      &
                 mixdims%j_start:mixdims%j_end,                       &
                 mixdims%k_start:mixdims%k_end)

REAL ::  mix_cf2   (mixdims%i_start:mixdims%i_end,                    &
                 mixdims%j_start:mixdims%j_end,                       &
                 mixdims%k_start:mixdims%k_end)
 
REAL ::  mix_rain  (mixdims%i_start:mixdims%i_end,                    &
                 mixdims%j_start:mixdims%j_end,                       &
                 mixdims%k_start:mixdims%k_end)

REAL ::  mix_graup  (mixdims%i_start:mixdims%i_end,                   &
                 mixdims%j_start:mixdims%j_end,                       &
                 mixdims%k_start:mixdims%k_end)

REAL :: conv                                                        

REAL :: moist (mixdims%i_start:mixdims%i_end,                         &
               mixdims%j_start:mixdims%j_end,                         &
               mixdims%k_start:mixdims%k_end)

INTEGER :: i, j, k



! ---------------------------------------------------------------------
! Section 1. convert mix_v, mix_cl,mix_cf to q, qcl,qcf
! ---------------------------------------------------------------------

IF (lhook) CALL dr_hook('MIX_TO_Q',zhook_in,zhook_handle)
DO k=mixdims%k_start,mixdims%k_end
  DO j=mixdims%j_start,mixdims%j_end
    DO i=mixdims%i_start,mixdims%i_end
       moist(i,j,k)=1.0+ mix_v (i,j,k)+mix_cl(i,j,k)+mix_cf(i,j,k)
      END DO
    END DO
  END DO

IF (l_mcr_qcf2) THEN
DO k=mixdims%k_start,mixdims%k_end
  DO j=mixdims%j_start,mixdims%j_end
    DO i=mixdims%i_start,mixdims%i_end
       moist(i,j,k)= moist(i,j,k)+mix_cf2(i,j,k)
      END DO
    END DO
  END DO
END IF
IF (l_mcr_qrain) THEN
DO k=mixdims%k_start,mixdims%k_end
  DO j=mixdims%j_start,mixdims%j_end
    DO i=mixdims%i_start,mixdims%i_end
       moist(i,j,k)= moist(i,j,k)+mix_rain(i,j,k)
      END DO
    END DO
  END DO
END IF
IF (l_mcr_qgraup) THEN
DO k=mixdims%k_start,mixdims%k_end
  DO j=mixdims%j_start,mixdims%j_end
    DO i=mixdims%i_start,mixdims%i_end
       moist(i,j,k)= moist(i,j,k)+mix_graup(i,j,k)
      END DO
    END DO
  END DO
END IF

DO k=mixdims%k_start,mixdims%k_end
  DO j=mixdims%j_start,mixdims%j_end
    DO i=mixdims%i_start,mixdims%i_end
        conv= 1./ moist(i,j,k)
        q  (i,j,k) = mix_v (i,j,k)*conv
        qcl(i,j,k) = mix_cl(i,j,k)*conv
        qcf(i,j,k) = mix_cf(i,j,k)*conv
      END DO
    END DO
  END DO

IF (l_mcr_qcf2) THEN
DO k=mixdims%k_start,mixdims%k_end
  DO j=mixdims%j_start,mixdims%j_end
    DO i=mixdims%i_start,mixdims%i_end
        conv= 1./moist(i,j,k)
        qcf2(i,j,k)  = mix_cf2(i,j,k) * conv
      END DO
    END DO
  END DO
END IF
IF (l_mcr_qrain) THEN
DO k=mixdims%k_start,mixdims%k_end
  DO j=mixdims%j_start,mixdims%j_end
    DO i=mixdims%i_start,mixdims%i_end
        conv= 1./moist(i,j,k)
        qrain(i,j,k)  = mix_rain(i,j,k) * conv
      END DO
    END DO
  END DO
END IF
IF (l_mcr_qgraup) THEN
DO k=mixdims%k_start,mixdims%k_end
  DO j=mixdims%j_start,mixdims%j_end
    DO i=mixdims%i_start,mixdims%i_end
        conv= 1./moist(i,j,k)
        qgraup(i,j,k)  = mix_graup(i,j,k) * conv
      END DO
    END DO
  END DO
END IF

! end of routine


IF (lhook) CALL dr_hook('EG_MIX_TO_Q',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_mix_to_q
END MODULE eg_mix_to_q_mod
