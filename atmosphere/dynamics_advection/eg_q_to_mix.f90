! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_q_to_mix_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_q_to_mix                                                &
                  (qdims,mixdims,                                     &
                   q, qcl, qcf,                                       &
                   qcf2, qrain, qgraup,                               &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,             &
                   mix_v, mix_cl, mix_cf,                             &
                   mix_cf2, mix_rain, mix_graup , swap_in             &
                   )


USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE Field_Types
USE proc_info_mod, ONLY : model_domain
USE domain_params
USE eg_swap_bounds_mod
USE atm_fields_bounds_mod, ONLY : array_dims

IMPLICIT NONE
!
! Description:
!          Convert from specific humidities to mixing ratios
!  
!
! Method:
!          See ENDGame Formulation version 3.03 Section 7.10
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

LOGICAL, OPTIONAL :: swap_in
LOGICAL :: swap

INTEGER                                                               &
  row_length                                                          &
                  ! number of points on a row    
, rows                                                                
                  ! number of rows of data

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
                  mixdims%j_start:mixdims%j_end,                      &
                  mixdims%k_start:mixdims%k_end)

REAL ::  mix_cf2   (mixdims%i_start:mixdims%i_end,                    &
                    mixdims%j_start:mixdims%j_end,                    &
                    mixdims%k_start:mixdims%k_end)
 
REAL ::  mix_rain  (mixdims%i_start:mixdims%i_end,                    &
                    mixdims%j_start:mixdims%j_end,                    &
                    mixdims%k_start:mixdims%k_end)

REAL ::  mix_graup  (mixdims%i_start:mixdims%i_end,                   &
                     mixdims%j_start:mixdims%j_end,                   &
                     mixdims%k_start:mixdims%k_end)

! local variables

REAL :: moist_conv (mixdims%i_start:mixdims%i_end,                    &
               mixdims%j_start:mixdims%j_end,                         &
               mixdims%k_start:mixdims%k_end)

INTEGER :: i, j, k


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_Q_TO_MIX',zhook_in,zhook_handle)

IF (PRESENT(swap_in)) THEN
  swap = swap_in
ELSE
  swap = .TRUE.
END IF

row_length = mixdims%i_end-mixdims%i_start+1-2*mixdims%halo_i
rows       = mixdims%j_end-mixdims%j_start+1-2*mixdims%halo_j

! ----------------------------------------------------------------------
! Section 1. convert q, qcl,qcf to mix_v, mix_cl,mix_cf
! ----------------------------------------------------------------------
DO k=mixdims%k_start,mixdims%k_end
  DO j=1,rows
    DO i=1,row_length

      moist_conv(i,j,k) = 1 - q(i,j,k) - qcl(i,j,k) - qcf(i,j,k)

    END DO
  END DO
END DO

IF (l_mcr_qcf2) THEN
  DO k=mixdims%k_start,mixdims%k_end
    DO j=1,rows
      DO i=1,row_length

        moist_conv(i,j,k)      = moist_conv(i,j,k) - qcf2(i,j,k)

      END DO
    END DO
  END DO
END IF

IF (l_mcr_qrain) THEN
  DO k=mixdims%k_start,mixdims%k_end
    DO j=1,rows
      DO i=1,row_length

        moist_conv(i,j,k)      = moist_conv(i,j,k) - qrain(i,j,k)

      END DO
    END DO
  END DO
END IF

IF (l_mcr_qgraup) THEN
  DO k=mixdims%k_start,mixdims%k_end
    DO j=1,rows
      DO i=1,row_length

        moist_conv(i,j,k)      = moist_conv(i,j,k) - qgraup(i,j,k)

      END DO
    END DO
  END DO
END IF

DO k=mixdims%k_start,mixdims%k_end
  DO j=1,rows
    DO i=1,row_length

        moist_conv(i,j,k)      = 1./moist_conv(i,j,k)

    END DO
  END DO
END DO


DO k = mixdims%k_start,mixdims%k_end
  DO j=1,rows
    DO i=1,row_length
      mix_v (i,j,k) = q  (i,j,k)*moist_conv(i,j,k)
      mix_cl(i,j,k) = qcl(i,j,k)*moist_conv(i,j,k)
      mix_cf(i,j,k) = qcf(i,j,k)*moist_conv(i,j,k)
    END DO
  END DO
END DO

IF(     l_mcr_qcf2.AND..NOT.l_mcr_qrain.AND..NOT.l_mcr_qgraup) CALL eq_conv1(   &
                                          qdims,mixdims,moist_conv,mix_cf2,qcf2)
IF(.NOT.l_mcr_qcf2.AND.     l_mcr_qrain.AND..NOT.l_mcr_qgraup) CALL eq_conv1(   &
                                        qdims,mixdims,moist_conv,mix_rain,qrain)
IF(.NOT.l_mcr_qcf2.AND..NOT.l_mcr_qrain.AND.     l_mcr_qgraup) CALL eq_conv1(   &
                                      qdims,mixdims,moist_conv,mix_graup,qgraup)

IF(     l_mcr_qcf2.AND.     l_mcr_qrain.AND..NOT.l_mcr_qgraup) CALL eq_conv2(   &
                           qdims,mixdims,moist_conv,mix_cf2,qcf2,mix_rain,qrain)
IF(.NOT.l_mcr_qcf2.AND.     l_mcr_qrain.AND.     l_mcr_qgraup) CALL eq_conv2(   &
                       qdims,mixdims,moist_conv,mix_rain,qrain,mix_graup,qgraup)
IF(     l_mcr_qcf2.AND..NOT.l_mcr_qrain.AND.     l_mcr_qgraup) CALL eq_conv2(   &
                         qdims,mixdims,moist_conv,mix_graup,qgraup,mix_cf2,qcf2)

IF(     l_mcr_qcf2.AND.     l_mcr_qrain.AND.     l_mcr_qgraup) THEN

  DO k = mixdims%k_start,mixdims%k_end
    DO j=1,rows
      DO i=1,row_length

        mix_cf2(i,j,k)   = qcf2(i,j,k)   * moist_conv(i,j,k)
        mix_rain(i,j,k)  = qrain(i,j,k)  * moist_conv(i,j,k)
        mix_graup(i,j,k) = qgraup(i,j,k) * moist_conv(i,j,k)

      END DO
    END DO
  END DO

END IF


IF (swap) THEN
  CALL eg_swap_bounds(mix_v  ,mixdims,fld_type_p,.FALSE.)
  CALL eg_swap_bounds(mix_cl ,mixdims,fld_type_p,.FALSE.)
  CALL eg_swap_bounds(mix_cf ,mixdims,fld_type_p,.FALSE.)

  IF (l_mcr_qcf2   ) CALL eg_swap_bounds(mix_cf2   ,mixdims,fld_type_p,.FALSE.)
  IF (l_mcr_qrain  ) CALL eg_swap_bounds(mix_rain  ,mixdims,fld_type_p,.FALSE.)
  IF (l_mcr_qgraup ) CALL eg_swap_bounds(mix_graup ,mixdims,fld_type_p,.FALSE.)

END IF

! end of routine

IF (lhook) CALL dr_hook('EG_Q_TO_MIX',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_q_to_mix

SUBROUTINE eq_conv1 (qdims,mixdims,moist_conv,mix,q)

USE atm_fields_bounds_mod, ONLY : array_dims

IMPLICIT NONE
!
! Description:
!          Convert from specific humidities to mixing ratios
!  
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

TYPE (array_dims) , INTENT(IN) :: qdims,mixdims

REAL :: q (qdims%i_start:qdims%i_end,                                 &
           qdims%j_start:qdims%j_end,                                 &
           qdims%k_start:qdims%k_end)                                 

REAL ::  mix    (mixdims%i_start:mixdims%i_end,                       &
                 mixdims%j_start:mixdims%j_end,                       &
                 mixdims%k_start:mixdims%k_end)

REAL ::  moist_conv (mixdims%i_start:mixdims%i_end,                   &
                     mixdims%j_start:mixdims%j_end,                   &
                     mixdims%k_start:mixdims%k_end)

INTEGER :: i,j,k

  DO k = mixdims%k_start,mixdims%k_end
    DO j= mixdims%j_start+mixdims%halo_j,mixdims%j_end-mixdims%halo_j
      DO i= mixdims%i_start+mixdims%halo_i,mixdims%i_end-mixdims%halo_i

        mix(i,j,k)   = q(i,j,k)   * moist_conv(i,j,k)

      END DO
    END DO
  END DO

END SUBROUTINE


SUBROUTINE eq_conv2 (qdims,mixdims,moist_conv,mix1,q1,mix2,q2)

USE atm_fields_bounds_mod, ONLY : array_dims

IMPLICIT NONE
!
! Description:
!          Convert from specific humidities to mixing ratios
!  
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

TYPE (array_dims) , INTENT(IN) :: qdims,mixdims

REAL :: q1(qdims%i_start:qdims%i_end,                                 &
           qdims%j_start:qdims%j_end,                                 &
           qdims%k_start:qdims%k_end)                                 

REAL ::  mix1   (mixdims%i_start:mixdims%i_end,                       &
                 mixdims%j_start:mixdims%j_end,                       &
                 mixdims%k_start:mixdims%k_end)

REAL :: q2(qdims%i_start:qdims%i_end,                                 &
           qdims%j_start:qdims%j_end,                                 &
           qdims%k_start:qdims%k_end)                                 

REAL ::  mix2   (mixdims%i_start:mixdims%i_end,                       &
                 mixdims%j_start:mixdims%j_end,                       &
                 mixdims%k_start:mixdims%k_end)

REAL ::  moist_conv (mixdims%i_start:mixdims%i_end,                   &
                     mixdims%j_start:mixdims%j_end,                   &
                     mixdims%k_start:mixdims%k_end)

INTEGER i,j,k

  DO k = mixdims%k_start,mixdims%k_end
    DO j= mixdims%j_start+mixdims%halo_j,mixdims%j_end-mixdims%halo_j
      DO i= mixdims%i_start+mixdims%halo_i,mixdims%i_end-mixdims%halo_i

        mix1(i,j,k)   = q1(i,j,k)   * moist_conv(i,j,k)
        mix2(i,j,k)   = q2(i,j,k)   * moist_conv(i,j,k)

      END DO
    END DO
  END DO

END SUBROUTINE

END MODULE eg_q_to_mix_mod
