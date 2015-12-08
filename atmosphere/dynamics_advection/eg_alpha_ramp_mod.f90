! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_alpha_ramp_mod

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

INTEGER, SAVE :: alpha_relax_type = imdi
INTEGER, SAVE :: alpha_relax_int  = 1

INTEGER, PARAMETER :: alpha_relax_type_constant = 1
INTEGER, PARAMETER :: alpha_relax_type_sqrt     = 2
INTEGER, PARAMETER :: alpha_relax_type_coslaw   = 3
INTEGER, PARAMETER :: alpha_relax_type_step     = 4


CONTAINS

SUBROUTINE eg_alpha_ramp(chain_number)

USE eg_alpha_mod
USE timestep_mod,    ONLY : timestep_number
USE proc_info_mod,   ONLY : me
USE yomhook,         ONLY : lhook, dr_hook
USE parkind1,        ONLY : jprb, jpim
USE conversions_mod, ONLY : pi
USE PrintStatus_mod
USE ereport_mod,     ONLY : ereport

IMPLICIT NONE
!
! Description: 
!        
!
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


LOGICAL, SAVE :: first_call = .TRUE.
REAL,    SAVE :: alpha_star,beta_star

INTEGER chain_number

CHARACTER(len=256)            :: Cmessage
CHARACTER(len=15)             :: Routine = 'eg_ramp_alpha'
INTEGER                       :: ICODE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('EG_RAMP_ALPHA',zhook_in,zhook_handle)

IF(first_call) THEN

    alpha_star = alpha_u
    beta_star  = 1.-alpha_u
    first_call = .FALSE.

END IF


SELECT CASE (alpha_relax_type)

  CASE (alpha_relax_type_sqrt)

    ALPHA_U     = beta_star*1./DBLE(timestep_number)**.5 + alpha_star

    ALPHA_V     = ALPHA_U
    ALPHA_W     = ALPHA_U
    ALPHA_RHO   = ALPHA_U
    ALPHA_P     = ALPHA_U

    alpha_changed = .TRUE.

  CASE (alpha_relax_type_constant)

!   nothing to be done here, we simply use the alphas from the namelist

  CASE (alpha_relax_type_coslaw)

    IF(timestep_number.lt.alpha_relax_int.and.chain_number.eq.0) THEN

      ALPHA_U   = beta_star*cos(pi*DBLE(timestep_number)/(2.*DBLE(alpha_relax_int)))**2 + alpha_star

    ELSE

      ALPHA_U   = alpha_star

    END IF

    ALPHA_V     = ALPHA_U
    ALPHA_W     = ALPHA_U
    ALPHA_RHO   = ALPHA_U
    ALPHA_P     = ALPHA_U
    alpha_changed = .TRUE.

  CASE (alpha_relax_type_step)

    IF (timestep_number.le.alpha_relax_int.and.chain_number.eq.0) THEN

      ALPHA_U   = 1.0

    ELSE

      ALPHA_U   = alpha_star

    END IF

    ALPHA_V     = ALPHA_U
    ALPHA_W     = ALPHA_U
    ALPHA_RHO   = ALPHA_U
    ALPHA_P     = ALPHA_U

    alpha_changed = .TRUE.

  CASE DEFAULT

    icode = 1
    cmessage = 'unknown alpha relax type'
    CALL Ereport(Routine,ICODE,CMESSAGE)

END SELECT


IF( me == 0 .AND. PrintStatus > PrStatus_Normal)     &
       WRITE(6,fmt='(A,E22.15)') 'ALPHA:', ALPHA_U

IF (lhook) CALL dr_hook('EG_RAMP_ALPHA',zhook_out,zhook_handle)

END SUBROUTINE
END MODULE
