! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate the effect of convection upon the large-scale cloud fractions.
!
MODULE pc2_environ_mod

IMPLICIT NONE

! Description:
! Calculate the effect of convection upon the large-scale cloud fractions.
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.

CONTAINS

! Subroutine Interface:
SUBROUTINE pc2_environ (k, npnts,                                         &
                    qclek, qclekp1, qcfek, qcfekp1,                       &
                    cflek, cflekp1,  cffek,  cffekp1,                     &
                    bcfek,  bcfekp1,                                      &
                    qclpk, qclpkp1, qcfpk, qcfpkp1,                       &
                    dqclek, dqcfek,                                       &
                    dqclekp1, dqcfekp1,                                   &
                    l_q_interact,                                         &
                    bterm,                                                &
                    ! Out
                    dcflek, dcffek, dbcfek,                               &
                    dcflekp1, dcffekp1, dbcfekp1)
                    

USE water_constants_mod, ONLY: lc, lf
USE cv_derived_constants_mod, ONLY: ls
USE atmos_constants_mod, ONLY: r, cp

USE cv_run_mod, ONLY: eff_dcfl, eff_dcff

USE cloud_inputs_mod, ONLY: i_pc2_conv_coupling,                        &
                            starticeTKelvin, alliceTdegC

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

INTEGER,INTENT(IN) :: k             ! present model layer
INTEGER,INTENT(IN) :: npnts         ! Number of points

REAL,INTENT(IN) :: qclek(npnts)     ! Env. qcl in layer k (kg/kg)
REAL,INTENT(IN) :: qclekp1(npnts)   ! Env. qcl in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qcfek(npnts)     ! Env. qcf in layer k (kg/kg)
REAL,INTENT(IN) :: qcfekp1(npnts)   ! Env. qcf in layer k+1 (kg/kg)
REAL,INTENT(IN) :: cflek(npnts)     ! Env. liquid cloud volume fraction 
                                    ! in layer k
REAL,INTENT(IN) :: cflekp1(npnts)   ! Env. liquid cloud volume fraction 
                                    ! in layer k+1
REAL,INTENT(IN) :: cffek(npnts)     ! Env. frozen cloud volume fraction 
                                    ! in layer k
REAL,INTENT(IN) :: cffekp1(npnts)   ! Env. frozen cloud volume fraction 
                                    ! in layer k+1
REAL,INTENT(IN) :: bcfek(npnts)     ! Env. total cloud volume fraction
                                    ! in layer k
REAL,INTENT(IN) :: bcfekp1(npnts)   ! Env. total cloud volume fraction
                                    ! in layer k+1
REAL,INTENT(IN) :: qclpk(npnts)     ! Par. qcl in layer k (kg/kg)
REAL,INTENT(IN) :: qclpkp1(npnts)   ! Par. qcl in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qcfpk(npnts)     ! Par. qcf in layer k (kg/kg)
REAL,INTENT(IN) :: qcfpkp1(npnts)   ! Par. qcf in layer k+1 (kg/kg)
REAL,INTENT(IN) :: dqclek(npnts)    ! Increment to qcl in layer k (kg/kg/s)
REAL,INTENT(IN) :: dqcfek(npnts)    ! Increment to qcf in layer k (kg/kg/s)
REAL,INTENT(IN) :: dqclekp1(npnts)  ! Increment to qcl in layer k+1 (kg/kg/s)
REAL,INTENT(IN) :: dqcfekp1(npnts)  ! Increment to qcf in layer k+1 (kg/kg/s)

LOGICAL,INTENT(IN) :: l_q_interact  ! True if PC2 is switched on
LOGICAL,INTENT(IN) :: bterm(npnts)  ! Mask for parcels which terminate 
                                    ! in layer k+1

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
! Convection increments to model fields at level k
REAL,INTENT(OUT) :: dcflek(npnts)   ! Increment to liquid cloud volume fraction
                                    ! in layer k
REAL,INTENT(OUT) :: dcffek(npnts)   ! Increment to frozen cloud volume fraction
                                    ! in layer k
REAL,INTENT(OUT) :: dbcfek(npnts)   ! Increment to total cloud volume fraction
                                    ! in layer k
! Convection increments to model fields at level k+1
REAL,INTENT(OUT) :: dcflekp1(npnts) ! Increment to liquid cloud volume fraction
                                    ! in layer k+1
REAL,INTENT(OUT) :: dcffekp1(npnts) ! Increment to frozen cloud volume fraction
                                    ! in layer k+1
REAL,INTENT(OUT) :: dbcfekp1(npnts) ! Increment to total cloud volume fraction
                                    ! in layer k+1

!-----------------------------------------------------------------------
! Variables that are defined locally
!-----------------------------------------------------------------------

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER :: i              ! loop counter 

REAL :: denom             ! Denominator in cloud increment calculation

! Parameters
REAL, PARAMETER :: ls0        = 5.0e-5  ! Minimum value for ls - l

IF (lhook) CALL dr_hook('PC2_ENVIRON',zhook_in,zhook_handle)

IF (l_q_interact) THEN

  DO i=1,npnts
!-----------------------------------------------------------------------
!   Increment to Liquid Cloud Volume Fraction
!-----------------------------------------------------------------------
    ! Increment at level k
    denom = MAX( (qclpk(i)-qclek(i) ), ls0)

    dcflek(i) = eff_dcfl * dqclek(i)*(1.0 - cflek(i)) / denom

    ! Increment at level k+1
    ! For efficiency only calc at termination. This test could be removed as
    ! the k+1 increments will always be overwritten by the k'th increments
    ! on the next iteration of the level loop except at termination.
    IF (bterm(i)) THEN
      denom = MAX( (qclpkp1(i)-qclekp1(i) ), ls0)

      dcflekp1(i) = eff_dcfl * dqclekp1(i)*(1.0 - cflekp1(i)) / denom
    ELSE
      dcflekp1(i) = 0.0
    END IF
    
!-----------------------------------------------------------------------
!   Increment to Frozen Cloud Volume Fraction
!-----------------------------------------------------------------------
    ! Increment at level k
    denom = MAX( (qcfpk(i)-qcfek(i) ), ls0)

    dcffek(i) = eff_dcff * dqcfek(i)*(1.0 - cffek(i)) / denom

    ! Increment at level k+1
    IF (bterm(i)) THEN
      denom = MAX( (qcfpkp1(i)-qcfekp1(i) ), ls0)

      dcffekp1(i) = eff_dcff * dqcfekp1(i)*(1.0 - cffekp1(i)) / denom
    ELSE
      dcffekp1(i) = 0.0
    END IF

!
!-----------------------------------------------------------------------
!   Increment to Total Cloud Volume Fraction
!-----------------------------------------------------------------------
    ! Increment at level k
    dbcfek(i) = dcflek(i) + dcffek(i)

    ! Increment at level k+1
    IF (bterm(i)) THEN
      dbcfekp1(i) = dcflekp1(i) + dcffekp1(i)
    ELSE
      dbcfekp1(i) = 0.0
    END IF
    
  END DO  !loop over npnts

ELSE
!-----------------------------------------------------------------------
!   Not PC2 so set all cloud fraction increments to zero
!-----------------------------------------------------------------------
  DO i=1,npnts
    dcflek(i)     = 0.0
    dcflekp1(i)   = 0.0
    dcffek(i)     = 0.0
    dcffekp1(i)   = 0.0
    dbcfek(i)     = 0.0
    dbcfekp1(i)   = 0.0
  END DO
  
END IF !l_q_interact

IF (lhook) CALL dr_hook('PC2_ENVIRON',zhook_out,zhook_handle)
RETURN
END SUBROUTINE pc2_environ
END MODULE pc2_environ_mod
