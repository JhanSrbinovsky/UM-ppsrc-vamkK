! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Lifts Parcel from layer k to k+1

MODULE lift_par_6a_mod

IMPLICIT NONE

!
! Description:
!   Lifts Parcel from layer k to k+1 taking account of 
!   entrainment, detrainment, phase changes and condensation
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

SUBROUTINE lift_par_6a (npnts, thek, thekp1,                         &
                        qek, qekp1, qclek, qcfek,                    &
                        qclekp1, qcfekp1,                            &
                        pk, pkp1, exkp1,                             &
                        thpk, qpk, qclpk, qcfpk,                     &
                        ekp14, ekp34,                                &
                        l_q_interact, bwk, bwkp1,                    &
                        !Out
                        bgmkp1, thpkp1, qpkp1,                       &
                        qclpkp1, qcfpkp1,                            &
                        Qlkp1, Qfkp1, Frezkp1)

USE water_constants_mod, ONLY: lc, lf
USE cv_derived_constants_mod, ONLY: ls, lfrcp
USE atmos_constants_mod, ONLY: cp, repsilon, kappa, rv
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Subroutine arguments
! Arguments with INTENT(IN):
INTEGER, INTENT(IN) :: npnts    ! Number of points

! Properties of the cloud environment
REAL,INTENT(IN) :: thek(npnts)      ! Env. pot. temperature in layer k (K)
REAL,INTENT(IN) :: thekp1(npnts)    ! Env. pot. temperature in layer k+1 (K)
REAL,INTENT(IN) :: qek(npnts)       ! Env. mixing ratio of in layer k (kg/kg)
REAL,INTENT(IN) :: qekp1(npnts)     ! Env. mixing ratio of in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qclek(npnts)     ! Env. liquid condensate mixing ratio 
                                    ! in layer k (kg/kg)
REAL,INTENT(IN) :: qcfek(npnts)     ! Env. frozen condensate mixing ratio 
                                    ! in layer k (kg/kg)
REAL,INTENT(IN) :: qclekp1(npnts)   ! Env. liquid condensate mixing ratio 
                                    ! in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qcfekp1(npnts)   ! Env. frozen condensate mixing ratio 
                                    ! in layer k+1 (kg/kg)
REAL,INTENT(IN) :: pk(npnts)        ! Pressure at level k (Pa)
REAL,INTENT(IN) :: pkp1(npnts)      ! Pressure at level k+1 (Pa)
REAL,INTENT(IN) :: exkp1(npnts)     ! Exner ratio at mid-point of layer k+1

! Properties of the parcel at layer k
REAL,INTENT(IN) :: thpk(npnts)      ! Par. pot. temperature in layer k (K)
REAL,INTENT(IN) :: qpk(npnts)       ! Par. mixing ratio of in layer k (kg/kg)
REAL,INTENT(IN) :: qclpk(npnts)     ! Par. liquid condensate mixing ratio 
                                    ! in layer k (kg/kg)
REAL,INTENT(IN) :: qcfpk(npnts)     ! Par. frozen condensate mixing ratio 
                                    ! in layer k (kg/kg)

!Entrainment and detrainment rates
REAL,INTENT(IN) :: ekp14(npnts)     ! entrainment coefficient at level k+1/4 
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: ekp34(npnts)     ! entrainment coefficient at level k+3/4 
                                    ! multiplied by appropriate layer thickness

LOGICAL,INTENT(IN) :: l_q_interact  ! True if PC2 is switched on

! Array  arguments with INTENT(IN):
LOGICAL,INTENT(IN) :: bwk(npnts)    ! Mask for whether condensate is liquid
                                    ! in layer k
LOGICAL,INTENT(IN) :: bwkp1(npnts)  ! Mask for whether condensate is liquid 
                                    ! in layer k+1

! Array  arguments with INTENT(OUT):
! Properties of the parcel at layer k+1 after entrainment, phase changes 
! and condensation
LOGICAL,INTENT(OUT) :: bgmkp1(npnts)! Is Parcel saturated in layer k+1?

REAL,INTENT(OUT) :: thpkp1(npnts)   ! Par. pot. temperature in layer k+1 (K)
REAL,INTENT(OUT) :: qpkp1(npnts)    ! Par. mixing ratio of in layer k+1 (kg/kg)
REAL,INTENT(OUT) :: qclpkp1(npnts)  ! Par. liquid condensate mixing ratio 
                                    ! in layer k+1 (kg/kg)
REAL,INTENT(OUT) :: qcfpkp1(npnts)  ! Par. frozen condensate mixing ratio 
                                    ! in layer k+1 (kg/kg)
! Amount of moist processes
REAL,INTENT(OUT) :: Qlkp1(npnts)    ! Amount of condensation to liquid water 
                                    ! multiplied by the layer thickness (kg/kg)
REAL,INTENT(OUT) :: Qfkp1(npnts)    ! Amount of deposition to ice water 
                                    ! multiplied by the layer thickness (kg/kg)
REAL,INTENT(OUT) :: Frezkp1(npnts)  ! Amount of freezing multiplied 
                                    ! by the layer thickness (kg/kg)

!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------
INTEGER :: i, j, ktra     ! Loop counters

REAL :: tpkp1(npnts)     ! Par. temperature in layer k+1 (K)
REAL :: thpkp1dry(npnts) ! Par. pot. temperature in layer k+1 
                         ! after dry ascent (K)
REAL :: qpkp1dry(npnts)  ! Par. mixing ratio of in layer k+1
                         ! after dry ascent (kg/kg)
REAL :: qspkp1(npnts)    ! Saturation Mixing ratio of parcel
                         ! after dry ascent (kg/kg)
REAL :: dqsdth           ! Rate of change of qsat with potential temperature
REAL :: el               ! Latent heat of gas-to-whatever-condenses PC2 defn
                         ! (J/kg)
REAL :: Factor(npnts)    ! Factor used in update calculation

REAL :: lbycpexner    !L/(cp*exner)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('LIFT_PAR_6A',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!   Initial 'dry' ascent
! ----------------------------------------------------------------------
DO i=1,npnts  
  Factor(i)  = 1.0/((1.0+ekp14(i))*(1.0+ekp34(i)))
  thpkp1(i)  = ( thpk(i)                                               &
               + ekp14(i)*thek(i)                                      &
               + (1.0+ekp14(i))*ekp34(i)*thekp1(i) )                   &
               * Factor(i)

  qpkp1(i)   = ( qpk(i)                                                &
               + ekp14(i)*qek(i)                                       &
               + (1.0+ekp14(i))*ekp34(i)*qekp1(i) )                    &
               * Factor(i)
END DO   ! i

IF (l_q_interact) THEN
  DO i=1,npnts
  !PC2 so entrainment from the environment
    qclpkp1(i) = ( qclpk(i)                                            &
               + ekp14(i)*qclek(i)                                     &
               + (1.0+ekp14(i))*ekp34(i)*qclekp1(i) )                  &
               * Factor(i)
               
    qcfpkp1(i) = ( qcfpk(i)                                            &
               + ekp14(i)*qcfek(i)                                     &
               + (1.0+ekp14(i))*ekp34(i)*qcfekp1(i) )                  &
               * Factor(i)
  END DO  ! i
ELSE
  DO i=1,npnts
  !Not PC2 so no entrainment from the environment
    qclpkp1(i) =  qclpk(i) * Factor(i)

    qcfpkp1(i) =  qcfpk(i) * Factor(i)
  END DO  ! i
END IF ! l_q_interact

! ----------------------------------------------------------------------
!       Currently mixed phase parcel is forbidden. Melt or freeze the
!       entrained layer cloud and adjust parcel temperature accordingly.
! ----------------------------------------------------------------------
DO i=1,npnts
  IF (bwkp1(i)) THEN
    Frezkp1(i) = -qcfpkp1(i)
    qcfpkp1(i) = 0.0
    qclpkp1(i) = qclpkp1(i) - Frezkp1(i)
  ELSE
    Frezkp1(i) = qclpkp1(i)
    qclpkp1(i) = 0.0
    qcfpkp1(i) = qcfpkp1(i) + Frezkp1(i)
  END IF
  thpkp1(i)  = thpkp1(i) + Frezkp1(i) * lfrcp / exkp1(i)

END DO   ! i

DO i=1,npnts
  tpkp1(i)     = thpkp1(i) * exkp1(i)
! Save dry ascent values
  thpkp1dry(i) = thpkp1(i)
  qpkp1dry(i)  = qpkp1(i)
END DO

! DEPENDS ON: qsat
CALL qsat (qspkp1,tpkp1,pkp1,npnts)

! ----------------------------------------------------------------------
!       Calculate theta and q if the parcel was brought to saturation
! ----------------------------------------------------------------------
  
DO j=1,3  !Three iterations should be sufficient

  DO i=1, npnts
      IF ( bwkp1(i) ) THEN
        el=lc
      ELSE
        el=ls
      END IF

      dqsdth = el * qspkp1(i) / ( rv * exkp1(i) * thpkp1(i) * thpkp1(i) )
      lbycpexner=el/(cp*exkp1(i))

      !Calculate the next estimate of the parcel's p.temp after condensation
      thpkp1(i) = ( thpkp1dry(i) + lbycpexner*(qpkp1dry(i) - qspkp1(i)  &
                   + thpkp1(i)*dqsdth) ) /                              &
                   (1.0 + lbycpexner*dqsdth)
      
      !Calculate the next estimate of the parcel's temp after condensation
      tpkp1(i) = thpkp1(i) * exkp1(i)
      
  END DO  !npnts

! Calculate qsat at the next estimate of the parcel's temp after condensation
! DEPENDS ON: qsat
  CALL qsat (qspkp1,tpkp1,pkp1,npnts)

END DO !End iterations


! ----------------------------------------------------------------------
! Update parcel properties to take account of condensation 
! or evaporation
! ----------------------------------------------------------------------

DO i=1,npnts

  IF ( bwkp1(i) ) THEN
    Qlkp1(i) = qpkp1dry(i) - qspkp1(i)
    Qfkp1(i) = 0.0
  ELSE
    Qfkp1(i) = qpkp1dry(i) - qspkp1(i)
    Qlkp1(i) = 0.0
  END IF

! Ql and Qf can be negative if evaporation rather than condensation.
! In this case, limit evaporation of condensate to the amount of 
! available condensate.
  Qlkp1(i)   = MAX(Qlkp1(i), -qclpkp1(i))
  Qfkp1(i)   = MAX(Qfkp1(i), -qcfpkp1(i))

! To prevent evaporation comment out previous two lines and uncomment the
! following two lines
!  Qlkp1(i)   = MAX(Qlkp1(i), 0.0)
!  Qfkp1(i)   = MAX(Qfkp1(i), 0.0)

  thpkp1(i)  = thpkp1dry(i) + (lc*Qlkp1(i) + ls*Qfkp1(i))  &
               / (cp * exkp1(i))
  tpkp1(i)   = thpkp1(i) * exkp1(i)

  qclpkp1(i) = qclpkp1(i) + Qlkp1(i)
  qcfpkp1(i) = qcfpkp1(i) + Qfkp1(i)

  qpkp1(i)   = qpkp1dry(i) - Qlkp1(i) - Qfkp1(i)

END DO  !points

! Calculate qsat at the parcel's temp after condensation
! or evaporation
! DEPENDS ON: qsat
  CALL qsat (qspkp1,tpkp1,pkp1,npnts)

DO i=1,npnts
  bgmkp1(i) = ( qpkp1(i) > 0.9999 * qspkp1(i) )
END DO


IF (lhook) CALL dr_hook('LIFT_PAR_6A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE lift_par_6a
END MODULE lift_par_6a_mod
