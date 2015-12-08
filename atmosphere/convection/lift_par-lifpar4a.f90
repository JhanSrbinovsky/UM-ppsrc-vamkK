! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Lifts Parcel from layer k to k+1
!
SUBROUTINE lift_par_4a (npnts,np_full,                               &
                        thpkp1,qpkp1,xsqkp1,bgmkp1,bwkp1,bwk,        &
                        thpk,qpk,xpk,thekp1,qekp1,thek,qek,qsekp1,   &
                        qclpkp1, qclpk, qclekp1, qclek, l_q_interact,&
                        qcfpkp1, qcfpk, qcfekp1, qcfek,              &
                        pk,pkp1,exkp1,ekp14,ekp34,l_mom_gk,upkp1,    &
                        vpkp1,upk,vpk,uek,uekp1,vek,vekp1,l_tracer,  &
                        ntra,trapkp1,trapk,traekp1,traek)

USE water_constants_mod, ONLY: lc, lf
USE atmos_constants_mod, ONLY: cp, repsilon, kappa
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE rv_mod, ONLY: rv
IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Lifts Parcel from layer k to k+1
!   Taking environment and moist processes into account
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------

! Subroutine arguments
INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,np_full              & ! Full vector length
 ,ntra                   ! Number of tracer variables

LOGICAL, INTENT(IN) :: & 
  l_tracer             & ! Switch for inclusion of tracers
 ,l_mom_gk             & ! Switch for inclusion of (Gregory-Kershaw) CMT
 ,l_q_interact           ! Switch for overwriting of parcel variables if 
                         ! calculating condensate increments

LOGICAL, INTENT(IN) :: & 
  bwkp1(npnts)         & ! Mask for whether condensate is liquid in layer k+1
 ,bwk(npnts)             ! Mask for whether condensate is liquid in layer k

REAL, INTENT(IN) :: &
  thek(npnts)       & ! Potential temperature of cloud environment 
                      ! in layer k (K)
 ,thekp1(npnts)     & ! Potential temperature of cloud environment
                      ! in layer k+1 (K)
 ,qek(npnts)        & ! Mixing ratio of cloud environment in layer k (kg/kg)
 ,qekp1(npnts)      & ! Mixing ratio of cloud environment in layer k+1 (kg/kg)
 ,qclek(npnts)      & ! Liquid condensate Mixing ratio of cloud
                      ! environment in layer k (kg/kg)
 ,qclekp1(npnts)    & ! Liquid condensate Mixing ratio of cloud
                      ! environment in layer k+1 (kg/kg)
 ,qcfek(npnts)      & ! Frozen condensate Mixing ratio of cloud
                      ! environment in layer k (kg/kg)
 ,qcfekp1(npnts)    & ! Frozen condensate Mixing ratio of cloud
                      ! environment in layer k+1 (kg/kg)
 ,uek(npnts)        & ! U of environment in layer k (m/s)
 ,uekp1(npnts)      & ! U of environment in layer k+1 (m/s)
 ,vek(npnts)        & ! V of environment in layer k (m/s)
 ,vekp1(npnts)        ! V of environment in layer k+1 (m/s)

REAL, INTENT(IN) ::    &
  traek(np_full,ntra)  & ! Tracer content of cloud environment in 
                         ! layer k (kg/kg)
 ,traekp1(np_full,ntra)& ! Tracer content of cloud environment in 
                         !layer k+1 (kg/kg)                                            &
 ,trapk(np_full,ntra)    ! Parcel Tracer content in layer k (kg/kg)

REAL, INTENT(IN) :: &
  qsekp1(npnts)     & ! Saturation Mixing ratio of cloud environment in
                      ! layer k+1 (kg/kg)
 ,thpk(npnts)       & ! Parcel Potential temperature in layer k  (K)
 ,qpk(npnts)        & ! Parcel Mixing ratio in layer k (kg/kg)
 ,xpk(npnts)        & ! Parcel condensate in layer k (kg/kg)
 ,qclpk(npnts)      & ! Parcel Liquid condensate Mixing ratio in layer k (kg/kg)
 ,qcfpk(npnts)      & ! Parcel Frozen condensate Mixing ratio in layer k (kg/kg)
 ,upk(npnts)        & ! Parcel U in layer k (m/s)
 ,vpk(npnts)        & ! Parcel V in layer k (m/s)
 ,pkp1(npnts)       & ! Pressure at level k+1 (Pa)
 ,pk(npnts)         & ! Pressure at centre of layer k
 ,exkp1(npnts)      & ! Exner ratio at mid-point of layer k+1
 ,ekp14(npnts)      & ! entrainment coefficient at level k+1/4 multiplied
                      ! by appropriate layer thickness
 ,ekp34(npnts)        ! entrainment coefficient at level k+3/4 multiplied
                      ! by appropriate layer thickness 


LOGICAL, INTENT(INOUT) :: & 
  bgmkp1(npnts)             ! IN Is Parcel saturated in layer k ?
                            ! OUT Is Parcel saturated in layer k+1 ?

REAL, INTENT(OUT) :: &
  thpkp1(npnts)      & ! Parcel Potential temperature in layer k+1 after
                       ! entrainment and latent heating (K)
 ,qpkp1(npnts)       & ! Parcel Mixing ratio in layer k+1 after
                       ! entrainment and latent heating (kg/kg)
 ,qclpkp1(npnts)     & ! Parcel Liquid condensate Mixing ratio
                       ! in layer k+1 after entrainment
 ,qcfpkp1(npnts)     & ! Parcel Frozen condensate Mixing ratio
                       ! in layer k+1 after entrainment
 ,upkp1(npnts)       & ! Parcel U in layer k+1 after entrainment (m/s)
 ,vpkp1(npnts)       & ! Parcel V in layer k+1 after entrainment (m/s)
 ,xsqkp1(npnts)        ! Excess Parcel water after lifting from layer k
                       ! to k+1(kg/kg)

REAL, INTENT(OUT) ::    &
  trapkp1(np_full,ntra)   ! Parcel tracer conent in layer k+1 after 
                          ! entrainment. (kg/kg)

!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------
INTEGER ::   &
  i,ktra         ! Loop counters

REAL ::           &
  thpkp1t(npnts)  & ! Initial estimate of Parcel temperature
                    ! in layer k+1 after entrainment (K)
 ,tt(npnts)       & ! Temporary temperature used in calculation
                    ! of saturation Mixing ratio (K)
 ,pt(npnts)       & ! Temporary pressure used in calculation
                    ! of saturation Mixing ratio (K)
 ,qspkp1(npnts)     ! saturation Mixing ratio of parcel
                    ! after dry ascent (kg/kg)

REAL ::     &
  dqsdt     & ! Rate of change of qsat with temperature
 ,el        & ! Latent heat of gas-to-whatever-condenses PC2 defn (J/kg)
 ,repss     & ! 1 / EPSS(K)
 ,lfbycp


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

PARAMETER ( lfbycp = lf / cp )
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('LIFT_PAR_4A',zhook_in,zhook_handle)
DO i=1,npnts

  ! ----------------------------------------------------------------------
  !   Lift parcel Mixing ratio, potential temperature, U, V and tracer
  !   to the next level
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  !   Initial 'dry' ascent
  !
  !   UM Documentation paper P27
  !   section (3), equations (11B), (12B)
  ! ----------------------------------------------------------------------

  repss = 1.0 / ((1.0+ekp14(i))*(1.0+ekp34(i)))

  thpkp1(i) = ( thpk(i)+ ekp14(i)*thek(i) + ekp34(i)*(1.0+ekp14(i))*thekp1(i)) &
               / ((1.0+ekp14(i))*(1.0+ekp34(i)))

  qpkp1(i) = (  qpk(i) + ekp14(i)*qek(i) + ekp34(i)*(1.0+ekp14(i))*qekp1(i) )  &
               / ((1.0+ekp14(i))*(1.0+ekp34(i)))

  IF (l_q_interact) THEN

    qclpkp1(i)= ( qclpk(i)                                                   &
                  + ekp14(i)*qclek(i) + ekp34(i)*(1.+ekp14(i))*qclekp1(i)    &
                    ) * repss

    qcfpkp1(i)= ( qcfpk(i)                                                   &
                  + ekp14(i)*qcfek(i) + ekp34(i)*(1.+ekp14(i))*qcfekp1(i)    &
                    ) * repss

    ! ----------------------------------------------------------------------
    ! IMPORTANT NOTE: The condensate term has the same form as those for the
    ! vapour and temperature in the initial ascent. If we had a clear idea
    ! on how to parametrize detrainment and mixing of condensate, it could
    ! be updated alongside the vapour process-by-process within convection.
    ! HOWEVER, we opt instead simply to STORE the humidity excess over Qsat
    ! in XSQKP1 as the vapour is updated and apply moisture conservation
    ! just before the increments to the environment variables are calculated
    ! (ie. XSQKP1 is added to QCLPKP1 or QCFPKP1 in CLOUDW1A).
    ! ----------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    !       Currently mixed phase parcel is forbidden. Melt or freeze the
    !       entrained layer cloud and adjust parcel temperature accordingly.
    ! ----------------------------------------------------------------------

    IF (bwkp1(i)) THEN
      qclpkp1(i) = qclpkp1(i) + qcfpkp1(i)
      thpkp1(i)  = ( thpkp1(i) - ( MAX( 0.0, qcfpkp1(i) ) * lf /  &
                                          (cp * exkp1(i)) ) )
      qcfpkp1(i) = 0.0
    ELSE

      qcfpkp1(i) = qclpkp1(i) + qcfpkp1(i)
      thpkp1(i)  = ( thpkp1(i) + ( MAX( 0.0, qclpkp1(i) ) * lf /  &
                                          (cp * exkp1(i)) ) )
      qclpkp1(i) = 0.0
    END IF

  ELSE ! L_q_interact_if1


    ! Allow for release of latent heat of fusion if crossing the freezing
    ! level.
    IF ( bwkp1(i) .AND. .NOT. bwk(i) ) THEN
      thpkp1(i) = thpkp1(i) - lfbycp * xpk(i) / exkp1(i)
    ELSE IF ( bwk(i) .AND. .NOT. bwkp1(i) ) THEN
      thpkp1(i) = thpkp1(i) + lfbycp * xpk(i) / exkp1(i)
    END IF

  END IF  ! L_q_interact_if1

END DO   ! i

IF (l_mom_gk) THEN   ! l_mom_gk set according to type of CMT required
                     ! This code does Gregory-Kershaw CMT

  DO i=1,npnts

    ! Calculations for cumulus convection done after termination

    upkp1(i) = (upk(i)+                                             &
                ekp14(i)*uek(i) + ekp34(i)*(1.0+ekp14(i))*uekp1(i)  &
                  ) / ((1.0+ekp14(i))*(1.0+ekp34(i)))

    vpkp1(i) = (vpk(i)+                                             &
                ekp14(i)*vek(i) + ekp34(i)*(1.0+ekp14(i))*vekp1(i)  &
                  ) / ((1.0+ekp14(i))*(1.0+ekp34(i)))
    upkp1(i) = upkp1(i) - (0.7*(uek(i)-uekp1(i))/(1.0+ekp34(i)))
    vpkp1(i) = vpkp1(i) - (0.7*(vek(i)-vekp1(i))/(1.0+ekp34(i)))

  END DO

END IF       ! l_mom_gk test

IF (l_tracer) THEN

  DO ktra = 1,ntra
    DO i = 1,npnts

      trapkp1(i,ktra) = ( trapk(i,ktra)+ ekp14(i)*traek(i,ktra)              &
                                   + ekp34(i)*(1.0+ekp14(i))*traekp1(i,ktra) &
                              ) / ((1.0+ekp14(i))*(1.0+ekp34(i)))
  
    END DO
  END DO

END IF

!-----------------------------------------------------------------------
!   Calculate where the parcel is supersaturated (i.e.where gamma(k+1)=1
!   See DCTN 29 Page 123)
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Convert potential temperature to temperature and calculate
! pressure of layer k for calculation of saturated mixing ratio
!-----------------------------------------------------------------------

DO i=1,npnts

  ! For Parcel saturated in layer k, use Exner pressure in k+1 with
  ! theta for layer k (gives partly compensating errors)
  ! For unsaturated Parcel in layer k, use the k+1 to test whether
  ! saturation is occurring there

  IF ( bgmkp1(i) ) THEN
    tt(i) = thpk(i) * exkp1(i)
  ELSE
    tt(i) = thpkp1(i) * exkp1(i)
  END IF
END DO

! DEPENDS ON: qsat
CALL qsat (qspkp1,tt,pkp1,npnts)

DO i=1, npnts

  bgmkp1(i) = bgmkp1(i) .OR. ( qspkp1(i)  <   qpkp1(i) )

  IF ( bgmkp1(i) ) THEN
    ! ----------------------------------------------------------------------
    !  Adjust parcel temperature to account for latent heating
    ! ----------------------------------------------------------------------
    IF ( bwkp1(i) ) THEN
      el = lc
    ELSE
      el = lc + lf
    END IF

    dqsdt = el * qspkp1(i) / ( rv * tt(i)*tt(i) )

    ! Consistent with QSAT calculation above:
    thpkp1t(i) = thpkp1(i) + ( pkp1(i) - pk(i) ) * kappa * tt(i)  &
                * dqsdt  * ( cp * tt(i) - repsilon * el ) /        &
                ( repsilon * exkp1(i)*( cp + el * dqsdt ) * pkp1(i) )

    ! ----------------------------------------------------------------------
    !  Calculate T from theta
    ! ----------------------------------------------------------------------
    tt(i) = thpkp1t(i) * exkp1(i)
  END IF
END DO

! ----------------------------------------------------------------------
!    Calculate a more accurate parcel saturated mixing ratio
! ----------------------------------------------------------------------

! DEPENDS ON: qsat
CALL qsat (qspkp1, tt, pkp1, npnts)

DO i=1,npnts
  IF ( bgmkp1(i) ) THEN
    IF ( bwkp1(i) ) THEN
      el = lc
    ELSE
      el = lc + lf
    END IF
    dqsdt = el * qspkp1(i) / ( rv * tt(i)*tt(i) )
    thpkp1t(i) = ( thpkp1(i) * cp   +                                     &
            ((qpkp1(i) - qspkp1(i))/exkp1(i) + thpkp1t(i) * dqsdt ) * el) &
              /( cp  + dqsdt * el )


    xsqkp1(i) = ( thpkp1t(i) - thpkp1(i) ) * cp * exkp1(i) / el

    bgmkp1(i) = ( xsqkp1(i)  >   0.0 ) .AND. (                            &
                ( qpkp1(i) - xsqkp1(i) )  >   0.85 * qspkp1(i))
  END IF
END DO
DO i=1, npnts
  IF ( bgmkp1(i) ) THEN
    thpkp1(i) = thpkp1t(i)
    qpkp1(i)  = qpkp1(i) - xsqkp1(i)
  ELSE
    xsqkp1(i) = 0.0
  END IF
END DO

IF (lhook) CALL dr_hook('LIFT_PAR_4A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE lift_par_4a
