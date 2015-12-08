! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate the effect of the downdraught on the large-scale atmosphere
!
! Subroutine Interface:  
!
SUBROUTINE dd_env(npnts, np_full, ntra                                     &
                  ,l_tracer, b_dd_end, bdd_start, bdd_on                   &
                  ,thdd_k, thdd_km1, qdd_k, qdd_km1, the_k, the_km1        &
                  ,qe_k, qe_km1, flx_dd_k, flx_dd_km1, delpk, delpkm1      &
                  ,deltd, delqd, amdetk, ekm14                             &
                  ,tradd_k, tradd_km1, trae_k, trae_km1, deltrad           &
                  ,dthbydt_k, dthbydt_km1, dqbydt_k, dqbydt_km1            &
                  ,dtrabydt_k, dtrabydt_km1 )

USE yomhook, ONLY: lhook, dr_hook  
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! 
! Description: Downdraught routine
!              Calculate the effect of the downdraught on the 
!              large-scale atmosphere
!
! Method: UM documentation paper 27
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------
INTEGER, INTENT(IN) :: &
  npnts                &  ! Vector length
 ,np_full              &  ! Full vector length
 ,ntra                    ! Number of tracers

LOGICAL, INTENT(IN) :: &
  l_tracer                ! Switch for tracers

LOGICAL, INTENT(IN) :: &
  b_dd_end(npnts)      & ! Mask for those points where downdraught is 
                         ! terminating
 ,bdd_start(npnts)     & ! Mask for those points where downdraught is 
                         ! starting
 ,bdd_on(npnts)          ! Mask for those points where downdraught is on

REAL, INTENT(IN) :: &
  thdd_k(npnts)     & ! Potential temperature of downdraught in layer k (K)
 ,thdd_km1(npnts)   & ! Potential temperature of downdraught in layer k-1 (K)
 ,qdd_k(npnts)      & ! Mixing ratio of downdraught in layer k (kg/kg)  
 ,qdd_km1(npnts)    & ! Mixing ratio of downdraught in layer k-1 (kg/kg)
 ,the_k(npnts)      & ! Potential temperature of environment in layer k (K)
 ,the_km1(npnts)    & ! Potential temperature of environment in layer k-1 (K)
 ,qe_k(npnts)       & ! Mixing ratio of environment in layer k  (kg/kg) 
 ,qe_km1(npnts)     & ! Mixing ratio of environment in layer k-1 (kg/kg)   
 ,flx_dd_k(npnts)   & ! Downdraught mass flux of layer k (Pa/s)
 ,flx_dd_km1(npnts) & ! Downdraught mass flux of layer k-1 (Pa/s)
 ,delpk(npnts)      & ! Change in pressure across layer k (Pa)
 ,delpkm1(npnts)    & ! Change in pressure across layer k-1  (Pa)
 ,deltd(npnts)      & ! Cooling necessary to achieve saturation (K)
 ,delqd(npnts)      & ! moistening necessary to achieve saturation (kg/kg)
 ,amdetk(npnts)     & ! Mixing detrainment rate
 ,ekm14(npnts)        ! Exner ratio at layer k-1/4

REAL, INTENT(IN) ::       &
  tradd_k(np_full,ntra)   & ! Downdraught tracer content of layer k (kg/kg)
 ,tradd_km1(npnts,ntra)   & ! Downdraught tracer content of layer k-1 (kg/kg)
 ,trae_k(np_full,ntra)    & ! Environment tracer content of layer k (kg/kg)
 ,trae_km1(np_full,ntra)  & ! Environment tracer content of layer k-1(kg/kg)
 ,deltrad(npnts,ntra)       ! Depletion of environment tracer due to 
                            ! downdraught formation (kg/kg)

REAL, INTENT(INOUT) :: &
  dthbydt_k(npnts)     & ! In  Increment to potential temperature of layer k
                         ! Out Updated increment potential temperature layer k
                         !           (K/s)
 ,dthbydt_km1(npnts)   & ! In  Increment to potential temperature of layer k-1
                         ! Out Updated increment potential temperature layer k-1
                         !           (K/s)
 ,dqbydt_k(npnts)      & ! In  Increment to mixing ratio of layer k
                         ! Out Updated increment mixing ratio layer k (kg/kg)
 ,dqbydt_km1(npnts)      ! In  Increment to mixing ratio  of layer k-1
                         ! Out Updated increment mixing ratio layer k-1 (kg/kg)

REAL, INTENT(INOUT) ::       &
  dtrabydt_k(np_full,ntra)   & ! In  Increment to tracer content of layer k
                               ! Out Updated increment tracer content layer k 
                               ! (kg/kg/s)
 ,dtrabydt_km1(np_full,ntra)   ! In  Increment to tracer content of layer k-1
                               ! Out Updated increment tracer content layer k-1
                               ! (kg/kg/s)

!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER ::       &
  i, ktra          ! Loop counters

REAL ::          &
  tempry           ! Used in calculations of the effect on the environment

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Calculate the effect on the environment in layer k
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DD_ENV',zhook_in,zhook_handle)

DO i=1,npnts
  IF (bdd_on(i)) THEN

  !-----------------------------------------------------------------------
  ! Subtract the energy used to form the downdraught
  !-----------------------------------------------------------------------

    tempry = flx_dd_k(i)/delpk(i)
    IF (bdd_start(i)) THEN
      dthbydt_k(i) = dthbydt_k(i)-tempry*deltd(i)
      dqbydt_k(i)  = dqbydt_k(i) -tempry*delqd(i)
    END IF

  !-----------------------------------------------------------------------
  ! Effect of convection and downdarught upon potential temperature of 
  ! layer k
  !-----------------------------------------------------------------------

    dthbydt_k(i) = dthbydt_k(i) + tempry * (                                 &

                        ! compensating subsidence term       
                  (1.0+ekm14(i)) * (1.0-amdetk(i)) * (the_km1(i)-the_k(i))   & 

                        ! Mixing detrainment term
                  +  amdetk(i)* (thdd_k(i)-the_k(i)) )

  !-----------------------------------------------------------------------
  ! Effect of convection and downdraught upon mixing ratio of
  ! layer k
  !-----------------------------------------------------------------------

    dqbydt_k(i) = dqbydt_k(i) + tempry * (                                  &

                        ! compensating subsidence term       
                  (1.0+ekm14(i)) * (1.0-amdetk(i)) * (qe_km1(i)-qe_k(i))    &

                        ! Mixing detrainment term
                  +  amdetk(i)* (qdd_k(i)-qe_k(i)) )


  !-----------------------------------------------------------------------
  ! Terminal detrainment and subsidence in terminal layer
  !-----------------------------------------------------------------------

    IF (b_dd_end(i)) THEN
      tempry         = flx_dd_km1(i)/delpkm1(i)

      dthbydt_km1(i) = dthbydt_km1(i)+tempry* (thdd_km1(i)-the_km1(i))  

      dqbydt_km1(i)  = dqbydt_km1(i)+tempry*(qdd_km1(i)-qe_km1(i))
    END IF

  END IF    ! test on bdd_on
END DO      ! loop over npnts

!-----------------------------------------------------------------------
! Effect of convection and downdraught upon tracer content of  layer k
!-----------------------------------------------------------------------

IF(l_tracer)THEN

  DO ktra=1,ntra
    DO i=1,npnts
      IF(bdd_on(i))THEN

        tempry = flx_dd_k(i)/delpk(i)
        IF(bdd_start(i))THEN
          dtrabydt_k(i,ktra) = dtrabydt_k(i,ktra)-tempry*deltrad(i,ktra)
        END IF
        dtrabydt_k(i,ktra) = dtrabydt_k(i,ktra) + tempry * (                  &

                        ! compensating subsidence term       
         (1.0+ekm14(i)) * (1.0-amdetk(i)) * (trae_km1(i,ktra)-trae_k(i,ktra)) &

                        ! Mixing detrainment term
            + amdetk(i)* (tradd_k(i,ktra)-trae_k(i,ktra))  )

        !--------------------------------------------------------------------
        ! Terminal Detrainment of tracer
        !--------------------------------------------------------------------

        IF(b_dd_end(i))THEN
          tempry = flx_dd_km1(i)/delpkm1(i)
          dtrabydt_km1(i,ktra)=dtrabydt_km1(i,ktra)+tempry*                  &
                                   (tradd_km1(i,ktra)-trae_km1(i,ktra))
        END IF

      END IF
    END DO
  END DO

END IF      ! l_tracer

IF (lhook) CALL dr_hook('DD_ENV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE dd_env

