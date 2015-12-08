! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to initialise the downdraught
!
! Subroutine Interface: (argument list does not yet obey coding standards)
SUBROUTINE dd_init(npnts, np_full                                           &
                  ,th_ud_k, q_ud_k, the_k, qe_k, qse_k, pk, exk             &
                  ,thdd_k, qdd_k, deltd, delqd                              &
                  ,bdd_start, K, bddi, bdd_on                               &
                  ,l_tracer                                                 &
                  ,ntra, tra_ud_k, trae_k, tradd_k, deltrad)


USE atmos_constants_mod, ONLY: c_virtual
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE cv_run_mod,  ONLY:                                            &
          i_convection_vn,                                              &
          i_convection_vn_4a,                                           &
          i_convection_vn_5a,                                           &
          i_convection_vn_6a

USE ereport_mod, ONLY : ereport

IMPLICIT NONE

! 
! Description: Routine to initialise the downdraught
!
! Method: UM documentataion paper 27
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
 ,ntra                 &  ! Number of tracers
 ,k                       ! Present model layer

LOGICAL, INTENT(IN) :: &
  l_tracer                ! Switch for tracers

LOGICAL, INTENT(IN) :: &
  bddi(npnts)          & ! Mask for points where downdraught may initiate
 ,bdd_on(npnts)          ! mask for those points where downdraught is on

REAL, INTENT(IN) ::  &
  th_ud_k(npnts)     & ! Parcel potential temperature of updraught, layer k (K)
 ,q_ud_k(npnts)      & ! Parcel mixing ratio of updraught, layer k (kg/kg)
 ,the_k(npnts)       & ! Potential temperature of environment in layer k (K)
 ,qe_k(npnts)        & ! Mixing ratio of environment in layer k (kg/kg)
 ,qse_k(npnts)       & ! Mixing ratio of environment qsat in layer k (kg/kg)
 ,pk(npnts)          & ! Pressure of layer k (Pa)
 ,exk(npnts)           ! Exner ratio layer k

REAL, INTENT(IN) ::     &
  trae_k(np_full,ntra)  & ! Tracer content of environment in layer k  (kg/kg)
 ,tra_ud_k(np_full,ntra)  ! Parcel tracer content of updraught, layer k (kg/kg)

LOGICAL, INTENT(INOUT) ::  &
  bdd_start(npnts)           ! input mask for those points where DD may start ?
                             ! output set to true if DD started on level k   

REAL, INTENT(OUT) ::     &
  thdd_k(npnts)          & ! Downdraught potential temperature of layer k (K)
 ,qdd_k(npnts)           & ! Downdraught mixing ratio of layer k (kg/kg)
 ,tradd_k(np_full,ntra)  & ! Downdraught tracer content of layer k (kg/kg)
 ,deltd(npnts)           & ! Cooling necessary to achieve saturation
 ,delqd(npnts)           & ! Moistening necessary to achieve saturation
 ,deltrad(np_full,ntra)    ! Depletion of environment tracer due to formation 
                           ! of Downdraught 
!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER ::       &
  i, ktra          ! Loop counters
  
INTEGER ::       & 
  errorstatus
  
CHARACTER (LEN=80) ::  cmessage
  
REAL ::               &
  th_mean(npnts)      & ! Mean potential temperature used in calculation of 
                        ! saturated downdraught potential temperature  in 
                        ! layer k (K)
 ,q_mean(npnts)       & ! Mean mixing ratio used in calculation of 
                        ! saturated downdraught mixing ratio in layer k (kg/kg)
 ,t_mean(npnts)       & ! Mean temperature used in calculation of 
                        ! saturated downdraught potential temperature in 
                        ! layer k (kg/kg)
 ,tra_mean(npnts,ntra)& ! Mean tracer used as initial tracer content of
                        ! downdraught in layer k (kg/kg) 
 ,thdds(npnts)        & ! Saturated downdraught potential temperature in 
                        ! layer k (K) 
 ,qdds(npnts)         & ! Saturated downdraught mixing ratio in 
                        ! layer k (kg/kg) 
 ,buoy(npnts)           ! Buoyancy of parcel in layer k

REAL ::             &
  thdd_v            &  ! Virtual potential temperature of parcel in layer k (K)          
 ,the_v                ! Virtual potential temperature of environment in
                       ! layer k (K)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
! Calculate mean temperature, mixing ratio, U, V  and tracer.
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DD_INIT',zhook_in,zhook_handle)

DO i=1,npnts
  th_mean(i) = (the_k(i)+th_ud_k(i))*0.5
  q_mean(i)  = (qe_k(i)+q_ud_k(i))*0.5
  t_mean(i)  = th_mean(i)*exk(i)
END DO


! For tracers

IF(l_tracer)THEN

   DO ktra=1,ntra
     DO i=1,npnts
       tra_mean(i,ktra) = (trae_k(i,ktra)+tra_ud_k(i,ktra))*0.5
     END DO
   END DO

END IF
 
!-----------------------------------------------------------------------
! Calculate saturated downdraught potential temperature for layer k
!-----------------------------------------------------------------------

SELECT CASE ( i_convection_vn )
  CASE ( i_convection_vn_4a )

! DEPENDS ON: satcal_4a
    CALL satcal_4a(npnts,k,t_mean,th_mean,pk,exk,q_mean,the_k,qse_k,qdds,thdds)
  
  CASE ( i_convection_vn_5a )
  
! DEPENDS ON: satcal
    CALL satcal(npnts,k,t_mean,th_mean,pk,exk,q_mean,the_k,qse_k,qdds,thdds)
 
  CASE ( i_convection_vn_6a )
  
! DEPENDS ON: satcal
    CALL satcal(npnts,k,t_mean,th_mean,pk,exk,q_mean,the_k,qse_k,qdds,thdds) 
 
  CASE DEFAULT ! i_convection_vn

     errorstatus = 10
     WRITE (cmessage,'(A)') 'Convection scheme version value not valid'
     WRITE (cmessage,'(A,I6)') '   i_convection_vn = ',i_convection_vn
     CALL Ereport ( 'DD_INIT', errorstatus, cmessage)

END SELECT ! i_convection_vn      

!-----------------------------------------------------------------------
! Is saturated parcel negatively buoyant compared to environment?
!-----------------------------------------------------------------------

DO I=1,npnts
  IF (.NOT. bdd_on(i) .AND. bddi(i) ) THEN
    thdd_v = thdds(i)*(1.0+c_virtual*qdds(i))
    the_v  = the_k(i)*(1.0+c_virtual*qe_k(i))
    buoy(i) = thdd_v - THE_V
 
    IF (buoy(i)  <   0.5 ) THEN

    !-----------------------------------------------------------------------
    ! Initiate downdraught
    !-----------------------------------------------------------------------

       thdd_k(i) = thdds(i)
       qdd_k(i) = qdds(i)
       bdd_start(i) = .TRUE.
 
    !-----------------------------------------------------------------------
    ! Calculate cooling and moistening to achieve saturation
    !-----------------------------------------------------------------------

       deltd(i) = thdds(i)-the_k(i)
       delqd(i) = qdds(i)-qe_k(i)
     END IF
   END IF
END DO


IF(l_tracer)THEN

  DO ktra=1,ntra
    DO i=1,npnts
      IF(.NOT.bdd_on(i).AND.bddi(i).AND.K >= 4)THEN
        IF(buoy(i) <  0.5)THEN
          tradd_k(i,ktra) = tra_mean(i,ktra)
          deltrad(i,ktra) = tradd_k(i,ktra)-trae_k(i,ktra)
        END IF
      END IF
    END DO
  END DO

END IF

IF (lhook) CALL dr_hook('DD_INIT',zhook_out,zhook_handle)
RETURN

END SUBROUTINE dd_init

