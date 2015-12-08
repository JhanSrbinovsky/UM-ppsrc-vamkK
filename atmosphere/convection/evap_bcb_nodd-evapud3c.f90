! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate convective precipitation reaching the surface
!
! Subroutine Interface:
SUBROUTINE evap_bcb_nodd(npnts, n_nodd, kct, iccb, index1,                &
                         bwater, timestep,                                &
                         pstar, p_layer_centres, p_layer_boundaries,      &
                         exner_layer_centres, exner_layer_boundaries,     &
                         the, qe, qse, cca,                               &
                         precip, dthbydt, dqbydt,                         &
                         rain, snow, rain_3d, snow_3d)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE cv_run_mod,  ONLY:                                            &
          i_convection_vn,                                              &
          i_convection_vn_4a,                                           &
          i_convection_vn_5a,                                           &
          i_convection_vn_6a

USE ereport_mod, ONLY : ereport

IMPLICIT NONE

! Description : To calculate the convective precipitation reaching the
!               surface.
!
!     Method : The evaporation below cloud base follows that done in
!              the downdraught routine for the environmental part of the
!              column. The points which are gathered here are those
!              points which have an updraught, but no downdraught.
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
 ,n_nodd               &  ! Compressed vector length for calculation 
 ,kct                     ! Convective cloud top layer

INTEGER, INTENT(IN) :: &
  iccb(npnts)       & ! Cloud base level
 ,index1(npnts)          ! Index of downdraughts not possible

LOGICAL, INTENT(IN) :: &
  bwater(npnts,2:kct+1)  ! mask for those points at which condensate is 
                         ! water in layer k

REAL, INTENT(IN) ::    &
  timestep               ! Timestep

REAL, INTENT(IN) ::                      &
  pstar(npnts)                        & ! Surface pressure (Pa)

 ,p_layer_centres(npnts,0:kct+1)      & ! pressure at layer centre (Pa)

 ,p_layer_boundaries(npnts,0:kct+1)   & ! pressure at layer boundaries (Pa)

 ,exner_layer_centres(npnts,0:kct+1)  & ! Exner function at layer centres
                                           ! starting at level k-1/2
 ,exner_layer_boundaries(npnts,0:kct+1)&! Exner function at layer boundaries
                                           ! starting at level k-1/2
 ,the(npnts,kct+1)                    & ! Model environmental potential
                                           ! temperature (K)
 ,qe(npnts,kct+1)                     & ! Environmental mixing ratio (kg/kg)
 ,qse(npnts,kct+1)                    & ! Environmental qsat (kg/kg)
 ,cca(npnts)                            ! Convective cloud amount

REAL, INTENT(INOUT) ::    &
  precip(npnts,kct+1)     & ! Precipitation added when descending from layer
                            ! K to K-1(kg/m**2/s)
 ,dthbydt(npnts,kct+1) & ! IN Increment to model potential temperature (K/s)
                            ! OUT Updated increment to model potential
                            !     temperature (K/s)
 ,dqbydt(npnts,kct+1)    ! IN Increment to model mixing ratio (kg/kg/s)
                            ! OUT Updated increment to model
                            !     mixing ratio (kg/kg/s)

REAL, INTENT(OUT) ::      &
  rain(npnts)          & ! Rainfall at surface (kg/m**2/s)
 ,snow(npnts)          & ! Snowfall at surface (kg/m**2/s)
 ,rain_3d(npnts,kct+1) & ! Rainfall flux (kg/m**2/s) 
 ,snow_3d(npnts,kct+1)   ! Snowfall flux  (kg/m**2/s) 

!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER ::       &
  i,k               ! Loop counters

INTEGER ::       &
  iccb_c(n_nodd)    ! Compressed cloud base level

LOGICAL ::       &
  bwater_k_c(n_nodd) ! Compressed mask for those points at which condensate
                     ! is water in layer k

REAL ::                &
  exner_km12_c(n_nodd) & ! Compressed Exner functions at layer k
 ,exner_kp12_c(n_nodd) & ! Compressed Exner functions at layer k+1
 ,exner_km32_c(n_nodd) & ! Compressed Exner functions at layer k-1
 ,pk(n_nodd)           & ! Pressure of layer k (Pa)
 ,pkm1_c(n_nodd)       & ! Pressure of layer k-1 (Pa)
 ,exk_c(n_nodd)        & ! Exner ratio for layer k
 ,exkm1_c(n_nodd)      & ! Exner ratio for layer k-1
 ,delpkm1_c(n_nodd)    & ! Pressure difference across layer k-1 (Pa)
 ,pstar_c(n_nodd)        ! Compressed surface pressure (Pa)

REAL ::                &
  precip_k_c(n_nodd)   & ! Compressed precipitation added when descending from
                         ! layer k to k-1 (kg/m**2/s)
 ,qe_k_c(n_nodd)       & ! Compressed parcel mixing ratio of layer k (kg/kg)

 ,the_k_c(n_nodd)      & ! Compressed parcel potential temperature of layer k
                         !  (K)
 ,qe_km1_c(n_nodd)     & ! Compressed parcel mixing ratio of layer k-1 (kg/kg)

 ,qse_km1_c(n_nodd)    & ! Compressed qsat environment of layer k-1 (kg/kg)
 ,the_km1_c(n_nodd)    & ! Compressed parcel potential temperature of layer 
                         !  k-1 (K)
 ,dthbydt_km1_c(n_nodd)& ! Compressed increment to model potential temperature
                         ! of layer k-1 (K/s)
 ,dqbydt_km1_c(n_nodd)   ! Compressed increment to model  mixing ratio
                         ! of layer k-1 (kg/kg/s)

REAL ::                &
  rain_c(n_nodd)       & ! Amount of rainfall passing through environment
                         ! (kg/m**2/s)
 ,snow_c(n_nodd)       & ! Amount of ssnowfall passing through environment
                         ! (kg/m**2/s)
 ,cca_c(n_nodd)          ! Compressed convective cloud amount


INTEGER ::       & 
  errorstatus
  
CHARACTER (LEN=80) ::  cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('EVAP_BCB_NODD',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Loop over levels moving downwards from cloud top 
!-----------------------------------------------------------------------
DO k = kct+1,2,-1

  DO i=1,n_nodd
    the_k_c(i)   = the(index1(i),k)
    the_km1_c(i) = the(index1(i),k-1)
    qe_k_c(i)    = qe(index1(i),k)
    qe_km1_c(i)  = qe(index1(i),k-1)
    qse_km1_c(i)  = qse(index1(i),k-1)
    dthbydt_km1_c(i) = dthbydt(index1(i),k-1)
    dqbydt_km1_c(i)  = dqbydt(index1(i),k-1)
    exner_km12_c(i) = exner_layer_boundaries(index1(i),k-1)
    exner_kp12_c(i) = exner_layer_boundaries(index1(i),k)
    exner_km32_c(i) = exner_layer_boundaries(index1(i),k-2)
    precip_k_c(i)   = precip(index1(i),k)
  END DO
  IF (k == kct+1) THEN
    DO i=1,n_nodd
      pstar_c(i) = pstar(index1(i))
      rain_c(i) = 0.0
      snow_c(i) = 0.0
      iccb_c(i) = iccb(index1(i))
      cca_c(i) = cca(index1(i))
    END DO
  END IF
  IF (k == kct+1) THEN
    DO i=1,n_nodd
      exk_c(i) = exner_layer_centres(index1(i),k)
    END DO
  ELSE
    DO i=1,n_nodd
      exk_c(i) = exkm1_c(i)
    END DO
  END IF
  DO i=1,n_nodd
    pkm1_c(i)    = p_layer_centres(index1(i),k-1)
    delpkm1_c(i) = p_layer_boundaries(index1(i),k-2) -                      &
                                         p_layer_boundaries(index1(i),k-1)
    exkm1_c(i)   = exner_layer_centres(index1(i),k-1)
  END DO

  DO i=1,n_nodd
    IF (bwater(index1(i),k)) THEN
      rain_c(i) = rain_c(i) + precip(index1(i),k)
    ELSE
      snow_c(i) = snow_c(i) + precip(index1(i),k)
    END IF
  END DO

!----------------------------------------------------------------------
! Carry out change of phase calculation for precipitation falling
! through environment
!----------------------------------------------------------------------

! DEPENDS ON: chg_phse
  CALL chg_phse (n_nodd,k,rain_c,snow_c,dthbydt_km1_c,            &
                 exk_c,exkm1_c,delpkm1_c,the_k_c,the_km1_c,       &
                 timestep,cca_c)

!----------------------------------------------------------------------
! Reset precipitation falling through environment if downdraught 
! terminates
!----------------------------------------------------------------------

  SELECT CASE ( i_convection_vn )
    CASE ( i_convection_vn_4a )

! DEPENDS ON: pevp_bcb_4a
  CALL pevp_bcb_4a (n_nodd,k-1,iccb_c,the_km1_c,pkm1_c,qe_km1_c,qse_km1_c,  &
                 delpkm1_c,rain_c,snow_c,dthbydt_km1_c,                  &
                 dqbydt_km1_c,exkm1_c,timestep,cca_c)

    CASE ( i_convection_vn_5a )
    
! DEPENDS ON: pevp_bcb
  CALL pevp_bcb (n_nodd,k-1,iccb_c,the_km1_c,pkm1_c,qe_km1_c,qse_km1_c,  &
                 delpkm1_c,rain_c,snow_c,dthbydt_km1_c,                  &
                 dqbydt_km1_c,exkm1_c,timestep,cca_c)    

    CASE ( i_convection_vn_6a )
    
! DEPENDS ON: pevp_bcb
  CALL pevp_bcb (n_nodd,k-1,iccb_c,the_km1_c,pkm1_c,qe_km1_c,qse_km1_c,  &
                 delpkm1_c,rain_c,snow_c,dthbydt_km1_c,                  &
                 dqbydt_km1_c,exkm1_c,timestep,cca_c)    

    CASE DEFAULT ! i_convection_vn

     errorstatus = 10
     WRITE (cmessage,'(A)') 'Convection scheme version value not valid'
     WRITE (cmessage,'(A,I6)') '   i_convection_vn = ',i_convection_vn
     CALL Ereport ( 'DD_INIT', errorstatus, cmessage)
     
  END SELECT ! i_convection_vn  
                   
  DO i=1,n_nodd
    dthbydt(index1(i),k-1) = dthbydt_km1_c(i)
    dqbydt(index1(i),k-1)  = dqbydt_km1_c(i)

    ! Zero precipitation, as is (slyly) done in downd3c
    precip(index1(i),k) = 0.0
  END DO

  IF (k == 2) THEN
    DO i=1,n_nodd
      rain(index1(i)) = rain(index1(i)) + rain_c(i)
      snow(index1(i)) = snow(index1(i)) + snow_c(i)
    END DO
  END IF

  DO i=1,n_nodd
    rain_3d(index1(i), k-1) = rain_3d(index1(i), k-1) + rain_c(i)
    snow_3d(index1(i), k-1) = snow_3d(index1(i), k-1) + snow_c(i)
  END DO


END DO      !  MAIN LOOP OVER LEVELS

!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('EVAP_BCB_NODD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE evap_bcb_nodd
