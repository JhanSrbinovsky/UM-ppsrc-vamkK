! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate convective precipitation reaching the surface
!
! Subroutine Interface:
SUBROUTINE evap_bcb_nodd_all (npnts,n_nodd,klev,kterm             &
                      ,iccb, index1, bwater                       &
                      ,exner_layer_centres,exner_layer_boundaries &
                      ,p_layer_centres, p_layer_boundaries,pstar  &
                      ,timestep , cca, the, qe, qse, precip       &
                      ,dthbydt, dqbydt, rain, snow                &
                      ,rain_3d, snow_3d)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE cv_run_mod,  ONLY:                                            &
          i_convection_vn,                                              &
          i_convection_vn_4a,                                           &
          i_convection_vn_5a,                                           &
          i_convection_vn_6a

USE ereport_mod, ONLY : ereport

IMPLICIT NONE

!  Description : To calculate the convective precipitation reaching the
!               surface.
!
!     Method : the evaporation below cloud base follows that done in
!              the downdraught routine for the environmental part of the
!              column. the points which are gathered here are those
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
 ,klev                    ! Number of levels (may be model levels-1 or a 
                          ! reduced set required here)

INTEGER, INTENT(IN) :: &
  kterm(npnts)         &  ! Convective cloud top layer
 ,iccb(npnts)          &  ! Cloud base level
 ,index1(npnts)           ! Index of points where downdraughts not possible

LOGICAL, INTENT(IN) :: &
  bwater(npnts,2:klev+1)   ! Mask for points at which condensate is liquid
 
REAL, INTENT(IN) ::                      &
  exner_layer_centres(npnts,0:klev+1)    & ! Exner function at layer centres
                                           ! starting at level k-1/2
 ,exner_layer_boundaries(npnts,0:klev+1) & ! Exner function at layer boundaries
                                           ! starting at level k-1/2
 ,p_layer_centres(npnts,0:klev+1)        & ! pressure at layer centre (Pa)

 ,p_layer_boundaries(npnts,0:klev+1)     & ! pressure at layer boundaries (Pa)

 ,pstar(npnts)                             ! Surface pressure (Pa)

REAL, INTENT(IN) ::    &
  timestep               ! timestep

REAL, INTENT(IN) ::    &
  cca(npnts)           & ! 2d convective cloud amount
 ,the(npnts,klev+1)    & ! Model enviromental potential temperature (K)
 ,qe(npnts,klev+1)     & ! Model enviromental mixing ratio (kg/kg)
 ,qse(npnts,klev+1)      ! Model enviromental qsat mixing ratio (kg/kg)

REAL, INTENT(INOUT) ::  &
  precip(npnts,klev+1)  & ! precipitation added when descending from layer
                          ! k to k-1 (kg/m**2/s)
 ,dthbydt(npnts,klev+1) & ! increment to model potential temperature (K/s)

 ,dqbydt(npnts,klev+1)  & ! increment to model mixing ratio (kg/kg/s)

 ,rain(npnts)           & ! rainfall at surface (kg/m**2/s)

 ,snow(npnts)           & ! snowfall at surface (kg/m**2/s)

 ,rain_3d(npnts,klev+1) & ! rainfall profile  (kg/m**2/s)

 ,snow_3d(npnts,klev+1)   ! snowfall profile  (kg/m**2/s) 

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

INTEGER ::        &
  i,k             & ! Loop counters
 ,nsofar            ! points in level

INTEGER ::        &
  index2(n_nodd)  & ! compress index for each level
 ,iccb_c(n_nodd)    ! Compressed cloud base level

REAL ::                &
  exner_km12_c(n_nodd) & ! Compressed exner function at layer K
 ,exner_kp12_c(n_nodd) & ! Compressed exner function at layer K+1
 ,exner_km32_c(n_nodd) & ! Compressed exner function at layer K-1
 ,pkm1_c(n_nodd)       & ! pressure of layer k-1 (Pa)
 ,exk_c(n_nodd)        & ! exner ratio for layer k
 ,exkm1_c(n_nodd)      & ! exner ratio for layer k-1
 ,delpkm1_c(n_nodd)    & ! Pressure difference across layer K-1 (Pa)
 ,pstar_c(n_nodd)        ! Compressed surface pressure (Pa)


REAL ::                &
  precip_k_c(n_nodd)   & ! Compressed precipitation added when descending from 
                         ! layer k to k-1 (kg/m**2/s)
 ,qe_k_c(n_nodd)       & ! Compressed parcel mixing ratio of layer K (kg/kg)
 ,qe_km1_c(n_nodd)     & ! Compressed parcel mixing ratio of layer K-1 (kg/kg)
 ,qse_km1_c(n_nodd)    & ! Compressed qsat environment of layer K-1 (kg/kg)
 ,the_k_c(n_nodd)      & ! Compressed parcel potential temperature of
                         ! layer k (K)
 ,the_km1_c(n_nodd)    & ! Compressed parcel potential temperature of 
                         ! layer k-1 (K)
 ,dthbydt_km1_c(n_nodd)& ! Compressed increment to model potential temperature 
                         ! of layer k-1 (K/s)
 ,dqbydt_km1_c(n_nodd) & ! Compressed increment to model mixing ratio of
                         ! layer k-1 (kg/kg/s)
 ,rain_c(n_nodd)       & ! Amount of rainfall passing through environment 
                         ! (kg/m**2/s)
 ,snow_c(n_nodd)       & ! Amount of snowfall passing through environment
                         ! (kg/m**2/s)
 ,cca_c(n_nodd)          ! Compressed convective cloud amount
                                    
                                    
                                                                 !
INTEGER ::       & 
  errorstatus
  
CHARACTER (LEN=80) ::  cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!----------------------------------------------------------------------
! Loop over levels working downwards from maximum termination level + 1
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('EVAP_BCB_NODD_ALL',zhook_in,zhook_handle)

DO k = klev+1,2,-1

! How many points terminated at or above this layer ?
! Only work on these points. For the bottom levels nsofar should
! equal n_nodd

  nsofar = 0
  DO i=1,n_nodd
    IF (kterm(index1(i))+1 >= k) THEN
      nsofar = nsofar + 1
      index2(nsofar) = index1(i)
    END IF
  END DO

! Compress to required points

  DO i = 1, nsofar
    the_k_c(i)   = the(index2(i),k)
    the_km1_c(i) = the(index2(i),k-1)
    qe_k_c(i)    = qe(index2(i),k)
    qe_km1_c(i)  = qe(index2(i),k-1)
    qse_km1_c(i)  = qse(index2(i),k-1)
    dthbydt_km1_c(i)  = dthbydt(index2(i),k-1)
    dqbydt_km1_c(i)   = dqbydt(index2(i),k-1)
    exner_km12_c(i)   = exner_layer_boundaries(index2(i),k-1)
    exner_kp12_c(i)   = exner_layer_boundaries(index2(i),k)
    exner_km32_c(i)   = exner_layer_boundaries(index2(i),k-2)
    precip_k_c(i)   = precip(index2(i),k)
    pstar_c(i) = pstar(index2(i))
    rain_c(i)  = rain(index2(i))
    snow_c(i)  = snow(index2(i))
    iccb_c(i)  = iccb(index2(i))
    cca_c(i)   = cca(index2(i))
    pkm1_c(i)    = p_layer_centres(index2(i),k-1)
    delpkm1_c(i) = p_layer_boundaries(index2(i),k-2) -            &
                          p_layer_boundaries(index2(i),k-1)
    exk_c(i)   = exner_layer_centres(index2(i),k)
    exkm1_c(i) = exner_layer_centres(index2(i),k-1)
  END DO

  DO i = 1,nsofar
    IF (bwater(index2(i),k)) THEN
      rain_c(i) = rain_c(i) + precip(index2(i),k)
    ELSE
      snow_c(i) = snow_c(i) + precip(index2(i),k)
    END IF
  END DO

!----------------------------------------------------------------------
! Carry out change of phase calculation for precipitation falling
! through environment
!----------------------------------------------------------------------

! DEPENDS ON: chg_phse
  CALL chg_phse (nsofar,k,rain_c,snow_c,dthbydt_km1_c,            &
                 exk_c,exkm1_c,delpkm1_c,the_k_c,the_km1_c,       &
                 timestep,cca_c)

!----------------------------------------------------------------------
! Reset precipitation falling through environment if downdraught
! terminates
!----------------------------------------------------------------------

  SELECT CASE ( i_convection_vn )
    CASE ( i_convection_vn_4a )

! DEPENDS ON: pevp_bcb_4a
  CALL pevp_bcb_4a (nsofar,k-1,iccb_c,the_km1_c,pkm1_c,qe_km1_c,qse_km1_c,  &
                 delpkm1_c,rain_c,snow_c,dthbydt_km1_c,                  &
                 dqbydt_km1_c,exkm1_c,timestep,cca_c)

    CASE ( i_convection_vn_5a )

! DEPENDS ON: pevp_bcb
  CALL pevp_bcb (nsofar,k-1,iccb_c,the_km1_c,pkm1_c,qe_km1_c,qse_km1_c,  &
                 delpkm1_c,rain_c,snow_c,dthbydt_km1_c,                  &
                 dqbydt_km1_c,exkm1_c,timestep,cca_c)

    CASE ( i_convection_vn_6a )

! DEPENDS ON: pevp_bcb
  CALL pevp_bcb (nsofar,k-1,iccb_c,the_km1_c,pkm1_c,qe_km1_c,qse_km1_c,  &
                 delpkm1_c,rain_c,snow_c,dthbydt_km1_c,                  &
                 dqbydt_km1_c,exkm1_c,timestep,cca_c)

    CASE DEFAULT ! i_convection_vn

     errorstatus = 10
     WRITE (cmessage,'(A)') 'Convection scheme version value not valid'
     WRITE (cmessage,'(A,I6)') '   i_convection_vn = ',i_convection_vn
     CALL Ereport ( 'DD_INIT', errorstatus, cmessage)
     
  END SELECT ! i_convection_vn  

! Expand output values back into full arrays

  DO i=1,nsofar
    dthbydt(index2(i),k-1) = dthbydt_km1_c(i)
    dqbydt(index2(i),k-1)  = dqbydt_km1_c(i)

! Zero precipitation, as is (slyly) done in downd3c

    precip(index2(i),k) = 0.0
    rain(index2(i)) = rain_c(i)
    snow(index2(i)) = snow_c(i)
  END DO

! Capture 3d rain/snow profiles
  DO i=1,nsofar
    rain_3d(index2(i),k-1) = rain_3d(index2(i),k-1) + rain_c(i)
    snow_3d(index2(i),k-1) = snow_3d(index2(i),k-1) + snow_c(i)
  END DO

END DO      ! main loop over levels
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('EVAP_BCB_NODD_ALL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE evap_bcb_nodd_all
