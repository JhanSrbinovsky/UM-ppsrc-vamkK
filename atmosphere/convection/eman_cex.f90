! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


!   Compress expand around Emanuel Downdraught code

SUBROUTINE eman_cex(npnts,nmid,kmax_term,nlev,trlev,ntra          &
,                      kterm,l_mid,l_tracer                       &
,                      exner_layer_centres,exner_layer_boundaries &
,                      p,ph,timestep,th,q,qse,tracer,precip       &
,                      dthbydt, dqbydt, dtrabydt                  &
,                      rain, snow ,dwn_flux                       &
                   )

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim 

IMPLICIT NONE

!  Description : Interface to Emanuel downdraught routine
! 
!  Method :
!   First compress to points with mid-level convection.
!   Then call  Emanuel Down draughts.
!   Expand results back to full grid.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards 8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

INTEGER, INTENT(IN) :: &
  npnts                & ! No. of deep convection points
, nmid                 & ! No. of mid level convection points
, nlev                 & ! No. of model layers
, trlev                & ! No. of tracer levels
, ntra                 & ! No. of tracer fields
, kmax_term            & ! highest level reached by convection
, kterm(npnts)           ! level reached by convection


LOGICAL, INTENT(IN) :: &
  l_tracer             & ! true in tracers present
, l_mid(npnts)           ! true  if mid level convection

REAL, INTENT(IN)    ::                 &
  exner_layer_centres(npnts,0:nlev)    & ! Exner
, exner_layer_boundaries(npnts,0:nlev) & ! Exner at half level above
, p(npnts,0:nlev)                      & ! Pressure  (Pa)
, ph(npnts,0:nlev)                     & ! Pressure at half level (Pa)
, timestep                             & ! timestep for model physics in seconds
, qse(npnts,nlev)                      & ! Saturation specific humidity of cloud
                                         ! environment (kg/kg)
, q (npnts,nlev)                       & ! specific humidity of cloud 
                                         ! environment(kg/kg)
, th (npnts,nlev)                      & ! theta of cloud environment(K)
, precip(npnts,nlev)                   & ! Precip from updraught G-R scheme
                                         ! (kg/m2/s) not units for wdtrain
, tracer(npnts,trlev,ntra)               !  Tracer on model levels  (kg/kg)


REAL, INTENT(INOUT)    :: &
  rain(npnts)             & ! rainfall at surface (kg/m**2/s)
, snow(npnts)               ! snowfall at surface (kg/m**2/s)

! increments
REAL, INTENT(INOUT)    ::   &
  dqbydt(npnts,nlev)        & ! increments to q (kg/kg/s)
, dthbydt(npnts,nlev)       & ! increments to potential temperature(K/s)
, dtrabydt(npnts,nlev,ntra)   ! increments to model tracers(kg/kg/s)


! Arguments with intent OUT:

REAL, INTENT(OUT)    ::    &
  dwn_flux(npnts,nlev)       ! Downdraught mass flux (Pa/s) diagnostic

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
! compressed version of most input arrays

INTEGER ::        &
  kterm_c(nmid)      ! compressed termination level

REAL ::                                   &
  rain_c(nmid)                            & ! compressed rain
, snow_c(nmid)                            &
, dwn_flux_c(nmid,nlev)                   &
, dqbydt_c(nmid,nlev)                     &
, dthbydt_c(nmid,nlev)                    &
, dtrabydt_c(nmid,nlev,ntra)              &
, exner_layer_centres_c(nmid,0:nlev)      &
, exner_layer_boundaries_c(nmid,0:nlev)   &
, p_c(nmid,0:nlev)                        &
, ph_c(nmid,0:nlev)                       &
, qse_c(nmid,nlev)                        &
, q_c(nmid,nlev)                          &
, th_c(nmid,nlev)                         &
, precip_c(nmid,nlev)                     &
, tracer_c(nmid,trlev,ntra)

INTEGER ::    &
 dmi(nmid)      ! index of convecting points

INTEGER ::    &
 i ,j ,k        ! loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Compress input arrays
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('EMAN_CEX',zhook_in,zhook_handle)

j=0
DO i = 1,npnts
  IF (l_mid(i)) THEN
    j = j+1
    dmi(j) = i
  END IF
END DO

DO j = 1, nmid
  kterm_c(j) = kterm(dmi(j))
  rain_c(j)  = rain(dmi(j))
  snow_c(j)  = snow(dmi(j))

END DO

DO k = 0, nlev
  DO j = 1, nmid
    exner_layer_centres_c(j,k)   = exner_layer_centres(dmi(j),k)
    exner_layer_boundaries_c(j,k)=exner_layer_boundaries(dmi(j),k)
    p_c(j,k)  = p(dmi(j),k)
    ph_c(j,k) = ph(dmi(j),k)
  END DO
END DO

DO k = 1, nlev
  DO j = 1, nmid
    q_c(j,k)   = q(dmi(j),k)
    th_c(j,k)  = th(dmi(j),k)
    qse_c(j,k) = qse(dmi(j),k)

    precip_c(j,k) = precip(dmi(j),k)

    dqbydt_c(j,k)  = dqbydt(dmi(j),k)
    dthbydt_c(j,k) = dthbydt(dmi(j),k)

  END DO
END DO
IF (l_tracer) THEN
  DO i = 1, ntra
    DO k = 1, trlev
      DO j = 1, nmid
        tracer_c(j,k,i) = tracer(dmi(j),k,i)
      END DO
    END DO
    DO k = 1, nlev
      DO j = 1, nmid
        dtrabydt_c(j,k,i) = dtrabydt(dmi(j),k,i)
      END DO
    END DO
  END DO
END IF

! ----------------------------------------------------------------------
! Call Emanuel scheme for just mid level convecting points
! ----------------------------------------------------------------------
! DEPENDS ON: eman_dd
  CALL eman_dd (nmid,kmax_term,nlev,trlev,ntra,kterm_c,l_tracer   &
,               exner_layer_centres_c,exner_layer_boundaries_c    &
,               p_c,ph_c,timestep,th_c,q_c,qse_c,tracer_c         &
,               precip_c,dthbydt_c,dqbydt_c,dtrabydt_c            &
,               rain_c, snow_c ,dwn_flux_c)

! ----------------------------------------------------------------------
! Expand results
! ----------------------------------------------------------------------
DO i = 1,nmid
  rain(dmi(i)) = rain_c(i)
  snow(dmi(i)) = snow_c(i)
END DO
DO k = 1, nlev
  DO i = 1, nmid
    dwn_flux(dmi(i),k) = dwn_flux_c(i,k)
    dqbydt(dmi(i),k)   = dqbydt_c(i,k)
    dthbydt(dmi(i),k)  = dthbydt_c(i,k)
  END DO
END DO

IF (l_tracer) THEN
  DO j = 1, ntra
    DO k = 1, nlev
      DO i = 1, nmid
        dtrabydt(dmi(i),k,j) = dtrabydt_c(i,k,j)
      END DO
    END DO
  END DO
END IF
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('EMAN_CEX',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eman_cex

