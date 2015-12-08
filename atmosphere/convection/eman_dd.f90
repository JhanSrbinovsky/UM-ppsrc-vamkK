! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


!   Emanuel Downdraught code

SUBROUTINE eman_dd(n_dp,kmax_term,nlev,trlev,ntra,kterm,l_tracer  &
,                      exner_layer_centres,exner_layer_boundaries &
,                      p, ph, timestep                            &
,                      th, q, qse, tracer, precip                 &
,                      dthbydt, dqbydt, dtrabydt                  &
,                      rain, snow ,down_flux                      &
    )

USE earth_constants_mod, ONLY: g

USE atmos_constants_mod, ONLY: r, cp, c_virtual
USE water_constants_mod, ONLY: lc, lf
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Description ;
!   Calculation of unsaturated Down draughts using Emanuel code.
!   Note at present the 4A convectin code deals in specific humidity
!   etc but the Emanuel code says it requires mixing ratios.
!   William Ingram appears to have taken the 4.3 version of the Emanuel
!   code and imported this into the old dynamics.
!   William Ingram chose to set various constants in a way consistent
!   with the UM. These choices mean that Emanuel's original choice
!   of a varying latent heat lv=lc+(cpv-cl)*(T-273.15)is reduced
!   to lv=lc.
!
!   Called by DEEP_CONV.
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
  n_dp                 & ! No. of deep convection points
, nlev                 & ! No. of model layers
, trlev                & ! No. of tracer levels
, ntra                 & ! No. of tracer fields
, kmax_term            & ! highest level reached by convection
, kterm(n_dp)            ! level reached by convection


LOGICAL, INTENT(IN) :: &
    & l_tracer           ! true in tracers present

REAL, INTENT(IN)    ::                &
  exner_layer_centres(n_dp,0:nlev)    & ! Exner
, exner_layer_boundaries(n_dp,0:nlev) & ! Exner at half level above
, p(n_dp,0:nlev)                      & ! Pressure  (Pa)
, ph(n_dp,0:nlev)                     & ! Pressure at half level (Pa)
, timestep                            & ! timestep for model physics in seconds
, qse(n_dp,nlev)                      & ! Saturation specific humidity of cloud
                                        ! environment (kg/kg)
, q (n_dp,nlev)                       & ! specific humidity of cloud 
                                        ! environment(kg/kg)
, th (n_dp,nlev)                      & ! theta of cloud environment(K)
, precip(n_dp,nlev)                   & ! Precip from updraught G-R scheme
                                        ! (kg/m2/s) not units for wdtrain
, tracer(n_dp,trlev,ntra)               !  Tracer on model levels  (kg/kg)


REAL, INTENT(INOUT)    :: &
  rain(n_dp)              & ! rainfall at surface (kg/m**2/s)
, snow(n_dp)                ! snowfall at surface (kg/m**2/s)

! increments
REAL, INTENT(INOUT)    ::  &
  dqbydt(n_dp,nlev)        & ! increments to q (kg/kg/s)
, dthbydt(n_dp,nlev)       & ! increments to potential temperature (K/s)
, dtrabydt(n_dp,nlev,ntra)   ! increments to model tracers (kg/kg/s)


REAL, INTENT(OUT)    ::    &
 down_flux(n_dp,nlev)        ! Downdraught mass flux (Pa/s) diagnostic



!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

REAL ::                    &
 wdtrain(n_dp)             & ! detrained precip in level k (kg/m/s3)
,wt(n_dp,nlev)             & ! fall speed for precipitation in (kg/m/s3)
,evap(n_dp,nlev)           & ! evaporation from Downdraught (/s)
,water(n_dp,nlev)          & ! rainwater in Downdraught (lp in paper)
,t(n_dp,nlev)              & ! temperature on model levels (K)
,h(n_dp,nlev)              & ! static energy (J ?)
,dpinv(n_dp)               & ! 1/dp
,gz(n_dp,nlev)             & ! gz a form of height measured from level 1
,dqbydt_dd1(n_dp,nlev)     & ! increments to q (kg/kg/s)
,dqbydt_dd2(n_dp,nlev)     & ! increments to q (kg/kg/s)
,dqbydt_dd(n_dp,nlev)      & ! increments to q (kg/kg/s)
,md(n_dp,nlev)             & ! Downdraught mass flux Emanuel scheme
                             ! (kg/m2/s not Pa/s)
,qp_down (n_dp,nlev)       & ! specific humidity of downdraught (kg/kg)
,trap_down(n_dp,nlev,ntra) & ! tracer values of downdraught (kg/kg)
,lv (n_dp,nlev)              ! latent heat release on evaporation

REAL ::       &
 qsm          & ! mean q
,afac         &
,b6,c6        &
,revap        &
,dhdp         & ! dz/dp
,rat          &
,ginv, lvcp   &
,tvx,tvy, dp  &
,coeff,fac,qstm

REAL   ::          &
 qsum(n_dp)        & ! qsum  conservation check
,qsum2(n_dp)       & ! qsum
,qsump(n_dp)       & ! qsum positive increments
,qsumn(n_dp)       & ! qsum negative increments
,qloss(n_dp)       & ! q lossed from atmosphere
,cfactor(n_dp)       ! correction factor

INTEGER ::         &
 i,j,k             & ! loop counters
, jtt(n_dp)        & ! location of 0.949*p(1)
, kphase(n_dp)       ! level of phase change from snow to rain

! parameters for Emanuel scheme

REAL, PARAMETER :: &
  coeffr = 1.0     & ! coefficient govering rate of evaporation of rain
, coeffs = 0.8     & ! coefficient govering rate of evaporation of snow
, omtrain = 50.0   & ! the assumed fall speed of rain (kg/m/s3 = Pa/s)
, omtsnow =  5.5   & ! the assumed fall speed of snow (kg/m/s3 = Pa/s)
, sigd = 0.05      & ! fractional area covered by unsaturated DD
, sigs = 0.12        ! fraction of precipitation falling outside cloud

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Model constants


!-----------------------------------------------------------------------
! 1.0 Initialisation of arrays
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('EMAN_DD',zhook_in,zhook_handle)

lvcp = lc/cp
ginv = 1./g


DO i= 1,n_dp
  qsum(i) = 0.0
  qsum2(i) = 0.0
  qsump(i) = 0.0
  qsumn(i) = 0.0
END DO

DO k=1,nlev
  DO i= 1,n_dp

! initialise arrays

    md(i,k) = 0.0
    wt(i,k) = 0.0
    evap(i,k) = 0.0
    water(i,k) = 0.0

    dqbydt_dd(i,k)  = 0.0
    dqbydt_dd1(i,k)  = 0.0
    dqbydt_dd2(i,k)  = 0.0

    jtt(i)      = 2
    kphase(i)   = 0

! convert from theta to T

    t(i,k)  = th(i,k) *exner_layer_centres(i,k)

    IF (t(i,k) >  273.15) THEN  ! inconsistent T freeze
      lv(i,k) = lc
    ELSE
      lv(i,k) = lc+lf
    END IF
  END DO
END DO

! Appears to be deriving some form of height using

!      p=rhoRT and dp/dz=-rhog    =>   gdz = -dp RT/p

!  equations here are using Tv  not T
!  Taking lowest model T level as zero height.

k = 1
  DO i= 1,n_dp
    gz(i,k)  = 0.0
    h(i,k)  = cp*t(i,k)+gz(i,k)
    qp_down(i,k) = q(i,k)
  END DO
DO k=2,nlev
  DO i= 1,n_dp
    tvx=t(i,k)  *(1.+c_virtual*q(i,k))
    tvy=t(i,k-1)*(1.+c_virtual*q(i,k-1))

    gz(i,k)=gz(i,k-1)+0.5*r*(tvx+tvy)*(p(i,k-1)-p(i,k))/ph(i,k-1)

    h(i,k)  = cp*t(i,k)+gz(i,k)
    qp_down(i,k) = q(i,k-1)
  END DO
END DO


! initialise tracer in down draught

IF (l_tracer) THEN
  DO j=1,ntra
    k=1
      DO i=1,n_dp
        trap_down(i,k,j) = tracer(i,k,j)
      END DO
    DO k=2,nlev
      DO i=1,n_dp
        trap_down(i,k,j) = tracer(i,k-1,j)
      END DO
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! 2.0 Main Down draught calculation
!-----------------------------------------------------------------------
! Assuming the array precip on model levels from G-R scheme corresponds
! to the Emanuel variable wdtrain - detrained precipitation.
!  Note In Emanuel scheme
!  wdtrain(i,k)=G*EP(I,k)*M(I,k)*CLW(I,k)
!  ep(i,k) - precip efficiency
!  M(i,k)  - up draught mass flux
!  clw(i,k) - condensed cloud water
!-----------------------------------------------------------------------

DO k=kmax_term+1,1,-1          ! level loop working Downwards

  DO i=1,n_dp           ! grid point loop

! start at level above termination
    IF (k <= kterm(i)+1.AND.kterm(i) /= 0) THEN
!-----------------------------------------------------------------------
!  Does Down draughts if precipitating efficiency >= 0.0001 at
!  top if cloud (Emanuel) => precip>0.0
!-----------------------------------------------------------------------
      IF (precip(i,kterm(i)) >  0.0) THEN  ! test for Downdraught

! detrained precipitation for level k


        wdtrain(i) = precip(i,k)*g

! Find rain water and evaporation using provisional estimates of qp

! Solution of equation (9) In Emanuel's paper to find lp(k) rain water
! at level k. Note as working Down from highest level know lp(k+1)

        ! Value of terminal velocity and coeffecient of evaporation for rain

        IF (t(i,k) >  273.15) THEN
          coeff=coeffr
          wt(i,k)=omtrain

        ELSE

        ! Value of terminal velocity and coeffecient of evaporation for snow

          coeff=coeffs
          wt(i,k)=omtsnow

        END IF

        IF (t(i,k) >  273.15.AND.t(i,k+1) <= 273.15) THEN
          kphase(i) = k       ! melting level
        END IF

! qsm - temporary estimate for water vapour of parcel at level k

        qsm=0.5*(q(i,k)+qp_down(i,k+1))


! Looks like expression for evaporation without sqrt (lp(k))

        afac=coeff*ph(i,k)*0.01*(qse(i,k)-qsm)                    &
                    /(1.0e4+2.0e3*0.01*ph(i,k)*qse(i,k))
        afac=MAX(afac,0.0)
        b6=(ph(i,k)-ph(i,k+1))*sigs*afac/wt(i,k)

        c6=(water(i,k+1)*wt(i,k+1)+wdtrain(i)/sigd)/wt(i,k)

        IF (c6 == 0.0) THEN
          revap =0.0     ! set because of numeric problems
        ELSE
          revap=0.5*(-b6+SQRT(b6*b6+4.*c6))
        END IF
        IF (revap <  0.0) THEN
          revap =0.0        ! reset
        END IF

        evap(i,k)=sigs*afac*revap

! water - rain water (lp in Emanuel paper)

        water(i,k)=revap*revap


! Calculate precipitating Downdraught mass flux under hydrostatic approx

        IF(k >  1) THEN
          dhdp=(h(i,k)-h(i,k-1))/(p(i,k-1)-p(i,k))
          dhdp=MAX(dhdp,0.1)      ! correct units
          md(i,k)=ginv*lv(i,k)*sigd*evap(i,k)/dhdp
          md(i,k)=MAX(md(i,k),0.0)


! Add a small amount of inertia to Downdraught
! (not referred to in paper?)

          fac=20.0*100./(ph(i,k-2)-ph(i,k-1))

          md(i,k)=(fac*md(i,k+1)+md(i,k))/(1.+fac)


! Force Downdraught mass flux to decrease linearly to zero between about
! 950mb and the surface.
! Actually coded as decreasing linearly to zero between
! 0.949*lowest level pressure and model T level (Is this what we really
! want it to do ?) (Did the original Emanuel code assume level 1 was
! the surface ?)


          IF(p(i,k) >  (0.949*p(i,1)))THEN
            jtt(i)=MAX(jtt(i),k)
            md(i,k)=md(i,jtt(i))*(p(i,1)-p(i,k))/(p(i,1)-p(i,jtt(i)))

          END IF

        END IF   ! k>1

! Find mixing ratio of precipitation Downdraught

        IF(k == 1)THEN
          qstm=qse(i,k)
        ELSE
          qstm=qse(i,k-1)
        END IF
        IF(md(i,k) >  md(i,k+1))THEN
          rat=md(i,k+1)/md(i,k)
          qp_down(i,k)=qp_down(i,k+1)*rat+q(i,k)*(1.0-rat)+ginv*   &
                       sigd*(ph(i,k-1)-ph(i,k))*(evap(i,k)/md(i,k))

        ELSE
          IF(md(i,k+1) >  0.0)THEN
            qp_down(i,k)=(gz(i,k+1)-gz(i,k)+qp_down(i,k+1)*lv(i,k) &
                             +cp*(t(i,k+1)-t(i,k)))/lv(i,k)
          END IF
        END IF
!  k+1 md(i,k)= 0.0 therefore qp(i,1) remains unchanged
        qp_down(i,k)=MIN(qp_down(i,k),qstm)
        qp_down(i,k)=MAX(qp_down(i,k),0.0)

      END IF
    END IF
  END DO   ! loop over points
END DO     ! loop over levels


! calculate tracers in downdraught plume

IF (l_tracer) THEN
  DO j=1,ntra
    DO k=kmax_term+1,1,-1   ! level loop working Downwards
      DO i=1,n_dp           ! grid point loop

        IF (k <= kterm(i)+1.AND.kterm(i) /= 0) THEN

          IF(md(i,k) >  md(i,k+1))THEN
            rat=md(i,k+1)/md(i,k)

            trap_down(i,k,j)=trap_down(i,k+1,j)*rat+tracer(i,k,j)*(1.-rat)
          ELSE
            IF(md(i,k+1) >  0.0)THEN
              trap_down(i,k,j)=trap_down(i,k+1,j)
            END IF
          END IF

        END IF

      END DO   ! loop over points
    END DO     ! loop over levels
  END DO     ! loop over tracers

END IF      ! test on tracers

! water corresponds to lp in paper
! Emanuel scheme assumes all surface precip rain. Note wt depends on
! whether rain or snow should we use this information?

!  now in mm/s = kg/m2/s  (wt in Pa/s)

DO i=1,n_dp
  IF (wt(i,1) == omtrain) THEN      ! falling as rain
    rain(i) = rain(i) + wt(i,1)*sigd*water(i,1)/g
  ELSE     ! falling as snow
    snow(i) = snow(i) + wt(i,1)*sigd*water(i,1)/g
  END IF
END DO


! ----------------------------------------------------------------------
! Calculation of increments
! ----------------------------------------------------------------------

!  level 1

k=1
DO i=1,n_dp
  IF (kterm(i) /= 0) THEN
    dpinv(i) = 1.0/(ph(i,0)-ph(i,1))


! increments to temperature    dtheta = dt/exner
! evaporation of precipitation

    dthbydt(i,k) = dthbydt(i,k) - (lv(i,k)/cp)*sigd*evap(i,1)       &
                      /exner_layer_centres(i,1)

! increments to q
! change in q due to transport by Down draught

    dqbydt_dd2(i,k)=dqbydt_dd2(i,k)+g*md(i,2)*(qp_down(i,2)-q(i,1))*dpinv(i)


! increase in q due to evaporation of precip

    dqbydt_dd1(i,k) = dqbydt_dd1(i,k)+sigd*evap(i,1)

! total

    dqbydt_dd(i,k) = dqbydt_dd1(i,k)+dqbydt_dd2(i,k)

  END IF
END DO      ! loop over gridpoints


!  level where phase changes add lf*water falling

DO i=1,n_dp

  IF (kphase(i) /= 0) THEN

! evaporation of precipitation

    dthbydt(i,kphase(i)) = dthbydt(i,kphase(i))                            &
                           - lf*water(i,kphase(i))/cp/timestep             &
                               /exner_layer_centres(i,kphase(i))
  END IF
END DO      ! loop over gridpoints

!  Levels above k=1

DO k=2,kmax_term+1        ! level loop
  DO i=1,n_dp            ! grid point loop

    IF (k <= kterm(i)+1.AND.kterm(i) /= 0) THEN

      dpinv(i)=1.0/(ph(i,k-1)-ph(i,k))

! Temperature increments - evaporation of precipitation

      dthbydt(i,k) = dthbydt(i,k) - (lv(i,k)/cp)*sigd*evap(i,k)          &
                      /exner_layer_centres(i,k)

! change in q due to transport by Down draught

      dqbydt_dd2(i,k) = dqbydt_dd2(i,k) +                                &
                        g*dpinv(i)*(md(i,k+1)*(qp_down(i,k+1)-q(i,k))    &
                                   -md(i,k)*(qp_down(i,k)-q(i,k-1)))


! increase in q due to evaporation of precip

      dqbydt_dd1(i,k) = dqbydt_dd1(i,k) + sigd*evap(i,k)
      dqbydt_dd(i,k) = dqbydt_dd1(i,k) +  dqbydt_dd2(i,k)


    END IF           ! test on whether level in convecting column
  END DO             ! gridpoint loop
END DO               ! level loop

IF (l_tracer) THEN
  DO j = 1, ntra
    k = 1
      DO i=1,n_dp
        IF (kterm(i) /= 0) THEN
         dpinv(i) = 1.0/(ph(i,1)-ph(i,2))
         dtrabydt(i,k,j) = dtrabydt(i,k,j)                                    &
                     +  g*dpinv(i)*md(i,k+1)*(trap_down(i,k+1,j)-tracer(i,k,j))
        END IF
      END DO
  END DO

  DO j = 1, ntra
    DO k= 2, kmax_term+1
      DO i=1,n_dp
        IF (k >= kterm(i).AND.kterm(i) /= 0) THEN

          dpinv(i)=1.0/(ph(i,k)-ph(i,k+1))
          dtrabydt(i,k,j) = dtrabydt(i,k,j)                                   &
                   +g*dpinv(i)*(md(i,k+1)*(trap_down(i,k+1,j)-tracer(i,k,j))- &
                                md(i,k)*(trap_down(i,k,j)-tracer(i,k-1,j)))

        END IF       ! test on whether level in convecting column
      END DO        ! gridpoint loop
    END DO           ! level loop
  END DO             ! tracer loop

END IF

! Conservation checks

DO k= 1, nlev
  DO i=1,n_dp
     dp=(ph(i,k-1)-ph(i,k))       ! correct for layer k
     qsum(i)=qsum(i) + dqbydt_dd(i,k)*dp/g
     qsum2(i)=qsum2(i) + precip(i,k)
     IF (dqbydt_dd(i,k) >  0.0) THEN
       qsump(i) = qsump(i) + dqbydt_dd(i,k)*dp/g
     ELSE
       qsumn(i) = qsumn(i) + dqbydt_dd(i,k)*dp/g
     END IF

! change units of Emanuel Downdraught mass flux for diagnostic output

     down_flux(i,k) = g*md(i,k)

  END DO           ! gridpoint loop
END DO             ! level

! Need to correct dq increments as moisture not conserved and non
! conservation is systematic leading to loss of q from atmosphere.

DO i = 1, n_dp
  qloss(i) = qsum2(i)-rain(i)-qsum(i)-snow(i)
  IF (qloss(i)  == 0.0) THEN
    cfactor(i) = 1.0
  ELSE
    IF (qsump(i) /= 0.0) THEN
      cfactor(i) = 1.0 + qloss(i)/qsump(i)
    ELSE
      IF (qsumn(i) == 0.0) THEN
! No increments to q but no moisture balance ie qsum2 not = rain/snow
! Cannot correct increments therefore correct rain
        cfactor(i) = 1.0
        IF (rain(i) >  0.0) THEN
          rain(i) = qsum2(i)
        ELSE  ! assume precip snow
          snow(i) = qsum2(i)
        END IF
      ELSE
        cfactor(i) = 1.0 - qloss(i)/qsumn(i)
      END IF
    END IF
  END IF
END DO

! Scale positive increments to increase them

DO k= 1, nlev
  DO i=1,n_dp

    IF (qsump(i) /= 0.0) THEN
      IF (dqbydt_dd(i,k) >  0.0) THEN
        dqbydt_dd(i,k) = dqbydt_dd(i,k)*cfactor(i)
      END IF
    ELSE       ! reduce negative increments
      IF (dqbydt_dd(i,k) <  0.0) THEN
        dqbydt_dd(i,k) = dqbydt_dd(i,k)*cfactor(i)
      END IF
    END IF

    dqbydt(i,k) = dqbydt(i,k) + dqbydt_dd(i,k)

  END DO           ! gridpoint loop
END DO             ! level

! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('EMAN_DD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eman_dd

